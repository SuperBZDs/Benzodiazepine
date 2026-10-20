# -*- coding: utf-8 -*-
"""一键跑通：读数据 → 差异分析 → 出图 → 写结果表。

用法：
  python src/run_all.py --demo
  python src/run_all.py --input data/GSE10072_series_matrix.txt.gz
  python src/run_all.py --input data/expr.csv --groups data/groups.csv
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import geo_io
import plots
import stats


def parse_args():
    p = argparse.ArgumentParser(description='GEO 数据差异表达分析')
    p.add_argument('--input', help='series matrix 文件（.txt 或 .txt.gz）')
    p.add_argument('--groups', help='自定义分组表（两列：sample, group），配合 CSV 输入')
    p.add_argument('--annot', help='可选 GPL 注释表（探针 → 基因符号）', default=None)
    p.add_argument('--demo', action='store_true', help='用模拟数据验证流程')
    p.add_argument('--outdir', default='results')
    p.add_argument('--fdr', type=float, default=0.05)
    p.add_argument('--lfc', type=float, default=1.0)
    p.add_argument('--top', type=int, default=30)
    p.add_argument('--collapse', choices=['max_mean', 'none'], default='max_mean',
                   help='同一基因多个探针时的合并方式：取平均表达最高的探针（默认）或不合并')
    return p.parse_args()


def collapse_to_gene(res, expr, how='max_mean'):
    """把探针层面的结果合并到基因层面。

    默认取「平均表达最高的那个探针」——GEO 里同一基因常有多个探针，
    直接取一个代表探针是最常见也最容易解释的做法；同时在 n_probes 列记下原有几个探针。
    """
    means = []
    for probe in res['probe']:
        means.append(float(np.nanmean(expr.loc[probe].values)) if probe in expr.index else np.nan)
    res = res.copy()
    res['mean_expr'] = means
    counts = res.groupby('symbol')['probe'].transform('count')
    res['n_probes'] = counts
    if how == 'max_mean':
        res = (res.sort_values(['mean_expr', 'p'], ascending=[False, True])
                  .drop_duplicates(subset='symbol', keep='first'))
    return res.sort_values('p').reset_index(drop=True)


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    if args.demo:
        expr, meta = geo_io.demo_matrix()
        source = 'demo（模拟数据，仅用于验证流程）'
        label = 'DEMO'
        converted = False
    elif args.input and args.groups:
        expr, meta = geo_io.read_csv_input(args.input, args.groups)
        source = os.path.basename(args.input)
        label = os.path.basename(args.input).split('_')[0].split('.')[0]
        expr, converted = geo_io.maybe_log2(expr)
    elif args.input:
        expr, meta = geo_io.read_series_matrix(args.input)
        source = os.path.basename(args.input)
        label = os.path.basename(args.input).split('_')[0].split('.')[0]
        expr, converted = geo_io.maybe_log2(expr)
    else:
        raise SystemExit('给一个 --input，或者用 --demo 先验证流程')

    meta = meta[meta['group'].isin(['tumor', 'normal'])]
    expr = expr[[c for c in expr.columns if c in meta.index]]
    n_t = int((meta['group'] == 'tumor').sum())
    n_n = int((meta['group'] == 'normal').sum())
    if n_t < 2 or n_n < 2:
        raise SystemExit('某一组样本少于 2 个：tumor=%d, normal=%d，检查分组关键词' % (n_t, n_n))

    keep = expr.notna().mean(axis=1) >= 0.8
    expr = expr.loc[keep]
    print('样本：tumor %d，normal %d；探针 %d' % (n_t, n_n, expr.shape[0]))

    tumor = expr[meta.index[meta['group'] == 'tumor']].values
    normal = expr[meta.index[meta['group'] == 'normal']].values

    rows = []
    for i, probe in enumerate(expr.index):
        t, df, p = stats.welch_ttest(tumor[i], normal[i])
        rows.append((probe,
                     float(np.nanmean(tumor[i]) - np.nanmean(normal[i])),
                     float(t), float(df), float(p)))
    res = pd.DataFrame(rows, columns=['probe', 'log2fc', 't', 'df', 'p'])
    res['fdr'] = stats.bh_fdr(res['p'].values)

    ann = geo_io.load_annotation(args.annot) if args.annot else None
    if ann is not None:
        res['symbol'] = res['probe'].map(ann)
        # 注释表里没有的探针退回用探针号做名字，避免空基因名被合并到一起
        res['symbol'] = res['symbol'].fillna(res['probe'])
    else:
        res['symbol'] = res['probe']

    res['direction'] = np.where(res['log2fc'] > 0, 'up', 'down')
    n_probe_rows = len(res)
    if args.collapse == 'max_mean' and ann is not None:
        res = collapse_to_gene(res, expr, how='max_mean')
        print('基因层面合并：%d 探针 → %d 个基因（每个基因取平均表达最高的探针）' % (n_probe_rows, len(res)))
    res = res.sort_values('p').reset_index(drop=True)
    res.to_csv(os.path.join(args.outdir, 'deg_results.csv'), index=False, float_format='%.6g')

    sig = res[(res['fdr'] < args.fdr) & (res['log2fc'].abs() >= args.lfc)]
    n_up = int((sig['direction'] == 'up').sum())
    n_down = int((sig['direction'] == 'down').sum())

    plots.volcano(res, os.path.join(args.outdir, 'volcano.png'),
                  fdr=args.fdr, lfc=args.lfc,
                  title='%s   up %d / down %d' % (label, n_up, n_down))

    pick = sig['probe'].head(args.top).tolist() or res['probe'].head(args.top).tolist()
    labels = res.set_index('probe').loc[pick, 'symbol'].astype(str).tolist()
    sub = expr.loc[pick]
    sub.index = labels
    plots.heatmap(sub, meta, labels, os.path.join(args.outdir, 'heatmap.png'),
                  title='Top %d genes by FDR' % len(pick))

    lines = [
        '输入：%s' % source,
        'log2 转换：%s' % ('是' if converted else '否（数据看起来已经是 log 尺度）'),
        '样本：tumor %d，normal %d；参与检验的探针 %d，合并后基因 %d' % (n_t, n_n, expr.shape[0], len(res)),
        '阈值：FDR < %.2f，|log2FC| >= %.2f' % (args.fdr, args.lfc),
        '显著：共 %d（上调 %d，下调 %d）' % (len(sig), n_up, n_down),
        '合并方式：%s' % ('每个基因取平均表达最高的探针' if args.collapse == 'max_mean' and ann is not None else '未合并（探针层面）'),
        '',
        'FDR 最小的 10 个：',
    ]
    for _, r in res.head(10).iterrows():
        lines.append('  %-14s %-10s log2FC=%+.2f  FDR=%.3g' %
                     (r['probe'], str(r['symbol'])[:10], r['log2fc'], r['fdr']))
    text = '\n'.join(lines)
    with open(os.path.join(args.outdir, 'summary.txt'), 'w', encoding='utf-8') as fh:
        fh.write(text + '\n')
    print(text)
    print('\n输出目录：%s（deg_results.csv / volcano.png / heatmap.png / summary.txt）' % args.outdir)


if __name__ == '__main__':
    main()
