# -*- coding: utf-8 -*-
"""读取 GEO series matrix 文件，或直接读 CSV；并给样本分组。"""
import gzip
import io
import os
import re

import numpy as np
import pandas as pd

# 分组关键词。标题或来源里命中就归到对应组，大小写不敏感。
TUMOR_KEYS = ('tumor', 'tumour', 'cancer', 'carcinoma', 'adenocarcinoma', 'malignant')
NORMAL_KEYS = ('normal', 'non-tumor', 'nontumor', 'adjacent', 'control', 'healthy')


def _open_text(path):
    if path.endswith('.gz'):
        return gzip.open(path, 'rt', encoding='utf-8', errors='replace')
    return open(path, 'r', encoding='utf-8', errors='replace')


def read_series_matrix(path):
    """返回 (表达矩阵 DataFrame, 样本信息 DataFrame)。

    表达矩阵：行是探针，列是 GSM 号。
    样本信息：index 是 GSM 号，列至少含 title、source、group。
    """
    titles, sources, chars, gsms = [], [], [], []
    table_lines = []
    in_table = False

    with _open_text(path) as fh:
        for line in fh:
            line = line.rstrip('\n')
            if line.startswith('!series_matrix_table_begin'):
                in_table = True
                continue
            if line.startswith('!series_matrix_table_end'):
                in_table = False
                continue
            if in_table:
                table_lines.append(line)
                continue
            if line.startswith('!Sample_geo_accession'):
                gsms = re.findall(r'"([^"]*)"', line)
            elif line.startswith('!Sample_title'):
                titles = re.findall(r'"([^"]*)"', line)
            elif line.startswith('!Sample_source_name_ch1'):
                sources = re.findall(r'"([^"]*)"', line)
            elif line.startswith('!Sample_characteristics_ch1'):
                chars = re.findall(r'"([^"]*)"', line)

    if not table_lines:
        raise ValueError('没找到 !series_matrix_table_begin，确认下载的是 series matrix 文件')

    expr = pd.read_csv(io.StringIO('\n'.join(table_lines)), sep='\t', quotechar='"')
    expr = expr.rename(columns={expr.columns[0]: 'ID_REF'}).set_index('ID_REF')
    expr = expr.apply(pd.to_numeric, errors='coerce')

    if not gsms:
        gsms = list(expr.columns)
    meta = pd.DataFrame({'title': titles[:len(gsms)] if titles else [''] * len(gsms),
                         'source': sources[:len(gsms)] if sources else [''] * len(gsms),
                         'chars': chars[:len(gsms)] if chars else [''] * len(gsms)},
                        index=gsms)
    meta['group'] = [guess_group(row) for _, row in meta.iterrows()]
    return expr, meta


def guess_group(row):
    text = ' '.join(str(v).lower() for v in row.values)
    if any(k in text for k in NORMAL_KEYS):
        return 'normal'
    if any(k in text for k in TUMOR_KEYS):
        return 'tumor'
    return 'unknown'


def read_csv_input(expr_path, group_path):
    """读自定义 CSV：表达矩阵（行=基因，列=样本）+ 分组成员表（两列：sample, group）。"""
    expr = pd.read_csv(expr_path, index_col=0)
    grp = pd.read_csv(group_path)
    grp.columns = ['sample', 'group']
    grp = grp.set_index('sample')
    expr = expr[[c for c in expr.columns if c in grp.index]]
    meta = grp.loc[expr.columns].copy()
    meta['title'] = meta.index
    meta['source'] = ''
    meta['chars'] = ''
    return expr, meta


def maybe_log2(expr):
    """判断是否需要 log2 转换：最大值很大、且存在非整数时认为还是线性值。"""
    safe = expr.where(expr > 0)
    if np.nanmax(safe.values) > 50 and (safe.stack() % 1 != 0).mean() > 0.5:
        return np.log2(safe + 1.0), True
    return expr, False


def load_annotation(path):
    """可选：GPL 注释表（两列：探针 ID、基因符号），返回 Series。"""
    if not path or not os.path.exists(path):
        return None
    ann = pd.read_csv(path, sep='\t', dtype=str)
    cols = list(ann.columns)
    probe_col = cols[0]
    sym_col = next((c for c in cols if c.lower() in ('gene symbol', 'symbol', 'gene_symbol', 'genesymbol')), None)
    if sym_col is None:
        return None
    ann = ann[[probe_col, sym_col]].dropna()
    ann = ann[ann[sym_col].str.strip() != '']
    ann = ann.drop_duplicates(subset=[probe_col])
    return ann.set_index(probe_col)[sym_col]


def demo_matrix(n_genes=2000, n_per_group=20, n_de=200, seed=7):
    """生成一份模拟数据，用来验证流程能不能跑通（结果不能当真实结论用）。"""
    rng = np.random.default_rng(seed)
    samples = ['DEMO_%02d' % (i + 1) for i in range(n_per_group * 2)]
    genes = ['GENE%04d' % (i + 1) for i in range(n_genes)]
    base = rng.normal(8, 1.5, size=(n_genes, 1))
    noise = rng.normal(0, 0.45, size=(n_genes, n_per_group * 2))
    expr = base + noise
    de_idx = rng.choice(n_genes, n_de, replace=False)
    direction = rng.choice([-1, 1], n_de)
    expr[de_idx, n_per_group:] += (direction * 1.6)[:, None]
    expr = pd.DataFrame(np.clip(expr, 0, None), index=genes, columns=samples)
    meta = pd.DataFrame({'title': samples,
                         'source': ['demo'] * len(samples),
                         'chars': [''] * len(samples),
                         'group': ['normal'] * n_per_group + ['tumor'] * n_per_group},
                        index=samples)
    return expr, meta
