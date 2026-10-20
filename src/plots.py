# -*- coding: utf-8 -*-
"""出图：火山图与热图。坐标轴用英文，避免不同系统缺中文字体。"""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

NAVY = '#1F3864'
GREY = '#BDBDBD'
LIGHT = '#9DC3E6'


def volcano(df, out_path, fdr=0.05, lfc=1.0, title=''):
    fig, ax = plt.subplots(figsize=(6.2, 5.0), dpi=300)
    sig = (df['fdr'] < fdr) & (df['log2fc'].abs() >= lfc)
    ax.scatter(df.loc[~sig, 'log2fc'], -np.log10(df.loc[~sig, 'fdr'].clip(lower=1e-300)),
               s=6, c=GREY, alpha=0.55, linewidths=0)
    ax.scatter(df.loc[sig, 'log2fc'], -np.log10(df.loc[sig, 'fdr'].clip(lower=1e-300)),
               s=8, c=NAVY, alpha=0.85, linewidths=0)
    ax.axhline(-np.log10(fdr), color=LIGHT, lw=0.8, ls='--')
    ax.axvline(-lfc, color=LIGHT, lw=0.8, ls='--')
    ax.axvline(lfc, color=LIGHT, lw=0.8, ls='--')
    ax.set_xlabel('log2 fold change (tumor / normal)', fontsize=9)
    ax.set_ylabel('-log10 FDR', fontsize=9)
    ax.set_title(title or 'Differential expression', fontsize=10, loc='left')
    for s in ('top', 'right'):
        ax.spines[s].set_visible(False)
    ax.tick_params(labelsize=8)
    fig.tight_layout()
    fig.savefig(out_path)
    plt.close(fig)


def heatmap(expr, meta, genes, out_path, title=''):
    """genes：要画的基因（行），按 z-score 画。"""
    sub = expr.loc[genes]
    z = sub.sub(sub.mean(axis=1), axis=0).div(sub.std(axis=1).replace(0, np.nan), axis=0)
    z = z.fillna(0)
    order = meta.loc[z.columns].sort_values('group').index
    z = z[order]
    fig, ax = plt.subplots(figsize=(7.2, 4.6), dpi=300)
    im = ax.imshow(z.values, aspect='auto', cmap='RdBu_r', vmin=-2, vmax=2,
                   interpolation='nearest')
    ax.set_yticks(range(len(genes)))
    ax.set_yticklabels(genes, fontsize=6.5)
    ax.set_xticks([])
    groups = meta.loc[order, 'group'].values
    for i in range(1, len(groups)):
        if groups[i] != groups[i - 1]:
            ax.axvline(i - 0.5, color='black', lw=1.2)
    ax.set_xlabel('samples (left: %s, right: %s)' % (groups[0], groups[-1]), fontsize=8)
    ax.set_title(title or 'Top genes (row z-score)', fontsize=10, loc='left')
    cb = fig.colorbar(im, ax=ax, fraction=0.025, pad=0.01)
    cb.ax.tick_params(labelsize=7)
    fig.tight_layout()
    fig.savefig(out_path)
    plt.close(fig)
