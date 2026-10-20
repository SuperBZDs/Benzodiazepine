# -*- coding: utf-8 -*-
"""差异表达统计：Welch t 检验 + Benjamini-Hochberg 校正。

只依赖 numpy，不依赖 scipy，避免装不上科学计算包时跑不动。
"""
import math
import numpy as np


def _betacf(a, b, x, itmax=200, eps=3e-12):
    """不完全贝塔函数的连分式（Lentz 算法）。"""
    qab, qap, qam = a + b, a + 1.0, a - 1.0
    c = 1.0
    d = 1.0 - qab * x / qap
    if abs(d) < 1e-30:
        d = 1e-30
    d = 1.0 / d
    h = d
    for m in range(1, itmax + 1):
        m2 = 2 * m
        aa = m * (b - m) * x / ((qam + m2) * (a + m2))
        d = 1.0 + aa * d
        if abs(d) < 1e-30:
            d = 1e-30
        c = 1.0 + aa / c
        if abs(c) < 1e-30:
            c = 1e-30
        d = 1.0 / d
        h *= d * c
        aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2))
        d = 1.0 + aa * d
        if abs(d) < 1e-30:
            d = 1e-30
        c = 1.0 + aa / c
        if abs(c) < 1e-30:
            c = 1e-30
        d = 1.0 / d
        delta = d * c
        h *= delta
        if abs(delta - 1.0) < eps:
            break
    return h


def _betai(a, b, x):
    """正则化不完全贝塔函数 I_x(a, b)。"""
    if x <= 0.0:
        return 0.0
    if x >= 1.0:
        return 1.0
    ln_bt = (math.lgamma(a + b) - math.lgamma(a) - math.lgamma(b)
             + a * math.log(x) + b * math.log(1.0 - x))
    bt = math.exp(ln_bt)
    if x < (a + 1.0) / (a + b + 2.0):
        return bt * _betacf(a, b, x) / a
    return 1.0 - bt * _betacf(b, a, 1.0 - x) / b


def t_two_sided_p(t, df):
    """t 分布双侧 p 值。"""
    if df <= 0 or not np.isfinite(t):
        return 1.0
    return _betai(df / 2.0, 0.5, df / (df + t * t))


def welch_ttest(x, y):
    """两样本 Welch t 检验，返回 (t, df, p)。

    x、y 是两个分组的表达值，长度可以不等，允许缺失（由调用方先处理）。
    """
    nx, ny = len(x), len(y)
    if nx < 2 or ny < 2:
        return np.nan, np.nan, 1.0
    vx, vy = np.var(x, ddof=1), np.var(y, ddof=1)
    se2 = vx / nx + vy / ny
    if se2 <= 0:
        return np.nan, np.nan, 1.0
    t = (np.mean(x) - np.mean(y)) / math.sqrt(se2)
    df = se2 ** 2 / ((vx / nx) ** 2 / (nx - 1) + (vy / ny) ** 2 / (ny - 1))
    return t, df, t_two_sided_p(t, df)


def bh_fdr(p):
    """Benjamini-Hochberg 校正，返回与输入同序的 q 值。"""
    p = np.asarray(p, dtype=float)
    n = len(p)
    order = np.argsort(p)
    ranked = p[order]
    q = ranked * n / (np.arange(n) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    q = np.clip(q, 0, 1)
    out = np.empty(n)
    out[order] = q
    return out
