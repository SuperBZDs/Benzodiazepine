# -*- coding: utf-8 -*-
"""最小自检：确认自己实现的 t 检验和 BH 校正没有算错。

跑法：python tests/test_stats.py
"""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'src'))

from stats import t_two_sided_p, welch_ttest, bh_fdr


def close(a, b, tol=1e-4):
    return abs(a - b) < tol


def main():
    # t 分布双侧 p 值的已知值（查表或由标准库计算得到）
    assert close(t_two_sided_p(1.0, 1), 0.5), 'df=1, t=1 应为 0.5'
    assert close(t_two_sided_p(2.0, 10), 0.07339, 1e-5), 'df=10, t=2 应为 0.07339'
    assert close(t_two_sided_p(1.96, 100), 0.05278, 1e-5), 'df=100, t=1.96 应为 0.05278'

    # 两组各相差 1 的等距数据，Welch t 应为 -1，自由度 8
    x = [1.0, 2, 3, 4, 5]
    y = [2.0, 3, 4, 5, 6]
    t, df, p = welch_ttest(x, y)
    assert close(t, -1.0, 1e-9) and close(df, 8.0, 1e-9), 'Welch 统计量不对'
    assert close(p, 0.3466, 1e-4), 'Welch p 值不对'

    # BH：所有 p 相同时，校正后 q 应等于最大 p
    q = bh_fdr([0.01, 0.02, 0.03, 0.04])
    assert all(close(v, 0.04, 1e-9) for v in q), 'BH 校正不对'

    print('全部通过')


if __name__ == '__main__':
    main()
