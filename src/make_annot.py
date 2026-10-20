# -*- coding: utf-8 -*-
"""把 GEO 的 GPL 注释文件（.annot.gz）转成两列 TSV：探针 ID、基因符号。

用法：python src/make_annot.py data/GPL96.annot.gz data/GPL96_annot.tsv
"""
import csv
import gzip
import sys


def convert(src, dst):
    rows = 0
    with gzip.open(src, 'rt', encoding='utf-8', errors='replace') as fh, \
            open(dst, 'w', encoding='utf-8', newline='') as out:
        w = csv.writer(out, delimiter='\t')
        w.writerow(['ID', 'Gene Symbol'])
        header = None
        for line in fh:
            if header is None:
                if line.startswith('ID\t'):
                    header = line.rstrip('\n').split('\t')
                continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) < len(header):
                parts += [''] * (len(header) - len(parts))
            rec = dict(zip(header, parts))
            pid = (rec.get('ID') or '').strip()
            sym = (rec.get('Gene symbol') or '').strip()
            if pid and sym not in ('', '---', 'NA'):
                w.writerow([pid, sym])
                rows += 1
    print('写出 %d 条映射：%s' % (rows, dst))


if __name__ == '__main__':
    if len(sys.argv) < 3:
        raise SystemExit('用法：python src/make_annot.py data/GPL96.annot.gz data/GPL96_annot.tsv')
    convert(sys.argv[1], sys.argv[2])
