#!/usr/bin/env python3
"""Filter pixy output files for NA in key columns and gzip them.

Usage: filter_and_gzip_pixy.py <outdir> <prefix>
Example: filter_and_gzip_pixy.py results/pixy/sitelevel/Bv1_1_10000 pixy_site

This script looks for files matching <outdir>/<prefix>_*.txt and for each file
selects an appropriate column to filter on (e.g. avg_dxy, avg_pi, avg_hudson_fst,
tajima_d, avg_watterson_theta). Rows where that column is 'NA' or empty are removed.
The filtered output is written to the same path with a .gz suffix and the original
text file is removed.
"""
import sys
import os
import gzip
import glob

COL_MAP = {
    '_dxy.txt': 'avg_dxy',
    '_fst.txt': 'avg_hudson_fst',
    '_pi.txt': 'avg_pi',
    '_tajima_d.txt': 'tajima_d',
    '_watterson_theta.txt': 'avg_watterson_theta',
}


def choose_col_from_header(header_fields, filename):
    # Prefer explicit mapping by filename suffix
    for suffix, col in COL_MAP.items():
        if filename.endswith(suffix):
            if col in header_fields:
                return col
    # Fallback: pick the first column from COL_MAP values that appears in header
    for col in COL_MAP.values():
        if col in header_fields:
            return col
    return None


def filter_file(path):
    with open(path, 'r') as inf:
        header = inf.readline().rstrip('\n').split('\t')
        col = choose_col_from_header(header, os.path.basename(path))
        if col is None:
            print(f"[skip] no known filter column in {path}")
            return
        idx = header.index(col)
        gzpath = path + '.gz'
        with gzip.open(gzpath, 'wt') as outf:
            outf.write('\t'.join(header) + '\n')
            for line in inf:
                parts = line.rstrip('\n').split('\t')
                val = parts[idx] if idx < len(parts) else ''
                if val != 'NA' and val != '':
                    outf.write('\t'.join(parts) + '\n')
    # remove original file
    try:
        os.remove(path)
    except Exception as e:
        print(f"warning: could not remove {path}: {e}")
    print(f"filtered+gzipped: {gzpath}")


def main():
    if len(sys.argv) != 3:
        print(__doc__)
        sys.exit(2)
    outdir = sys.argv[1]
    prefix = sys.argv[2]
    pattern = os.path.join(outdir, f"{prefix}_*.txt")
    files = glob.glob(pattern)
    if not files:
        print(f"no files matching {pattern}")
        return
    for f in files:
        try:
            filter_file(f)
        except Exception as e:
            print(f"error filtering {f}: {e}")


if __name__ == '__main__':
    main()
