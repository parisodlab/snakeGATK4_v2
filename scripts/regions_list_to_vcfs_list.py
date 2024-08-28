#!/usr/bin/env python
import sys
import os


FN = sys.argv[1]
ZEROPAD = 6


with open(FN,'r') as fh:
    for i,l in enumerate(fh):
        f = l.strip().split()
        print(f"{VCFS_DIR}/pfilt-{i:06}.{f[0]}:{f[1]}-{f[2]}.vcf.gz")