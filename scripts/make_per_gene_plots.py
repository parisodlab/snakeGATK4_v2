#!/usr/bin/env python3
import sys
import os
import csv
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from itertools import chain

def plot_per_gene_per_pop(populations, out_png):
    flattened = list(chain.from_iterable(populations))
    # each entry is (fixed_hetero, total_sites)
    var_sites = [g[1] for g in flattened]
    fixed_hetero = [g[0] for g in flattened]
    plt.figure(figsize=(6,6))
    plt.scatter(var_sites, fixed_hetero, s=6, alpha=0.6)
    plt.xlabel('Number of variable sites')
    plt.ylabel('Number of fixed heterozygotes')
    plt.tight_layout()
    plt.savefig(out_png)
    plt.close()

def main():
    if len(sys.argv) < 4:
        print('Usage: make_per_gene_plots.py <het_tsv> <out_dir> <prefix>')
        sys.exit(2)
    het_tsv = sys.argv[1]
    out_dir = sys.argv[2]
    prefix = sys.argv[3]
    os.makedirs(out_dir, exist_ok=True)

    with open(het_tsv) as fh:
        reader = csv.reader(fh, delimiter='\t')
        header = next(reader)
        pops = header[3:]
        data = {p: [] for p in pops}
        for row in reader:
            for i, p in enumerate(pops, start=3):
                val = row[i]
                try:
                    total, het = val.split(',')
                    total = int(total)
                    het = int(het)
                except Exception:
                    total, het = 0, 0
                data[p].append((het, total))

    out_png = os.path.join(out_dir, f"{prefix}_per_gene_per_pop.png")
    plot_per_gene_per_pop(list(data.values()), out_png)
    print('wrote', out_png)

if __name__ == '__main__':
    main()
