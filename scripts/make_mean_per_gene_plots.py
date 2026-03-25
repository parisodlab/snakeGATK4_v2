#!/usr/bin/env python3
import sys
import os
import csv
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from itertools import chain

def plot_mean_per_gene(populations, out_png):
    populations_num = len(populations)
    result = [[gene[0] / gene[1] if gene[1] != 0 else 0 for gene in population] for population in populations]
    result = [sum(values) / populations_num for values in zip(*result)]
    plt.figure(figsize=(6,4))
    plt.hist(result, bins=40)
    plt.xlabel('Mean fraction of fixed heterozygotes per gene')
    plt.ylabel('Number of genes')
    plt.tight_layout()
    plt.savefig(out_png)
    plt.close()

def main():
    if len(sys.argv) < 4:
        print('Usage: make_mean_per_gene_plots.py <het_tsv> <out_dir> <prefix>')
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

    out_png = os.path.join(out_dir, f"{prefix}_mean_per_gene.png")
    plot_mean_per_gene(list(data.values()), out_png)
    print('wrote', out_png)

if __name__ == '__main__':
    main()
