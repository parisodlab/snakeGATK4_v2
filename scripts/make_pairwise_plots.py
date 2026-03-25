#!/usr/bin/env python3
import sys
import os
import csv
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def plot_pairwise(pop1, pop2, out_png):
    # pop1/pop2 are lists of (het_count, total_sites)
    vals1 = [ (g[0] / g[1]) if g[1] != 0 else 0 for g in pop1 ]
    vals2 = [ (g[0] / g[1]) if g[1] != 0 else 0 for g in pop2 ]
    plt.figure(figsize=(6,6))
    plt.scatter(vals1, vals2, s=6, alpha=0.6)
    plt.xlabel('Mean fraction of fixed heterozygotes per gene')
    plt.ylabel('Mean fraction of fixed heterozygotes per gene')
    plt.tight_layout()
    plt.savefig(out_png)
    plt.close()

def main():
    if len(sys.argv) < 4:
        print('Usage: make_pairwise_plots.py <het_tsv> <out_dir> <prefix>')
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

    for i in range(len(pops)):
        for j in range(i+1, len(pops)):
            a = pops[i]
            b = pops[j]
            out_png = os.path.join(out_dir, f"{prefix}_pairwise_{a}_vs_{b}.png")
            plot_pairwise(data[a], data[b], out_png)
    print('plots written to', out_dir)

if __name__ == '__main__':
    main()
