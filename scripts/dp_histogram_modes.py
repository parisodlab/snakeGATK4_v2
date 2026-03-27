import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages
from scipy.signal import find_peaks
import os

# Access inputs and outputs from Snakemake
input_file = snakemake.input[0]
output_geno_pdf = snakemake.output[0]
output_sites_pdf = snakemake.output[1]


def write_placeholder_pdf(path, title, reason):
    with PdfPages(path) as pdf:
        fig, ax = plt.subplots(figsize=(10, 5))
        ax.axis('off')
        ax.text(0.5, 0.60, title, ha='center', va='center', fontsize=16)
        ax.text(0.5, 0.42, reason, ha='center', va='center', fontsize=11)
        pdf.savefig(fig)
        plt.close(fig)

# Load CSV (tab-separated)
df = pd.read_csv(input_file, sep='\t')

# Filter bins
df = df[(df['Bin'] >= 4) & (df['Bin'] <= 200)]

if df.empty:
    write_placeholder_pdf(output_geno_pdf, "NumGeno Distribution", "No depth bins between 4 and 200 were available.")
    write_placeholder_pdf(output_sites_pdf, "NumSites Distribution", "No depth bins between 4 and 200 were available.")
    raise SystemExit(0)

# Unique samples
samples = df['Sample'].unique()

if len(samples) == 0:
    write_placeholder_pdf(output_geno_pdf, "NumGeno Distribution", "No samples remained after filtering the DP table.")
    write_placeholder_pdf(output_sites_pdf, "NumSites Distribution", "No samples remained after filtering the DP table.")
    raise SystemExit(0)

# Plot NumGeno vs Bin
with PdfPages(output_geno_pdf) as pdf:
    geno_peaks = []
    for sample in samples:
        df_s = df[df['Sample'] == sample].sort_values('Bin')
        plt.figure(figsize=(10,5))
        x = np.arange(len(df_s))
        heights = df_s['NumGeno'].values
        plt.bar(x, heights, color='skyblue', align='center', edgecolor='black', linewidth=0.2)
        # Find peaks on the heights array
        counts = heights
        peaks, _ = find_peaks(counts)
        # Select top 2 peaks by height
        if len(peaks) > 0:
            top_peaks = peaks[counts[peaks].argsort()[::-1][:2]]
            plt.scatter(x[top_peaks], counts[top_peaks], color='red', s=60, zorder=10, label='Mode')
            # record peak info
            for rank, pk in enumerate(top_peaks, start=1):
                geno_peaks.append({
                    'Sample': sample,
                    'Metric': 'NumGeno',
                    'Bin': int(df_s['Bin'].iloc[pk]),
                    'Value': float(counts[pk]),
                    'Rank': rank
                })
        plt.title(f"NumGeno Distribution - {sample}")
        plt.xlabel("Bin")
        plt.ylabel("NumGeno")
        # set xticks to bin labels; reduce number of labels if too many
        if len(x) > 40:
            step = int(np.ceil(len(x) / 40))
        else:
            step = 1
        plt.xticks(x[::step], df_s['Bin'].astype(str).values[::step], rotation=90)
        plt.legend()
        plt.tight_layout()
        pdf.savefig()
        plt.close()
    # save geno peaks table next to PDF
    try:
        geno_csv = os.path.splitext(output_geno_pdf)[0] + '_peaks.csv'
        pd.DataFrame(geno_peaks).to_csv(geno_csv, index=False)
    except Exception:
        pass

# Plot NumSites vs Bin
with PdfPages(output_sites_pdf) as pdf:
    sites_peaks = []
    for sample in samples:
        df_s = df[df['Sample'] == sample].sort_values('Bin')
        plt.figure(figsize=(10,5))
        x = np.arange(len(df_s))
        heights = df_s['NumSites'].values
        plt.bar(x, heights, color='lightgreen', align='center', edgecolor='black', linewidth=0.2)
        # Find peaks on the heights array
        counts = heights
        peaks, _ = find_peaks(counts)
        # Select top 2 peaks by height
        if len(peaks) > 0:
            top_peaks = peaks[counts[peaks].argsort()[::-1][:2]]
            plt.scatter(x[top_peaks], counts[top_peaks], color='red', s=60, zorder=10, label='Mode')
            # record peak info
            for rank, pk in enumerate(top_peaks, start=1):
                sites_peaks.append({
                    'Sample': sample,
                    'Metric': 'NumSites',
                    'Bin': int(df_s['Bin'].iloc[pk]),
                    'Value': float(counts[pk]),
                    'Rank': rank
                })
        plt.title(f"NumSites Distribution - {sample}")
        plt.xlabel("Bin")
        plt.ylabel("NumSites")
        if len(x) > 40:
            step = int(np.ceil(len(x) / 40))
        else:
            step = 1
        plt.xticks(x[::step], df_s['Bin'].astype(str).values[::step], rotation=90)
        plt.legend()
        plt.tight_layout()
        pdf.savefig()
        plt.close()
    # save sites peaks table next to PDF
    try:
        sites_csv = os.path.splitext(output_sites_pdf)[0] + '_peaks.csv'
        pd.DataFrame(sites_peaks).to_csv(sites_csv, index=False)
    except Exception:
        pass
