#!/usr/bin/env python3
import argparse
import gzip
import pandas as pd
import re
from tqdm import tqdm

# ----------------------------
# Functions
# ----------------------------

def parse_gt(gt):
    """Convert genotype string to dosage (0-1). Handles any ploidy."""
    if "." in gt or gt == "":
        return None
    alleles = re.split(r"[\/|]", gt)
    alleles = [a for a in alleles if a != ""]
    if not alleles:
        return None
    try:
        alleles = list(map(int, alleles))
    except ValueError:
        return None
    ploidy = len(alleles)
    alt_count = sum(1 for a in alleles if a == 1)
    return alt_count / ploidy

def is_fixed_hetero(dosages, tol=1e-6):
    """Strict fixed heterozygote check for any ploidy."""
    dosages = [d for d in dosages if d is not None]
    if not dosages:
        return False
    if not all(0 < d < 1 for d in dosages):
        return False
    first = dosages[0]
    return all(abs(d - first) < tol for d in dosages)

def is_variable(dosages):
    """Check if a site is variable (some alleles differ across individuals)."""
    dosages = [d for d in dosages if d is not None]
    if not dosages:
        return False
    return any(d > 0 for d in dosages) and any(d < 1 for d in dosages)

def load_samples(samples_file):
    """Load sample-population mapping."""
    df = pd.read_csv(samples_file, sep="\t", header=None, names=["sample","pop"])
    pop_map = df.groupby("pop")["sample"].apply(list).to_dict()
    return df, pop_map

def load_genes_from_gff3(gff_file):
    """
    Parse MAKER GFF3 and return a DataFrame with gene coordinates.
    Only extracts 'gene' features.
    """
    genes = []
    open_func = gzip.open if gff_file.endswith(".gz") else open
    with open_func(gff_file, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 9:
                continue
            chrom, feature_type, start, end, attr = parts[0], parts[2], int(parts[3]), int(parts[4]), parts[8]
            if feature_type.lower() == "gene":
                # Extract gene ID from attributes (ID=gene:XYZ)
                match = re.search(r"ID=gene:([^;]+)", attr)
                gene_id = match.group(1) if match else f"{chrom}_{start}_{end}"
                genes.append([chrom, start, end, gene_id])
    return pd.DataFrame(genes, columns=["chrom","start","end","gene_id"])


def process_vcf(vcf_file, samples_df, pop_map, genes_df, missing_threshold=0.33):
    """Stream VCF and compute per-gene, per-population variable and fixed heterozygote counts.

    This avoids loading the entire genotype matrix into memory. For each variant, the code
    finds overlapping genes (from genes_df) and updates running counts per gene and population.
    """
    import bisect
    open_func = gzip.open if vcf_file.endswith(".gz") else open

    # Build genes lookup by chromosome: sorted starts, ends, gids
    genes_by_chrom = {}
    gene_ids = []
    for _, g in genes_df.iterrows():
        chrom = g.chrom
        genes_by_chrom.setdefault(chrom, {"starts": [], "ends": [], "gids": []})
        genes_by_chrom[chrom]["starts"].append(int(g.start))
        genes_by_chrom[chrom]["ends"].append(int(g.end))
        genes_by_chrom[chrom]["gids"].append(g.gene_id)
        gene_ids.append(g.gene_id)

    # Ensure lists are sorted by start (they should be already if genes_df is ordered, but safe)
    for chrom, d in genes_by_chrom.items():
        ordered = sorted(zip(d["starts"], d["ends"], d["gids"]))
        d["starts"], d["ends"], d["gids"] = map(list, zip(*ordered)) if ordered else ([], [], [])

    # Initialize counters: gene_id -> pop -> {'variable': int, 'fixed': int, 'sites': int}
    gene_counts = {gid: {pop: {'variable': 0, 'fixed': 0, 'sites': 0} for pop in pop_map} for gid in gene_ids}

    header_samples = None
    # Precompute per-pop sample column indices after reading header
    pop_sample_indices = {}

    # Use tqdm to show progress while streaming the VCF
    with open_func(vcf_file, "rt") as f:
        # iterate with tqdm to show activity
        for line in tqdm(f, desc="Processing VCF", unit="lines"):
            if line.startswith("#CHROM"):
                header = line.strip().split("\t")
                header_samples = header[9:]
                # map sample name to column index in parts
                sample_to_col = {s: i + 9 for i, s in enumerate(header_samples)}
                # build per-pop indices (only for samples present in VCF)
                for pop, samples in pop_map.items():
                    pop_sample_indices[pop] = [sample_to_col[s] for s in samples if s in sample_to_col]
                # Sanity: print counts once header parsed
                import sys
                print(f"Found {len(header_samples)} samples in VCF; populations: {len(pop_sample_indices)}", flush=True)
                continue
            if line.startswith("#"):
                continue

            parts = line.strip().split("\t")
            if len(parts) < 2:
                continue
            chrom = parts[0]
            try:
                pos = int(parts[1])
            except ValueError:
                continue

            if chrom not in genes_by_chrom:
                continue

            d = genes_by_chrom[chrom]
            # bisect to find all genes with start <= pos
            idx = bisect.bisect_right(d["starts"], pos)
            # iterate backwards through candidates until end < pos
            for i in range(idx - 1, -1, -1):
                if d["ends"][i] < pos:
                    break
                gid = d["gids"][i]

                # For each population, collect dosages for the samples in that pop
                for pop, indices in pop_sample_indices.items():
                    if not indices:
                        continue
                    dosages = []
                    for col_idx in indices:
                        if col_idx >= len(parts):
                            # missing sample column
                            dosages.append(None)
                            continue
                        sample_field = parts[col_idx]
                        gt = sample_field.split(":")[0] if sample_field != "." else "."
                        dosages.append(parse_gt(gt))

                    # If all samples missing, skip
                    if not dosages:
                        continue
                    missing = sum(1 for dval in dosages if dval is None)
                    if missing / len(dosages) >= missing_threshold:
                        continue

                    # Update per-gene-per-pop counters
                    gene_counts[gid][pop]['sites'] += 1
                    if is_variable(dosages):
                        gene_counts[gid][pop]['variable'] += 1
                    if is_fixed_hetero(dosages):
                        gene_counts[gid][pop]['fixed'] += 1

    # Build output DataFrame in same order as genes_df
    out_rows = []
    for _, g in genes_df.iterrows():
        gid = g.gene_id
        chrom = g.chrom
        start = int(g.start)
        end = int(g.end)
        row = [chrom, start, end]
        for pop in pop_map:
            counts = gene_counts.get(gid, {}).get(pop, {'variable': 0, 'fixed': 0})
            row.append(f"{counts['variable']},{counts['fixed']}")
        out_rows.append(row)

    header = ["contig", "start", "end"] + list(pop_map.keys())
    out_df = pd.DataFrame(out_rows, columns=header)
    return out_df

# ----------------------------
# Main
# ----------------------------
def main():
    parser = argparse.ArgumentParser(description="Per-gene fixed heterozygote and variable site counting with mixed ploidy (MAKER GFF3)")
    parser.add_argument("-v", "--vcf", required=True, help="Input VCF (can be gzipped)")
    parser.add_argument("-s", "--samples", required=True, help="Sample-to-population file (tab-delimited: sample pop)")
    parser.add_argument("-a", "--annotation", required=True, help="MAKER GFF3 annotation file")
    parser.add_argument("-m", "--missing", default=0.33, type=float, help="Missingness threshold per population")
    parser.add_argument("-o", "--output", required=True, help="Output file")
    args = parser.parse_args()

    samples_df, pop_map = load_samples(args.samples)
    genes_df = load_genes_from_gff3(args.annotation)
    out_df = process_vcf(args.vcf, samples_df, pop_map, genes_df, args.missing)
    out_df.to_csv(args.output, sep="\t", index=False)

if __name__ == "__main__":
    main()
