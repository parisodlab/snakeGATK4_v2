#!/usr/bin/env python3
"""
scripts/scantools_run.py
Usage: scantools_run.py <vcf.gz> <sample_map.txt> <scaffolds.list> <outdir> <window>

Produces per-chromosome windowed TSVs with columns: chrom, start, end, n_sites, mean_depth, mean_af, mean_pi, het_count, missing_prop
Tries to use ScanTools if importable; otherwise falls back to cyvcf2 implementation.
"""
import sys
import os
import gzip
import shutil
import subprocess
import math
import datetime
import time
import glob
from collections import defaultdict
import pandas as pd
from pathlib import Path


def read_sample_map(path):
    sm = {}
    with open(path) as f:
        for line in f:
            if line.strip() == "":
                continue
            parts = line.strip().split("\t")
            # expect sample\tploidy or at least sample in first column
            sample = parts[0]
            ploidy = int(parts[1]) if len(parts) > 1 and parts[1].isdigit() else None
            sm[sample] = ploidy
    return sm


def read_scaffolds(path):
    scaff = []
    with open(path) as f:
        for line in f:
            if line.strip() == "":
                continue
            parts = line.strip().split()
            if len(parts) == 1:
                name = parts[0]
                length = None
            else:
                name, length = parts[0], int(parts[1])
            scaff.append((name, length))
    return scaff


def fallback_cyvcf2(vcf_path, samples_map, scaffolds, outdir, window):
    try:
        from cyvcf2 import VCF
        import numpy as np
    except Exception as e:
        raise RuntimeError("cyvcf2 is required for the fallback path: " + str(e))

    vcf = VCF(vcf_path)
    samples = vcf.samples

    # pre-create per-chrom bins
    bins = {}
    for chrom, length in scaffolds:
        if length is None:
            # will infer from sites
            continue
        edges = list(range(1, length + window, window))
        bins[chrom] = edges

    # accumulate per-bin stats: dict chrom->bin_start->stats
    stats = defaultdict(
        lambda: defaultdict(
            lambda: {
                "n_sites": 0,
                "depth": 0,
                "af_sum": 0.0,
                "pi_sum": 0.0,
                "het": 0,
                "missing": 0,
            }
        )
    )

    for rec in vcf:
        chrom = rec.CHROM
        pos = rec.POS
        if chrom not in bins:
            # skip sites on scaffolds with unknown length (or create single bin)
            # create a single bin keyed by pos rounded down to window
            bstart = ((pos - 1) // window) * window + 1
        else:
            edges = bins[chrom]
            bstart = ((pos - 1) // window) * window + 1

        # compute per-site totals using AD if available
        try:
            AD = rec.format("AD")
        except Exception:
            AD = None

        ref_sum = 0
        alt_sum = 0
        missing = 0
        het_count = 0
        for i, s in enumerate(samples):
            refc = None
            altc = None
            if AD is not None:
                try:
                    a = AD[i]
                    if a is None:
                        missing += 1
                        continue
                    # AD may be bytes or array
                    if len(a) >= 2:
                        refc = int(a[0])
                        altc = int(a[1])
                    else:
                        # fallback
                        refc = 0
                        altc = 0
                except Exception:
                    missing += 1
                    continue
            else:
                gt = rec.genotypes[i]
                alleles = gt[: len(gt) - 1] if isinstance(gt, list) else gt
                if any([a == -1 for a in alleles]):
                    missing += 1
                    continue
                refc = sum(1 for a in alleles if a == 0)
                altc = sum(1 for a in alleles if a == 1)

            ref_sum += refc
            alt_sum += altc
            if refc > 0 and altc > 0:
                het_count += 1

        depth = ref_sum + alt_sum
        if depth == 0:
            af = float("nan")
            pi = float("nan")
        else:
            af = alt_sum / depth
            pi = 2 * af * (1 - af)

        st = stats[chrom][bstart]
        st["n_sites"] += 1
        st["depth"] += depth
        st["af_sum"] += 0 if (af != af) else af
        st["pi_sum"] += 0 if (pi != pi) else pi
        st["het"] += het_count
        st["missing"] += missing

    # write per-chrom output files
    os.makedirs(outdir, exist_ok=True)
    combined = []
    for chrom, bins_dict in stats.items():
        rows = []
        for bstart in sorted(bins_dict.keys()):
            st = bins_dict[bstart]
            n = st["n_sites"]
            if n == 0:
                continue
            mean_depth = st["depth"] / n
            mean_af = st["af_sum"] / n
            mean_pi = st["pi_sum"] / n
            het = st["het"]
            missing_prop = st["missing"] / (n * len(samples))
            end = bstart + window - 1
            rows.append(
                (chrom, bstart, end, n, mean_depth, mean_af, mean_pi, het, missing_prop)
            )
            combined.append(
                (chrom, bstart, end, n, mean_depth, mean_af, mean_pi, het, missing_prop)
            )

        outf = os.path.join(outdir, f"{chrom}_windows.tsv")
        with open(outf, "w") as fh:
            fh.write(
                "chrom\tstart\tend\tn_sites\tmean_depth\tmean_af\tmean_pi\thet_count\tmissing_prop\n"
            )
            for r in rows:
                fh.write("\t".join(map(str, r)) + "\n")

    # combined
    comb = os.path.join(outdir, "combined_windows.tsv")
    with open(comb, "w") as fh:
        fh.write(
            "chrom\tstart\tend\tn_sites\tmean_depth\tmean_af\tmean_pi\thet_count\tmissing_prop\n"
        )
        for r in combined:
            fh.write("\t".join(map(str, r)) + "\n")

    print("Wrote scantools-like windowed outputs to:", outdir)


def main():
    if len(sys.argv) < 6:
        print(__doc__)
        sys.exit(1)
    # positional args
    vcf = sys.argv[1]
    sample_map = sys.argv[2]
    scaffolds = sys.argv[3]
    outdir = sys.argv[4]
    window = int(sys.argv[5])

    # optional args: scantools path and mode
    scantools_path = None
    mode = "local"  # local, pbs, slurm
    pops_arg = None
    if len(sys.argv) >= 7:
        scantools_path = sys.argv[6]
    if len(sys.argv) >= 8:
        mode = sys.argv[7]
    if len(sys.argv) >= 9:
        pops_arg = sys.argv[8]

    samples = read_sample_map(sample_map)
    scafs = read_scaffolds(scaffolds)

    # If a ScanTools path was provided, try to use it (without installing):
    if scantools_path:
        scantools_path = os.path.abspath(scantools_path)
        if not os.path.exists(scantools_path):
            raise FileNotFoundError(f"ScanTools path {scantools_path} not found")

        # create a patched copy for SLURM if requested
        patched_path = scantools_path
        if mode.lower() == "slurm":
            patched_path = os.path.join(scantools_path, "ScanTools_slurm")
            if os.path.exists(patched_path):
                print(f"Patched ScanTools already exists at {patched_path}")
            else:
                print(f"Creating SLURM-patched copy of ScanTools at {patched_path}")
                shutil.copytree(scantools_path, patched_path)
                # simple textual replacements: PBS -> SBATCH and qsub -> sbatch
                for pyfile in glob.glob(
                    os.path.join(patched_path, "**", "*.py"), recursive=True
                ):
                    try:
                        with open(pyfile, "r") as fh:
                            txt = fh.read()
                        newtxt = txt.replace("#PBS", "#SBATCH")
                        newtxt = newtxt.replace("qsub", "sbatch")
                        # replace common PBS resource flags with SBATCH equivalents heuristically
                        newtxt = newtxt.replace("walltime=", "--time=")
                        newtxt = newtxt.replace("mem=", "--mem=")
                        with open(pyfile, "w") as fh:
                            fh.write(newtxt)
                    except Exception:
                        continue

        # ensure we can import ScanTools from the patched path
        sys.path.insert(0, patched_path)
        try:
            import ScanTools as ST
        except Exception as e:
            print("Failed to import ScanTools from", patched_path, "error:", e)
            print("Proceeding without ScanTools import; fall back will be used")
            ST = None

        # prepare vcf directory inside ScanTools repo
        vcf_dir_name = "polyplChapter"
        vcf_dir = os.path.join(patched_path, vcf_dir_name)
        os.makedirs(vcf_dir, exist_ok=True)
        # copy vcf into that folder (preserve index)
        vcf_basename = os.path.basename(vcf)
        dest_vcf = os.path.join(vcf_dir, vcf_basename)
        print(f"Copying {vcf} -> {dest_vcf}")
        shutil.copy2(vcf, dest_vcf)
        if os.path.exists(vcf + ".tbi"):
            shutil.copy2(vcf + ".tbi", dest_vcf + ".tbi")

        # prepare PopKey
        popkey_path = os.path.join(patched_path, "PartA", "PopKey.txt")
        os.makedirs(os.path.dirname(popkey_path), exist_ok=True)
        # if pops_arg is a path to a popkey file, use it
        if pops_arg and os.path.exists(pops_arg):
            shutil.copy2(pops_arg, popkey_path)
            print("Copied provided PopKey to", popkey_path)
        else:
            # try to create PopKey from sample_map: expect a third column 'pop' or second column as pop
            df = pd.read_csv(sample_map, sep="\t", header=None)
            # heuristics: if header present
            try:
                header = pd.read_csv(sample_map, sep="\t", nrows=0).columns.tolist()
                dfh = pd.read_csv(sample_map, sep="\t")
                if "pop" in dfh.columns:
                    popcol = "pop"
                elif "population" in dfh.columns:
                    popcol = "population"
                elif "popcode" in dfh.columns:
                    popcol = "popcode"
                else:
                    # fallback: if second column non-numeric, treat as pop
                    df2 = pd.read_csv(sample_map, sep="\t", header=None)
                    if df2.shape[1] >= 2 and not df2[1].dtype == int:
                        popcol = 1
                    else:
                        popcol = None
            except Exception:
                # fallback simplistic parse
                df2 = pd.read_csv(sample_map, sep="\t", header=None)
                if df2.shape[1] >= 2:
                    popcol = 1
                else:
                    popcol = None

            if popcol is None:
                raise ValueError(
                    "Could not infer population column from sample_map. Provide a PopKey file or include a pop column in sample_map."
                )

            # write PopKey as two columns: sample\tPOP
            with open(popkey_path, "w") as pk:
                if isinstance(popcol, int):
                    df2 = pd.read_csv(sample_map, sep="\t", header=None)
                    for _, row in df2.iterrows():
                        sample = row[0]
                        pop = row[popcol]
                        pk.write(f"{sample}\t{pop}\n")
                else:
                    dfh = pd.read_csv(sample_map, sep="\t")
                    for _, row in dfh.iterrows():
                        sample = row.iloc[0]
                        pop = row[popcol]
                        pk.write(f"{sample}\t{pop}\n")
            print("Wrote PopKey to", popkey_path)

        # initialize ScanTools instance if imported
        if ST is not None:
            print("Initializing ScanTools object...")
            test = ST.scantools(patched_path)

            # now run splitVCFsNorepol as per instructions
            print("Running splitVCFsNorepol ...")
            # default parameter choices, adapt as needed
            test.splitVCFsNorepol(
                vcf_dir=vcf_dir_name,
                min_dp="8",
                mffg="0.2",
                mem="16",
                time_scratch="02:00:00",
                ncpu="12",
                overwrite=True,
                scratch_gb="1",
                keep_intermediates=False,
                use_scratch=True,
                scratch_path=os.environ.get("SCRATCHDIR", "/scratch"),
                pops=None,
                print1=False,
            )

            # create pops list from PopKey
            pops = []
            with open(popkey_path) as pk:
                for line in pk:
                    s = line.strip().split("\t")
                    if len(s) >= 2:
                        pop = s[1]
                        if pop not in pops:
                            pops.append(pop)

            # run within-population metrics
            print("Running calcwpm ...")
            test.calcwpm(
                recode_dir=f"VCF_{vcf_dir_name}_DP8.M0.2",
                window_size=window,
                min_snps=50,
                pops=pops,
                mem=1,
                ncpu=1,
                scratch_gb=1,
                use_repol=False,
                time_scratch="1:20:00",
                overwrite=True,
                sampind=7,
                print1=False,
            )

            # run pairwise metrics
            print("Running calcPairwisebpm ...")
            test.calcPairwisebpm(
                recode_dir=f"VCF_{vcf_dir_name}_DP8.M0.2",
                pops=pops,
                window_size=window,
                min_snps=50,
                mem=1,
                ncpu=1,
                use_repol=False,
                keep_intermediates=False,
                time_scratch="0:40:00",
                scratch_gb=1,
                print1=False,
            )

        else:
            print(
                "ScanTools not available for API calls; the script will NOT run ScanTools steps. It did copy VCF and created PopKey at",
                patched_path,
            )

        # after ScanTools steps we still run the fallback generator in case user wants a plain combined_windows.tsv
        fallback_cyvcf2(vcf, samples, scafs, outdir, window)

    else:
        # no ScanTools path provided: fallback only
        fallback_cyvcf2(vcf, samples, scafs, outdir, window)


if __name__ == "__main__":
    main()
