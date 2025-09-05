#!/usr/bin/env python3
"""
Generate per-ploidy pixy populations files from a metadata TSV.
Usage: generate_ploidy_populations.py <METAFILE> <OUTDIR>
Writes files like OUTDIR/populations_<ploidy>.txt containing sample<tab>pop lines.
"""
import sys
import os
import pandas as pd

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: generate_ploidy_populations.py METAFILE OUTDIR", file=sys.stderr)
        sys.exit(2)
    metafile = sys.argv[1]
    outdir = sys.argv[2]
    os.makedirs(outdir, exist_ok=True)
    df = pd.read_csv(metafile, sep="\t", dtype=str)
    if (
        "ploidy" not in df.columns
        or "sample" not in df.columns
        or "pop" not in df.columns
    ):
        raise ValueError("METAFILE must contain `sample`, `pop`, and `ploidy` columns")
    for pl in df["ploidy"].unique().tolist():
        sub = df.loc[df["ploidy"] == str(pl), ["sample", "pop"]]
        out = os.path.join(outdir, f"populations_{pl}.txt")
        if sub.empty:
            open(out, "w").close()
        else:
            sub.to_csv(out, sep="\t", header=False, index=False)
    print("Wrote ploidy population files to", outdir)
