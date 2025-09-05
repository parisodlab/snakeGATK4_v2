#!/usr/bin/env python3
import pandas as pd
import sys
from pathlib import Path

p = Path(sys.argv[1])
if not p.exists():
    print(f"File {p} not found", file=sys.stderr)
    sys.exit(1)

df = pd.read_csv(p, sep="\t", dtype=str)
if "sample" not in df.columns:
    print('No "sample" column found in file', file=sys.stderr)
    sys.exit(1)

# Create or replace 'pop' column with first 3 characters of 'sample'
df["pop"] = df["sample"].astype(str).str.slice(0, 3)

# Preserve original column order; ensure 'pop' is the second column if possible
cols = list(df.columns)
if cols[0] != "sample":
    # move sample to front
    cols.remove("sample")
    cols.insert(0, "sample")
if "pop" in cols:
    cols.remove("pop")
    cols.insert(1, "pop")

# Reorder and write back
df = df[cols]
# Write without index, preserving tab-separated format
df.to_csv(p, sep="\t", index=False)
print(f"Updated {p} (set pop = first 3 chars of sample)")
