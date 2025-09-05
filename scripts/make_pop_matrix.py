#!/usr/bin/env python3
import sys
from pathlib import Path

if len(sys.argv) < 3:
    print("Usage: make_pop_matrix.py <pop_map.txt> <out.tsv>", file=sys.stderr)
    sys.exit(1)

pop_map = Path(sys.argv[1])
out_file = Path(sys.argv[2])
if not pop_map.exists():
    print(f"File {pop_map} not found", file=sys.stderr)
    sys.exit(1)

samples = []
pop_of = {}
with pop_map.open() as fh:
    for line in fh:
        line = line.strip()
        if not line:
            continue
        parts = line.split()
        sample = parts[0]
        pop = sample[:3]
        samples.append(sample)
        pop_of[sample] = pop

unique_pops = sorted(sorted(set(pop_of.values())))

# write matrix
with out_file.open("w") as out:
    out.write("Samples")
    for p in unique_pops:
        out.write("\t" + p)
    out.write("\n")
    for s in samples:
        out.write(s)
        for p in unique_pops:
            out.write("\t" + ("1" if pop_of[s] == p else "0"))
        out.write("\n")
print(f"Wrote {out_file} with {len(samples)} samples and {len(unique_pops)} pops")
