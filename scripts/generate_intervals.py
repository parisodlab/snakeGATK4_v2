#!/usr/bin/env python3
"""Generate tiled intervals TSV from a VCF header's contig lines.
Usage: generate_intervals.py <vcf> <out_tsv> <window_size> <step>
Writes a TSV with columns: chrom\tstart\tend (1-based inclusive).
"""
"""
Generate tiled intervals TSV from a VCF header's contig lines.
Usage: generate_intervals.py <vcf> <out_tsv> <window_size> <step>

If <window_size> is 0, the script will emit a single interval per contig
covering the whole chromosome (start=1, end=length). Otherwise it tiles
each contig using <window_size> and <step> as before.

Writes a TSV with columns: chrom\tstart\tend (1-based inclusive).
"""
import sys
import subprocess
import re
import os


def read_contigs_from_vcf(vcf_path):
    # call bcftools view -h and parse ##contig lines
    cmd = ["bcftools", "view", "-h", vcf_path]
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    out, err = p.communicate()
    if p.returncode != 0:
        raise SystemExit(f"bcftools failed: {err.strip()}")
    contigs = []
    rx = re.compile(r"ID=([^,>]+).*length=([0-9]+)")
    for line in out.splitlines():
        if line.startswith("##contig"):
            m = rx.search(line)
            if m:
                contigs.append((m.group(1), int(m.group(2))))
    return contigs


def tile_contig(chrom, length, win, step):
    # if win <= 0 treat as whole-contig interval
    if win <= 0:
        yield chrom, 1, length
        return

    start = 1
    while start <= length:
        end = min(start + win - 1, length)
        yield chrom, start, end
        start += step


def main(argv):
    if len(argv) != 5:
        print(__doc__)
        sys.exit(2)
    _, vcf, out_tsv, win_s, step_s = argv
    win = int(win_s)
    step = int(step_s)
    os.makedirs(os.path.dirname(out_tsv) or ".", exist_ok=True)
    contigs = read_contigs_from_vcf(vcf)
    with open(out_tsv, "w") as out:
        for chrom, length in contigs:
            for c, s, e in tile_contig(chrom, length, win, step):
                out.write("\t".join([c, str(s), str(e)]) + "\n")


if __name__ == "__main__":
    main(sys.argv)
