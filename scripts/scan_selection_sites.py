#!/usr/bin/env python
"""Stream-process piawka per-site output for two populations.

The mixed per-site piawka output can be extremely large (tens of GB). This script
streams through it and extracts only the rows needed to build a compact per-site
DAS vs DAG (or any two pops) table.

Input example (from this repo):
  results/piawka/mixed/persite/piawka_tajima_persite.tsv.gz

Outputs:
- sites_<HIGH>_vs_<LOW>.tsv.gz : one row per site with per-pop pi, log10(pi ratio), and
    between-pop Dxy (and Fst if present)
- genomewide_summary_sites.json : mean + (approx) 2.5%/97.5% from reservoir samples
- top_sites.<type>.tsv : top-N sites for simple signals (optional)

Notes
- We assume input is grouped by locus (VCF order), so we can flush on locus changes.
- Per-site file locus often looks like: <vcf_path>_Bv1_1030837
  We parse chrom/pos from the last 2 underscore-separated tokens.
- Tajima's D is a window-based statistic; it is not meaningful per single site.
    The per-site piawka output therefore typically does NOT contain TajD.

Run under a conda env with Python (no pandas required):
  conda run -n python_plot python scripts/scan_selection_sites.py ...
"""

from __future__ import annotations

import argparse
import gzip
import json
import math
import random
import sys
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple


def parse_site_locus(locus: str) -> Tuple[str, int]:
    """Parse per-site locus into (chrom, pos).

    Handles strings like:
      /path/to/file.vcf.gz_Bv1_1030837
    and falls back to:
      Bv1_1030837
    """
    # Windowed locus looks like Chrom:start-end; do not accept here
    if ":" in locus and "-" in locus:
        raise ValueError(f"Locus looks windowed, not per-site: {locus}")

    # Split from the right to tolerate underscores in path
    parts = locus.rsplit("_", 2)
    if len(parts) < 3:
        raise ValueError(f"Unrecognized per-site locus format: {locus!r}")
    chrom = parts[-2]
    pos = int(parts[-1])
    return chrom, pos


@dataclass
class Reservoir:
    k: int
    values: List[float]
    n_seen: int = 0

    def add(self, x: float, rng: random.Random):
        if x is None or math.isnan(x) or math.isinf(x):
            return
        self.n_seen += 1
        if len(self.values) < self.k:
            self.values.append(float(x))
            return
        j = rng.randint(1, self.n_seen)
        if j <= self.k:
            self.values[j - 1] = float(x)


def _summary_from_reservoir(values: List[float]):
    if not values:
        return {"mean": float("nan"), "q025": float("nan"), "q975": float("nan")}
    v = sorted(values)
    n = len(v)
    mean = sum(v) / n

    def q(p: float) -> float:
        if n == 1:
            return v[0]
        idx = p * (n - 1)
        lo = int(math.floor(idx))
        hi = int(math.ceil(idx))
        if lo == hi:
            return v[lo]
        frac = idx - lo
        return v[lo] * (1 - frac) + v[hi] * frac

    return {"mean": float(mean), "q025": float(q(0.025)), "q975": float(q(0.975))}


def write_tsv_row(out_fh, fields: List[object]):
    out_fh.write("\t".join("" if f is None else str(f) for f in fields) + "\n")


def main(argv: Optional[List[str]] = None) -> int:
    ap = argparse.ArgumentParser(description="Process piawka per-site output for two pops")
    ap.add_argument("--piawka-persite", type=Path, required=True, help="piawka_tajima_persite.tsv.gz")
    ap.add_argument("--high-pop", required=True, help="High elevation pop label (e.g. DAS)")
    ap.add_argument("--low-pop", required=True, help="Low elevation pop label (e.g. DAG)")
    ap.add_argument("--outdir", type=Path, required=True)
    ap.add_argument("--topn", type=int, default=1000, help="Top N sites to track for simple signals (default: 1000)")
    ap.add_argument(
        "--reservoir-size",
        type=int,
        default=200000,
        help="Reservoir sample size for approximate genome-wide 95%% intervals (default: 200000)",
    )
    ap.add_argument(
        "--max-sites",
        type=int,
        default=0,
        help="Debug: stop after writing this many sites (0 = no limit)",
    )
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument(
        "--pi-pseudocount",
        type=float,
        default=0.0,
        help="Add this pseudocount to BOTH pi_high and pi_low when computing log10(piH/piL) (default: 0, i.e. require pi>0)",
    )
    ap.add_argument(
        "--progress-every",
        type=int,
        default=1_000_000,
        help="Write a progress line to stderr every N input lines (default: 1,000,000)",
    )

    args = ap.parse_args(argv)

    pi_pseudocount = float(args.pi_pseudocount or 0.0)
    if pi_pseudocount < 0:
        raise ValueError("--pi-pseudocount must be >= 0")

    args.outdir.mkdir(parents=True, exist_ok=True)

    out_sites = args.outdir / f"sites_{args.high_pop}_vs_{args.low_pop}.tsv.gz"

    rng = random.Random(args.seed)
    res_fst = Reservoir(args.reservoir_size, [])
    res_dxy = Reservoir(args.reservoir_size, [])
    res_pi_h = Reservoir(args.reservoir_size, [])
    res_pi_l = Reservoir(args.reservoir_size, [])
    res_lr = Reservoir(args.reservoir_size, [])

    # Keep top sites by simple scores without needing z-scores
    # store tuples (score, chrom, pos, ...)
    import heapq

    top_abs_lr: List[Tuple[float, str, int]] = []
    top_dxy: List[Tuple[float, str, int]] = []
    top_fst: List[Tuple[float, str, int]] = []

    def push_top(heap: List[Tuple[float, str, int]], score: float, chrom: str, pos: int, topn: int):
        if math.isnan(score) or math.isinf(score):
            return
        item = (float(score), chrom, int(pos))
        if len(heap) < topn:
            heapq.heappush(heap, item)
        else:
            if item[0] > heap[0][0]:
                heapq.heapreplace(heap, item)

    # Streaming state
    current_locus: Optional[str] = None
    current_chrom: Optional[str] = None
    current_pos: Optional[int] = None
    pi_h = None
    pi_l = None
    dxy = None
    fst_hud = None

    def flush():
        nonlocal pi_h, pi_l, dxy, fst_hud
        nonlocal current_locus, current_chrom, current_pos
        nonlocal written

        if current_locus is None:
            return
        # Only output if we have something relevant
        if all(v is None for v in [pi_h, pi_l, dxy, fst_hud]):
            pi_h = pi_l = dxy = fst_hud = None
            return

        log10_pi_ratio = None
        if pi_h is not None and pi_l is not None:
            if pi_pseudocount > 0:
                a = pi_h + pi_pseudocount
                b = pi_l + pi_pseudocount
                if a > 0 and b > 0 and math.isfinite(a) and math.isfinite(b):
                    log10_pi_ratio = math.log10(a / b)
            else:
                if pi_h > 0 and pi_l > 0:
                    log10_pi_ratio = math.log10(pi_h / pi_l)

        write_tsv_row(
            out_fh,
            [
                current_chrom,
                current_pos,
                current_locus,
                dxy,
                fst_hud,
                pi_h,
                pi_l,
                log10_pi_ratio,
            ],
        )

        # update reservoirs
        if fst_hud is not None:
            res_fst.add(fst_hud, rng)
        if dxy is not None:
            res_dxy.add(dxy, rng)
        if pi_h is not None:
            res_pi_h.add(pi_h, rng)
        if pi_l is not None:
            res_pi_l.add(pi_l, rng)
        if log10_pi_ratio is not None:
            res_lr.add(log10_pi_ratio, rng)

        # update top lists
        if log10_pi_ratio is not None:
            push_top(top_abs_lr, abs(log10_pi_ratio), current_chrom, current_pos, args.topn)
        if dxy is not None:
            push_top(top_dxy, dxy, current_chrom, current_pos, args.topn)
        if fst_hud is not None:
            push_top(top_fst, fst_hud, current_chrom, current_pos, args.topn)

        written += 1
        # reset
        pi_h = pi_l = dxy = fst_hud = None

    written = 0
    lines_read = 0
    t0 = time.time()
    progress_every = max(1, int(args.progress_every))
    sys.stderr.write(
        f"[scan_selection_sites] START in={args.piawka_persite} high={args.high_pop} low={args.low_pop} out={args.outdir}\n"
    )

    with gzip.open(args.piawka_persite, "rt") as in_fh, gzip.open(out_sites, "wt") as out_fh:
        header = in_fh.readline().rstrip("\n").split("\t")
        idx = {name: i for i, name in enumerate(header)}
        required = ["locus", "pop1", "pop2", "metric", "value"]
        for r in required:
            if r not in idx:
                raise ValueError(f"Missing required column {r!r} in header")

        # output header
        write_tsv_row(
            out_fh,
            [
                "chrom",
                "pos",
                "locus",
                f"dxy_{args.high_pop}_vs_{args.low_pop}",
                "fst",
                f"pi_{args.high_pop}",
                f"pi_{args.low_pop}",
                f"log10_pi_ratio_{args.high_pop}_over_{args.low_pop}",
            ],
        )

        for line in in_fh:
            lines_read += 1
            if lines_read % progress_every == 0:
                dt = max(1e-9, time.time() - t0)
                rate = lines_read / dt
                sys.stderr.write(
                    f"[scan_selection_sites] lines={lines_read:,} sites_written={written:,} rate={rate:,.0f} lines/s\n"
                )
            parts = line.rstrip("\n").split("\t")
            if len(parts) < len(header):
                continue

            locus = parts[idx["locus"]]
            metric = parts[idx["metric"]]

            # fast reject most rows
            if metric not in ("pi", "Dxy", "Fst_HUD"):
                continue

            # Detect locus change
            if current_locus is None:
                current_locus = locus
                current_chrom, current_pos = parse_site_locus(locus)
                current_chrom = str(current_chrom)
            elif locus != current_locus:
                flush()
                if args.max_sites and written >= args.max_sites:
                    break
                current_locus = locus
                current_chrom, current_pos = parse_site_locus(locus)
                current_chrom = str(current_chrom)

            pop1 = parts[idx["pop1"]]
            pop2 = parts[idx["pop2"]]
            try:
                value = float(parts[idx["value"]])
            except Exception:
                continue

            if metric == "pi":
                if pop2 == ".":
                    if pop1 == args.high_pop:
                        pi_h = value
                    elif pop1 == args.low_pop:
                        pi_l = value
            elif metric == "Dxy":
                if (pop1 == args.high_pop and pop2 == args.low_pop) or (pop1 == args.low_pop and pop2 == args.high_pop):
                    dxy = value
            elif metric == "Fst_HUD":
                if (pop1 == args.high_pop and pop2 == args.low_pop) or (pop1 == args.low_pop and pop2 == args.high_pop):
                    fst_hud = value

        # flush last
        flush()

    dt = max(1e-9, time.time() - t0)
    sys.stderr.write(
        f"[scan_selection_sites] DONE lines={lines_read:,} sites_written={written:,} elapsed={dt/3600:.2f}h\n"
    )

    # write summary + tops
    summary = {
        "pop_high": args.high_pop,
        "pop_low": args.low_pop,
        "pi_pseudocount_for_ratio": float(pi_pseudocount),
        "reservoir_size": args.reservoir_size,
        "metrics": {
            "dxy": _summary_from_reservoir(res_dxy.values),
            "fst": _summary_from_reservoir(res_fst.values),
            "pi_high": _summary_from_reservoir(res_pi_h.values),
            "pi_low": _summary_from_reservoir(res_pi_l.values),
            "log10_pi_ratio": _summary_from_reservoir(res_lr.values),
        },
        "n_sites_written": written,
    }

    with open(args.outdir / "genomewide_summary_sites.json", "w") as fh:
        json.dump(summary, fh, indent=2)

    def write_top(path: Path, heap: List[Tuple[float, str, int]], score_name: str):
        rows = sorted(heap, key=lambda x: x[0], reverse=True)
        with open(path, "w") as fh:
            fh.write("score\tchrom\tpos\n")
            for s, c, p in rows:
                fh.write(f"{s}\t{c}\t{p}\n")

    write_top(args.outdir / "top_sites.abs_log10_pi_ratio.tsv", top_abs_lr, "abs_log10_pi_ratio")
    write_top(args.outdir / "top_sites.dxy.tsv", top_dxy, "dxy")
    write_top(args.outdir / "top_sites.fst.tsv", top_fst, "fst")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
