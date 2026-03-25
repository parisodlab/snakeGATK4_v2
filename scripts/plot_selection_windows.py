#!/usr/bin/env python
"""Plot stacked selection-scan panels around top windows.

Inputs are produced by scripts/scan_selection_windows.py:
- windows_{HIGH}_vs_{LOW}.tsv.gz
- genomewide_summary.json
- top_windows.<selection_type>.tsv

Produces per-window plots as PNG + PDF with 4 stacked panels:
1) Fst
2) pi (high + low) with genome-wide mean + 95% interval
3) log10(pi_high/pi_low) with genome-wide 95% interval
4) Tajima's D (high + low) with DH-DL (secondary axis) + 95% interval

Run with a conda env that has pandas/numpy/matplotlib:
  conda run -n python_plot python scripts/plot_selection_windows.py ...
"""

from __future__ import annotations

import argparse
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple


@dataclass(frozen=True)
class MetricSummary:
    mean: float
    q025: float
    q975: float


def _parse_gff_attributes(attr: str) -> Dict[str, str]:
    out: Dict[str, str] = {}
    for part in (attr or "").split(";"):
        part = part.strip()
        if not part:
            continue
        if "=" in part:
            k, v = part.split("=", 1)
            out[k] = v
        elif " " in part:
            # tolerate GTF-ish format key "value"
            k, v = part.split(" ", 1)
            out[k] = v.strip().strip('"')
    return out


def load_genes_from_gff(gff_path: Path, feature_types: Tuple[str, ...] = ("gene",)):
    """Load gene intervals from a GFF into a per-chrom list.

    Returns: dict chrom -> list of (start, end, name)
    Coordinates are 1-based inclusive as in GFF.
    """
    genes: Dict[str, List[Tuple[int, int, str]]] = {}
    if gff_path is None or not gff_path.exists():
        return genes
    with open(gff_path, "rt") as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            chrom, _src, ftype, start_s, end_s, _score, _strand, _phase, attrs = parts
            if ftype not in feature_types:
                continue
            try:
                start = int(start_s)
                end = int(end_s)
            except Exception:
                continue
            if end < start:
                start, end = end, start
            d = _parse_gff_attributes(attrs)
            name = d.get("Name") or d.get("gene") or d.get("gene_name") or d.get("ID") or "gene"
            genes.setdefault(chrom, []).append((start, end, name))

    # sort per chrom
    for c in genes:
        genes[c].sort(key=lambda x: (x[0], x[1], x[2]))
    return genes


def genes_overlapping(genes_by_chrom: Dict[str, List[Tuple[int, int, str]]], chrom: str, start: int, end: int):
    """Return list of (start,end,name) overlapping [start,end] (1-based)."""
    if chrom not in genes_by_chrom:
        return []
    hits = []
    for s, e, n in genes_by_chrom[chrom]:
        if e < start:
            continue
        if s > end:
            break
        hits.append((s, e, n))
    return hits


def _padded_limits(lo: float, hi: float, pad_frac: float = 0.25, clamp: Optional[Tuple[float, float]] = None):
    import math

    if lo is None or hi is None:
        return None
    if any(math.isnan(x) for x in [lo, hi]):
        return None
    if hi < lo:
        lo, hi = hi, lo
    rng = hi - lo
    if rng == 0:
        pad = max(1e-6, abs(lo) * 0.1)
    else:
        pad = rng * pad_frac
    ymin = lo - pad
    ymax = hi + pad
    if clamp is not None:
        c0, c1 = clamp
        ymin = max(c0, ymin)
        ymax = min(c1, ymax)
    return ymin, ymax


def _load_summary(path: Path) -> Dict[str, MetricSummary]:
    with open(path) as fh:
        raw = json.load(fh)
    out: Dict[str, MetricSummary] = {}
    for k, v in raw.get("metrics", {}).items():
        out[k] = MetricSummary(mean=v.get("mean"), q025=v.get("q025"), q975=v.get("q975"))
    return out


def _read_table(path: Path):
    import pandas as pd

    return pd.read_csv(path, sep="\t", compression="infer", low_memory=False)


def _finite_minmax(values):
    import numpy as np

    arr = np.asarray(values, dtype=float)
    arr = arr[np.isfinite(arr)]
    if arr.size == 0:
        return None
    return float(arr.min()), float(arr.max())


def compute_global_ylims(windows) -> Dict[str, Tuple[float, float]]:
    """Compute common y-limits across all plots from the full windows table.

    This avoids clipping outlier windows in individual plots.
    """
    ylims: Dict[str, Tuple[float, float]] = {}

    if "fst" in windows.columns:
        mm = _finite_minmax(windows["fst"].to_numpy())
        if mm:
            lo, hi = mm
            lo = max(0.0, lo)
            hi = min(1.0, hi)
            lim = _padded_limits(lo, hi, pad_frac=0.1, clamp=(0.0, 1.0))
            if lim:
                ylims["fst"] = lim

    if "dxy" in windows.columns:
        mm = _finite_minmax(windows["dxy"].to_numpy())
        if mm:
            lo, hi = mm
            lo = max(0.0, lo)
            lim = _padded_limits(lo, hi, pad_frac=0.1, clamp=(0.0, float("inf")))
            if lim:
                ylims["dxy"] = lim

    if "pi_high" in windows.columns and "pi_low" in windows.columns:
        mm_h = _finite_minmax(windows["pi_high"].to_numpy())
        mm_l = _finite_minmax(windows["pi_low"].to_numpy())
        if mm_h or mm_l:
            lo = min([x for x in [mm_h[0] if mm_h else None, mm_l[0] if mm_l else None] if x is not None])
            hi = max([x for x in [mm_h[1] if mm_h else None, mm_l[1] if mm_l else None] if x is not None])
            lo = max(0.0, lo)
            lim = _padded_limits(lo, hi, pad_frac=0.1, clamp=(0.0, float("inf")))
            if lim:
                ylims["pi"] = lim

    if "log10_pi_ratio" in windows.columns:
        mm = _finite_minmax(windows["log10_pi_ratio"].to_numpy())
        if mm:
            lim = _padded_limits(mm[0], mm[1], pad_frac=0.1)
            if lim:
                ylims["log10_pi_ratio"] = lim

    if "tajima_high" in windows.columns and "tajima_low" in windows.columns:
        mm_h = _finite_minmax(windows["tajima_high"].to_numpy())
        mm_l = _finite_minmax(windows["tajima_low"].to_numpy())
        if mm_h or mm_l:
            lo = min([x for x in [mm_h[0] if mm_h else None, mm_l[0] if mm_l else None] if x is not None])
            hi = max([x for x in [mm_h[1] if mm_h else None, mm_l[1] if mm_l else None] if x is not None])
            lim = _padded_limits(lo, hi, pad_frac=0.1)
            if lim:
                ylims["tajima"] = lim

    if "tajima_diff" in windows.columns:
        mm = _finite_minmax(windows["tajima_diff"].to_numpy())
        if mm:
            lim = _padded_limits(mm[0], mm[1], pad_frac=0.1)
            if lim:
                ylims["tajima_diff"] = lim

    return ylims


def _region_subset(windows, chrom: str, start: int, end: int, flank_bp: int):
    lo = max(0, start - flank_bp)
    hi = end + flank_bp
    sub = windows[(windows["chrom"] == chrom) & (windows["mid"] >= lo) & (windows["mid"] <= hi)].copy()
    return sub.sort_values("mid")


def _format_kb(x_bp: float) -> float:
    return float(x_bp) / 1000.0


def plot_one_region(
    sub,
    summaries: Dict[str, MetricSummary],
    ylims: Dict[str, Tuple[float, float]],
    pop_high: str,
    pop_low: str,
    chrom: str,
    w_start: int,
    w_end: int,
    selection_type: str,
    out_png: Path,
    out_pdf: Path,
    title: str,
    genes_by_chrom: Optional[Dict[str, List[Tuple[int, int, str]]]] = None,
    gene_pad_bp: int = 5000,
    max_gene_labels: int = 8,
):
    import numpy as np
    import matplotlib.pyplot as plt

    x_kb = sub["mid"].astype(float).apply(_format_kb).to_numpy()
    # identify the exact selected window within the plotted region
    sel = (sub["start"].astype(int) == int(w_start)) & (sub["end"].astype(int) == int(w_end))

    fig, axes = plt.subplots(
        nrows=7,
        ncols=1,
        sharex=True,
        figsize=(9, 10.2),
        gridspec_kw={"height_ratios": [1.0, 1.0, 1.1, 1.0, 1.1, 1.0, 0.45], "hspace": 0.05},
    )

    # Common highlight span
    span_lo = _format_kb(w_start)
    span_hi = _format_kb(w_end)
    for ax in axes:
        ax.axvspan(span_lo, span_hi, color="#f5e642", alpha=0.35, lw=0)

    # Panel 1: Fst
    ax = axes[0]
    if "fst" in sub.columns:
        y = sub["fst"].astype(float).to_numpy()
        ax.scatter(x_kb, y, color="black", s=10, alpha=0.75)
        if ("fst" in selection_type) or (selection_type in ("high_fst", "sweep_high", "sweep_low", "combined_abs")):
            ax.scatter(x_kb[sel.to_numpy()], y[sel.to_numpy()], color="#d62728", s=45, zorder=5)
    s = summaries.get("fst")
    if s:
        ax.axhspan(s.q025, s.q975, color="0.85", alpha=0.7, zorder=0)
    if "fst" in ylims:
        ax.set_ylim(*ylims["fst"])
    ax.set_ylabel("Fst")

    # Panel 2: Dxy
    ax = axes[1]
    if "dxy" in sub.columns:
        y = sub["dxy"].astype(float).to_numpy()
        ax.scatter(x_kb, y, color="black", s=10, alpha=0.75)
        if ("dxy" in selection_type):
            ax.scatter(x_kb[sel.to_numpy()], y[sel.to_numpy()], color="#d62728", s=45, zorder=5)
    s = summaries.get("dxy")
    if s:
        ax.axhspan(s.q025, s.q975, color="0.85", alpha=0.7, zorder=0)
    if "dxy" in ylims:
        ax.set_ylim(*ylims["dxy"])
    ax.set_ylabel("Dxy")

    # Panel 3: Pi
    ax = axes[2]
    y_high = sub["pi_high"].astype(float).to_numpy()
    y_low = sub["pi_low"].astype(float).to_numpy()
    ax.scatter(x_kb, y_high, color="#1f77b4", s=10, alpha=0.75, label=f"{pop_high} (high)")
    ax.scatter(x_kb, y_low, color="#d62728", s=10, alpha=0.75, label=f"{pop_low} (low)")
    if selection_type in ("sweep_high", "sweep_low"):
        ax.scatter(x_kb[sel.to_numpy()], y_high[sel.to_numpy()], color="#d62728", s=45, zorder=5)
        ax.scatter(x_kb[sel.to_numpy()], y_low[sel.to_numpy()], color="#d62728", s=45, zorder=5)

    # genome-wide mean + CI (single fixed band + mean line)
    s_high = summaries.get("pi_high")
    s_low = summaries.get("pi_low")
    if s_high and s_low:
        lo = min(s_high.q025, s_low.q025)
        hi = max(s_high.q975, s_low.q975)
        ax.axhspan(lo, hi, color="0.85", alpha=0.7, zorder=0)
        pooled_mean = (s_high.mean + s_low.mean) / 2.0
        ax.axhline(pooled_mean, color="0.25", ls="--", lw=1.0, label="Genome avg π")
    if "pi" in ylims:
        ax.set_ylim(*ylims["pi"])

    ax.set_ylabel("π")
    ax.legend(loc="upper right", frameon=False, fontsize=8)

    # Panel 4: log10 pi ratio
    ax = axes[3]
    y = sub["log10_pi_ratio"].astype(float).to_numpy()
    ax.scatter(x_kb, y, color="#2ca02c", s=10, alpha=0.8)
    if selection_type in ("pi_ratio_abs", "combined_abs"):
        ax.scatter(x_kb[sel.to_numpy()], y[sel.to_numpy()], color="#d62728", s=45, zorder=5)
    ax.axhline(0.0, color="0.4", lw=0.8, ls=":")
    s = summaries.get("log10_pi_ratio")
    if s:
        ax.axhspan(s.q025, s.q975, color="0.85", alpha=0.7, zorder=0)
    if "log10_pi_ratio" in ylims:
        ax.set_ylim(*ylims["log10_pi_ratio"])
    ax.set_ylabel("log10(πH/πL)")

    # Panel 5: Tajima D per-pop
    ax = axes[4]
    y_high = sub["tajima_high"].astype(float).to_numpy()
    y_low = sub["tajima_low"].astype(float).to_numpy()
    ax.scatter(x_kb, y_high, color="#1f77b4", s=10, alpha=0.55, label=f"{pop_high} (high)")
    ax.scatter(x_kb, y_low, color="#d62728", s=10, alpha=0.55, label=f"{pop_low} (low)")
    ax.axhline(0.0, color="0.4", lw=0.8, ls=":")
    ax.set_ylabel("Tajima D")

    # fixed genome-wide CI band for Tajima D (pooled across pops)
    s_th = summaries.get("tajima_high")
    s_tl = summaries.get("tajima_low")
    if s_th and s_tl:
        lo = min(s_th.q025, s_tl.q025)
        hi = max(s_th.q975, s_tl.q975)
        ax.axhspan(lo, hi, color="0.85", alpha=0.7, zorder=0)
    if "tajima" in ylims:
        ax.set_ylim(*ylims["tajima"])

    # Panel 6: DH-DL (Tajima diff)
    ax = axes[5]
    ydiff = sub["tajima_diff"].astype(float).to_numpy()
    ax.scatter(x_kb, ydiff, color="#2ca02c", s=10, alpha=0.8)
    if selection_type in ("tajima_diff_abs", "sweep_high", "sweep_low", "combined_abs"):
        ax.scatter(x_kb[sel.to_numpy()], ydiff[sel.to_numpy()], color="#d62728", s=45, zorder=5)
    s = summaries.get("tajima_diff")
    if s:
        ax.axhspan(s.q025, s.q975, color="0.85", alpha=0.7, zorder=0)
    if "tajima_diff" in ylims:
        ax.set_ylim(*ylims["tajima_diff"])
    ax.axhline(0.0, color="0.4", lw=0.8, ls=":")
    ax.set_ylabel("DH-DL")

    # Gene panel: show genes overlapping window +/- pad
    gene_ax = axes[6]
    gene_ax.set_ylabel("Genes")
    gene_ax.set_yticks([])
    gene_ax.spines["left"].set_visible(False)
    gene_ax.spines["right"].set_visible(False)
    gene_ax.spines["top"].set_visible(False)
    gene_ax.set_ylim(0.0, 1.0)
    if genes_by_chrom is not None:
        g_start = max(1, int(w_start) - int(gene_pad_bp))
        g_end = int(w_end) + int(gene_pad_bp)
        hits = genes_overlapping(genes_by_chrom, chrom=chrom, start=g_start, end=g_end)
        if hits:
            yrows = [0.25, 0.5, 0.75]
            for j, (gs, ge, name) in enumerate(hits[:max_gene_labels]):
                y = yrows[j % len(yrows)]
                x0 = _format_kb(gs)
                x1 = _format_kb(ge)
                xm = _format_kb((gs + ge) / 2.0)
                gene_ax.plot([x0, x1], [y, y], color="black", lw=4, solid_capstyle="butt")
                gene_ax.text(xm, y + 0.08, name, ha="center", va="bottom", fontsize=7, color="black", rotation=0)
        else:
            gene_ax.text(0.01, 0.5, "(no genes within window±pad)", transform=gene_ax.transAxes, ha="left", va="center", fontsize=7, color="0.35")
    else:
        gene_ax.text(0.01, 0.5, "(no GFF loaded)", transform=gene_ax.transAxes, ha="left", va="center", fontsize=7, color="0.35")

    axes[-1].set_xlabel(f"Position on {chrom} (kb)")
    fig.suptitle(title, y=0.99, fontsize=11)

    for ax in axes:
        ax.spines["top"].set_visible(False)

    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=200, bbox_inches="tight")
    fig.savefig(out_pdf, bbox_inches="tight")
    plt.close(fig)


def main(argv: Optional[List[str]] = None) -> int:
    p = argparse.ArgumentParser(description="Plot stacked panels for top selection windows")
    p.add_argument("--windows", type=Path, required=True, help="windows_{HIGH}_vs_{LOW}.tsv.gz")
    p.add_argument("--summary", type=Path, required=True, help="genomewide_summary.json")
    p.add_argument("--tops", type=Path, default=None, help="top_windows.<type>.tsv")
    p.add_argument("--outdir", type=Path, required=True, help="Output dir for plots")
    p.add_argument("--high-pop", required=True)
    p.add_argument("--low-pop", required=True)
    p.add_argument("--flank-bp", type=int, default=500_000, help="Flank around the window in bp (default: 500000 ~ 1Mbp view)")
    p.add_argument("--max-plots", type=int, default=10, help="Max windows to plot (0 = no limit; default: 10)")
    p.add_argument(
        "--plot-fst-outliers",
        action="store_true",
        help="Plot all Fst outlier windows (fst >= genomewide fst q975 unless overridden)",
    )
    p.add_argument(
        "--fst-outlier-threshold",
        type=float,
        default=None,
        help="Override Fst outlier threshold (default: genomewide summary fst q975)",
    )
    p.add_argument(
        "--plot-dxy-outliers",
        action="store_true",
        help="Plot all Dxy outlier windows (dxy >= genomewide dxy q975 unless overridden)",
    )
    p.add_argument(
        "--dxy-outlier-threshold",
        type=float,
        default=None,
        help="Override Dxy outlier threshold (default: genomewide summary dxy q975)",
    )
    p.add_argument(
        "--gff",
        type=Path,
        default=Path("data/reference/Varia.all.maker.renamed.gff"),
        help="GFF for gene overlap annotation (default: data/reference/Varia.all.maker.renamed.gff)",
    )
    p.add_argument("--gene-pad-bp", type=int, default=5000, help="Pad around window when checking gene overlaps (default: 5000)")
    p.add_argument("--max-gene-labels", type=int, default=8, help="Max gene labels to annotate per plot (default: 8)")

    args = p.parse_args(argv)

    windows = _read_table(args.windows)
    summaries = _load_summary(args.summary)
    ylims = compute_global_ylims(windows)

    genes_by_chrom = None
    if args.gff is not None and args.gff.exists():
        genes_by_chrom = load_genes_from_gff(args.gff, feature_types=("gene",))

    args.outdir.mkdir(parents=True, exist_ok=True)

    if args.plot_fst_outliers and args.plot_dxy_outliers:
        raise SystemExit("Use only one of --plot-fst-outliers or --plot-dxy-outliers per run")

    if args.plot_fst_outliers:
        thr = args.fst_outlier_threshold
        if thr is None:
            thr = summaries.get("fst").q975 if summaries.get("fst") else None
        if thr is None:
            raise SystemExit("Cannot determine Fst outlier threshold; summary is missing fst.q975")
        sel_df = windows.copy()
        sel_df["fst"] = sel_df["fst"].astype(float)
        sel_df = sel_df[sel_df["fst"].notna()]
        sel_df = sel_df[sel_df["fst"] >= float(thr)]
        sel_df = sel_df.sort_values("fst", ascending=False)
        # build a tops-like table
        tops = sel_df[["chrom", "start", "end", "mid", "locus", "fst"]].copy()
        tops = tops.rename(columns={"fst": "score"})
        selection_type = "fst_outliers"
    elif args.plot_dxy_outliers:
        thr = args.dxy_outlier_threshold
        if thr is None:
            thr = summaries.get("dxy").q975 if summaries.get("dxy") else None
        if thr is None:
            raise SystemExit("Cannot determine Dxy outlier threshold; summary is missing dxy.q975")
        sel_df = windows.copy()
        if "dxy" not in sel_df.columns:
            raise SystemExit("windows file has no 'dxy' column; regenerate windows table with scan_selection_windows.py")
        sel_df["dxy"] = sel_df["dxy"].astype(float)
        sel_df = sel_df[sel_df["dxy"].notna()]
        sel_df = sel_df[sel_df["dxy"] >= float(thr)]
        sel_df = sel_df.sort_values("dxy", ascending=False)
        tops = sel_df[["chrom", "start", "end", "mid", "locus", "dxy"]].copy()
        tops = tops.rename(columns={"dxy": "score"})
        selection_type = "dxy_outliers"
    else:
        if args.tops is None:
            raise SystemExit("--tops is required unless --plot-fst-outliers is set")
        tops = _read_table(args.tops)
        # infer selection type from tops filename if available
        selection_type = args.tops.name
        if selection_type.startswith("top_windows.") and selection_type.endswith(".tsv"):
            selection_type = selection_type[len("top_windows.") : -len(".tsv")]

    if args.max_plots > 0:
        tops = tops.head(args.max_plots)
    for i, row in enumerate(tops.itertuples(index=False), start=1):
        chrom = getattr(row, "chrom")
        start = int(getattr(row, "start"))
        end = int(getattr(row, "end"))
        locus = getattr(row, "locus")
        score = getattr(row, "score") if hasattr(row, "score") else None

        sub = _region_subset(windows, chrom=chrom, start=start, end=end, flank_bp=args.flank_bp)

        label_score = "" if score is None else f" | score={score:.3g}"
        title = f"{args.high_pop} vs {args.low_pop} | {selection_type} | {locus}{label_score}"

        stem = f"{selection_type}.top{i:02d}.{chrom}_{start}_{end}"
        out_png = args.outdir / f"{stem}.png"
        out_pdf = args.outdir / f"{stem}.pdf"

        plot_one_region(
            sub=sub,
            summaries=summaries,
            ylims=ylims,
            pop_high=args.high_pop,
            pop_low=args.low_pop,
            chrom=chrom,
            w_start=start,
            w_end=end,
            selection_type=selection_type,
            genes_by_chrom=genes_by_chrom,
            gene_pad_bp=args.gene_pad_bp,
            max_gene_labels=args.max_gene_labels,
            out_png=out_png,
            out_pdf=out_pdf,
            title=title,
        )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
