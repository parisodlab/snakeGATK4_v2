#!/usr/bin/env python
"""Manhattan-style genome-wide plots for windowed selection metrics.

Reads the merged windows table produced by scripts/scan_selection_windows.py
(e.g. results/selection_windows/DAS_vs_DAG/windows_DAS_vs_DAG.tsv.gz) and
creates a multi-panel figure with one panel per metric.

Example:
  conda run -n python_plot python scripts/plot_selection_manhattan.py \
    --windows results/selection_windows/DAS_vs_DAG/windows_DAS_vs_DAG.tsv.gz \
    --summary results/selection_windows/DAS_vs_DAG/genomewide_summary.json \
    --out-prefix results/selection_windows/DAS_vs_DAG/manhattan/DAS_vs_DAG
"""

from __future__ import annotations

import argparse
import json
import math
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple


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


def load_genes_from_gff(gff_path: Optional[Path], feature_types: Tuple[str, ...] = ("gene",)):
    """Load gene intervals from a GFF into a per-chrom list.

    Returns: dict chrom -> list of (start, end, name)
    Coordinates are 1-based inclusive as in GFF.
    """
    genes: Dict[str, List[Tuple[int, int, str]]] = {}
    if gff_path is None or not Path(gff_path).exists():
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


@dataclass(frozen=True)
class MetricSummary:
    mean: float
    q025: float
    q975: float


def _load_summary(path: Optional[Path]) -> Dict[str, MetricSummary]:
    if path is None:
        return {}
    with open(path) as fh:
        raw = json.load(fh)
    out: Dict[str, MetricSummary] = {}
    for k, v in raw.get("metrics", {}).items():
        out[k] = MetricSummary(mean=v.get("mean"), q025=v.get("q025"), q975=v.get("q975"))
    return out


def _read_table(path: Path):
    import pandas as pd

    return pd.read_csv(path, sep="\t", compression="infer", low_memory=False)


def _chrom_sort_key(chrom: str):
    """Natural-ish sorting: keep numeric suffixes ordered."""
    if chrom is None:
        return (2, "")
    s = str(chrom)
    m = re.match(r"^([A-Za-z_\-\.]+?)(\d+)$", s)
    if m:
        prefix, num = m.group(1), int(m.group(2))
        return (0, prefix, num)
    m = re.match(r"^(chr|chrom|scaffold)[_\-]?(\d+)$", s, flags=re.IGNORECASE)
    if m:
        return (0, m.group(1).lower(), int(m.group(2)))
    return (1, s)


def build_genome_axis(df, gap_bp: int = 0):
    """Add genome-wide x coordinate by laying chromosomes end-to-end."""
    import pandas as pd

    chrom_lengths = (
        df.groupby("chrom", as_index=False)["end"].max().rename(columns={"end": "chrom_len"})
    )
    chrom_lengths = chrom_lengths.sort_values("chrom", key=lambda s: s.map(_chrom_sort_key))

    offsets: Dict[str, int] = {}
    x = 0
    for row in chrom_lengths.itertuples(index=False):
        chrom = getattr(row, "chrom")
        offsets[str(chrom)] = int(x)
        x += int(getattr(row, "chrom_len")) + int(gap_bp)

    df = df.copy()
    df["chrom"] = df["chrom"].astype(str)
    df["offset"] = df["chrom"].map(offsets).astype("int64")
    df["genome_x"] = df["mid"].astype(float) + df["offset"].astype(float)

    # chromosome centers for ticks
    chrom_centers: List[Tuple[str, float]] = []
    for row in chrom_lengths.itertuples(index=False):
        chrom = str(getattr(row, "chrom"))
        clen = float(getattr(row, "chrom_len"))
        chrom_centers.append((chrom, offsets[chrom] + clen / 2.0))

    # chromosome boundaries for vertical lines
    chrom_bounds: List[Tuple[str, float]] = []
    for row in chrom_lengths.itertuples(index=False):
        chrom = str(getattr(row, "chrom"))
        chrom_bounds.append((chrom, float(offsets[chrom])))

    return df, chrom_centers, chrom_bounds


def filter_chromosomes(
    df,
    keep_chroms: Optional[List[str]] = None,
    chrom_prefix: Optional[str] = None,
    chrom_min: Optional[int] = None,
    chrom_max: Optional[int] = None,
):
    """Filter dataframe to a subset of chromosomes.

    - keep_chroms: explicit list (exact match after stringification)
    - chrom_prefix + chrom_min/chrom_max: generate chrom names like f"{prefix}{i}"
    """
    if keep_chroms is None and chrom_prefix is None:
        return df

    d = df.copy()
    d["chrom"] = d["chrom"].astype(str)

    if keep_chroms is not None:
        ks = {str(c) for c in keep_chroms}
        return d[d["chrom"].isin(ks)].copy()

    if chrom_prefix is not None:
        if chrom_min is None or chrom_max is None:
            raise SystemExit("--chrom-prefix requires --chrom-min and --chrom-max")
        ks = {f"{chrom_prefix}{i}" for i in range(int(chrom_min), int(chrom_max) + 1)}
        return d[d["chrom"].isin(ks)].copy()

    return d


def _finite_minmax(arr) -> Optional[Tuple[float, float]]:
    import numpy as np

    a = np.asarray(arr, dtype=float)
    a = a[np.isfinite(a)]
    if a.size == 0:
        return None
    return float(a.min()), float(a.max())


def _padded_limits(lo: float, hi: float, pad_frac: float = 0.08, clamp: Optional[Tuple[float, float]] = None):
    if lo is None or hi is None:
        return None
    if any((isinstance(x, float) and math.isnan(x)) for x in [lo, hi]):
        return None
    if hi < lo:
        lo, hi = hi, lo
    rng = hi - lo
    pad = (rng * pad_frac) if rng > 0 else max(1e-6, abs(lo) * 0.1)
    ymin = lo - pad
    ymax = hi + pad
    if clamp is not None:
        c0, c1 = clamp
        ymin = max(c0, ymin)
        ymax = min(c1, ymax)
    return ymin, ymax


def _metric_ylim(df, col: str, clamp: Optional[Tuple[float, float]] = None, quantile_clip: Optional[Tuple[float, float]] = None):
    import numpy as np

    v = df[col].to_numpy(dtype=float)
    v = v[np.isfinite(v)]
    if v.size == 0:
        return None
    if quantile_clip is not None:
        qlo, qhi = quantile_clip
        lo = float(np.quantile(v, qlo))
        hi = float(np.quantile(v, qhi))
    else:
        lo = float(v.min())
        hi = float(v.max())
    return _padded_limits(lo, hi, pad_frac=0.08, clamp=clamp)


def _finite_quantiles_from_df(df, col: str) -> Optional[Tuple[float, float]]:
    import numpy as np

    if col not in df.columns:
        return None
    v = df[col].to_numpy(dtype=float)
    v = v[np.isfinite(v)]
    if v.size == 0:
        return None
    return float(np.quantile(v, 0.025)), float(np.quantile(v, 0.975))


def _get_q025_q975(df, col: str, summaries: Dict[str, MetricSummary]) -> Optional[Tuple[float, float]]:
    s = summaries.get(col)
    if s and all(isinstance(getattr(s, k), (int, float)) and math.isfinite(getattr(s, k)) for k in ["q025", "q975"]):
        return float(s.q025), float(s.q975)
    return _finite_quantiles_from_df(df, col)


def _other_panel_outlier_count(df, summaries: Dict[str, MetricSummary]):
    """Count (per-row) how many *non-Fst* panels are outliers.

    For Manhattan plots, "outlier" is defined to match the grey CI band:
    a value is an outlier if it falls outside [q025, q975] (two-sided).
    """
    import numpy as np

    specs = [
        "dxy",
        "pi_high",
        "pi_low",
        "log10_pi_ratio",
        "tajima_high",
        "tajima_low",
        "tajima_diff",
    ]

    counts = np.zeros(len(df), dtype=int)
    for col in specs:
        if col not in df.columns:
            continue
        qs = _get_q025_q975(df, col=col, summaries=summaries)
        if qs is None:
            continue
        q025, q975 = qs
        if not (math.isfinite(q025) and math.isfinite(q975)):
            continue

        v = df[col].to_numpy(dtype=float)
        m = np.isfinite(v)
        o = m & ((v <= q025) | (v >= q975))
        counts += o.astype(int)
    return counts


def plot_manhattan(
    df,
    summaries: Dict[str, MetricSummary],
    out_png: Path,
    out_pdf: Path,
    title: str,
    point_size: float,
    alpha: float,
    gap_bp: int,
    quantile_clip: Optional[Tuple[float, float]],
    highlight_fst_outliers: bool,
    fst_outlier_threshold: Optional[float],
    fst_highlight_min_other_outliers: int,
    fst_outlier_size_mult: float,
    fst_outlier_color: str,
    genes_by_chrom: Optional[Dict[str, List[Tuple[int, int, str]]]] = None,
    gene_pad_bp: int = 5000,
):
    import numpy as np
    import matplotlib.pyplot as plt

    df, chrom_centers, chrom_bounds = build_genome_axis(df, gap_bp=gap_bp)

    metrics = [
        ("fst", "Fst", (0.0, 1.0)),
        ("dxy", "Dxy", (0.0, float("inf"))),
        ("pi_high", "π (high)", (0.0, float("inf"))),
        ("pi_low", "π (low)", (0.0, float("inf"))),
        ("log10_pi_ratio", "log10(πH/πL)", None),
        ("tajima_high", "Tajima D (high)", None),
        ("tajima_low", "Tajima D (low)", None),
        ("tajima_diff", "DH−DL", None),
    ]

    n = len(metrics)
    fig, axes = plt.subplots(
        nrows=n,
        ncols=1,
        sharex=True,
        figsize=(14, 2.0 + 1.25 * n),
        gridspec_kw={"hspace": 0.05},
    )

    # alternating colors by chrom
    chroms = [c for c, _ in chrom_centers]
    chrom_to_color = {c: ("#4C72B0" if i % 2 == 0 else "#55A868") for i, c in enumerate(chroms)}
    colors = df["chrom"].map(chrom_to_color).fillna("#4C72B0").to_numpy()

    x = df["genome_x"].to_numpy(dtype=float)

    other_outlier_count = None
    if highlight_fst_outliers and fst_highlight_min_other_outliers > 0:
        other_outlier_count = _other_panel_outlier_count(df, summaries=summaries)

    # Precompute per-metric outlier masks (outside grey CI band) for panels that can support an Fst highlight.
    supporting_metrics = [
        "dxy",
        "pi_high",
        "pi_low",
        "log10_pi_ratio",
        "tajima_high",
        "tajima_low",
        "tajima_diff",
    ]
    metric_outlier_masks: Dict[str, object] = {}
    for col in supporting_metrics:
        if col not in df.columns:
            continue
        qs = _get_q025_q975(df, col=col, summaries=summaries)
        if qs is None:
            continue
        q025, q975 = qs
        if not (math.isfinite(q025) and math.isfinite(q975)):
            continue
        yy = df[col].to_numpy(dtype=float)
        metric_outlier_masks[col] = (np.isfinite(yy) & ((yy <= q025) | (yy >= q975)))

    # Determine fst outlier threshold
    fst_thr = None
    if highlight_fst_outliers:
        if fst_outlier_threshold is not None:
            fst_thr = float(fst_outlier_threshold)
        else:
            s = summaries.get("fst")
            if s and isinstance(s.q975, (int, float)) and math.isfinite(s.q975):
                fst_thr = float(s.q975)
            elif "fst" in df.columns:
                import numpy as np

                v = df["fst"].to_numpy(dtype=float)
                v = v[np.isfinite(v)]
                if v.size:
                    fst_thr = float(np.quantile(v, 0.975))

    # Determine which windows are *actually highlighted* in the Fst panel
    highlight_mask = None
    if highlight_fst_outliers and fst_thr is not None and "fst" in df.columns:
        yfst = df["fst"].to_numpy(dtype=float)
        base = np.isfinite(x) & np.isfinite(yfst) & (yfst >= float(fst_thr))
        if other_outlier_count is not None:
            base = base & (other_outlier_count >= int(fst_highlight_min_other_outliers))
        highlight_mask = base

    # Styling for supporting highlights (other panels)
    support_size_mult = max(4.0, float(fst_outlier_size_mult) * 0.55)

    for ax, (col, ylabel, clamp) in zip(axes, metrics):
        if col not in df.columns:
            ax.text(0.01, 0.5, f"(missing column: {col})", transform=ax.transAxes, ha="left", va="center", color="0.4")
            ax.set_ylabel(ylabel)
            continue

        y = df[col].to_numpy(dtype=float)
        m = np.isfinite(x) & np.isfinite(y)
        ax.scatter(x[m], y[m], c=colors[m], s=point_size, alpha=alpha, linewidths=0, rasterized=True)

        # Overlay: highlight supporting outliers in this panel for the same windows highlighted in Fst
        if col != "fst" and highlight_mask is not None and col in metric_outlier_masks:
            mo = m & highlight_mask & metric_outlier_masks[col]
            if mo.any():
                ax.scatter(
                    x[mo],
                    y[mo],
                    c=fst_outlier_color,
                    s=point_size * support_size_mult,
                    alpha=min(1.0, alpha + 0.15),
                    edgecolors="k",
                    linewidths=0.20,
                    rasterized=True,
                    zorder=5,
                )

        # Overlay: Fst outliers highlighted with bigger dots
        if col == "fst" and highlight_fst_outliers and fst_thr is not None:
            mo = highlight_mask if highlight_mask is not None else (m & (y >= fst_thr))
            if mo.any():
                ax.scatter(
                    x[mo],
                    y[mo],
                    c=fst_outlier_color,
                    s=point_size * float(fst_outlier_size_mult),
                    alpha=min(1.0, alpha + 0.15),
                    edgecolors="k",
                    linewidths=0.25,
                    rasterized=True,
                    zorder=5,
                )
                ax.axhline(fst_thr, color=fst_outlier_color, lw=0.9, ls=":")

                # Gene flags for highlighted outliers (if GFF provided and overlap exists)
                idxs = np.flatnonzero(mo)
                for j, i in enumerate(idxs.tolist()):
                    label = "intergenic"
                    if genes_by_chrom:
                        try:
                            chrom = str(df.iloc[i]["chrom"])
                            wstart = int(df.iloc[i]["start"])
                            wend = int(df.iloc[i]["end"])
                            lo = max(1, wstart - int(gene_pad_bp))
                            hi = wend + int(gene_pad_bp)
                        except Exception:
                            chrom = None
                            lo = None
                            hi = None
                        if chrom is not None and lo is not None and hi is not None:
                            hits = genes_overlapping(genes_by_chrom, chrom=chrom, start=lo, end=hi)
                            names = []
                            for _s, _e, n in hits:
                                if n not in names:
                                    names.append(n)
                            if names:
                                # Keep label compact but informative
                                if len(names) > 3:
                                    label = ",".join(names[:3]) + f",+{len(names)-3}"
                                else:
                                    label = ",".join(names)

                    # Stagger offsets a bit to reduce label collisions
                    yoff = 16 + 7 * (j % 7)
                    ax.annotate(
                        label,
                        xy=(x[i], y[i]),
                        xytext=(0, yoff),
                        textcoords="offset points",
                        ha="center",
                        va="bottom",
                        fontsize=7,
                        rotation=90,
                        color="0.15",
                        arrowprops={"arrowstyle": "-", "lw": 0.7, "color": "0.15"},
                        zorder=6,
                        annotation_clip=False,
                    )

        s = summaries.get(col)
        if s and all(isinstance(getattr(s, k), (int, float)) and math.isfinite(getattr(s, k)) for k in ["q025", "q975"]):
            ax.axhspan(s.q025, s.q975, color="0.85", alpha=0.8, zorder=0)
        if s and isinstance(s.mean, (int, float)) and math.isfinite(s.mean):
            ax.axhline(s.mean, color="0.35", lw=0.9, ls="--")

        ylim = _metric_ylim(df, col=col, clamp=clamp, quantile_clip=quantile_clip)
        if ylim:
            ax.set_ylim(*ylim)

        ax.set_ylabel(ylabel)
        ax.spines["top"].set_visible(False)

        # chromosome boundary lines
        for _c, xb in chrom_bounds:
            ax.axvline(xb, color="0.9", lw=0.6, zorder=0)

    # X ticks as chromosome names
    tick_pos = [p for _c, p in chrom_centers]
    tick_lab = [c for c, _p in chrom_centers]
    axes[-1].set_xticks(tick_pos)
    axes[-1].set_xticklabels(tick_lab, rotation=0, fontsize=8)
    axes[-1].set_xlabel("Genome position (chromosomes concatenated)")

    fig.suptitle(title, y=0.995, fontsize=12)

    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=220, bbox_inches="tight")
    fig.savefig(out_pdf, bbox_inches="tight")
    plt.close(fig)


def main(argv: Optional[List[str]] = None) -> int:
    p = argparse.ArgumentParser(description="Create multi-panel Manhattan plots for windowed selection metrics")
    p.add_argument("--windows", type=Path, required=True, help="windows_{HIGH}_vs_{LOW}.tsv.gz")
    p.add_argument("--summary", type=Path, default=None, help="genomewide_summary.json (optional; adds mean/CI bands)")
    p.add_argument("--out-prefix", type=Path, required=True, help="Output prefix (writes .png and .pdf)")
    p.add_argument("--title", default=None, help="Figure title (default: derived from filename)")
    p.add_argument("--point-size", type=float, default=2.0, help="Scatter marker size (default: 2.0)")
    p.add_argument("--alpha", type=float, default=0.7, help="Point alpha (default: 0.7)")
    p.add_argument("--gap-bp", type=int, default=0, help="Gap between chromosomes in bp (default: 0)")
    p.add_argument(
        "--chrom",
        default=None,
        help="Comma-separated chromosome list to keep (exact match), e.g. 'Bv1,Bv2,...,Bv9'",
    )
    # Back-compat alias
    p.add_argument(
        "--keep-chroms",
        default=None,
        help="Alias for --chrom (comma-separated chromosome list)",
    )
    p.add_argument(
        "--chrom-prefix",
        default=None,
        help="Keep chromosomes named like PREFIX+N using --chrom-min/--chrom-max (e.g. prefix 'Bv')",
    )
    p.add_argument("--chrom-min", type=int, default=None, help="Min chromosome number for --chrom-prefix")
    p.add_argument("--chrom-max", type=int, default=None, help="Max chromosome number for --chrom-prefix")
    p.add_argument("--bv1-bv9", action="store_true", help="Shortcut for --chrom-prefix Bv --chrom-min 1 --chrom-max 9")

    p.add_argument(
        "--highlight-fst-outliers",
        action="store_true",
        help="Highlight Fst outlier windows in the Fst panel with bigger dots",
    )
    p.add_argument(
        "--fst-highlight-min-other-outliers",
        type=int,
        default=None,
        help=(
            "When highlighting Fst outliers, require the same window to also be an outlier in at least N other panels "
            "(default: 3 when chromosome filtering is used; otherwise 0)"
        ),
    )
    p.add_argument(
        "--fst-outlier-threshold",
        type=float,
        default=None,
        help="Override Fst outlier threshold (default: summary fst q975, else empirical 97.5th percentile)",
    )
    p.add_argument(
        "--fst-outlier-size-mult",
        type=float,
        default=8.0,
        help="Size multiplier for highlighted Fst outliers (default: 8.0)",
    )
    p.add_argument(
        "--fst-outlier-color",
        default="#d62728",
        help="Color for highlighted Fst outliers (default: #d62728)",
    )
    p.add_argument(
        "--gff",
        type=Path,
        default=Path("data/reference/Varia.all.maker.renamed.gff"),
        help="GFF for gene overlap annotation of highlighted outliers (default: data/reference/Varia.all.maker.renamed.gff)",
    )
    p.add_argument(
        "--gene-pad-bp",
        type=int,
        default=5000,
        help="Pad around window when checking gene overlaps (default: 5000)",
    )
    p.add_argument(
        "--ylim-quantiles",
        type=str,
        default=None,
        help="Optional y-limit clipping quantiles like '0.001,0.999' (default: no clipping; show extremes)",
    )

    args = p.parse_args(argv)

    df = _read_table(args.windows)
    summaries = _load_summary(args.summary)

    keep_chroms = None
    if args.bv1_bv9:
        args.chrom_prefix = "Bv"
        args.chrom_min = 1
        args.chrom_max = 9

    chrom_arg = args.chrom if args.chrom is not None else args.keep_chroms
    if chrom_arg:
        # preserve order but deduplicate
        seen = set()
        keep_chroms = []
        for c in [x.strip() for x in str(chrom_arg).split(",") if x.strip()]:
            if c not in seen:
                keep_chroms.append(c)
                seen.add(c)
    df = filter_chromosomes(
        df,
        keep_chroms=keep_chroms,
        chrom_prefix=args.chrom_prefix,
        chrom_min=args.chrom_min,
        chrom_max=args.chrom_max,
    )
    if len(df) == 0:
        raise SystemExit("No rows remain after chromosome filtering")

    if args.title is None:
        args.title = f"Manhattan metrics: {args.windows.name}"

    quantile_clip = None
    if args.ylim_quantiles:
        parts = [p.strip() for p in args.ylim_quantiles.split(",") if p.strip()]
        if len(parts) != 2:
            raise SystemExit("--ylim-quantiles must be like '0.001,0.999'")
        qlo, qhi = float(parts[0]), float(parts[1])
        if not (0.0 <= qlo < qhi <= 1.0):
            raise SystemExit("--ylim-quantiles must satisfy 0 <= qlo < qhi <= 1")
        quantile_clip = (qlo, qhi)

    out_png = Path(str(args.out_prefix) + ".png")
    out_pdf = Path(str(args.out_prefix) + ".pdf")

    # Default behavior requested: with chromosome filtering, only highlight Fst outliers that are also outliers in >=3 other panels.
    # Keep backward-compatible behavior (highlight all Fst outliers) when not filtering unless explicitly overridden.
    chrom_filtering_used = bool(keep_chroms) or bool(args.chrom_prefix) or bool(args.bv1_bv9)
    if args.fst_highlight_min_other_outliers is None:
        fst_min_other = 3 if chrom_filtering_used else 0
    else:
        fst_min_other = int(args.fst_highlight_min_other_outliers)
        if fst_min_other < 0:
            raise SystemExit("--fst-highlight-min-other-outliers must be >= 0")

    genes_by_chrom = load_genes_from_gff(args.gff)

    plot_manhattan(
        df=df,
        summaries=summaries,
        out_png=out_png,
        out_pdf=out_pdf,
        title=args.title,
        point_size=args.point_size,
        alpha=args.alpha,
        gap_bp=args.gap_bp,
        quantile_clip=quantile_clip,
        highlight_fst_outliers=bool(args.highlight_fst_outliers),
        fst_outlier_threshold=args.fst_outlier_threshold,
        fst_highlight_min_other_outliers=fst_min_other,
        fst_outlier_size_mult=args.fst_outlier_size_mult,
        fst_outlier_color=str(args.fst_outlier_color),
        genes_by_chrom=genes_by_chrom,
        gene_pad_bp=int(args.gene_pad_bp),
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
