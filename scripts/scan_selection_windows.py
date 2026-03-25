#!/usr/bin/env python
"""Scan windowed selection signals between two populations.

This script consumes piawka windowed outputs and produces a single per-window
TSV with Fst, pi, Tajima's D, and derived signals between two populations.

It also ranks top windows for a handful of selection-signal "types" and writes
per-type Top-N tables with empirical (genome-wide) tail probabilities.

Expected inputs (from this repo's workflow):
- results/piawka/mixed/window/piawka_main.tsv.gz   (Fst_HUD, ...)
- results/piawka/mixed/window/piawka_tajima.tsv.gz (TajD, pi, Dxy)

Note: the Tajima file does not contain Fst. If you want an Fst-like statistic
derived from Tajima-file stats, use --fst-source nei, which computes
    fst_nei_raw = 1 - mean(pi_high, pi_low) / Dxy
and uses a clipped version (0..1) as the primary 'fst' column.

Run with the repo's conda env that has pandas/numpy:
  conda run -n python_plot python scripts/scan_selection_windows.py ...
"""

from __future__ import annotations

import argparse
import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple


def parse_locus(locus: str) -> Tuple[str, int, int]:
    """Parse piawka locus string like 'Bv1:240001-250000'."""
    if ":" not in locus or "-" not in locus:
        raise ValueError(f"Unrecognized locus format: {locus!r}")
    chrom, rest = locus.split(":", 1)
    start_s, end_s = rest.split("-", 1)
    start = int(start_s)
    end = int(end_s)
    if end < start:
        start, end = end, start
    return chrom, start, end


def zscore(series):
    mean = series.mean(skipna=True)
    std = series.std(skipna=True, ddof=0)
    if std == 0 or math.isnan(std):
        return (series * 0)  # all zeros
    return (series - mean) / std


def empirical_p_high(values, x: float) -> float:
    """One-sided empirical p for high tail: P(X >= x)."""
    values = values.dropna()
    if len(values) == 0 or math.isnan(x):
        return float("nan")
    return float((values >= x).sum() / len(values))


def empirical_p_low(values, x: float) -> float:
    """One-sided empirical p for low tail: P(X <= x)."""
    values = values.dropna()
    if len(values) == 0 or math.isnan(x):
        return float("nan")
    return float((values <= x).sum() / len(values))


def empirical_p_two_sided(values, x: float) -> float:
    """Two-sided empirical p based on absolute deviation from median via ECDF."""
    values = values.dropna()
    if len(values) == 0 or math.isnan(x):
        return float("nan")
    # two-sided around median using tails
    p_hi = float((values >= x).sum() / len(values))
    p_lo = float((values <= x).sum() / len(values))
    return float(2.0 * min(p_hi, p_lo))


@dataclass(frozen=True)
class GenomeSummary:
    mean: float
    q025: float
    q975: float


def summarize(values) -> GenomeSummary:
    values = values.dropna()
    if len(values) == 0:
        return GenomeSummary(mean=float("nan"), q025=float("nan"), q975=float("nan"))
    return GenomeSummary(
        mean=float(values.mean()),
        q025=float(values.quantile(0.025)),
        q975=float(values.quantile(0.975)),
    )


def _read_piawka_main(path: Path):
    import pandas as pd

    usecols = ["locus", "nSites", "pop1", "pop2", "nUsed", "metric", "value"]
    return pd.read_csv(
        path,
        sep="\t",
        usecols=usecols,
        dtype={"locus": "string", "pop1": "string", "pop2": "string", "metric": "string"},
        compression="infer",
        low_memory=False,
    )


def _read_piawka_tajima_filtered(path: Path, pops: Iterable[str], pop_high: str, pop_low: str):
    import pandas as pd

    pops_set = set(pops)
    usecols = ["locus", "nSites", "pop1", "pop2", "nUsed", "metric", "value"]

    parts = []
    # this file can be large; chunk to reduce memory
    for chunk in pd.read_csv(
        path,
        sep="\t",
        usecols=usecols,
        dtype={"locus": "string", "pop1": "string", "pop2": "string", "metric": "string"},
        compression="infer",
        low_memory=False,
        chunksize=500_000,
    ):
        # Keep: TajD + pi for each pop, and Dxy between the two pops.
        keep_taj = (chunk["metric"] == "TajD") & (chunk["pop1"].isin(pops_set))
        keep_pi = (chunk["metric"] == "pi") & (chunk["pop1"].isin(pops_set))
        keep_dxy = (chunk["metric"] == "Dxy") & (
            ((chunk["pop1"] == pop_high) & (chunk["pop2"] == pop_low))
            | ((chunk["pop1"] == pop_low) & (chunk["pop2"] == pop_high))
        )
        chunk = chunk[keep_taj | keep_pi | keep_dxy]
        if len(chunk) > 0:
            parts.append(chunk)
    if parts:
        return pd.concat(parts, ignore_index=True)
    return pd.DataFrame(columns=usecols)


def build_windows(
    piawka_main: Path,
    piawka_tajima: Path,
    pop_high: str,
    pop_low: str,
    min_used_sites: int = 0,
    pi_pseudocount: float = 0.0,
    fst_source: str = "hudson",
):
    import numpy as np
    import pandas as pd

    fst_source = str(fst_source or "hudson").lower()
    if fst_source not in {"hudson", "nei"}:
        raise ValueError("fst_source must be one of: hudson, nei")

    main = _read_piawka_main(piawka_main)

    # Hudson Fst between the two pops (from piawka_main)
    fst_hud = main[main["metric"] == "Fst_HUD"].copy()
    fst_hud = fst_hud[
        ((fst_hud["pop1"] == pop_high) & (fst_hud["pop2"] == pop_low))
        | ((fst_hud["pop1"] == pop_low) & (fst_hud["pop2"] == pop_high))
    ]
    if min_used_sites > 0:
        fst_hud = fst_hud[fst_hud["nUsed"] >= min_used_sites]
    fst_hud = (
        fst_hud[["locus", "value"]]
        .groupby("locus", as_index=False)["value"]
        .mean()
        .rename(columns={"value": "fst_hudson"})
    )

    # Use the Tajima file for pi and Dxy as well (keeps consistency with TajD windows).
    # NOTE: this file does NOT contain Fst; if fst_source=='nei', we derive a proxy Fst from pi and Dxy.
    taj = _read_piawka_tajima_filtered(piawka_tajima, pops=[pop_high, pop_low], pop_high=pop_high, pop_low=pop_low)
    if min_used_sites > 0:
        taj = taj[taj["nUsed"] >= min_used_sites]

    loci = pd.DataFrame({"locus": taj["locus"].dropna().drop_duplicates().astype("string")})

    # Dxy between the two pops (symmetric; average if both orientations exist)
    dxy = taj[taj["metric"] == "Dxy"][["locus", "value"]].copy()
    if len(dxy) > 0:
        dxy["value"] = dxy["value"].astype(float)
        dxy = dxy.groupby("locus", as_index=False)["value"].mean().rename(columns={"value": "dxy"})
    else:
        dxy = pd.DataFrame(columns=["locus", "dxy"])

    # pi per pop
    pi = taj[(taj["metric"] == "pi") & (taj["pop1"].isin([pop_high, pop_low]))].copy()
    pi_high = (
        pi[pi["pop1"] == pop_high][["locus", "value"]]
        .rename(columns={"value": "pi_high"})
        .drop_duplicates(subset=["locus"], keep="first")
    )
    pi_low = (
        pi[pi["pop1"] == pop_low][["locus", "value"]]
        .rename(columns={"value": "pi_low"})
        .drop_duplicates(subset=["locus"], keep="first")
    )

    tajd = taj[taj["metric"] == "TajD"].copy()

    taj_high = (
        tajd[tajd["pop1"] == pop_high][["locus", "value"]]
        .rename(columns={"value": "tajima_high"})
        .drop_duplicates(subset=["locus"], keep="first")
    )
    taj_low = (
        tajd[tajd["pop1"] == pop_low][["locus", "value"]]
        .rename(columns={"value": "tajima_low"})
        .drop_duplicates(subset=["locus"], keep="first")
    )

    windows = loci.merge(dxy, on="locus", how="left")
    windows = windows.merge(pi_high, on="locus", how="left").merge(pi_low, on="locus", how="left")
    windows = windows.merge(taj_high, on="locus", how="left").merge(taj_low, on="locus", how="left")
    windows = windows.merge(fst_hud, on="locus", how="left")

    # parse locus coordinates
    chroms: List[str] = []
    starts: List[int] = []
    ends: List[int] = []
    mids: List[float] = []
    for locus in windows["locus"].astype(str).tolist():
        c, s, e = parse_locus(locus)
        chroms.append(c)
        starts.append(s)
        ends.append(e)
        mids.append((s + e) / 2.0)
    windows.insert(0, "chrom", pd.Series(chroms, dtype="string"))
    windows.insert(1, "start", pd.Series(starts, dtype="int64"))
    windows.insert(2, "end", pd.Series(ends, dtype="int64"))
    windows.insert(3, "mid", pd.Series(mids, dtype="float64"))

    # Fst used for downstream analysis
    # - fst_source='hudson': use fst_hudson from piawka_main
    # - fst_source='nei': derive a proxy from Tajima-file stats: fst = 1 - mean(pi_high,pi_low) / dxy
    fst = np.full(len(windows), np.nan, dtype=float)
    fst_nei = np.full(len(windows), np.nan, dtype=float)
    if fst_source == "hudson":
        if "fst_hudson" in windows.columns:
            fst = windows["fst_hudson"].to_numpy(dtype=float)
    else:
        pi_h = windows["pi_high"].to_numpy(dtype=float)
        pi_l = windows["pi_low"].to_numpy(dtype=float)
        dxy_arr = windows["dxy"].to_numpy(dtype=float)
        pi_within = (pi_h + pi_l) / 2.0
        valid = np.isfinite(pi_within) & np.isfinite(dxy_arr) & (dxy_arr > 0)
        fst_nei[valid] = 1.0 - (pi_within[valid] / dxy_arr[valid])
        # clip to [0,1] for stability/interpretability
        fst = np.clip(fst_nei, 0.0, 1.0)

    windows["fst_nei_raw"] = fst_nei
    windows["fst"] = fst

    # derived signals
    # By default we require pi>0 in both pops; optionally add a small pseudocount to keep 0-PI windows
    # while preventing +/-inf ratios.
    pi_high_arr = windows["pi_high"].to_numpy(dtype=float)
    pi_low_arr = windows["pi_low"].to_numpy(dtype=float)

    pc = float(pi_pseudocount or 0.0)
    if pc < 0:
        raise ValueError("pi_pseudocount must be >= 0")

    ratio = np.full(len(windows), np.nan, dtype=float)
    if pc > 0:
        # symmetric pseudocount: log10((piH+pc)/(piL+pc)); 0/0 -> 0, and very small denominators are stabilized
        a = pi_high_arr + pc
        b = pi_low_arr + pc
        valid = np.isfinite(a) & np.isfinite(b) & (b > 0) & (a > 0)
        np.divide(a, b, out=ratio, where=valid)
    else:
        valid = np.isfinite(pi_high_arr) & np.isfinite(pi_low_arr) & (pi_high_arr > 0) & (pi_low_arr > 0)
        np.divide(pi_high_arr, pi_low_arr, out=ratio, where=valid)

    logratio = np.full(len(windows), np.nan, dtype=float)
    m = (ratio > 0) & np.isfinite(ratio)
    logratio[m] = np.log10(ratio[m])
    windows["log10_pi_ratio"] = logratio
    windows["tajima_diff"] = windows["tajima_high"] - windows["tajima_low"]

    # sort
    windows = windows.sort_values(["chrom", "start"]).reset_index(drop=True)
    return windows


def recompute_log10_pi_ratio_inplace(windows, pi_pseudocount: float):
    """Recompute log10(pi_high/pi_low) in-place using an optional pseudocount."""
    import numpy as np

    pi_high_arr = windows["pi_high"].to_numpy(dtype=float)
    pi_low_arr = windows["pi_low"].to_numpy(dtype=float)
    pc = float(pi_pseudocount or 0.0)
    if pc < 0:
        raise ValueError("pi_pseudocount must be >= 0")

    ratio = np.full(len(windows), np.nan, dtype=float)
    if pc > 0:
        a = pi_high_arr + pc
        b = pi_low_arr + pc
        valid = np.isfinite(a) & np.isfinite(b) & (b > 0) & (a > 0)
        np.divide(a, b, out=ratio, where=valid)
    else:
        valid = np.isfinite(pi_high_arr) & np.isfinite(pi_low_arr) & (pi_high_arr > 0) & (pi_low_arr > 0)
        np.divide(pi_high_arr, pi_low_arr, out=ratio, where=valid)

    logratio = np.full(len(windows), np.nan, dtype=float)
    m = (ratio > 0) & np.isfinite(ratio)
    logratio[m] = np.log10(ratio[m])
    windows["log10_pi_ratio"] = logratio


def rank_top_windows(windows, topn: int) -> Dict[str, "object"]:
    import numpy as np
    import pandas as pd

    df = windows.copy()

    # precompute z-scores for composite signals
    df["z_fst"] = zscore(df["fst"])
    df["z_pi_high_low"] = zscore(-df["pi_high"])  # higher => lower pi_high
    df["z_pi_low_low"] = zscore(-df["pi_low"])  # higher => lower pi_low
    df["z_taj_high_low"] = zscore(-df["tajima_high"])  # higher => more negative
    df["z_taj_low_low"] = zscore(-df["tajima_low"])  # higher => more negative
    df["z_abs_logratio"] = zscore(np.abs(df["log10_pi_ratio"]))
    df["z_abs_tajdiff"] = zscore(np.abs(df["tajima_diff"]))

    df["score_high_fst"] = df["fst"]
    df["score_pi_ratio_abs"] = np.abs(df["log10_pi_ratio"])
    df["score_tajima_diff_abs"] = np.abs(df["tajima_diff"])

    df["score_sweep_high"] = df[["z_fst", "z_pi_high_low", "z_taj_high_low"]].sum(axis=1, skipna=False)
    df["score_sweep_low"] = df[["z_fst", "z_pi_low_low", "z_taj_low_low"]].sum(axis=1, skipna=False)

    df["score_combined_abs"] = df[["z_fst", "z_abs_logratio", "z_abs_tajdiff"]].sum(axis=1, skipna=False)

    # metric distributions for empirical p-values
    dist_fst = df["fst"].dropna()
    dist_pi_high = df["pi_high"].dropna()
    dist_pi_low = df["pi_low"].dropna()
    dist_logratio = df["log10_pi_ratio"].dropna()
    dist_taj_high = df["tajima_high"].dropna()
    dist_taj_low = df["tajima_low"].dropna()
    dist_tajdiff = df["tajima_diff"].dropna()

    rankings: Dict[str, Tuple[str, bool]] = {
        "high_fst": ("score_high_fst", True),
        "pi_ratio_abs": ("score_pi_ratio_abs", True),
        "tajima_diff_abs": ("score_tajima_diff_abs", True),
        "sweep_high": ("score_sweep_high", True),
        "sweep_low": ("score_sweep_low", True),
        "combined_abs": ("score_combined_abs", True),
    }

    out: Dict[str, object] = {}
    out["rankings"] = rankings
    out["df_scored"] = df

    top_tables: Dict[str, pd.DataFrame] = {}
    for name, (col, desc) in rankings.items():
        t = df.dropna(subset=[col]).sort_values(col, ascending=not desc).head(topn).copy()
        t.insert(0, "selection_type", name)
        t.insert(1, "score", t[col])
        # empirical p on score (high-tail)
        scores = df[col].dropna()
        t["p_empirical_high_tail"] = t[col].apply(lambda x: empirical_p_high(scores, float(x)))

        # empirical p-values for underlying metrics (where available)
        t["p_fst_high_tail"] = t["fst"].apply(lambda x: empirical_p_high(dist_fst, float(x)) if pd.notna(x) else float("nan"))
        t["p_pi_high_low_tail"] = t["pi_high"].apply(lambda x: empirical_p_low(dist_pi_high, float(x)) if pd.notna(x) else float("nan"))
        t["p_pi_low_low_tail"] = t["pi_low"].apply(lambda x: empirical_p_low(dist_pi_low, float(x)) if pd.notna(x) else float("nan"))
        t["p_log10_pi_ratio_two_sided"] = t["log10_pi_ratio"].apply(
            lambda x: empirical_p_two_sided(dist_logratio, float(x)) if pd.notna(x) else float("nan")
        )
        t["p_tajima_high_low_tail"] = t["tajima_high"].apply(lambda x: empirical_p_low(dist_taj_high, float(x)) if pd.notna(x) else float("nan"))
        t["p_tajima_low_low_tail"] = t["tajima_low"].apply(lambda x: empirical_p_low(dist_taj_low, float(x)) if pd.notna(x) else float("nan"))
        t["p_tajima_diff_two_sided"] = t["tajima_diff"].apply(
            lambda x: empirical_p_two_sided(dist_tajdiff, float(x)) if pd.notna(x) else float("nan")
        )
        top_tables[name] = t[
            [
                "selection_type",
                "chrom",
                "start",
                "end",
                "locus",
                "score",
                "p_empirical_high_tail",
                "p_fst_high_tail",
                "p_pi_high_low_tail",
                "p_pi_low_low_tail",
                "p_log10_pi_ratio_two_sided",
                "p_tajima_high_low_tail",
                "p_tajima_low_low_tail",
                "p_tajima_diff_two_sided",
                "fst",
                "pi_high",
                "pi_low",
                "log10_pi_ratio",
                "tajima_high",
                "tajima_low",
                "tajima_diff",
            ]
        ]

    out["top_tables"] = top_tables
    return out


def write_outputs(
    windows,
    scored,
    top_tables,
    outdir: Path,
    pop_high: str,
    pop_low: str,
    pi_pseudocount: float,
):
    import pandas as pd

    outdir.mkdir(parents=True, exist_ok=True)

    # genome-wide summaries (for plotting)
    summary = {
        "pop_high": pop_high,
        "pop_low": pop_low,
        "pi_pseudocount_for_ratio": float(pi_pseudocount or 0.0),
        "metrics": {
            "fst": summarize(windows["fst"]).__dict__,
            "fst_hudson": summarize(windows["fst_hudson"]).__dict__ if "fst_hudson" in windows.columns else summarize(windows["fst"]).__dict__,
            "fst_nei_raw": summarize(windows["fst_nei_raw"]).__dict__ if "fst_nei_raw" in windows.columns else summarize(windows["fst"]).__dict__,
            "dxy": summarize(windows["dxy"]).__dict__,
            "pi_high": summarize(windows["pi_high"]).__dict__,
            "pi_low": summarize(windows["pi_low"]).__dict__,
            "log10_pi_ratio": summarize(windows["log10_pi_ratio"]).__dict__,
            "tajima_high": summarize(windows["tajima_high"]).__dict__,
            "tajima_low": summarize(windows["tajima_low"]).__dict__,
            "tajima_diff": summarize(windows["tajima_diff"]).__dict__,
        },
    }

    # merged windows table
    windows_out = outdir / f"windows_{pop_high}_vs_{pop_low}.tsv.gz"
    windows.to_csv(windows_out, sep="\t", index=False, compression="gzip")

    scored_out = outdir / f"windows_{pop_high}_vs_{pop_low}.scored.tsv.gz"
    scored.to_csv(scored_out, sep="\t", index=False, compression="gzip")

    with open(outdir / "genomewide_summary.json", "w") as fh:
        json.dump(summary, fh, indent=2)

    # per-type tops
    all_tops = []
    for name, table in top_tables.items():
        p = outdir / f"top_windows.{name}.tsv"
        table.to_csv(p, sep="\t", index=False)
        all_tops.append(table)

    if all_tops:
        pd.concat(all_tops, ignore_index=True).to_csv(outdir / "top_windows.ALL.tsv", sep="\t", index=False)


def main(argv: Optional[List[str]] = None) -> int:
    p = argparse.ArgumentParser(description="Rank selection-signal windows between two populations using piawka outputs")
    p.add_argument("--piawka-main", type=Path, required=True, help="piawka_main.tsv(.gz) with Fst_HUD and pi")
    p.add_argument("--piawka-tajima", type=Path, required=True, help="piawka_tajima.tsv(.gz) with TajD")
    p.add_argument("--high-pop", required=True, help="High-elevation population label (e.g. DAS)")
    p.add_argument("--low-pop", required=True, help="Low-elevation population label (e.g. DAG)")
    p.add_argument("--min-used-sites", type=int, default=0, help="Min nUsed sites per row to keep (default: 0)")
    p.add_argument(
        "--fst-source",
        choices=["hudson", "nei"],
        default="hudson",
        help="Source for primary 'fst' column: hudson (Fst_HUD from piawka_main) or nei (derived from pi & Dxy in piawka_tajima)",
    )
    p.add_argument(
        "--pi-pseudocount",
        type=float,
        default=0.0,
        help="Add this pseudocount to BOTH pi_high and pi_low when computing log10(piH/piL) (default: 0, i.e. require pi>0)",
    )
    p.add_argument(
        "--pi-pseudocount-auto",
        action="store_true",
        help="Auto-choose a pseudocount as 0.5 * (1st percentile of positive pi across both pops); overrides --pi-pseudocount",
    )
    p.add_argument("--topn", type=int, default=10, help="Top N windows per selection type (default: 10)")
    p.add_argument("--outdir", type=Path, required=True, help="Output directory")

    args = p.parse_args(argv)

    windows = build_windows(
        piawka_main=args.piawka_main,
        piawka_tajima=args.piawka_tajima,
        pop_high=args.high_pop,
        pop_low=args.low_pop,
        min_used_sites=args.min_used_sites,
        pi_pseudocount=float(args.pi_pseudocount or 0.0),
        fst_source=str(args.fst_source),
    )

    # determine pseudocount used for ratio
    pi_pseudocount = float(args.pi_pseudocount or 0.0)
    if args.pi_pseudocount_auto:
        # Auto-choose pseudocount from already-loaded data: 0.5 * 1st percentile of positive pi
        import numpy as np
        pooled = windows[["pi_high", "pi_low"]].to_numpy(dtype=float).ravel()
        pooled = pooled[np.isfinite(pooled) & (pooled > 0)]
        if pooled.size == 0:
            pi_pseudocount = 0.0
        else:
            pi_pseudocount = 0.5 * float(np.quantile(pooled, 0.01))
        recompute_log10_pi_ratio_inplace(windows, pi_pseudocount=pi_pseudocount)

    ranked = rank_top_windows(windows, topn=args.topn)
    scored = ranked["df_scored"]
    top_tables = ranked["top_tables"]

    write_outputs(
        windows=windows,
        scored=scored,
        top_tables=top_tables,
        outdir=args.outdir,
        pop_high=args.high_pop,
        pop_low=args.low_pop,
        pi_pseudocount=pi_pseudocount,
    )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
