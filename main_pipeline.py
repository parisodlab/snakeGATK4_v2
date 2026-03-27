#!/usr/bin/env python3
"""Convenience launcher for the snakeSNPcalling Snakemake workflow.

This script is optional: you can always run Snakemake directly.

Why this exists:
- provides a single CLI entrypoint for the numbered rule files in workflow/
- standardizes logs and runtime reporting

Examples:
    python main_pipeline.py --stages mapping,calling_gatk
  python main_pipeline.py --stages qc,trim,mapping --profile slurm
    python main_pipeline.py --dry-run --stages mapping,calling_gatk -- --cores 1
"""

from __future__ import annotations

import argparse
import os
import subprocess
import sys
import time
from pathlib import Path

import yaml


STAGES: dict[str, str] = {
    "qc": "workflow/01-quality_control.rules",
    "trim": "workflow/02-trim.rules",
    "mapping": "workflow/03-mapping.rules",
    "calling_gatk": "workflow/04-calling_gatk4.rules",
    "calling_gatk4": "workflow/04-calling_gatk4.rules",
    "calling_freebayes": "workflow/04-calling_frebayes.rules",
    "filter_vcf_gatk": "workflow/05-filter_vcf_gatk4.rules",
    "filter_vcf_gatk4": "workflow/05-filter_vcf_gatk4.rules",
    "filter_vcf_freebayes": "workflow/05-filter_vcf_freebayes.rules",
    "vcf_annotation": "workflow/00-vcf_annotation.rules",
    "structure": "workflow/07-structure_analysis.rules",
    "pixy": "workflow/08-pixy.rules",
    "piawka": "workflow/09-piawka_mixed_ploidy.rules",
    "publication_ready": "workflow/10-publication_ready.rules",
    "popgen_investigate": "workflow/11-popgen_inversitgate.rules",
}

DEFAULT_STAGE_ORDER_BY_CALLER = {
    "gatk": ["mapping", "calling_gatk", "vcf_annotation", "filter_vcf_gatk"],
    "freebayes": ["mapping", "calling_freebayes", "vcf_annotation", "filter_vcf_freebayes"],
}

CALLING_STAGE_BY_CALLER = {
    "gatk": "calling_gatk",
    "freebayes": "calling_freebayes",
}

FILTER_STAGE_BY_CALLER = {
    "gatk": "filter_vcf_gatk",
    "freebayes": "filter_vcf_freebayes",
}


def _format_elapsed(seconds: float) -> str:
    total_seconds = int(round(seconds))
    hours = total_seconds // 3600
    total_seconds %= 3600
    minutes = total_seconds // 60
    total_seconds %= 60
    return f"{hours}:{minutes:02d}:{total_seconds:02d}"


def _tee_run(cmd: list[str], log_path: Path) -> int:
    log_path.parent.mkdir(parents=True, exist_ok=True)
    with log_path.open("wt", encoding="utf-8") as log_fh:
        log_fh.write("$ " + " ".join(cmd) + "\n")
        log_fh.flush()
        proc = subprocess.Popen(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
        )
        assert proc.stdout is not None
        for line in proc.stdout:
            sys.stdout.write(line)
            log_fh.write(line)
        return proc.wait()


def parse_args(argv: list[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Launch snakeSNPcalling Snakemake stages")
    parser.add_argument(
        "--caller",
        choices=["gatk", "freebayes"],
        default=None,
        help="Variant caller to use for default stage selection and the 'calling' stage alias (default: config 'variant_caller' or gatk)",
    )
    parser.add_argument(
        "--config",
        default="config_main.yaml",
        help="Path to the main config YAML (default: config_main.yaml)",
    )
    parser.add_argument(
        "--stages",
        default=None,
        help=f"Comma-separated stage list (known: {', '.join(STAGES)}, plus aliases: calling, filter_vcf)",
    )
    parser.add_argument(
        "--profile",
        default=None,
        help="Snakemake profile (e.g. 'slurm'); if unset, runs without --profile",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Pass -n to Snakemake (build DAG without executing)",
    )
    parser.add_argument(
        "--rerun-triggers",
        default="mtime",
        help="Value for --rerun-triggers (default: mtime)",
    )
    parser.add_argument(
        "snakemake_args",
        nargs=argparse.REMAINDER,
        help="Extra args passed to snakemake after '--'",
    )
    return parser.parse_args(argv)


def main(argv: list[str]) -> int:
    args = parse_args(argv)
    config_path = Path(args.config)
    if not config_path.exists():
        raise SystemExit(f"Config not found: {config_path}")

    with config_path.open("rt", encoding="utf-8") as fh:
        cfg = yaml.safe_load(fh) or {}

    project = str(cfg.get("PROJECT", "(unknown)"))
    log_folder = Path(str(cfg.get("log_folder", "logs")))
    log_folder.mkdir(parents=True, exist_ok=True)

    caller = str(args.caller or cfg.get("variant_caller", "gatk")).strip().lower()
    if caller not in DEFAULT_STAGE_ORDER_BY_CALLER:
        raise SystemExit(f"Unsupported caller '{caller}'. Expected one of: {', '.join(DEFAULT_STAGE_ORDER_BY_CALLER)}")

    requested_stage_str = args.stages or ",".join(DEFAULT_STAGE_ORDER_BY_CALLER[caller])
    requested = [s.strip() for s in str(requested_stage_str).split(",") if s.strip()]
    resolved_requested = []
    for stage in requested:
        if stage == "calling":
            resolved_requested.append(CALLING_STAGE_BY_CALLER[caller])
        elif stage == "filter_vcf":
            resolved_requested.append(FILTER_STAGE_BY_CALLER[caller])
        else:
            resolved_requested.append(stage)
    requested = resolved_requested
    unknown = [s for s in requested if s not in STAGES]
    if unknown:
        raise SystemExit(f"Unknown stage(s): {', '.join(unknown)}")

    runtime_log = log_folder / "log_running_time.txt"
    with runtime_log.open("at", encoding="utf-8") as rt:
        rt.write(f"\nProject name: {project}\n")
        rt.write(f"Start time: {time.ctime()}\n")

    extra = list(args.snakemake_args)
    if extra and extra[0] == "--":
        extra = extra[1:]

    for stage in requested:
        snakefile = Path(STAGES[stage])
        if not snakefile.exists():
            raise SystemExit(f"Snakefile for stage '{stage}' not found: {snakefile}")

        print(f"\n=== Stage: {stage} (project={project}) ===")
        stage_start = time.time()

        cmd = [
            "snakemake",
            "-p",
            "--rerun-triggers",
            str(args.rerun_triggers),
            "-s",
            str(snakefile),
            "--configfile",
            str(config_path),
        ]
        if args.dry_run:
            cmd.append("-n")
        if args.profile:
            cmd.extend(["--profile", str(args.profile)])
        cmd.extend(extra)

        rc = _tee_run(cmd, log_folder / f"log_{stage}.txt")
        stage_end = time.time()

        with runtime_log.open("at", encoding="utf-8") as rt:
            rt.write(f"Time of running {stage}: {_format_elapsed(stage_end - stage_start)}\n")

        if rc != 0:
            print(f"Stage '{stage}' failed with exit code {rc}")
            return rc

    with runtime_log.open("at", encoding="utf-8") as rt:
        rt.write(f"End time: {time.ctime()}\n")

    print("\nAll requested stages finished.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
