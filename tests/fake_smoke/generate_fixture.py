#!/usr/bin/env python3
"""Generate mixed-ploidy smoke-test data from a Biscutella genomic region.

This follows the requested practical scheme, adapted to a real Biscutella region:
- clip one genomic region from the Biscutella reference assembly
- mutate 0.3% of sites using the substitution model C<->G and A<->T
- generate one diploid draft and one tetraploid draft from the same region
- use tetraploid dosages 1/4 AAAa, 1/2 AAaa, 1/4 Aaaa
- record SNP positions and dosages per ploidy
- create VCFs with snp-sites from aligned haplotype FASTAs
- simulate NovaSeq paired reads with InSilicoSeq at configured coverages

The script also writes the metadata, clipped annotation, combined reference, and
config.fake.yaml required by the smoke workflow.
"""
from __future__ import annotations

import argparse
import gzip
import math
import random
import shutil
import subprocess
from pathlib import Path


import yaml
import csv


BASE_SWAP = {"A": "T", "T": "A", "C": "G", "G": "C"}
TETRAPLOID_ALT_COPIES = {"AAAa": 1, "AAaa": 2, "Aaaa": 3}
DEFAULT_CONFIG = Path(__file__).with_name("test_data_config.yaml")
LANE = "L001"


def parse_args() -> argparse.Namespace:
  parser = argparse.ArgumentParser(description="Generate mixed-ploidy smoke-test data")
  parser.add_argument(
    "--config",
    type=Path,
    default=DEFAULT_CONFIG,
    help=f"Path to the generator config YAML (default: {DEFAULT_CONFIG.name})",
  )
  parser.add_argument("--seed", type=int, default=20260327, help="Random seed for reproducible mutations")
  return parser.parse_args()


def load_generator_config(config_path: Path) -> dict:
  with config_path.open("rt", encoding="ascii") as handle:
    data = yaml.safe_load(handle) or {}
  for key in ("reference_fasta", "annotation_gff", "contig", "region_start", "region_length", "coverage_levels", "lanes_per_sample", "variant_panel_tsv"):
    if key not in data:
      raise ValueError(f"Generator config must define '{key}'")
  return data


def parse_variant_panel(panel_path: Path) -> dict:
  """Parse variant_panel.tsv and return dict: {(sample, ploidy): [variant_dicts]}"""
  panel = {}
  with panel_path.open("rt", encoding="utf-8") as f:
    reader = csv.DictReader(f, delimiter="\t")
    for row in reader:
      key = (row["sample"], int(row["ploidy"]))
      panel.setdefault(key, []).append(row)
  return panel


def wrap_fasta(sequence: str, width: int = 80) -> str:
  return "\n".join(sequence[i:i + width] for i in range(0, len(sequence), width))


def read_fasta_region(path: Path, contig: str, start_1based: int, length: int) -> str:
  sequence_parts: list[str] = []
  collecting = False
  with path.open("rt", encoding="ascii", errors="ignore") as handle:
    for line in handle:
      if line.startswith(">"):
        name = line[1:].strip().split()[0]
        collecting = name == contig
        if sequence_parts and not collecting:
          break
        continue
      if collecting:
        sequence_parts.append(line.strip())
  if not sequence_parts:
    raise ValueError(f"Contig {contig} not found in {path}")
  sequence = "".join(sequence_parts).upper()
  start = start_1based - 1
  end = start + length
  region = sequence[start:end]
  if len(region) != length:
    raise ValueError(f"Requested region {contig}:{start_1based}-{start_1based + length - 1} exceeds contig length")
  return region


def clip_gff(gff_path: Path, contig: str, start_1based: int, length: int, out_path: Path) -> int:
  region_end = start_1based + length - 1
  kept = 0
  with gff_path.open("rt", encoding="ascii", errors="ignore") as src, out_path.open("wt", encoding="ascii") as dst:
    dst.write("##gff-version 3\n")
    for line in src:
      if not line or line.startswith("#"):
        continue
      fields = line.rstrip("\n").split("\t")
      if len(fields) != 9 or fields[0] != contig:
        continue
      feature_start = int(fields[3])
      feature_end = int(fields[4])
      if feature_end < start_1based or feature_start > region_end:
        continue
      clipped_start = max(feature_start, start_1based) - start_1based + 1
      clipped_end = min(feature_end, region_end) - start_1based + 1
      fields[0] = "chrSynthetic"
      fields[3] = str(clipped_start)
      fields[4] = str(clipped_end)
      dst.write("\t".join(fields) + "\n")
      kept += 1
  return kept


def reset_output_dirs(paths: list[Path]) -> None:
  for path in paths:
    if path.exists():
      shutil.rmtree(path)
    path.mkdir(parents=True, exist_ok=True)


def clear_reference_artifacts(ref_dir: Path) -> None:
  for pattern in (
    "synthetic.fa",
    "synthetic.fa.*",
    "synthetic.dict",
    "synthetic.gff",
    "*_2x.fas",
    "*_2x_alignment.fas",
    "*_2x.vcf",
    "*_2x.vcf.gz",
    "*_2x.vcf.gz.tbi",
    "*_2x_snp_pos.txt",
    "*_4x.fas",
    "*_4x_alignment.fas",
    "*_4x.vcf",
    "*_4x.vcf.gz",
    "*_4x.vcf.gz.tbi",
    "*_4x_snp_pos.txt",
  ):
    for path in ref_dir.glob(pattern):
      path.unlink(missing_ok=True)


def weighted_choice(rng: random.Random, weighted_values: list[tuple[str, float]]) -> str:
  threshold = rng.random()
  running = 0.0
  for value, weight in weighted_values:
    running += weight
    if threshold <= running:
      return value
  return weighted_values[-1][0]


def choose_mutation_positions(sequence: str, mutation_rate: float, rng: random.Random) -> list[int]:
  mutable = [idx for idx, base in enumerate(sequence) if base in BASE_SWAP]
  count = max(1, round(len(sequence) * mutation_rate))
  if count > len(mutable):
    raise ValueError("Mutation rate selects more positions than mutable bases available")
  return sorted(rng.sample(mutable, count))


def build_diploid_haplotypes(reference: str, mutation_positions: list[int], rng: random.Random) -> tuple[list[tuple[str, str]], list[dict[str, object]]]:
  haplotypes = [list(reference) for _ in range(2)]
  records: list[dict[str, object]] = []
  for position in mutation_positions:
    alt_base = BASE_SWAP[reference[position]]
    alt_haplotype = rng.randrange(2)
    haplotypes[alt_haplotype][position] = alt_base
    records.append(
      {
        "position": position + 1,
        "ref": reference[position],
        "alt": alt_base,
        "dosage": "Aa",
        "alt_haplotypes": f"hap{alt_haplotype + 1}",
      }
    )
  return [(f"hap{idx + 1}", "".join(sequence)) for idx, sequence in enumerate(haplotypes)], records


def build_tetraploid_haplotypes(
  reference: str,
  mutation_positions: list[int],
  dosage_weights: list[tuple[str, float]],
  rng: random.Random,
) -> tuple[list[tuple[str, str]], list[dict[str, object]]]:
  haplotypes = [list(reference) for _ in range(4)]
  records: list[dict[str, object]] = []
  for position in mutation_positions:
    dosage = weighted_choice(rng, dosage_weights)
    alt_copies = TETRAPLOID_ALT_COPIES[dosage]
    alt_base = BASE_SWAP[reference[position]]
    chosen_haplotypes = sorted(rng.sample(range(4), alt_copies))
    for hap_index in chosen_haplotypes:
      haplotypes[hap_index][position] = alt_base
    records.append(
      {
        "position": position + 1,
        "ref": reference[position],
        "alt": alt_base,
        "dosage": dosage,
        "alt_haplotypes": ",".join(f"hap{hap + 1}" for hap in chosen_haplotypes),
      }
    )
  return [(f"hap{idx + 1}", "".join(sequence)) for idx, sequence in enumerate(haplotypes)], records


def write_fasta_records(path: Path, records: list[tuple[str, str]]) -> None:
  with path.open("wt", encoding="ascii") as handle:
    for name, sequence in records:
      handle.write(f">{name}\n{wrap_fasta(sequence)}\n")


def write_snp_positions(path: Path, records: list[dict[str, object]]) -> None:
  with path.open("wt", encoding="ascii") as handle:
    handle.write("position\tref\talt\tdosage\talt_haplotypes\n")
    for record in records:
      handle.write(
        f"{record['position']}\t{record['ref']}\t{record['alt']}\t{record['dosage']}\t{record['alt_haplotypes']}\n"
      )


def run_command(command: list[str], error_hint: str) -> None:
  try:
    subprocess.run(command, check=True)
  except FileNotFoundError as exc:
    raise RuntimeError(f"Required command not found while {error_hint}: {command[0]}") from exc
  except subprocess.CalledProcessError as exc:
    raise RuntimeError(f"Command failed while {error_hint}: {' '.join(command)}") from exc


def create_vcf_with_snp_sites(alignment_fasta: Path, output_vcf: Path) -> None:
  run_command(["snp-sites", "-v", "-o", str(output_vcf), str(alignment_fasta)], "running snp-sites")


def compute_n_reads(sequence_length: int, coverage: int, read_length: int) -> int:
  return max(2, math.ceil((coverage * sequence_length) / read_length))


def gzip_if_needed(path: Path) -> Path:
  if path.suffix == ".gz":
    return path
  gz_path = path.with_suffix(path.suffix + ".gz")
  with path.open("rb") as src, gzip.open(gz_path, "wb") as dst:
    shutil.copyfileobj(src, dst)
  path.unlink()
  return gz_path


def simulate_reads(
  draft_fasta: Path,
  output_prefix: Path,
  coverage: int,
  sequence_length: int,
  read_length: int,
  model: str,
  cpu_count: int,
) -> tuple[Path, Path]:
  n_reads = compute_n_reads(sequence_length, coverage, read_length)
  run_command(
    [
      "iss",
      "generate",
      "--n_reads",
      str(n_reads),
      "--draft",
      str(draft_fasta),
      "--model",
      model,
      "--cpus",
      str(cpu_count),
      "--output",
      str(output_prefix),
    ],
    f"simulating {coverage}x reads from {draft_fasta.name}",
  )
  r1 = output_prefix.with_name(output_prefix.name + "_R1.fastq")
  r2 = output_prefix.with_name(output_prefix.name + "_R2.fastq")
  if not r1.exists() or not r2.exists():
    raise RuntimeError(f"InSilicoSeq did not create expected FASTQs for prefix {output_prefix}")
  return gzip_if_needed(r1), gzip_if_needed(r2)


def build_pipeline_config(
  root: Path,
  ref_path: Path,
  gff_path: Path,
  metadata_path: Path,
  generator_cfg: dict,
) -> dict:
  return {
    "PROJECT": "synthetic_smoke",
    "variant_caller": "gatk",
    "QC": "no",
    "TRIMMING": "no",
    "SEQTYPE": "ddrad",
    "GENOME": str(ref_path),
    "ANNOTATION": str(gff_path),
    "SCAFFOLDS": "ALL",
    "METAFILE": str(metadata_path),
    "INPUT": str(root / "data" / "fastq"),
    "FINALOUTPUT": str(root / "results"),
    "log_folder": str(root / "logs"),
    "tmpdir": str(root / "tmp"),
    "scatter_mode": "chromosome",
    "haplotypecaller_scatter_mode": "chromosome",
    "jointgenotyping_scatter_mode": "chromosome",
    "chromosome_groups": "all",
    "haplotypecaller_chromosome_groups": "all",
    "jointgenotyping_chromosome_groups": "all",
    "intervals_dir": str(root / "intervals"),
    "scatter_count": 1,
    "haplotypecaller_scatter_count": 1,
    "jointgenotyping_scatter_count": 1,
    "HaplotypeCaller": {
      "minConfidenceForVariantCalling": "20",
      "minimum_mapping_quality": 10,
      "by_scatter_mode": {
        "chromosome": {
          "threads": 1,
          "mem_mb": 4000,
          "time": "1-0",
          "java_xms": "1G",
          "java_xmx": "2G",
        }
      },
      "min_depth": 1,
      "qual": 10.0,
      "strand_odds_ratio": 10.0,
      "quality_by_depth": 0.5,
      "fisherstrand": 100.0,
      "RMSMappingQuality": 20.0,
      "MappingQualityRankSumTest": -20.0,
      "ReadPosRankSum": -20.0,
      "high_missing": 1.0,
      "high_heterozygosity": 1.0,
    },
    "GenotypeGVCFs": {
      "minConfidenceForVariantCalling": "10",
      "by_scatter_mode": {
        "chromosome": {
          "threads": 1,
          "mem_mb": 4000,
          "time": "1-0",
          "java_heap_mb": 2000,
          "max_alternate_alleles": 12,
          "only_output_calls_starting_in_intervals": False,
        }
      },
    },
    "GenomicsDBImport": {
      "by_scatter_mode": {
        "chromosome": {
          "threads": 1,
          "mem_mb": 4000,
          "java_xms": "1G",
          "java_xmx": "2G",
          "batch_size": 50,
          "interval_padding": 0,
          "merge_input_intervals": False,
        }
      }
    },
    "GATK": {
      "min_depth": {"tetra": 1, "di": 1},
      "max_depth": {"tetra": 500, "di": 500},
      "INDEL": {"min_support": 1},
    },
    "FreeBayes": {
      "ploidy": 4,
      "chunks": 50000,
      "theta": 0.001,
      "max_complex_gap": -1,
      "use_best_n_alleles": 16,
      "min_base_quality": 0,
      "min_mapping_quality": 10,
    },
    "freebayes_params": {},
    "VariantFiltering": {
      "missingness_max": 1.0,
      "freebayes_site_filter": {
        "min_qual": 1.0,
        "min_saf": 1,
        "min_sar": 1,
        "min_rpr": 1,
        "min_rpl": 1,
      },
      "neutral_af_filter": {"enabled": False, "min_af": 0.0, "max_af": 1.0},
      "popgen_stats_use_af_filtered": False,
    },
    "SyntheticFixture": {
      "reference_fasta": generator_cfg["reference_fasta"],
      "annotation_gff": generator_cfg["annotation_gff"],
      "contig": generator_cfg["contig"],
      "region_start": generator_cfg["region_start"],
      "region_length": generator_cfg["region_length"],
      "coverage_levels": generator_cfg["coverage_levels"],
      "ploidies": [2, 4],
      "mutation_rate": generator_cfg["mutation_rate"],
      "substitution_model": "C<->G, A<->T",
      "diploid_dosage_scheme": {"Aa": 1.0},
      "tetraploid_dosage_scheme": {"AAAa": 0.25, "AAaa": 0.5, "Aaaa": 0.25},
      "iss_model": generator_cfg["iss_model"],
      "read_length": generator_cfg["read_length"],
    },
    "final_vcf_neutral": None,
    "final_vcf": None,
    "final_vcf_full": None,
  }


def main() -> None:

    args = parse_args()
    config_path = args.config.resolve()
    generator_cfg = load_generator_config(config_path)
    root = Path(__file__).resolve().parent

    data_dir = root / "data"
    ref_dir = data_dir / "reference"
    fastq_dir = data_dir / "fastq"
    results_dir = root / "results"
    logs_dir = root / "logs"
    tmp_dir = root / "tmp"
    intervals_dir = root / "intervals"

    ref_dir.mkdir(parents=True, exist_ok=True)
    clear_reference_artifacts(ref_dir)
    reset_output_dirs([fastq_dir, results_dir, logs_dir, tmp_dir, intervals_dir])

    rng = random.Random(args.seed)
    mutation_rate = float(generator_cfg.get("mutation_rate", 0.003))
    read_length = int(generator_cfg.get("read_length", 150))
    iss_model = str(generator_cfg.get("iss_model", "novaseq"))
    iss_cpus = int(generator_cfg.get("iss_cpus", 1))
    tetraploid_dosage_weights = [("AAAa", 0.25), ("AAaa", 0.5), ("Aaaa", 0.25)]
    lanes_per_sample = int(generator_cfg.get("lanes_per_sample", 1))
    variant_panel_path = Path(generator_cfg["variant_panel_tsv"])

    reference_fasta = Path(str(generator_cfg["reference_fasta"]))
    annotation_gff = Path(str(generator_cfg["annotation_gff"]))
    contig = str(generator_cfg["contig"])
    region_start = int(generator_cfg["region_start"])
    region_length = int(generator_cfg["region_length"])
    region_name = str(generator_cfg.get("region_name", f"{contig}_{region_start}_{region_start + region_length - 1}"))

    reference_sequence = read_fasta_region(reference_fasta, contig, region_start, region_length)
    ref_path = ref_dir / "synthetic.fa"
    ref_path.write_text(">chrSynthetic\n" + wrap_fasta(reference_sequence) + "\n", encoding="ascii")

    gff_path = ref_dir / "synthetic.gff"
    clipped_features = clip_gff(annotation_gff, contig, region_start, region_length, gff_path)
    if clipped_features == 0:
        raise ValueError("No annotation features overlapped the selected Biscutella region")

    # Parse variant panel
    panel = parse_variant_panel(variant_panel_path)

    metadata_lines = ["sample\tlane\tpop\tploidy\tfq1\tfq2"]
    for (sample, ploidy), variants in panel.items():
        label = f"{ploidy}x"
        # Build haplotypes for this sample using the variant panel
        hap_count = ploidy
        haplotypes = [list(reference_sequence) for _ in range(hap_count)]
        snp_records = []
        for v in variants:
            pos = int(v["position"]) - 1
            alt = v["alt"]
            dosage = int(v["dosage"])
            # Assign alt alleles to haplotypes (simple: first N)
            for h in range(dosage):
                haplotypes[h][pos] = alt
            snp_records.append({
                "position": pos + 1,
                "ref": reference_sequence[pos],
                "alt": alt,
                "dosage": dosage,
                "alt_haplotypes": ",".join(f"hap{h+1}" for h in range(dosage))
            })
        haplotype_records = [(f"{sample}_{label}_hap{h+1}", "".join(seq)) for h, seq in enumerate(haplotypes)]
        draft_path = ref_dir / f"{sample}_{label}.fas"
        alignment_path = ref_dir / f"{sample}_{label}_alignment.fas"
        snp_pos_path = ref_dir / f"{sample}_{label}_snp_pos.txt"
        vcf_path = ref_dir / f"{sample}_{label}.vcf"

        write_fasta_records(draft_path, haplotype_records)
        write_fasta_records(
            alignment_path,
            [("chrSynthetic", reference_sequence)] + haplotype_records,
        )
        write_snp_positions(snp_pos_path, snp_records)
        create_vcf_with_snp_sites(alignment_path, vcf_path)

        for coverage in generator_cfg["coverage_levels"]:
            for lane in range(1, lanes_per_sample + 1):
                lane_id = f"L{lane:03d}"
                fq_sample = f"{sample}_{label}_cov{coverage}x"
                output_prefix = fastq_dir / f"{fq_sample}_{lane_id}"
                fq1, fq2 = simulate_reads(
                    draft_fasta=draft_path,
                    output_prefix=output_prefix,
                    coverage=int(coverage),
                    sequence_length=len(reference_sequence),
                    read_length=read_length,
                    model=iss_model,
                    cpu_count=iss_cpus,
                )
                metadata_lines.append(f"{fq_sample}\t{lane_id}\t{label}\t{ploidy}\t{fq1}\t{fq2}")

    metadata_path = data_dir / "metadata.tsv"
    metadata_path.write_text("\n".join(metadata_lines) + "\n", encoding="ascii")

    config_fake_path = root / "config.fake.yaml"
    config_fake_path.write_text(
        yaml.safe_dump(build_pipeline_config(root, ref_path, gff_path, metadata_path, generator_cfg), sort_keys=False),
        encoding="ascii",
    )

    print(f"Wrote Biscutella region fixture under {root}")
    print(f"Region: {contig}:{region_start}-{region_start + region_length - 1}")
    print(f"Reference: {ref_path}")
    print(f"Annotation: {gff_path}")
    print(f"Metadata: {metadata_path}")
    print(f"Config: {config_fake_path}")


if __name__ == "__main__":
  main()
