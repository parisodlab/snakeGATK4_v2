# snakeSNPcalling

This repository runs a resequencing variant-calling workflow for mixed ploidy datasets, with support for:

- optional read QC and trimming
- mapping and lane-level BAM merging
- WGS and ddRAD handling
- GATK HaplotypeCaller gVCF generation
- joint genotyping with GenomicsDBImport + GenotypeGVCFs
- VCF filtering, depth masking, SNP/indel outputs, and downstream population-genetic analysis

The main user-controlled file is [config_main.yaml](config_main.yaml).

## Running The Pipeline

The pipeline entry point is [main_pipeline.py](main_pipeline.py).

```bash
python main_pipeline.py --stages mapping,calling_gatk,vcf_annotation,filter_vcf --profile slurm
python main_pipeline.py --caller freebayes --profile slurm
```

You can also run any stage directly with Snakemake (recommended for debugging):

```bash
snakemake -s workflow/03-mapping.rules --configfile config_main.yaml --profile slurm
snakemake -s workflow/04-calling_gatk4.rules --configfile config_main.yaml --profile slurm
snakemake -s workflow/04-calling_freebayes.rules --configfile config_main.yaml --profile slurm
snakemake -s workflow/05-filter_vcf_gatk4.rules --configfile config_main.yaml --profile slurm
snakemake -s workflow/05-filter_vcf_freebayes.rules --configfile config_main.yaml --profile slurm
snakemake -s workflow/00-vcf_annotation.rules --configfile config_main.yaml --profile slurm
```

### Stage Names

[main_pipeline.py](main_pipeline.py) runs one or more named stages. Each stage corresponds to a single rules file under `workflow/`:

- `qc` → `workflow/01-quality_control.rules`
- `trim` → `workflow/02-trim.rules`
- `mapping` → `workflow/03-mapping.rules`
- `calling_gatk` → `workflow/04-calling_gatk4.rules`
- `calling_gatk4` → `workflow/04-calling_gatk4.rules`
- `calling_freebayes` → `workflow/04-calling_freebayes.rules`
- `filter_vcf_gatk` → `workflow/05-filter_vcf_gatk4.rules`
- `filter_vcf_gatk4` → `workflow/05-filter_vcf_gatk4.rules`
- `filter_vcf_freebayes` → `workflow/05-filter_vcf_freebayes.rules`
- `vcf_annotation` → `workflow/00-vcf_annotation.rules`
- `structure` → `workflow/07-structure_analysis.rules`
- `pixy` → `workflow/08-pixy.rules`
- `piawka` → `workflow/09-piawka_mixed_ploidy.rules`
- `publication_ready` → `workflow/10-publication_ready.rules`
- `popgen_investigate` → `workflow/11-popgen_inversitgate.rules`

The launcher also accepts a caller alias:

- `--caller gatk` uses default stages `mapping,calling_gatk,vcf_annotation,filter_vcf`
- `--caller freebayes` uses default stages `mapping,calling_freebayes,vcf_annotation,filter_vcf`
- if you specify `--stages calling`, the launcher resolves that alias to the selected caller-specific stage
- if you specify `--stages filter_vcf`, the launcher resolves that alias to the selected caller-specific filtering stage

For compatibility, the older stage names `calling_gatk4` and `filter_vcf_gatk4` still work, but `calling_gatk` and `filter_vcf_gatk` remain valid launcher-facing names.
The numbered FreeBayes entrypoint is a thin wrapper around the canonical implementation file. Annotation is owned by `workflow/00-vcf_annotation.rules`.

Default caller is read from `variant_caller` in `config_main.yaml` and falls back to `gatk`.

### Dry Runs

To build the DAG and validate configuration without executing jobs:

```bash
python main_pipeline.py --dry-run --stages mapping,calling_gatk -- --cores 1
python main_pipeline.py --dry-run --caller freebayes -- --cores 1
```

## Required Inputs

You need the following before running:

- raw FASTQ files
- a reference FASTA
- a FASTA index, usually created with `samtools faidx ref.fasta`
- a BWA index, usually created with `bwa index ref.fasta`
- a sequence dictionary, usually created with `gatk CreateSequenceDictionary -R ref.fasta -O ref.dict`
- a metadata table referenced by `METAFILE`

## Metadata File

The file referenced by `METAFILE` should contain one row per lane or one row per sample, depending on your data organization.

Common columns used by this workflow:

- `sample`: biological sample ID
- `prefix`: legacy sample identifier; still supported
- `sample_lane`: preferred unique lane-level identifier if already available
- `sample_barcode`: optional, used to build `sample_lane` if present
- `lane`: sequencing lane
- `pop` or `population`: population label
- `ploidy`: numeric ploidy, typically `2` or `4`
- `fq1`: read 1 FASTQ path
- `fq2`: read 2 FASTQ path

The calling workflow can reconstruct `sample_lane` automatically from:

- `sample + lane + sample_barcode`
- `sample + lane`
- or `prefix`

## Config Guide

This section explains the main config keys and when to change them.

### Core Project Settings

- `PROJECT`
  - Used to create the output directory under `FINALOUTPUT`.
  - Change this when starting a new analysis or when you want a separate result tree.

- `FINALOUTPUT`
  - Root output directory.
  - Change this when moving results to a different filesystem or scratch location.

- `INPUT`
  - Input root directory for raw reads.
  - Change this when your FASTQ files live elsewhere.

- `GENOME`
  - Reference FASTA.
  - Change this whenever you change the reference assembly.
  - If this changes, all mapping and calling results should be considered reference-specific and usually need recomputation.

- `ANNOTATION`
  - GFF/GTF-style annotation used by downstream analyses.
  - Change this when switching species, assembly version, or gene annotation release.

- `METAFILE`
  - Sample sheet.
  - Change this whenever sample composition, population assignments, or ploidy assignments differ.

### Workflow Toggles

- `QC`
  - `yes` or `no`.
  - Set to `yes` only when you want to run the QC stage explicitly.
  - In the current launcher, QC is not automatically run unless you wire it back into [main_pipeline.py](main_pipeline.py).

- `TRIMMING`
  - `yes` or `no`.
  - Set to `yes` when reads need adapter/quality trimming.
  - Keep `no` if you already have cleaned reads and want to avoid duplicate work.

- `SEQTYPE`
  - `wgs` or `ddrad`.
  - `wgs`: MarkDuplicates is used and downstream calling uses BAMs from `dedup/`.
  - `ddrad`: duplicate marking is skipped for downstream calling and BAMs from `merged_samples/` are used.
  - Change this whenever the sequencing protocol changes.
  - This is one of the most important settings because it changes BAM selection in mapping and calling.

### Reference Scope

- `SCAFFOLDS`
  - Either `ALL` or a text file containing scaffold names.
  - Use `ALL` for full-genome analyses.
  - Use a scaffold list for testing, debugging, quick pilot runs, or targeted analyses.
  - Restricting scaffolds is one of the easiest ways to reduce runtime during development.

- `tmpdir`
  - Directory for temporary files.
  - Change this to a fast local scratch filesystem if available.
  - This matters especially for GATK and large joint-genotyping runs.

### Scatter Strategy

The calling workflow now supports separate scatter strategies for:

- HaplotypeCaller
- GenomicsDBImport + GenotypeGVCFs

This is useful because the best scatter scheme is often different for per-sample calling and joint genotyping.

#### Fallback Keys

- `scatter_mode`
- `scatter_count`
- `chromosome_groups`

These act as defaults if the stage-specific keys below are not set.

#### HaplotypeCaller Scatter Keys

- `haplotypecaller_scatter_mode`
  - `chromosome` or `interval`
  - Use `chromosome` when you want cleaner per-contig gVCF partitioning, especially on chromosome-scale assemblies.
  - Use `interval` when the chromosomes are extremely uneven in size and you want more balanced job lengths.

- `haplotypecaller_scatter_count`
  - Used only when `haplotypecaller_scatter_mode: interval`.
  - Increase this to make more, smaller jobs.
  - Decrease this to reduce DAG size and scheduling overhead.

- `haplotypecaller_chromosome_groups`
  - Used only when `haplotypecaller_scatter_mode: chromosome`.
  - Allowed values:
    - `all`
    - a pattern like `Bv1-9+rest`
  - Use `all` to treat every scaffold as its own unit.
  - Use a range like `Bv1-9+rest` when you want main chromosomes separate and all minor contigs bundled into one `REST` unit.
  - This is usually preferable on assemblies with many small contigs.

#### Joint Genotyping Scatter Keys

- `jointgenotyping_scatter_mode`
  - `chromosome` or `interval`
  - Use `interval` when joint genotyping becomes unbalanced or memory-heavy on whole chromosomes.
  - Use `chromosome` when the assembly is compact and chromosome-level splitting is sufficient.

- `jointgenotyping_scatter_count`
  - Used only when `jointgenotyping_scatter_mode: interval`.
  - Increase when GenomicsDBImport or GenotypeGVCFs jobs are too large.
  - Decrease if the DAG becomes too fragmented.

- `jointgenotyping_chromosome_groups`
  - Same semantics as the HaplotypeCaller version, but applied only to joint genotyping.

- `intervals_dir`
  - Base directory for interval scatter files.
  - The workflow now creates stage-specific subdirectories beneath this path.

#### Recommended Patterns

Examples of good combinations:

- Small test run:
  - `SCAFFOLDS` as a small list
  - both scatter modes set to `interval`
  - small `scatter_count`

- Chromosome-scale assembly with many minor contigs:
  - `haplotypecaller_scatter_mode: chromosome`
  - `haplotypecaller_chromosome_groups: Bv1-9+rest`
  - `jointgenotyping_scatter_mode: interval`

- Very small reference or tutorial dataset:
  - keep both stages on `interval`
  - moderate `scatter_count`

### HaplotypeCaller Filter Thresholds

The `HaplotypeCaller` block contains two kinds of settings:

- resource settings under `by_scatter_mode`
- variant filter thresholds used downstream in VCF filtration

Change resource settings when jobs are too slow, too memory-hungry, or too fragmented.

Change filter thresholds only when you have a strong reason, such as:

- a different ploidy regime
- lower or higher sequencing depth
- protocol-specific artifacts
- benchmarking against a validated truth or prior project standard

Important fields:

- `qual`
- `strand_odds_ratio`
- `quality_by_depth`
- `fisherstrand`
- `RMSMappingQuality`
- `MappingQualityRankSumTest`
- `ReadPosRankSum`
- `high_missing`
- `high_heterozygosity`

If you are unsure, keep these unchanged and tune depth thresholds first.

### GenomicsDBImport And GenotypeGVCFs Resources

Under:

- `GenomicsDBImport.by_scatter_mode`
- `GenotypeGVCFs.by_scatter_mode`

the main knobs are:

- `threads`
- `mem_mb`
- Java heap settings
- `batch_size`
- `interval_padding`
- `merge_input_intervals`
- `only_output_calls_starting_in_intervals`

Change these when:

- joint genotyping crashes due to memory
- interval mode produces too much overhead
- chromosome mode produces very imbalanced jobs
- temporary disk usage becomes too large

Practical guidance:

- reduce `mem_mb` and `java_heap_mb` only if you are constrained by scheduler limits
- increase `scatter_count` before aggressively increasing memory
- use `merge_input_intervals: false` in chromosome mode when you want more explicit region handling

### Depth Filtering

The `GATK` block controls genotype depth masking by ploidy.

- `GATK.min_depth.di`
- `GATK.max_depth.di`
- `GATK.min_depth.tetra`
- `GATK.max_depth.tetra`

Change these when:

- your diploids and tetraploids have different expected depth ranges
- coverage is much lower or higher than in the current project
- you see excessive genotype masking or too many suspicious high-depth calls

`GATK.INDEL.min_support` controls the minimum per-sample DP used in indel masking.

### FreeBayes Block

The `FreeBayes` section is only relevant if you run the FreeBayes workflow.

Do not change it for a standard GATK-only analysis.

Change these settings when:

- you intentionally switch to FreeBayes calling
- you work on problematic loci with many alternate alleles
- you need to constrain complexity or runtime in dense regions

Important fields:

- `variant_caller`
  - Used by `main_pipeline.py` to choose the default calling path.
  - Set to `freebayes` if you want `python main_pipeline.py` to default to FreeBayes instead of GATK.

- `FreeBayes.ploidy`
  - Used as the default/fallback ploidy for samples missing a metadata ploidy.
  - For mixed-ploidy datasets, the workflow generates a sample-level CNV map and passes it to FreeBayes with `--cnv-map`.

- `freebayes_params`
  - Optional extra CLI flags/arguments passed directly to FreeBayes.
  - This uses the same dictionary style as the standalone germline FreeBayes pipeline.
  - Boolean `true` values are emitted as flag-only parameters.

FreeBayes outputs are written under `{FINALOUTPUT}/{PROJECT}/freebayes/final_vcf/`.
Downstream rules now follow `variant_caller` automatically when `final_vcf`, `final_vcf_full`, and tool-specific `vcf` overrides are left unset.

### Downstream Analysis Blocks

These sections mostly point to downstream resources and analysis preferences:

- `final_vcf_neutral`
- `final_vcf`
- `final_vcf_full`
- `Structure`
- `PCA`
- `popgen_investigate`

Change these when:

- you produce a different final VCF target
- you want a different Structure run mode or `K` range
- you want different tree, outgroup, or Twisst settings

## Downstream Output Locations (Stages 05–11)

All downstream outputs are written beneath:

- `{FINALOUTPUT}/{PROJECT}` (referred to as `final_path` inside the rules)

Key directories:

- **Stage 05 (filter VCF)**
  - GATK outputs: `{final_path}/gatk4/final_vcf/`
  - FreeBayes outputs: `{final_path}/freebayes/final_vcf/`
  - Neutral SNP branch: `af_filtered.vcf.gz` → `4_fold_degenerate_filtered.vcf.gz` → `4_fold_degenerate_filtered.ld_pruned.vcf.gz`
  - Broad popgen-stat branch: `filtered.vcf.gz` by default, or `af_filtered.vcf.gz` if `VariantFiltering.popgen_stats_use_af_filtered: true`
  - Derived metadata (generated as tracked workflow outputs): `{final_path}/metadata/`
    - Shared caller-agnostic outputs: `sample_pop_unique.tsv`, `sample_ploidy_unique.tsv`
    - Additional GATK-specific QC and masking outputs remain under `{final_path}/gatk4/`

- **Stage 07 (STRUCTURE + PCA)**
  - Structure outputs: `{final_path}/structure/`
  - PCA outputs: `{final_path}/pca/`
  - Default neutral VCF follows `variant_caller` unless `final_vcf_neutral` or `final_vcf` is set.

- **Stage 08 (pixy)**
  - Outputs: `{final_path}/pixy/` (per-ploidy subdirectories)
  - Default VCF follows `variant_caller` unless `final_vcf_full` or `final_vcf` is set.
  - By default this uses the broad filtered VCF rather than the neutral-SNP branch.

- **Stage 09 (piawka)**
  - Outputs: `{final_path}/piawka/`
  - Default VCF follows `variant_caller` unless `final_vcf_full` or `final_vcf` is set.
  - By default this uses the broad filtered VCF rather than the neutral-SNP branch.

## AF Filtering Options

The shared filter behavior is configured under `VariantFiltering` in `config_main.yaml`:

- `missingness_max`
  - Maximum allowed `F_MISSING` fraction for the SNP-filtering branch.

- `neutral_af_filter.enabled`
  - Controls whether the neutral-SNP branch applies allele-frequency filtering before fourfold-site extraction.

- `neutral_af_filter.min_af` / `neutral_af_filter.max_af`
  - AF bounds for the neutral-SNP branch.

- `popgen_stats_use_af_filtered`
  - If `false` (recommended), pixy and piawka stay on the broader filtered VCF.
  - If `true`, their default caller-following VCF becomes `af_filtered.vcf.gz`.

- **Stage 10 (publication_ready)**
  - Outputs: `{final_path}/publication/`
  - Inputs expected to already exist:
    - pixy merged site-level stats from stage 08 (e.g. `{final_path}/pixy/{ploidy}/pixy_sitelevel_pi_merged.tsv.gz`)
    - 4-fold sites BED from stage 06 (`{final_path}/degeneracy/degeneracy-fourfold.bed`)

- **Stage 11 (popgen_investigate)**
  - Outputs: `{final_path}/popgen_investigate/`
  - Default VCF follows `variant_caller` unless `popgen_investigate.vcf` or `final_vcf` is set.

## New QC Outputs

The filtering workflow now writes per-sample mean DP summaries from two VCFs:

- `missingness_filtered.mean_dp_per_sample.tsv`
- `biallelic_snp_filtered.mean_dp_per_sample.tsv`
- `mean_dp_per_sample_qc.tsv`

Use these to compare depth behavior before and after SNP restriction.

## Which Settings Should You Usually Change First?

For a new project, the usual order is:

1. `PROJECT`
2. `GENOME`
3. `ANNOTATION`
4. `METAFILE`
5. `INPUT`
6. `SEQTYPE`
7. scaffold/split settings
8. depth thresholds if coverage differs strongly from earlier projects

## Which Settings Should You Usually Leave Alone?

Unless you are debugging or benchmarking, leave these unchanged at first:

- most `HaplotypeCaller` site-filter thresholds
- `GenomicsDBImport` batch and padding settings
- `GenotypeGVCFs` alternate-allele limits
- `FreeBayes` settings when not using FreeBayes

## Example Configurations

### Whole-Genome WGS

```yaml
SEQTYPE: "wgs"
SCAFFOLDS: "ALL"
haplotypecaller_scatter_mode: "chromosome"
haplotypecaller_chromosome_groups: "Bv1-9+rest"
jointgenotyping_scatter_mode: "interval"
jointgenotyping_scatter_count: 100
```

### ddRAD

```yaml
SEQTYPE: "ddrad"
TRIMMING: "yes"
haplotypecaller_scatter_mode: "interval"
jointgenotyping_scatter_mode: "interval"
```

### Small Debug Run

```yaml
SCAFFOLDS: "configs/test_scaffolds.txt"
haplotypecaller_scatter_mode: "interval"
haplotypecaller_scatter_count: 10
jointgenotyping_scatter_mode: "interval"
jointgenotyping_scatter_count: 10
```

## Notes

- The workflow supports mixed-ploidy calling by pulling per-sample ploidy from `METAFILE` (column `ploidy`). Ensure it is present and correct.
- The QC/trim/mapping stages accept metadata in either lane-resolved form (`sample_lane`) or can reconstruct lane IDs from `sample + lane [+ sample_barcode]`.
- The FreeBayes workflow uses the same BAM selection logic as GATK: `dedup/` for `SEQTYPE: wgs` and `merged_samples/` for `SEQTYPE: ddrad`.
- The FreeBayes workflow now supports mixed-ploidy samples via a generated sample-level copy-number map (`--cnv-map`).
