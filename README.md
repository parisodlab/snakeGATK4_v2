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
python main_pipeline.py
```

At present, [main_pipeline.py](main_pipeline.py) launches mapping and the GATK calling stage directly. Other workflow stages are driven from the Snakemake rules and downstream targets.

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

### Downstream Analysis Blocks

These sections mostly point to downstream resources and analysis preferences:

- `final_vcf`
- `final_vcf_full`
- `Structure`
- `PCA`
- `popgen_investigate`

Change these when:

- you produce a different final VCF target
- you want a different Structure run mode or `K` range
- you want different tree, outgroup, or Twisst settings

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

- The repository currently still contains some older rule naming in [main_pipeline.py](main_pipeline.py), while the workflow files are organized as numbered rule files under `workflow/`.
- If you want, the next useful cleanup would be to align [main_pipeline.py](main_pipeline.py) with the actual numbered rule filenames and document the expected execution order more explicitly.
