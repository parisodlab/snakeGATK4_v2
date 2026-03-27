# Fake Smoke Test Data Generation

This fixture uses the requested mutation and read-simulation scheme to simulate reads in a real Biscutella genomic region.

## Region Used

The generator currently clips this Biscutella varia region from the configured assembly:

- contig: `Bv3`
- start: `1`
- length: `10000000`

These values are defined in `tests/fake_smoke/test_data_config.yaml`.

## What Gets Simulated

The generator creates two datasets from the same Biscutella region:

- one diploid draft (`2x`)
- one tetraploid draft (`4x`)

Both use the same mutation scheme:

- `0.3%` polymorphic sites
- substitution model:
  - `C <-> G`
  - `A <-> T`

Ploidy-specific haplotype construction:

- diploid: all selected SNPs are simulated as `Aa`
- tetraploid dosage mix:
  - `1/4` `AAAa`
  - `1/2` `AAaa`
  - `1/4` `Aaaa`

For each ploidy, the generator writes:

- `*_2x.fas` or `*_4x.fas`
- `*_2x_snp_pos.txt` or `*_4x_snp_pos.txt`
- `*_2x.vcf` or `*_4x.vcf` generated with `snp-sites`

## Coverage Levels

The current configured coverages are:

- `20x`
- `5x`

These are applied to both diploid and tetraploid drafts.

## Read Simulation

Reads are simulated with `InSilicoSeq` using the NovaSeq model:

```bash
iss generate --n_reads <computed> --draft <draft.fas> --model novaseq --output <prefix>
```

The read count is computed as:

```text
coverage * sequence_length / read_length
```

The generator currently uses `150 bp` reads.

## Required Tools

Install these tools in the environment used to run the generator:

```bash
conda install -c bioconda insilicoseq snp-sites
```

## Generate Test Data

From the repository root:

```bash
conda run -n snake_env python tests/fake_smoke/generate_fixture.py --config tests/fake_smoke/test_data_config.yaml
```

This produces:

- clipped Biscutella reference and clipped GFF in `tests/fake_smoke/data/reference/`
- ploidy-specific draft FASTAs in `tests/fake_smoke/data/reference/`
- ploidy-specific SNP tables in `tests/fake_smoke/data/reference/`
- `snp-sites` VCFs in `tests/fake_smoke/data/reference/`
- simulated FASTQs in `tests/fake_smoke/data/fastq/`
- metadata in `tests/fake_smoke/data/metadata.tsv`
- workflow config in `tests/fake_smoke/config.fake.yaml`

## Run the Smoke Workflow

```bash
bash tests/fake_smoke/run_smoke.sh freebayes
```

The runner uses `tests/fake_smoke/test_data_config.yaml` by default.

To override the generator config:

```bash
FIXTURE_CONFIG=tests/fake_smoke/test_data_config.yaml bash tests/fake_smoke/run_smoke.sh gatk
```

## Notes

- The clipped region keeps the smoke test tied to a real Biscutella genomic context.
- The clipped GFF comes from the real Biscutella annotation, not a synthetic placeholder.
- Generated reads, logs, tmp files, intervals, indexes, and results are gitignored.
