#!/usr/bin/env bash
set -euo pipefail

ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
REPO=$(cd "${ROOT}/../.." && pwd)
CALLER=${1:-both}
CORES=${CORES:-2}
CONDA_ENV=${CONDA_ENV:-snake_env}
FIXTURE_CONFIG=${FIXTURE_CONFIG:-tests/fake_smoke/test_data_config.yaml}

cd "${REPO}"
conda run -n "${CONDA_ENV}" python tests/fake_smoke/generate_fixture.py --config "${FIXTURE_CONFIG}"

REF=tests/fake_smoke/data/reference/synthetic.fa
conda run -n "${CONDA_ENV}" samtools faidx "${REF}"
conda run -n "${CONDA_ENV}" gatk CreateSequenceDictionary -R "${REF}" -O tests/fake_smoke/data/reference/synthetic.dict
conda run -n "${CONDA_ENV}" bwa index "${REF}"

run_pipeline() {
  local caller_name=$1
  conda run -n "${CONDA_ENV}" python main_pipeline.py \
    --config tests/fake_smoke/config.fake.yaml \
    --caller "${caller_name}" \
    --stages mapping,calling,vcf_annotation,filter_vcf \
    -- --cores "${CORES}" --use-conda --conda-frontend conda
}

case "${CALLER}" in
  gatk)
    run_pipeline gatk
    ;;
  freebayes)
    run_pipeline freebayes
    ;;
  both)
    run_pipeline gatk
    run_pipeline freebayes
    ;;
  *)
    echo "Usage: $0 [gatk|freebayes|both]" >&2
    exit 2
    ;;
esac