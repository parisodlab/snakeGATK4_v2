#!/usr/bin/env bash
# Clone a ScanTools repository and produce a SLURM-patched copy for cluster use.
# Usage: ./scripts/clone_and_patch_scantools.sh <git_url> <clone_dir> [partition] [account]
# Example: ./scripts/clone_and_patch_scantools.sh https://github.com/pmonnahan/ScanTools ./scripts/ScanTools pibu_el8 rchoudhury

set -euo pipefail

GIT_URL=${1:-}
CLONE_DIR=${2:-}
PARTITION=${3:-}
ACCOUNT=${4:-}

if [[ -z "$GIT_URL" || -z "$CLONE_DIR" ]]; then
  echo "Usage: $0 <git_url> <clone_dir> [partition] [account]"
  exit 1
fi

if ! command -v git >/dev/null 2>&1; then
  echo "git not found in PATH. Please install git or run this script on a system with git." >&2
  exit 1
fi

echo "Cloning $GIT_URL -> $CLONE_DIR"
git clone "$GIT_URL" "$CLONE_DIR"

PATCHED_DIR="${CLONE_DIR}_slurm"
if [[ -d "$PATCHED_DIR" ]]; then
  echo "Patched dir $PATCHED_DIR already exists. Removing and recreating." 
  rm -rf "$PATCHED_DIR"
fi

cp -a "$CLONE_DIR" "$PATCHED_DIR"
echo "Created patched copy at $PATCHED_DIR"

echo "Applying heuristic PBS -> SLURM conversions in $PATCHED_DIR (this is best-effort; inspect results)"

# find files likely to contain cluster submission commands
find "$PATCHED_DIR" -type f \( -iname "*.sh" -o -iname "*.py" -o -iname "*.pbs" -o -iname "*.qsub" -o -iname "*.pl" \) -print0 | while IFS= read -r -d '' file; do
  # create a backup
  cp "$file" "$file.bak" || true
  # replace PBS directive header with SBATCH
  sed -e 's/^#\s*PBS/#SBATCH/g' \
      -e 's/\bqsub\b/sbatch/g' \
      -e 's/walltime=/--time=/g' \
      -e 's/mem=/--mem=/g' \
      -e 's/\bmodule add\b/module load/g' \
      "$file.bak" > "$file"

  # add a sample SBATCH header to any file that had a #PBS line originally
  if grep -q "#PBS" "$file.bak" 2>/dev/null; then
    echo "# NOTE: original PBS header lines were converted to SBATCH-like flags. Review and add cluster-specific parameters (partition, account)." >> "$file"
  fi
done

# Create a README describing manual checks
README_PATCH="$PATCHED_DIR/README_SLURM_PATCH.txt"
cat > "$README_PATCH" <<'EOF'
This directory was created by scripts/clone_and_patch_scantools.sh.

What was done (heuristic):
- The repository was copied to a new folder with suffix _slurm.
- Files with extensions .sh, .py, .pbs, .qsub, .pl were backed up with .bak and processed with simple text replacements:
  - '#PBS' -> '#SBATCH'
  - 'qsub' -> 'sbatch'
  - 'walltime=' -> '--time='
  - 'mem=' -> '--mem='
  - 'module add' -> 'module load'

What you must check manually:
- Inspect job headers in scripts under the patched folder and adapt them to your SLURM environment.
  Add lines such as: #SBATCH --partition=YOUR_PARTITION and #SBATCH --account=YOUR_ACCOUNT
- Verify module names (e.g., python module) and change 'module load' targets to modules available on your cluster.
- Check any hard-coded paths (e.g., to GATK3.7, fastsimcoal2) and update them to the cluster installation paths or module commands.
- Test small example runs interactively before submitting many batch jobs.

Suggested SBATCH header example to insert at top of scripts:
  #!/bin/bash
  #SBATCH --job-name=ScanTools_job
  #SBATCH --partition=YOUR_PARTITION
  #SBATCH --account=YOUR_ACCOUNT
  #SBATCH --time=02:00:00
  #SBATCH --mem=16000
  #SBATCH --cpus-per-task=12

EOF

echo "Created $README_PATCH with instructions."

echo "You should now inspect $PATCHED_DIR and adjust partition/account/module names and any site-specific paths."
echo "Example: add module load commands for GATK3.7 and fastsimcoal2 or ensure they are in PATH."

echo "Done. Patched copy at: $PATCHED_DIR"

echo "To run ScanTools with this patched copy, call the Python wrapper script with the patched path as 6th arg, e.g.:"
echo "python3 scripts/scantools_run.py /path/to/af_filtered.vcf.gz sample_map.txt scaffolds.list results/scantools 50000 $PATCHED_DIR slurm"

exit 0
