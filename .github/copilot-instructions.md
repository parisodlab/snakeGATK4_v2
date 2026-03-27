---
description: Bioinformatics & HPC workflow coding and project guidelines for Copilot. Load for any bioinformatics, genomics, or computational biology project in this workspace.
# applyTo: 'Bioinformatics, genomics, computational biology, Snakemake, SLURM, Python, R, shell scripts, HPC workflows'
---

## Project Context & Coding Guidelines

- Use for bioinformatics, genomics, and computational biology projects.
- Tools: Snakemake, Python, R, shell scripts, SLURM.
- Each project: data/, scripts/, results/, config/ folders.
- Prioritize reproducibility, modularity, and clear documentation.

### Coding Preferences
- Write readable, well-commented code with descriptive, domain-relevant names (e.g., sample_id, fasta_path).
- Scripts must include a docstring/header and usage example.
- Prefer relative paths; avoid hardcoded directories.

### Workflow & Pipeline Guidance
- Snakemake: Use config files, modular rules, clear input/output paths.
- R: Prefer tidyverse style, reproducible research practices.
- Shell: Ensure portability, error handling, resource checks.
- Output results to organized results/ or logs/ folders.

### HPC Execution Policy
- Never run compute-heavy jobs on the login node.
- Always use a persistent `screen` session before requesting a compute node.
- Before running any long job, get into the habit of checking the current node with `hostname` and ensuring you are on a compute node, not the login node. Login node is `login8.hpc.binf.unibe.ch`. Go to the compute node using `srun` or `sbatch` and then run your job from there. for example:
  shopt -s expand_aliases && source "$HOME/.bash_profile" && srun -p pibu_el8 -c 1 --mem=2G --time=00:15:00 --pty /bin/bash
- Estimate walltime, memory, and CPU based on task/input size and benchmarks; document assumptions.
- Each task must run in its own isolated session or scheduler allocation.

### Data Handling
- For large files (FASTA, VCF, BAM): use efficient file handling/streaming.
- Organize outputs/logs by project and analysis step.

### Collaboration
- Write code that is easy to share and rerun by others.
- Update README files and document new scripts/workflows.
- Use environment files (conda, requirements.txt) for reproducibility.

- Write code that is easy to share and rerun by others.
- Update README files and document any new scripts or workflows.
- Use environment files (e.g., conda, requirements.txt) for reproducibility.

---

## How to Use These Instructions

- Place these guidelines in your copilot-instructions.md file.
- Review them before starting new scripts or pipelines.
- Update as your workflow or tools evolve.
- When using Copilot or Copilot Chat, refer to these instructions to ensure code suggestions match your standards.