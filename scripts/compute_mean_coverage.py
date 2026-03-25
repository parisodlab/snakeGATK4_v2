"""
Compute mean site depth (INFO/DP) from a VCF and write to the provided output file.
This script is intended to be used from a Snakemake rule via `script:`. It reads the
wildcard/params from Snakemake's runtime.
"""

from snakemake.shell import shell
import subprocess

# snakemake provides `snakemake` object when executed as script
vcf = snakemake.params.get('coverage_vcf') if hasattr(snakemake, 'params') else None
out = snakemake.output[0]

if vcf is None:
    raise SystemExit('coverage_vcf parameter not provided to compute_mean_coverage.py')

cmd = f"bcftools query -f '%INFO/DP\n' -i 'INFO/DP>0' {vcf}"
proc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
vals = []
for line in proc.stdout:
    line = line.strip()
    if line == '' or line == '.':
        continue
    try:
        vals.append(float(line))
    except Exception:
        continue
proc.stdout.close()
proc.wait()
mean = sum(vals)/len(vals) if len(vals) > 0 else 0.0
with open(out, 'w') as fh:
    fh.write(str(mean) + '\n')
