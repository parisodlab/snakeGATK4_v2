"""
Filter indels with GATK-like thresholds, excess-depth exclusion and indel length limit.
This script is intended to be called from Snakemake via `script:`. It uses the `snakemake`
object to access inputs, outputs, params and threads.
"""
import os
import tempfile
import subprocess

# Inputs
dp_masked = snakemake.input.dp_masked
mean_cov = snakemake.input.mean_cov if isinstance(snakemake.input.mean_cov, str) else snakemake.input[1]

# Outputs
out_vcf = snakemake.output[0]
out_stats = snakemake.output[1]

# Params
params = snakemake.params
max_len = int(params.get('max_len', 20))
excess_mult = float(params.get('excess_mult', 2))
FS = params.get('FS', 40.0)
SOR = params.get('SOR', 3.0)
MQ = params.get('MQ', 30)
MQRankSum = params.get('MQRankSum', -12.5)
QD = params.get('QD', 2.0)
ReadPosRankSum = params.get('ReadPosRankSum', -4.0)

threads = snakemake.threads

# read mean coverage
mean = 0.0
with open(mean_cov) as m:
    txt = m.read().strip()
    mean = float(txt) if txt not in ('', '\n') else 0.0

max_depth = float(excess_mult) * mean

tmp = tempfile.mkdtemp()
tmp_v = os.path.join(tmp, 'indels_filtered_site.vcf.gz')

# Build bcftools include expression
expr = (
    f"INFO/FS <= {FS} && INFO/SOR <= {SOR} && INFO/MQ >= {MQ} && "
    f"INFO/MQRankSum >= {MQRankSum} && INFO/QD >= {QD} && INFO/ReadPosRankSum >= {ReadPosRankSum} && "
    f"INFO/DP <= {max_depth}"
)

# run bcftools view with the expression
subprocess.run([
    'bcftools', 'view', '-i', expr, '--threads', str(threads), '-Oz', '-o', tmp_v, dp_masked
], check=True)
subprocess.run(['bcftools', 'index', '--tbi', tmp_v], check=True)

# write header and body (filtered by indel absolute length)
hdr = os.path.join(tmp, 'hdr.vcf')
body = os.path.join(tmp, 'body.vcf')

subprocess.run(f"bcftools view -h {tmp_v} > {hdr}", shell=True, check=True)

# Filter by indel absolute length using Python: keep records where any ALT allele has
# abs(len(REF) - len(ALT)) <= max_len
with open(body, 'w') as bout:
    p = subprocess.Popen(["bcftools", "view", "-H", tmp_v], stdout=subprocess.PIPE, text=True)
    for line in p.stdout:
        line = line.rstrip('\n')
        if line == '':
            continue
        cols = line.split('\t')
        if len(cols) < 5:
            continue
        ref = cols[3]
        alts = cols[4].split(',') if cols[4] not in ('.', '') else []
        keep = False
        for alt in alts:
            try:
                d = abs(len(ref) - len(alt))
            except Exception:
                d = 0
            if d <= max_len:
                keep = True
                break
        if keep:
            bout.write(line + '\n')
    p.stdout.close()
    p.wait()

# concatenate header and body, bgzip and index
subprocess.run(f"cat {hdr} {body} | bgzip -c > {out_vcf}", shell=True, check=True)
subprocess.run(['bcftools', 'index', '--tbi', out_vcf], check=True)
subprocess.run(f"bcftools stats -s- --threads {threads} --verbose {out_vcf} > {out_stats}", shell=True, check=True)

# cleanup
subprocess.run(['rm', '-rf', tmp], check=True)
