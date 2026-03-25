#!/usr/bin/env python3
import csv
from pathlib import Path

mapping_text = '''sample	pop	elevation
A2-GA5	DAG	low
A2-GB13	DAG	low
A2-GB14	DAG	low
A2-GB16	DAG	low
A2-GB17	DAG	low
A2-GB9	DAG	low
A2-GC18	DAG	low
A2-GX	DAG	low
A2-Sch-12	DAS	high
A2-Sch-16	DAS	high
A2-Sch-17	DAS	high
A2-Sch-19	DAS	high
A2-Sch-21	DAS	high
A2-Sch-5	DAS	high
A2-Sch-8	DAS	high
A2-Sch-X	DAS	high
B2-AUR-X	DBM	high
B2-MV1	DBM	high
B2-MV10	DBM	high
B2-MV18	DBM	high
B2-MV2	DBM	high
B2-MV20	DBM	high
B2-MV3	DBM	high
B2-MV6	DBM	high
B2-MV9	DBM	high
BUC-10	TBU	low
BUC-15	TBU	low
BUC-16	TBU	low
BUC-18	TBU	low
BUC-3	TBU	low
BUC-5	TBU	low
BUC-7	TBU	low
BUC-9	TBU	low
FUL-1	TFU	low
FUL-10	TFU	low
FUL-11	TFU	low
FUL-14	TFU	low
FUL-3	TFU	low
FUL-4	TFU	low
FUL-7	TFU	low
FUL-9	TFU	low
K2-s2-12	DKS	low
K2-s2-3	DKS	low
K2-s2-4	DKS	low
K2-s2-9	DKS	low
K2-s2-A12	DKS	low
K2-s2-A17	DKS	low
K2-s2-A5	DKS	low
K2-s2-X	DKS	low
NAY-Aa3	TNA	high
NAY-Cb4	TNA	high
NAY-Da2	TNA	high
NAY-Fa18	TNA	high
NAY-Ha7	TNA	high
NAY-Ic2	TNA	high
NAY-Lc5	TNA	high
NAY-Z19	TNA	high
PRE-10	TPR	low
PRE-11	TPR	low
PRE-2	TPR	low
PRE-4	TPR	low
PRE-5	TPR	low
PRE-7	TPR	low
PRE-8	TPR	low
PRE-9	TPR	low
RAM-1	TRA	high
RAM-11	TRA	high
RAM-15	TRA	high
RAM-16	TRA	high
RAM-18	TRA	high
RAM-20	TRA	high
RAM-3	TRA	high
RAM-7	TRA	high
RCB-A16	DPR	high
RCB-B15	DPR	high
RCB-B17	DPR	high
RCB-B18	DPR	high
RCB-B19	DPR	high
RCB-C19	DPR	high
RCB-C8	DPR	high
RCB-X	DPR	high
s3-11	TDS	high
s3-12	TDS	high
s3-13	TDS	high
s3-3	TDS	high
s3-5	TDS	high
s3-6	TDS	high
s3-7	TDS	high
s3-8	TDS	high
V2-b11	DVB	low
V2-b18	DVB	low
V2-b19	DVB	low
V2-b20	DVB	low
V2-b3	DVB	low
V2-bA4	DVB	low
V2-bA6	DVB	low
V2-bX	DVB	low
'''

map_reader = csv.DictReader(mapping_text.strip().splitlines(), delimiter='\t')
mapping = {r['sample']: {'pop': r['pop'], 'elevation': r['elevation']} for r in map_reader}

infile = Path('configs/metadata_fixed.tsv')
outfile = infile

# read existing
with infile.open() as fh:
    reader = csv.reader(fh, delimiter='\t')
    rows = list(reader)

if not rows:
    raise SystemExit('Empty metadata file')

header = rows[0]
# ensure pop and elevation columns exist
if 'pop' not in header:
    header += ['pop']
if 'elevation' not in header:
    header += ['elevation']

# index of sample column
try:
    sample_idx = header.index('sample')
except ValueError:
    raise SystemExit('No sample column in metadata')

# build updated rows
new_rows = [header]
for r in rows[1:]:
    # pad row if shorter than header
    if len(r) < len(header):
        r = r + [''] * (len(header) - len(r))
    sname = r[sample_idx]
    if sname in mapping:
        r[header.index('pop')] = mapping[sname]['pop']
        r[header.index('elevation')] = mapping[sname]['elevation']
    new_rows.append(r)

# write back
with outfile.open('w') as fh:
    writer = csv.writer(fh, delimiter='\t', lineterminator='\n')
    writer.writerows(new_rows)

print('Updated', outfile)
