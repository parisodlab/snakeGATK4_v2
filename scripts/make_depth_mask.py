import argparse
import gzip
import numpy

parser = argparse.ArgumentParser(description='Creates depth mask for loci with depth >mean+2*sd')
parser.add_argument('-v', '--vcf')
parser.add_argument('-o', '--output')
#!/usr/bin/env python3
"""Create depth mask: compute per-sample mean/std and select loci with
depth above per-sample cutoff (mean+2*sd) in more than N samples.

Writes:
- per-locus counts TSV
- selected loci TSV (chrom\tpos)
- histogram TSV
- human-readable stats and TSV summary
- merged BED of consecutive selected positions
"""

import argparse
import gzip
import os
from collections import defaultdict
import numpy as np


def open_vcf(path):
    if path.endswith('.gz'):
        return gzip.open(path, 'rt')
    return open(path, 'r')


def main():
    parser = argparse.ArgumentParser(description='Creates depth mask for loci with depth > mean+2*sd')
    parser.add_argument('-v', '--vcf', required=True, help='Input VCF (gz allowed)')
    parser.add_argument('-o', '--output', required=True, help='Output list of passing loci (chrom\tpos per line)')
    parser.add_argument('-c', '--countout', required=True, help='Per-locus counts output (tsv chrom\tpos\tcount)')
    parser.add_argument('-p', '--plothist', required=True, help='Histogram output (tsv count\tfrequency)')
    parser.add_argument('-n', '--nsamples', required=True, help='Number of samples threshold')
    parser.add_argument('-l', '--locinumber', required=True, help='Estimated number of loci to pre-allocate')
    parser.add_argument('-m', '--min_depth', required=True, help='Minimum depth to include in mean/std')
    parser.add_argument('--stats', required=False, help='Human-readable stats output path')
    parser.add_argument('--summary', required=False, help='TSV summary (sample,mean,std,cutoff) path')
    parser.add_argument('--bed', required=False, help='Merged BED output path')

    args = parser.parse_args()

    vcf_input = args.vcf
    out_loci = args.output
    counts = args.countout
    histcounts = args.plothist
    nsamples = int(round(float(args.nsamples), 0))
    nloci = int(args.locinumber)
    min_depth = int(args.min_depth)
    stats_out = args.stats
    summary_tsv = args.summary
    bed_out = args.bed

    # ensure parent dirs exist for outputs
    for p in (out_loci, counts, histcounts, stats_out, summary_tsv, bed_out):
        if p:
            d = os.path.dirname(p)
            if d:
                os.makedirs(d, exist_ok=True)

    inp = open_vcf(vcf_input)

    loci_counter = 0
    sample_ids = []
    #!/usr/bin/env python3
    """Create depth mask: compute per-sample mean/std and select loci with
    depth above per-sample cutoff (mean+2*sd) in more than N samples.

    Writes:
    - per-locus counts TSV
    - selected loci TSV (chrom\tpos)
    - histogram TSV
    - human-readable stats and TSV summary
    - merged BED of consecutive selected positions
    """

    import argparse
    import gzip
    import os
    from collections import defaultdict
    import numpy as np


    def open_vcf(path):
        if path.endswith('.gz'):
            return gzip.open(path, 'rt')
        return open(path, 'r')


    def main():
        parser = argparse.ArgumentParser(description='Creates depth mask for loci with depth > mean+2*sd')
        parser.add_argument('-v', '--vcf', required=True, help='Input VCF (gz allowed)')
        parser.add_argument('-o', '--output', required=True, help='Output list of passing loci (chrom\tpos per line)')
        parser.add_argument('-c', '--countout', required=True, help='Per-locus counts output (tsv chrom\tpos\tcount)')
        parser.add_argument('-p', '--plothist', required=True, help='Histogram output (tsv count\tfrequency)')
        parser.add_argument('-n', '--nsamples', required=True, help='Number of samples threshold')
        parser.add_argument('-l', '--locinumber', required=True, help='Estimated number of loci to pre-allocate')
        parser.add_argument('-m', '--min_depth', required=True, help='Minimum depth to include in mean/std')
        parser.add_argument('--stats', required=False, help='Human-readable stats output path')
        parser.add_argument('--summary', required=False, help='TSV summary (sample,mean,std,cutoff) path')
        parser.add_argument('--bed', required=False, help='Merged BED output path')

        args = parser.parse_args()

        vcf_input = args.vcf
        out_loci = args.output
        counts = args.countout
        histcounts = args.plothist
        nsamples = int(round(float(args.nsamples), 0))
        nloci = int(args.locinumber)
        min_depth = int(args.min_depth)
        stats_out = args.stats
        summary_tsv = args.summary
        bed_out = args.bed

        # ensure parent dirs exist for outputs
        for p in (out_loci, counts, histcounts, stats_out, summary_tsv, bed_out):
            if p:
                d = os.path.dirname(p)
                if d:
                    os.makedirs(d, exist_ok=True)

        inp = open_vcf(vcf_input)

        loci_counter = 0
        sample_ids = []
        sample_length = 0
        data_array = None
        name_array = None

        for line in inp:
            if not line.strip():
                continue
            if line.startswith('##'):
                continue
            if line.startswith('#'):
                cols = line.strip().split('\t')
                sample_ids = cols[9:]
                sample_length = len(sample_ids)
                # preallocate
                data_array = np.zeros((nloci, sample_length), dtype=np.int32)
                name_array = np.full(nloci, '', dtype=object)
                continue

            # skip malformed lines
            if '\t' not in line:
                continue

            if loci_counter >= data_array.shape[0]:
                add_size = 100000
                data_array = np.vstack((data_array, np.zeros((add_size, sample_length), dtype=np.int32)))
                name_array = np.concatenate((name_array, np.full(add_size, '', dtype=object)))
                print('Warning: --locinumber was too small; grew arrays by', add_size)

            cols = line.strip().split('\t')
            if len(cols) < 9:
                continue
            chrom = cols[0]
            pos = cols[1]
            name_array[loci_counter] = chrom + '\t' + pos

            fmt = cols[8].split(':') if len(cols) > 8 else []
            dp_idx = fmt.index('DP') if 'DP' in fmt else None
            samples_fields = cols[9:]
            for i, sf in enumerate(samples_fields):
                if dp_idx is None:
                    data_array[loci_counter, i] = 0
                    continue
                fields = sf.split(':')
                if len(fields) <= dp_idx:
                    data_array[loci_counter, i] = 0
                    continue
                val = fields[dp_idx]
                if val == '.' or val == '':
                    data_array[loci_counter, i] = 0
                    continue
                try:
                    data_array[loci_counter, i] = int(val)
                except ValueError:
                    data_array[loci_counter, i] = 0

            loci_counter += 1

        inp.close()

        if loci_counter == 0:
            raise SystemExit('No loci parsed from VCF - exiting')

        # compute per-sample mean/std using values >= min_depth
        n_samples = len(sample_ids)
        data_mean = np.zeros(n_samples, dtype=float)
        data_std = np.zeros(n_samples, dtype=float)
        for j in range(n_samples):
            col = data_array[:loci_counter, j]
            valid = col[col >= min_depth]
            if valid.size > 0:
                data_mean[j] = float(np.mean(valid))
                data_std[j] = float(np.std(valid))
            else:
                data_mean[j] = 0.0
                data_std[j] = 0.0

        cutoff = data_mean + (2.0 * data_std)

        # default stats paths if not provided
        if stats_out is None:
            stats_out = counts + '.sample_depth_stats.txt'
        if summary_tsv is None:
            summary_tsv = counts + '.sample_mean_sd.tsv'
        if bed_out is None:
            # default bed path based on out_loci
            if out_loci.endswith('.tsv'):
                bed_out = out_loci.replace('.tsv', '.bed')
            else:
                bed_out = out_loci + '.bed'

        # write stats and summary
        with open(stats_out, 'w') as sf, open(summary_tsv, 'w') as st:
            st.write('sample\tmean\tstd\tcutoff\n')
            for i, sid in enumerate(sample_ids):
                sf.write('Sample %s: Mean depth: %s, Standard deviation: %s, Depth cutoff: %s\n' %
                         (sid, data_mean[i], data_std[i], cutoff[i]))
                st.write(f"{sid}\t{data_mean[i]}\t{data_std[i]}\t{cutoff[i]}\n")

        # histogram and select
        hist_dict = defaultdict(int)
        selected_sites = []
        with open(counts, 'w') as cntf, open(out_loci, 'w') as outf:
            for i in range(loci_counter):
                depths = data_array[i, :].astype(float)
                depth_count = int((depths > cutoff).sum())
                hist_dict[depth_count] += 1
                cntf.write(f"{name_array[i]}\t{depth_count}\n")
                if depth_count > nsamples:
                    outf.write(f"{name_array[i]}\n")
                    try:
                        chrom, pos_s = name_array[i].split('\t')
                        pos_i = int(pos_s)
                        selected_sites.append((chrom, pos_i))
                    except Exception:
                        pass

        # write histogram
        with open(histcounts, 'w') as hf:
            for k in sorted(hist_dict.keys()):
                hf.write(f"{k}\t{hist_dict[k]}\n")

        # produce merged BED from selected_sites (0-based start, end-exclusive)
        sites_by_chrom = defaultdict(list)
        for chrom, pos in selected_sites:
            sites_by_chrom[chrom].append(pos)

        with open(bed_out, 'w') as mb:
            for chrom in sorted(sites_by_chrom.keys()):
                positions = sorted(set(sites_by_chrom[chrom]))
                if not positions:
                    continue
                interval_start = positions[0]
                prev = interval_start
                for p in positions[1:]:
                    if p == prev + 1:
                        prev = p
                        continue
                    mb.write(f"{chrom}\t{interval_start - 1}\t{prev}\n")
                    interval_start = p
                    prev = p
                mb.write(f"{chrom}\t{interval_start - 1}\t{prev}\n")

        print('Wrote:', out_loci, counts, histcounts, stats_out, summary_tsv, bed_out)


    if __name__ == '__main__':
        main()
