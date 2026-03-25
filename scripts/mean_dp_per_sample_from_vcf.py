from collections import defaultdict
import gzip


def open_text(path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "rt")


vcf_path = str(snakemake.input.vcf)
output_path = str(snakemake.output.tsv)
label = str(snakemake.params.label)

sample_names = []
total_sites = 0
depth_sums = defaultdict(float)
depth_counts = defaultdict(int)

with open_text(vcf_path) as handle:
    for line in handle:
        if line.startswith("##"):
            continue
        if line.startswith("#CHROM"):
            parts = line.rstrip("\n").split("\t")
            sample_names = parts[9:]
            continue

        parts = line.rstrip("\n").split("\t")
        if len(parts) < 10:
            continue

        format_keys = parts[8].split(":")
        if "DP" not in format_keys:
            continue

        total_sites += 1
        dp_index = format_keys.index("DP")
        for sample_name, sample_field in zip(sample_names, parts[9:]):
            values = sample_field.split(":")
            if dp_index >= len(values):
                continue
            raw_dp = values[dp_index]
            if raw_dp in {".", ""}:
                continue
            try:
                depth_value = float(raw_dp)
            except ValueError:
                continue
            depth_sums[sample_name] += depth_value
            depth_counts[sample_name] += 1

with open(output_path, "wt") as out_handle:
    out_handle.write("sample\tmean_dp\tsites_with_dp\tfraction_sites_with_dp\tvcf_label\n")
    for sample_name in sample_names:
        site_count = depth_counts.get(sample_name, 0)
        if site_count == 0:
            mean_dp = "NA"
            fraction = "0"
        else:
            mean_dp = f"{depth_sums[sample_name] / site_count:.6f}"
            fraction = f"{site_count / total_sites:.6f}" if total_sites else "0"
        out_handle.write(
            f"{sample_name}\t{mean_dp}\t{site_count}\t{fraction}\t{label}\n"
        )
