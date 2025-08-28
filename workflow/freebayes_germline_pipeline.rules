##############################################################################
# freebayes_germline_pipeline.smk
##############################################################################
import os

# 1) Load the config file
configfile: "config.yaml"

# Retrieve config values
BAM_LIST_FILE     = config["bam_list_file"]
REFERENCE_GENOME  = config["reference_genome"]
SCATTER_MODE      = config["scatter_mode"]         # "chromosome" or "interval"
CHROMS            = config.get("chromosomes", [])
SCATTER_COUNT     = config.get("scatter_count", 25)
FREEBAYES_PARAMS  = config.get("freebayes_params", {})
RESULTS_DIR       = config["results_dir"]
LOGS_DIR          = config["logs_dir"]
INTERVALS_DIR     = config["intervals_dir"]
FINAL_MERGED_VCF  = config["final_merged_vcf"]

FREEBAYES_ENV     = config["freebayes_env"]  # e.g. "freebayes"
GATK_ENV          = config["gatk_env"]       # e.g. "gatk"
BCFTOOLS_ENV      = config["bcftools_env"]   # e.g. "bcftools"

# Make sure directories exist
os.makedirs(RESULTS_DIR, exist_ok=True)
os.makedirs(LOGS_DIR, exist_ok=True)
os.makedirs(INTERVALS_DIR, exist_ok=True)

##############################################################################
# (A) Additional subfolder for BAM list files
##############################################################################
BAM_LIST_DIR = os.path.join(RESULTS_DIR, "bam_lists")
os.makedirs(BAM_LIST_DIR, exist_ok=True)

##############################################################################
# 2) Read the BAM list file so Snakemake knows the inputs
##############################################################################
with open(BAM_LIST_FILE, "r") as f:
    BAM_FILES = [line.strip() for line in f if line.strip()]

print(f"DEBUG: Found {len(BAM_FILES)} BAM(s) in {BAM_LIST_FILE}")

##############################################################################
# 3) Function to convert freebayes_params dict to a command-line string
#    If the value is "true", treat it as a flag-only (no value appended).##############################################################################
def freebayes_param_str(params_dict):
    """
    Convert a dictionary of FreeBayes params to a command-line string.
    If the config says e.g. --standard-filters: "true", output just "--standard-filters".
    Otherwise, output "--key value".
    """
    parts = []
    for key, val in params_dict.items():
        if val.lower() == "true":
            parts.append(key)
        else:
            parts.append(f"{key} {val}")
    return " ".join(parts)

FREEBAYES_CMD_PARAMS = freebayes_param_str(FREEBAYES_PARAMS)

##############################################################################
# 4) Define naming scheme for GATK's scattered intervals
##############################################################################
def scattered_interval_name(i):
    """Return e.g. '0000-scattered.interval_list' for i=0."""
    return f"{i:04d}-scattered.interval_list"

def scattered_bed_name(i):
    """Return the corresponding bed name, e.g. '0000-scattered.bed'."""
    return f"{i:04d}-scattered.bed"

SCATTERED_INTERVAL_LISTS = [
    os.path.join(INTERVALS_DIR, scattered_interval_name(i))
    for i in range(SCATTER_COUNT)
]
SCATTERED_BED_FILES = [
    os.path.join(INTERVALS_DIR, scattered_bed_name(i))
    for i in range(SCATTER_COUNT)
]

##############################################################################
# 5) The final "all" rule - produce one merged VCF
##############################################################################
rule all:
    """
    Final pipeline target: merged VCF from scattered FreeBayes calls.
    """
    input:
        FINAL_MERGED_VCF

##############################################################################
# 6) Scatter intervals by Ns (if SCATTER_MODE == 'interval')
##############################################################################
rule scatter_intervals_by_ns:
    """
    Scatter intervals across the reference genome's Ns using GATK ScatterIntervalsByNs.
    If SCATTER_MODE != 'interval', just create an empty file.
    """
    input:
        reference=REFERENCE_GENOME
    output:
        intervals_list=os.path.join(INTERVALS_DIR, "scattered.interval_list")
    log:
        os.path.join(LOGS_DIR, "scatter_intervals_by_ns.log")
    conda:
        GATK_ENV
    shell:
        r"""
        if [ "{SCATTER_MODE}" = "interval" ]; then
            set -e
            echo "Starting ScatterIntervalsByNs" > {log}
            gatk ScatterIntervalsByNs \
                -R {input.reference} \
                -O {output.intervals_list} &>> {log}
            echo "Finished ScatterIntervalsByNs" >> {log}
        else
            # If not interval mode, produce an empty intervals file
            touch {output.intervals_list}
        fi
        """

##############################################################################
# 7) Split intervals => produce scattered files + bed
##############################################################################
rule split_intervals:
    """
    Split intervals into SCATTER_COUNT parts using GATK SplitIntervals,
    then convert each .interval_list => .bed
    """
    input:
        intervals_list=os.path.join(INTERVALS_DIR, "scattered.interval_list"),
        reference=REFERENCE_GENOME
    output:
        # Mark these scattered intervals and bed files as temp so they can be removed after usage
        interval_lists=temp(SCATTERED_INTERVAL_LISTS),
        bed_files=temp(SCATTERED_BED_FILES)
    params:
        scatter_count=SCATTER_COUNT
    log:
        os.path.join(LOGS_DIR, "split_intervals.log")
    conda:
        GATK_ENV
    shell:
        r"""
        if [ "{SCATTER_MODE}" = "interval" ]; then
            set -e
            echo "Starting SplitIntervals" > {log}
            gatk SplitIntervals \
                -R {input.reference} \
                -L {input.intervals_list} \
                --scatter-count {params.scatter_count} \
                -O {INTERVALS_DIR} &>> {log}
            echo "Finished SplitIntervals" >> {log}

            # Convert each .interval_list => .bed
            for f in {INTERVALS_DIR}/*-scattered.interval_list; do
                bed="${{f%.interval_list}}.bed"
                sed 's/:\|-/\t/g' "$f" \
                    | grep -v "^@" \
                    | cut -f1-3 \
                    > "$bed"
            done
        else
            # Not scattering => create empty placeholders
            for i in $(seq 0 $(( {{params.scatter_count}} - 1 ))); do
                touch {INTERVALS_DIR}/$(printf "%04d-scattered.interval_list" $i)
                touch {INTERVALS_DIR}/$(printf "%04d-scattered.bed" $i)
            done
        fi
        """

##############################################################################
# 8) Determine scatter units (chromosome or interval)
##############################################################################
def get_scatter_units():
    if SCATTER_MODE == "chromosome":
        return CHROMS
    elif SCATTER_MODE == "interval":
        return [f"{i:04d}-scattered" for i in range(SCATTER_COUNT)]
    else:
        return []

SCATTER_UNITS = get_scatter_units()

##############################################################################
# 9) Exponential memory usage + threads for freebayes_scatter
##############################################################################
def exponential_mem_mb(wildcards, attempt):
    base_mem = 16384  # 16 GB in MB
    return base_mem * (2 ** (attempt - 1))

def exponential_threads(wildcards, attempt):
    base_threads = 2
    return base_threads * (2 ** (attempt - 1))

##############################################################################
# 9) FreeBayes scatter rule
##############################################################################
rule freebayes_scatter:
    """
    Run FreeBayes on each scatter chunk (chromosome or interval),
    compress the VCF with bgzip, then index it with tabix.
    Doubling memory & threads each retry.
    Mark the scattered FreeBayes VCF as temp so it can be removed post-merge.
    """
    input:
        reference=REFERENCE_GENOME,
        bam_files=BAM_FILES,
        bed_file=lambda wc: os.path.join(INTERVALS_DIR, f"{wc.scatter_id}.bed")
    output:
        vcf=temp(f"{RESULTS_DIR}/freebayes_{{scatter_id}}.vcf.gz"),
        log=f"{LOGS_DIR}/freebayes_{{scatter_id}}.log"
    threads: exponential_threads
    resources:
        mem_mb = exponential_mem_mb
    conda:
        FREEBAYES_ENV
    shell:
        r"""
        set -e
        echo "Starting FreeBayes on {wildcards.scatter_id} (threads={threads}, mem={resources.mem_mb}MB)" >&2

        REGION_ARG=""
        if [ "{SCATTER_MODE}" = "chromosome" ]; then
            REGION_ARG="-r chr{wildcards.scatter_id}"
        elif [ "{SCATTER_MODE}" = "interval" ]; then
            REGION_ARG="--targets {input.bed_file}"
        fi

        # (A) Store the BAM list in a subfolder
        BAM_LIST_TXT="{BAM_LIST_DIR}/bam_list_{wildcards.scatter_id}.txt"
        rm -f "$BAM_LIST_TXT"
        for b in {input.bam_files}; do
            echo "$b" >> "$BAM_LIST_TXT"
        done

        # Run FreeBayes, pipe to bgzip, then index with tabix
        freebayes \
            -f {input.reference} \
            $REGION_ARG \
            -L "$BAM_LIST_TXT" \
            {FREEBAYES_CMD_PARAMS} \
            2> {output.log} \
        | bgzip -c > {output.vcf}

        tabix -p vcf {output.vcf}

        echo "Finished FreeBayes on {wildcards.scatter_id}" >&2
        """

##############################################################################
# 10) Merge and normalize scattered .vcf.gz => final VCF
##############################################################################
rule merge_vcfs:
    """
    Merge scattered FreeBayes VCFs into final merged VCF (bgzipped + indexed).
    """
    input:
        expand(
            os.path.join(RESULTS_DIR, "freebayes_{scatter_id}.vcf.gz"),
            scatter_id=SCATTER_UNITS
        )
    output:
        merged_vcf=FINAL_MERGED_VCF
    threads: 12
    resources:
        mem_mb = 16384
    log:
        os.path.join(LOGS_DIR, "merge_freebayes.log")
    conda:
        BCFTOOLS_ENV
    shell:
        r"""
        set -e
        echo "Merging scattered FreeBayes VCFs" > {log}

        # 1) Concatenate with naive mode (-n),
        # 2) normalize with additional arguments,
        # 3) fill tags,
        # 4) write final compressed VCF
        bcftools concat \
            -n \
            -O z \
            -- threads 4 \
            {input} 2>> {log} \
        | bcftools norm \
            -m-any --force -a --atom-overlaps . \
            -f {REFERENCE_GENOME} \
            --threads 4 \
            - 2>> {log} \
        | bcftools +fill-tags \
            -O z \
            -o {output.merged_vcf} \
            --threads 4 \
            -W=tbi - -- \
            2>> {log}

        echo "Done merging final VCF." >> {log}
        """
