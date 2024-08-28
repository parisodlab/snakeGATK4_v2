#!/usr/bin/env python
import pandas as pd
import yaml
import os
import time

with open("config_main.yaml") as yamlfile:
    config = yaml.full_load(yamlfile)

# Parameters to control the workflow
project = config["PROJECT"]

## Do you need to do quality control?
qc = config["QC"]
print("Is quality control required?\n", qc)

## Do you need to do trimming?
trim = config["TRIMMING"]
print("Is trimming required?\n", trim)


# Start the workflow
print("Start workflow on project: " + project)

## write the running time in a log file
file_log_time = open("logs/log_running_time.txt", "a+")
file_log_time.write("\nProject name: " + project + "\n")
file_log_time.write("Start time: " + time.ctime() + "\n")


def spend_time(start_time, end_time):
    seconds = end_time - start_time
    hours = seconds // 3600
    seconds %= 3600
    minutes = seconds // 60
    seconds %= 60

    return "%d:%02d:%02d" % (hours, minutes, seconds)


if qc:
    print("Start QC!")
    start_time = time.time()
    os.system(
        "snakemake --rerun-triggers mtime -s workflow/quality_control.rules --profile slurm 2>&1 | tee logs/log_quality_control.txt"
    )
    end_time = time.time()
    file_log_time.write(
        "Time of running QC: " + spend_time(start_time, end_time) + "\n"
    )
    print(
        "Quality control is done!\n Please check the report and decide whether trimming is needed\n Please remember to turn off the QC in the config file!"
    )
    os._exit(0)
else:
    if trim:
        print("Start Trimming!")
        start_time = time.time()
        os.system(
            "snakemake -s workflow/trim.rules --profile slurm 2>&1 | tee logs/log_trim.txt"
        )
        end_time = time.time()
        file_log_time.write(
            "Time of running trimming:" + spend_time(start_time, end_time) + "\n"
        )
        print("Trimming is done!")
    else:
        print("Trimming is not required")


print("Mapping is starting!")
start_time = time.time()
os.system(
    "snakemake --rerun-triggers mtime -s workflow/mapping.rules --profile slurm 2>&1 | tee logs/log_map.txt"
)
end_time = time.time()
file_log_time.write(
    "Time of running mapping:" + spend_time(start_time, end_time) + "\n"
)
print("Mapping is done!")

print("Variant calling using GATK4 is starting!")
start_time = time.time()
os.system(
    "snakemake --rerun-triggers mtime -s workflow/calling_gatk4.rules --profile slurm 2>&1 | tee logs/log_calling_gatk4.txt"
)
end_time = time.time()
file_log_time.write(
    "Time of running variant calling:" + spend_time(start_time, end_time) + "\n"
)
print("Variant calling is done!")

# print("Variant calling using freebayes is starting!")
# start_time = time.time()
# os.system(
#     "snakemake -np --rerun-triggers mtime -s workflow/calling_freebayes.rules --groups VariantCallingFreebayes=group_{sample} --group-components group_{sample}=2 --profile slurm 2>&1 | tee logs/log_calling_freebayes.txt"
# )
# end_time = time.time()
# file_log_time.write(
#     "Time of running variant calling:" + spend_time(start_time, end_time) + "\n"
# )
# print("Variant calling is done!")
