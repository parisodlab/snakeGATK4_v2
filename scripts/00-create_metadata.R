library(tidyverse)
library(data.table)
file <- "/data/users/rchoudhury/snakeSNPcalling/configs/samples.tsv"

data <- fread(file, sep = "\t", header = TRUE)

sdat <- data %>%
    mutate(
        sample_barcode = sample,
        lane = str_extract(forward_file_name, "_L[0-9]+"), # before it was just "L[0-9]+" so was not working with FUL or samples with L in name
        # Extract everything until R1 or R2 and remove the last underscore
        prefix = gsub("_R.*", "", forward_file_name),
        ploidy = ifelse(str_detect(library_name, "tetras"), 4, 2),
        fq1 = paste0("/data/users/rchoudhury/Biscutella_wgs_remapping_2025/data/raw/all/", forward_file_name),
        fq2 = paste0("/data/users/rchoudhury/Biscutella_wgs_remapping_2025/data/raw/all/", reverse_file_name)
    ) %>%
    select(
        sample,
        prefix,
        sample_barcode,
        lane,
        ploidy,
        fq1,
        fq2
    )

# Create a metadata file
write_tsv(sdat, "/data/users/rchoudhury/snakeSNPcalling/configs/metadata.tsv", col_names = TRUE)
