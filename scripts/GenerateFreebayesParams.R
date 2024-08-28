#!/usr/bin/env Rscript
log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

## Generate freebayes params ##
# 1) make bed for parallelising freebayes
library(dplyr)
library(data.table)
library(glue)


load_metadata <- function(metadata_path) {
    # Check the file extension and load metadata accordingly
    if (tools::file_ext(metadata_path) == "xlsx") {
        metadata <- readxl::read_excel(metadata_path)
    } else if (tools::file_ext(metadata_path) == "tsv") {
        metadata <- data.table::fread(metadata_path, sep = "\t")
    } else if (tools::file_ext(metadata_path) == "txt") {
        metadata <- data.table::fread(metadata_path, sep = "\t")
    } else if (tools::file_ext(metadata_path) == "csv") {
        metadata <- data.table::fread(metadata_path, sep = ",")
    } else {
        stop("Metadata file must be .xlsx, .tsv, or .csv")
    }
    return(metadata)
}

# read inputs
final_path <- snakemake@params[["final_path"]]
scaffolds <- snakemake@params[["scaffold"]]
fai <- fread(snakemake@input[["index"]])

# select scaffolds we want, and start, end columns
fai <- fai[fai$V1 %in% scaffolds, c(1, 2)]


# 2) Make bamlist and populations.tsv file

metadata <- load_metadata(snakemake@params[["metadata"]])
metadata$bams <- paste0(final_path, "/dedup/", metadata$sample, ".bam")

metadata %>%
    select(bams, pop) %>%
    unique() %>%
    fwrite(., snakemake@output[["pops"]], sep = "\t", row.names = FALSE, col.names = FALSE)

metadata %>%
    select(bams) %>%
    unique() %>%
    fwrite(., snakemake@output[["bamlist"]], sep = "\t", row.names = FALSE, col.names = FALSE)

sessionInfo()
