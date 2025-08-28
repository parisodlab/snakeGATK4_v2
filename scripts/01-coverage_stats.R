library(tidyverse)
library(data.table)
# Define the path
path <- "results/Biscutella_remapping_varia/bamcoverage/samtools"

# List all the files matching the pattern
files <- list.files(path = path, pattern = "_coverage.txt", full.names = TRUE)

coverage.table <- do.call(rbind, lapply(files, function(file) {
    # Read the file
    data <- fread(file)

    # Extract the sample name from the filename
    sample_name <- gsub("_coverage.txt", "", basename(file))

    # Add a new column with the sample name
    data$sample <- sample_name

    return(data)
}))


coverage.table %>% head()

# select when the first column is Bv1 to Bv9

coverate.stats = coverage.table %>%
    filter(V1 %in% paste0("Bv", 1:9)) %>%
    select(V1, V7, sample) %>%
    rename(Chromosome = V1, Coverage = V7) %>%
    mutate(Coverage = as.numeric(Coverage)) %>%
    # make Chromosomes 1 to 9 as columns
    pivot_wider(names_from = Chromosome, values_from = Coverage)

file_stat <- list.files(path = "results/Biscutella_remapping_varia/coverage", pattern = "_summaryStats.txt", full.names = TRUE)

summary.stat = do.call(rbind, (lapply(file_stat, function(file) {
    # Read the file
    data <- fread(file, header = FALSE)

    # Extract the sample name from the filename
    sample_name <- gsub("_summaryStats.txt", "", basename(file))

    # Add a new column with the sample name
    data$sample <- sample_name
    # rename the columns

    return(data)
})))

# rename the columns
colnames(summary.stat) <- c("stat", "value", "sample")
summary.stat <- summary.stat %>%
    mutate(stat = gsub(" ", "_", stat)) %>%
    pivot_wider(names_from = stat, values_from = value)
# rename the columns
colnames(summary.stat) <- c("sample", "Average_read_depth_across_genome_>0mapq", "%_genome_(breadth)_covered_>0x", "%_genome_(breadth)_covered_>10x")

# get read coverage from seqkit stats
file_stat <- list.files(path = "results/Biscutella_remapping_varia/data/seqkit", pattern = "_stats.txt", full.names = TRUE)

seqkit.stat <- do.call(rbind, lapply(file_stat, function(file) {
    # Read the file
    data <- fread(file, header = TRUE, sep = "\t")

    # Extract the sample name from the filename
    prefix_name <- gsub("_R.*", "", basename(file))

    # Add a new column with the sample name
    data$prefix <- prefix_name
    # rename the columns

    return(data)
}))
# get prefix to sample name conversion from metadata
metadata <- read.table("configs/metadata.tsv", header = TRUE, sep = "\t") %>%
    select(sample, prefix, ploidy)
seqkit.stat <- seqkit.stat %>%
    left_join(metadata, by = c("prefix")) %>%
    select(sample, sum_len, ploidy) %>%
    rename("Total_length" = sum_len) %>%
    group_by(sample) %>%
    reframe(Total_length = sum(Total_length), ploidy = ploidy) %>%
    unique() %>%
    mutate(estimated_read_depth = Total_length / 832241633)


# join the two tables
coverage.stats_full <- coverate.stats %>%
    left_join(summary.stat, by = "sample") %>%
    left_join(seqkit.stat, by = "sample") %>%
    select(sample, ploidy, Bv1:Bv9, "Total_length", "estimated_read_depth", "Average_read_depth_across_genome_>0mapq", "%_genome_(breadth)_covered_>0x", "%_genome_(breadth)_covered_>10x")


# write to a file
write.table(coverage.stats_full, file = "results/Biscutella_remapping_varia/coverage_stats.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
