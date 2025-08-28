library(tidyverse)
library(data.table)

# Load the scaffold lengths
scaffold_lengths <- fread("/data/users/rchoudhury/snakeSNPcalling/data/reference/Varia_review2.FINAL.Syri.fa.fai", sep = "\t", header = FALSE)

# Load the scaffolds that have gene annotations
scaffold_genes <- fread(cmd = "grep -v '^#' /data/users/rchoudhury/snakeSNPcalling/data/reference/Varia.all.maker.renamed.gff", sep = "\t", header = FALSE) %>%
    select(V1) %>%
    distinct() %>%
    rename(scaffold = V1)

# get the scaffold lengths of the scaffolds that have gene annotations
scaffold_genes_lengths <- scaffold_lengths %>%
    filter(V1 %in% scaffold_genes$scaffold) %>%
    select(scaffold = V1, length = V2)

# get the minimum length of the scaffolds that have gene annotations
min_length <- scaffold_genes_lengths %>%
    summarise(min_length = min(length)) %>%
    pull(min_length)


# there are 267 scaffolds and 217 have gene annotations, I decided to use all the scaffolds
