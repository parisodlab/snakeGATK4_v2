#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(SNPRelate)
  library(StAMPP)
  library(ape)
  library(dplyr)
  library(ggplot2)
  library(ggtree)
  library(tibble)
  library(adegenet)
  library(vcfR)
  library(dartR)
  library(pegas)
})

write_placeholder_plot <- function(out_plot, reason) {
  out_dir <- dirname(out_plot)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  pdf(out_plot, width = 10, height = 8)
  plot.new()
  text(0.5, 0.55, "NJ tree skipped", cex = 1.6)
  text(0.5, 0.45, reason, cex = 1.0)
  dev.off()
  writeLines(reason, file.path(out_dir, "nj_tree.skipped.txt"))
}

if (!requireNamespace("SNPfiltR", quietly = TRUE)) {
  install.packages("SNPfiltR", repos = "https://cloud.r-project.org")
}
library(SNPfiltR)

# install gplot if not already installed
if (!requireNamespace("gplots", quietly = TRUE)) {
  install.packages("gplots", repos = "https://cloud.r-project.org")
}

if (!requireNamespace("igraph", quietly = TRUE)) {
  install.packages("igraph", repos = "https://cloud.r-project.org")
}

if (!requireNamespace("phangorn", quietly = TRUE)) {
  install.packages("phangorn", repos = "https://cloud.r-project.org")
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: script.R <vcf_file> <metadata_file> <out_plot>")
}
# --- input args ---
# vcf_file      <- "/data/users/rchoudhury/snakeSNPcalling/results/Biscutella_remapping_varia/gatk4/final_vcf/4_fold_degenerate_filtered.ld_pruned.vcf.gz"
# metadata_file  <- "/data/users/rchoudhury/snakeSNPcalling/configs/metadata_fixed.tsv"
# out_plot      <- "/data/users/rchoudhury/snakeSNPcalling/results/Biscutella_remapping_varia/phylogeny/nj_tree.pdf"
vcf_file <- args[1]
metadata_file <- args[2]
out_plot <- args[3]

# --- read metadata ---
metadata <- read.table(metadata_file, header = TRUE, stringsAsFactors = FALSE)

# --- standardize metadata column names ---
col_names <- tolower(names(metadata))
map_col <- function(possible) {
  idx <- which(col_names %in% tolower(possible))
  if (length(idx) > 0) {
    return(names(metadata)[idx[1]])
  }
  return(NA)
}
sample_col <- map_col(c("sample", "sampleid", "sample_id", "ind", "id"))
pop_col <- map_col(c("pop", "population", "popid", "population_name"))
ploidy_col <- map_col(c("ploidy", "ploidy_level"))

if (is.na(sample_col) || is.na(pop_col) || is.na(ploidy_col)) {
  stop(sprintf(
    "Metadata must contain sample, population, ploidy. Found: sample=%s pop=%s ploidy=%s",
    sample_col, pop_col, ploidy_col
  ))
}

metadata <- metadata %>%
  select(!!sym(sample_col), !!sym(pop_col), !!sym(ploidy_col), dplyr::any_of("group")) %>%
  distinct()
names(metadata)[1:3] <- c("sampleID", "population", "ploidy")

## Convert VCF to genlight using vcfR and filter loci by call rate with dartR
vcf <- read.vcfR(vcf_file)
if (nrow(vcf@fix) == 0) {
  write_placeholder_plot(out_plot, "No fourfold SNPs were available after filtering and LD pruning.")
  quit(save = "no", status = 0)
}
vcf_filt_missing <- SNPfiltR::missing_by_snp(vcf, cutoff = .95)
vcf <- vcf_filt_missing
if (nrow(vcf@fix) == 0) {
  write_placeholder_plot(out_plot, "All fourfold SNPs were removed by missingness filtering before tree building.")
  quit(save = "no", status = 0)
}

helper_script <- "scripts/vcfR2genlight.tetra.R"
if (!file.exists(helper_script)) {
  stop(sprintf("Required helper script not found: '%s'. Make sure it exists and the working directory is correct.", helper_script))
}
# source safely and report errors
tryCatch(
  {
    source(helper_script)
  },
  error = function(e) {
    stop(sprintf("Failed to source helper '%s': %s", helper_script, conditionMessage(e)))
  }
)
# ensure the expected function was provided
if (!exists("vcfR2genlight", mode = "function")) {
  stop(sprintf("vcfR2genlight() not found after sourcing '%s'. Check the helper defines that function.", helper_script))
}

gl <- vcfR2genlight.tetra(vcf)
if (nLoc(gl) == 0 || nInd(gl) < 2) {
  write_placeholder_plot(out_plot, "Too few loci or samples remained to infer a neighbor-joining tree.")
  quit(save = "no", status = 0)
}

# sample IDs from VCF/genlight
vcf_samples <- indNames(gl)

# reorder genlight to match metadata order
keep_idx <- match(metadata$sampleID, indNames(gl))
if (any(is.na(keep_idx))) {
  missing_now <- metadata$sampleID[is.na(keep_idx)]
  stop(sprintf(
    "Sample ID mismatch: %d metadata sampleIDs not found in VCF: %s",
    length(missing_now), paste(head(missing_now, 10), collapse = ", ")
  ))
}
gl <- gl[keep_idx]
# quick checks of sample ID matching
mismatch <- which(indNames(gl) != metadata$sampleID)
if (length(mismatch)) {
  cat("Mismatches (showing up to 20):\n")
  print(data.frame(
    i = mismatch,
    indNames = indNames(gl)[mismatch],
    metadata = metadata$sampleID[mismatch],
    stringsAsFactors = FALSE
  )[1:20, ])
} else {
  cat("All sample IDs match and are in the same order.\n")
}

# attach population and ploidy from metadata
pop(gl) <- as.factor(metadata$population)

# --- ensure ploidy numeric (fallbacks for text labels) ---
parse_ploidy <- function(x) {
  x_chr <- as.character(x)
  num <- suppressWarnings(as.numeric(x_chr))
  is_num <- !is.na(num)
  res <- rep(NA_real_, length(x_chr))
  res[is_num] <- num[is_num]
  txt_idx <- which(is.na(res))
  if (length(txt_idx)) {
    txt <- tolower(x_chr[txt_idx])
    res[txt_idx[grepl("di", txt)]] <- 2
    res[txt_idx[grepl("tetra|4", txt)]] <- 4
    # default fallback to 2 (diploid) if still NA
    res[is.na(res)] <- 2
  }
  return(res)
}
ploidy_vec <- parse_ploidy(metadata$ploidy)
ploidy(gl) <- as.integer(ploidy_vec)


stampp_geno <- gl


# --- Nei’s distance ---
D <- stamppNeisD(stampp_geno, pop = FALSE) # Nei's 1972 distance between indivs


# get the folder of the output plot to save additional files
out_dir <- dirname(out_plot)

stamppPhylip(D, file = file.path(out_dir, "all_individuals_Neis_distance.phy.dst")) # write to phylip format

D.pop <- stamppNeisD(stampp_geno, pop = TRUE) # Nei's 1972 distance between pops
stamppPhylip(D.pop, file = file.path(out_dir, "all_pops_Neis_distance.phy.dst"))


############################################
# additional analyses using these distances
############################################


# --- NJ tree ---
tree <- nj(as.dist(D))

# --- prepare tip metadata ---
tip_df <- tibble(sampleID = tree$tip.label) %>%
  left_join(metadata, by = "sampleID") %>%
  # drop any pre-existing 'label' to avoid duplicate-name errors
  select(-dplyr::any_of("label")) %>%
  mutate(
    # add a 'label' column so %<+% can match tip labels
    label = sampleID,
    population = as.character(population),
    ploidy = ifelse(is.na(ploidy), "unknown", as.character(ploidy)),
    # normalize ploidy strings to match the scale_color_manual keys
    ploidy = dplyr::case_when(
      ploidy %in% c("2", "2.0", "2.00", "diploid") ~ "diploid",
      ploidy %in% c("4", "4.0", "4.00", "tetraploid") ~ "tetraploid"
    ),
    group = if ("group" %in% names(.)) as.character(group) else NA_character_
  )

# attach tip metadata to the tree data frame explicitly to avoid duplicate 'label' renaming
p <- ggtree(tree)

# left-join tip metadata: match tree$label to tip_df$sampleID
p$data <- p$data %>% dplyr::left_join(tip_df, by = c("label" = "sampleID"))

# use the tree's label column for tiplab and use ploidy for color
p <- p +
  geom_tiplab(aes(label = label, color = ploidy), size = 3) +
  scale_color_manual(values = c(diploid = "blue", tetraploid = "red", unknown = "gray")) +
  theme_tree2()

if ("group" %in% names(tip_df) && any(!is.na(tip_df$group))) {
  group_positions <- p$data %>%
    filter(isTip, !is.na(group)) %>%
    group_by(group) %>%
    summarise(x = mean(x), y = mean(y), .groups = "drop")
  p <- p + geom_label(
    data = group_positions,
    aes(x = x, y = y, label = group),
    inherit.aes = FALSE,
    fill = "white",
    alpha = 0.6,
    size = 3
  )
}

ggsave(out_plot, p, width = 10, height = 8)
write.tree(tree, file = file.path(out_dir, "all_individuals_Neis_distance.nwk")) # write to newick format
summary(D)
# write summary to a text file
summary_file <- file.path(out_dir, "Neis_distance_summary.txt")
sink(summary_file)
cat("Summary of Nei's Distance (individuals):\n")
print(summary(D))
sink()
summary(D.pop)
# write summary to a text file
summary_pop_file <- file.path(out_dir, "Neis_distance_pop_summary.txt")
sink(summary_pop_file)
cat("Summary of Nei's Distance (populations):\n")
print(summary(D.pop))
sink()

### heatmap of the pop distance matrix
pdf(file.path(out_dir, "Neis_distance_pop_heatmap.pdf"), width = 7, height = 7)
gplots::heatmap.2(D.pop, trace = "none", cexRow = 0.7, cexCol = 0.7)
dev.off()

### heatmap of the ind distance matrix
pdf(file.path(out_dir, "Neis_distance_ind_heatmap.pdf"), width = 7, height = 7)
gplots::heatmap.2(D, trace = "none", cexRow = 0.7, cexCol = 0.7)
dev.off()

### calculate AMOVA (differentiation among populations)
pops <- as.factor(metadata$population)
res <- pegas::amova(D ~ pops) # by default nperm=1000

# 3) inspect results (printed output usually contains the permutation p-value)
print(res)
str(res)

# 4) compute simple SSD proportion (effect size = SSD_ploidy / total_SSD)
tab <- res$tab
ssd_prop <- tab["ploidy_factor", "SSD"] / sum(tab[, "SSD"])
cat("SSD proportion (ploidy):", round(ssd_prop, 4), "\n")

# 5) save results to a file for record
outf <- file.path(out_dir, "AMOVA_pop_print.txt")
structure_file <- file.path(out_dir, "AMOVA_pop_structure.txt")
cat("AMOVA object print()\n\n", file = outf)
capture.output(print(res), file = outf, append = TRUE)
cat("\n\nAMOVA object str()\n\n", file = outf, append = TRUE)
capture.output(str(res), file = outf, append = TRUE)

cat("\n\nVariance components (sigma2, P.value) and percent of total variance\n\n", file = outf, append = TRUE)
vc <- res$varcomp$sigma2
pv <- res$varcomp$P.value
df_vc <- data.frame(
  component = rownames(res$varcomp), sigma2 = vc, P.value = pv,
  percent = round(vc / sum(vc) * 100, 3), row.names = NULL
)
write.table(df_vc, file = outf, sep = "\t", row.names = FALSE, quote = FALSE, append = TRUE)
# write also the varcoef
capture.output(res$varcomp, file = structure_file)


### calculate AMOVA (differentiation among ploidy levels)
ploidy_factor <- as.factor(tip_df$ploidy)
res_ploidy <- pegas::amova(D ~ ploidy_factor)
### write AMOVA results to a text file
# 3) inspect results (printed output usually contains the permutation p-value)
print(res_ploidy)
str(res_ploidy)

# 4) compute simple SSD proportion (effect size = SSD_ploidy / total_SSD)
tab <- res_ploidy$tab
ssd_prop <- tab["ploidy_factor", "SSD"] / sum(tab[, "SSD"])
cat("SSD proportion (ploidy):", round(ssd_prop, 4), "\n")

# 5) save results to a file for record
outf <- file.path(out_dir, "AMOVA_ploidy_print.txt")
structure_file <- file.path(out_dir, "AMOVA_ploidy_structure.txt")
cat("AMOVA object print()\n\n", file = outf)
capture.output(print(res_ploidy), file = outf, append = TRUE)
cat("\n\nAMOVA object str()\n\n", file = outf, append = TRUE)
capture.output(str(res_ploidy), file = outf, append = TRUE)

cat("\n\nVariance components (sigma2, P.value) and percent of total variance\n\n", file = outf, append = TRUE)
vc <- res_ploidy$varcomp$sigma2
pv <- res_ploidy$varcomp$P.value
df_vc <- data.frame(
  component = rownames(res_ploidy$varcomp), sigma2 = vc, P.value = pv,
  percent = round(vc / sum(vc) * 100, 3), row.names = NULL
)
write.table(df_vc, file = outf, sep = "\t", row.names = FALSE, quote = FALSE, append = TRUE)
# write also the varcoef
capture.output(res_ploidy$varcomp, file = structure_file)


### calculate AMOVA (differentiation among ploidy levels*population)
ploidy_factor <- as.factor(tip_df$ploidy)
pops <- as.factor(metadata$population)
res_all <- pegas::amova(D ~ ploidy_factor / pops)
### write AMOVA results to a text file
# 3) inspect results (printed output usually contains the permutation p-value)
print(res_all)
str(res_all)

# 4) compute simple SSD proportion (effect size = SSD_ploidy / total_SSD)
tab <- res_all$tab
ssd_prop <- tab["ploidy_factor", "SSD"] / sum(tab[, "SSD"])
cat("SSD proportion (ploidy):", round(ssd_prop, 4), "\n")

# 5) save results to a file for record
outf <- file.path(out_dir, "AMOVA_ploidy_and_population_print.txt")
structure_file <- file.path(out_dir, "AMOVA_ploidy_and_population_structure.txt")
cat("AMOVA object print()\n\n", file = outf)
capture.output(print(res_all), file = outf, append = TRUE)
cat("\n\nAMOVA object str()\n\n", file = outf, append = TRUE)
capture.output(str(res_all), file = outf, append = TRUE)

cat("\n\nVariance components (sigma2, P.value) and percent of total variance\n\n", file = outf, append = TRUE)
vc <- res_all$varcomp$sigma2
pv <- res_all$varcomp$P.value
df_vc <- data.frame(
  component = rownames(res_all$varcomp), sigma2 = vc, P.value = pv,
  percent = round(vc / sum(vc) * 100, 3), row.names = NULL
)
write.table(df_vc, file = outf, sep = "\t", row.names = FALSE, quote = FALSE, append = TRUE)
# write also the varcoef
capture.output(res_all$varcomp, file = structure_file)




mySum <- glSum(gl, alleleAsUnit = TRUE)
head(mySum)
png(file.path(out_dir, "AFS.png"), width = 1000, height = 500, units = "px")
barplot(table(mySum), col = "blue", space = 0, xlim = c(0, 80), xlab = "Allele counts", main = "Distribution of ALT allele counts in total dataset")
dev.off()

png(file.path(out_dir, "AFS.total.png"), width = 1000, height = 500, units = "px")
barplot(table(mySum), col = "blue", space = 0, xlab = "Allele counts", main = "Distribution of ALT allele counts in total dataset")
dev.off()



# Minimum Spanning Tree
net <- ape::mst(D)

pdf(file.path(out_dir, "MST_plot.pdf"), width = 8, height = 8)
plot(net, layout = igraph::layout_with_fr) # requires igraph
dev.off()

# neighborNet (prettier plot with coloured terminal branches)
nnet <- phangorn::neighborNet(as.dist(D), ord = NULL)
# Refine / fit the splits with least squares
# Calculates the expected pairwise distances between taxa implied by the splits (given some weights).
# Adjusts the split weights so that the expected distances match the observed distances as closely as possible.
# This is done by least squares: minimizing the sum of squared differences between observed and expected distances.
# The lambda and gamma arguments control penalty terms (regularization to avoid overfitting).
# Returns the optimized weights along with both:
# The unconstrained fit (pure least squares).
# The constrained fit (penalized).
# NeighborNet alone: gives you a network topology with heuristic weights.
# NeighborNet + splitsNetwork: gives you the same topology but with statistically optimized weights that better explain your original distance data.
# This is closer to what you’d do with trees:
# Neighbor-Joining → tree topology.
# Least-squares tree fitting → branch lengths.

fit <- phangorn::splitsNetwork(as.dist(D), splits = nnet$splits)
splitsnet <- phangorn::as.networx(fit)

# create colours for populations (fallback to ploidy if population missing)
tip_meta <- tip_df %>% tibble::deframe() # ensure tip_df is available
# match tip labels to metadata rows
tip_idx <- match(nnet$tip.label, tip_df$sampleID)
tip_pops <- as.character(tip_df$population[tip_idx])
tip_ploidy <- as.character(tip_df$ploidy[tip_idx])

# build a palette for populations (use ploidy as fallback group if population is NA)
groups <- ifelse(is.na(tip_pops) | tip_pops == "", tip_ploidy, tip_pops)
groups <- as.character(groups)
groups[is.na(groups)] <- "unknown"
uniq_groups <- unique(groups)

# diagnostics (helpful when tips end up grey)
cat("Number of tips:", length(groups), "\n")
cat("Unique groups:", length(uniq_groups), "\n")
print(table(groups)[1:20])
# pick colours robustly (use RColorBrewer if available, otherwise use rainbow)
if (requireNamespace("RColorBrewer", quietly = TRUE)) {
  # create an exact-length palette using colorRampPalette from a brewer base
  base_pal <- RColorBrewer::brewer.pal(min(8, max(3, length(uniq_groups))), "Set2")
  pal_cols <- grDevices::colorRampPalette(base_pal)(length(uniq_groups))
} else {
  pal_cols <- grDevices::rainbow(length(uniq_groups))
}
group_col_map <- setNames(pal_cols, uniq_groups)
tip_cols <- group_col_map[groups]
# final fallback for any NA colours
tip_cols[is.na(tip_cols)] <- "gray50"

# compute edge colours: terminal edges get the tip colour, internal edges are grey
ntip <- length(nnet$tip.label)
edge_mat <- nnet$edge
is_tip_edge <- edge_mat[, 2] <= ntip
edge_cols <- rep("gray60", nrow(edge_mat))
edge_cols[is_tip_edge] <- tip_cols[edge_mat[is_tip_edge, 2]]

pdf(file.path(out_dir, "neighbournet_plot_pop_colours.pdf"), width = 10, height = 8)
op <- par(mar = c(1, 1, 1, 1))
# plot should accept edge.color / edge.width; use smaller tip labels and no margins
plot(nnet,
  edge.color = edge_cols, edge.width = ifelse(is_tip_edge, 2, 1),
  tip.cex = 0.7, cex = 0.8, no.margin = TRUE
)
# add coloured tip labels (in case plot() didn't colour them)
tiplabels(pch = 19, col = tip_cols, cex = 0.8)
# legend for groups
legend("topright", legend = names(group_col_map), col = group_col_map, pch = 19, pt.cex = 1.0, bty = "n", ncol = 1)
par(op)
dev.off()


pdf(file.path(out_dir, "splitsNetwork_plot_pop_colours.pdf"), width = 8, height = 8)
op <- par(mar = c(1, 1, 1, 1))
# plot should accept edge.color / edge.width; use smaller tip labels and no margins
plot(splitsnet,
  edge.color = edge_cols, edge.width = ifelse(is_tip_edge, 2, 1),
  tip.cex = 0.7, cex = 0.8, no.margin = TRUE
)
# add coloured tip labels (in case plot() didn't colour them)
tiplabels(pch = 19, col = tip_cols, cex = 0.8)
# legend for groups
legend("topright", legend = names(group_col_map), col = group_col_map, pch = 19, pt.cex = 1.0, bty = "n", ncol = 1)
par(op)
dev.off()



# neighbourNet coloured by ploidy only (blue = diploid, red = tetraploid)
# match tip labels to metadata rows
tip_idx <- match(nnet$tip.label, tip_df$sampleID)
# take raw ploidy values from metadata (no silent coercion to "unknown")
tip_ploidy_raw <- as.character(tip_df$ploidy[tip_idx])

# basic cleaning: trim and lowercase
tip_ploidy_clean <- tolower(trimws(tip_ploidy_raw))

# detect any missing / empty entries and fail loudly (you said unknown is impossible)
missing_idx <- which(is.na(tip_ploidy_clean) | tip_ploidy_clean == "")
if (length(missing_idx) > 0) {
  stop(
    "Missing ploidy for samples (ploidy must be provided for all tips): ",
    paste0(nnet$tip.label[missing_idx], collapse = ", ")
  )
}

# normalize common representations to canonical labels
tip_ploidy <- tip_ploidy_clean
tip_ploidy[grepl("(^|[^0-9])2($|[^0-9])|^dip", tip_ploidy, perl = TRUE)] <- "diploid"
tip_ploidy[grepl("(^|[^0-9])4($|[^0-9])|tetra", tip_ploidy, perl = TRUE)] <- "tetraploid"

# after normalization, ensure all entries are one of the expected values; fail otherwise
bad_idx <- which(!tip_ploidy %in% c("diploid", "tetraploid"))
if (length(bad_idx) > 0) {
  stop(
    "Unrecognized ploidy values for samples: ",
    paste0(nnet$tip.label[bad_idx], " (raw='", tip_ploidy_raw[bad_idx], "')", collapse = "; ")
  )
}

# colours strictly by ploidy
ploidy_colors <- c(diploid = "blue", tetraploid = "red")
tip_cols <- ploidy_colors[tip_ploidy]

# compute edge colours: terminal edges get the tip colour, internal edges are grey
ntip <- length(nnet$tip.label)
edge_mat <- nnet$edge
is_tip_edge <- edge_mat[, 2] <= ntip
edge_cols <- rep("gray60", nrow(edge_mat))
edge_cols[is_tip_edge] <- tip_cols[edge_mat[is_tip_edge, 2]]

pdf(file.path(out_dir, "neighbournet_plot_ploidy_colours.pdf"), width = 10, height = 8)
op <- par(mar = c(1, 1, 1, 1))
plot(nnet,
  edge.color = edge_cols, edge.width = ifelse(is_tip_edge, 2, 1),
  tip.cex = 0.7, cex = 0.8, no.margin = TRUE
)
tiplabels(pch = 19, col = tip_cols, cex = 0.8)
tiplabels(text = nnet$tip.label, adj = c(0.5, -0.5), frame = "none", cex = 0.65)
legend("topright",
  legend = names(ploidy_colors), col = ploidy_colors, pch = 19,
  pt.cex = 1.0, bty = "n", ncol = 1
)
par(op)
dev.off()




pdf(file.path(out_dir, "splitsNetwork_plot_ploidy_colours.pdf"), width = 8, height = 8)
op <- par(mar = c(1, 1, 1, 1))
plot(splitsnet,
  edge.color = edge_cols, edge.width = ifelse(is_tip_edge, 2, 1),
  tip.cex = 0.7, cex = 0.8, no.margin = TRUE
)
tiplabels(pch = 19, col = tip_cols, cex = 0.8)
tiplabels(text = splitsnet$tip.label, adj = c(0.5, -0.5), frame = "none", cex = 0.65)
legend("topright",
  legend = names(ploidy_colors), col = ploidy_colors, pch = 19,
  pt.cex = 1.0, bty = "n", ncol = 1
)
par(op)
dev.off()
