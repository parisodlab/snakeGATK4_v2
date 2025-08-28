library(SNPRelate)
library(argparse)


# -------- INPUT VCF FILE --------
# get as parameters
vcf_file <- commandArgs(trailingOnly = TRUE)[1]
gds_file <- commandArgs(trailingOnly = TRUE)[2]
outdir <- commandArgs(trailingOnly = TRUE)[3]

# -------- CONVERT VCF TO GDS FORMAT --------
snpgdsVCF2GDS(
    vcf.fn = vcf_file,
    out.fn = gds_file,
    method = "biallelic.only"
)

# -------- OPEN GDS FILE --------
genofile <- snpgdsOpen(gds_file)

# -------- LD PRUNING --------
# Adjust LD threshold as needed (e.g., 0.2 = r² ≤ 0.2)
ld_thresh <- 0.2

# Perform pruning (this selects one SNP from each set of highly correlated SNPs)
snpset <- snpgdsLDpruning(genofile,
    method = "corr",
    ld.threshold = ld_thresh,
    autosome.only = FALSE,
    slide.max.bp = 50000, # optional window size
    maf = NaN, # disables MAF filter
    missing.rate = NaN
) # disables missing data filter

# -------- GET PRUNED SNP IDs --------
pruned_ids <- unlist(snpset)

# Get full SNP metadata
snp_ids <- read.gdsn(index.gdsn(genofile, "snp.id"))
snp_chrom <- read.gdsn(index.gdsn(genofile, "snp.chromosome"))
snp_pos <- read.gdsn(index.gdsn(genofile, "snp.position"))

# Subset metadata to pruned SNPs
keep_idx <- match(pruned_ids, snp_ids)
chroms <- snp_chrom[keep_idx]
positions <- snp_pos[keep_idx]

# Combine chromosome and position into GATK interval format (CHR:POS)
gatk_intervals <- paste0(chroms, ":", positions)

# Save to a file for GATK
writeLines(gatk_intervals, paste0(outdir, "/pruned_snps.intervals"))
# Close the GDS file
snpgdsClose(genofile)
