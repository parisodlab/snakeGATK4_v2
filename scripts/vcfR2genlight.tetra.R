# Define custom function for tetraploid VCF to genlight
vcfR2genlight.tetra <- function(x, n.cores = 1) {
    bi <- is.biallelic(x)
    if (sum(!bi) > 0) {
        msg <- paste("Found", sum(!bi), "loci with more than two alleles.")
        msg <- c(msg, "\n", paste("Objects of class genlight only support loci with two alleles."))
        msg <- c(msg, "\n", paste(sum(!bi), "loci will be omitted from the genlight object."))
        warning(msg)
        x <- x[bi, ]
    }
    x <- addID(x)
    CHROM <- x@fix[, "CHROM"]
    POS <- x@fix[, "POS"]
    ID <- x@fix[, "ID"]
    x <- extract.gt(x)
    x[x == "0|0"] <- 0
    x[x == "0|1"] <- 1
    x[x == "1|0"] <- 1
    x[x == "1|1"] <- 2
    x[x == "0/0"] <- 0
    x[x == "0/1"] <- 1
    x[x == "1/0"] <- 1
    x[x == "1/1"] <- 2
    x[x == "1/1/1/1"] <- 4
    x[x == "0/1/1/1"] <- 3
    x[x == "0/0/1/1"] <- 2
    x[x == "0/0/0/1"] <- 1
    x[x == "0/0/0/0"] <- 0
    if (requireNamespace("adegenet")) {
        x <- new("genlight", t(x), n.cores = n.cores)
    } else {
        warning("adegenet not installed")
    }
    adegenet::chromosome(x) <- CHROM
    adegenet::position(x) <- POS
    adegenet::locNames(x) <- ID
    return(x)
}
