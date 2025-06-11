install.packages("LEA")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("LEA")

library(LEA)

setwd("/Users/samjohnson/Desktop/")
vcf2geno("hard_variants_bial_noindels_q20_mindep8_maxdep75_maf30_miss95.recode.vcf", 
         output.file = "hard_variants_bial_noindels_q20_mindep8_maxdep75_maf30_miss95.geno")

geno_lines <- readLines("hard_variants_bial_noindels_q20_mindep8_maxdep75_maf30_miss95.geno")
geno_split <- strsplit(geno_lines, split = "")
geno_mat <- do.call(rbind, lapply(geno_split, as.integer))
geno_mat <- t(geno_mat)
str(geno_mat)
table(geno_mat)
dim(geno_mat)

if (!requireNamespace("VariantAnnotation", quietly = TRUE)) {
  BiocManager::install("VariantAnnotation")
}
library(VariantAnnotation)

vcf_path <- "hard_variants_bial_noindels_q20_mindep8_maxdep75_maf30_miss95.recode.vcf"
vcf <- readVcf(vcf_path)

sample_ids <- samples(header(vcf)) # length == 1184

rownames(geno_mat) <- sample_ids

geno_mat[geno_mat == 9] <- -9
table(geno_mat)
colnames(geno_mat) <- paste0("SNP", 1:ncol(geno_mat))

from_012 <- read.table(file = "hard_variants_bial_noindels_q20_mindep8_maxdep75_maf30_miss95.012_conv", 
                       header = FALSE, sep = "\t")
from_012 <- from_012[, -1]
ind <- read.table(file = "hard_variants_bial_noindels_q20_mindep8_maxdep75_maf30_miss95.012.indv", 
                  header = FALSE)
rownames(from_012) <- ind$V1
colnames(from_012) <- paste0("SNP", 1:ncol(from_012))

compare_geno_matrices <- function(mat1, mat2) {
  # Check if dimensions match
  if (!all(dim(mat1) == dim(mat2))) {
    cat("Dimensions differ:\n")
    cat("Matrix 1:", dim(mat1), "\n")
    cat("Matrix 2:", dim(mat2), "\n")
    return(FALSE)
  }
  
  # Check if matrices are identical
  if (all(mat1 == mat2, na.rm = TRUE)) {
    cat("✅ Matrices are identical in all non-NA cells.\n")
    return(TRUE)
  } else {
    cat("❌ Matrices differ.\n")
    
    # Find mismatching indices
    mismatches <- which(mat1 != mat2, arr.ind = TRUE)
    cat("Number of differing cells:", nrow(mismatches), "\n")
    cat("First few mismatches (row, col):\n")
    print(head(mismatches))
    
    # Optionally inspect values at first mismatch
    i <- mismatches[1,1]
    j <- mismatches[1,2]
    cat("Example difference at [", i, ",", j, "]: ",
        "mat1 =", mat1[i, j], ", mat2 =", mat2[i, j], "\n")
    
    return(FALSE)
  }
}

comp <- compare_geno_matrices(geno_mat, from_012)
# well... that's what you want i suppose!


