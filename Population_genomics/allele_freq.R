# load libraries
library(tidyverse)
library(patchwork)
library(vcfR)
library(adegenet)
library(RColorBrewer)

# Set working directory
setwd("/home/bferrette/sailfish_cl/population.analysis/sailfish/seascape/genomic_data")
# Load files
vcf <- read.vcfR("sailfish.filtered.vcf.gz")
pop.data <- read.delim("popdata.tsv", header = TRUE, sep = "\t")
all(colnames(vcf@gt)[-1] == pop.data$sample)
gl <- vcfR2genlight(vcf, n.cores = 24)
ploidy(gl) <- 2
pop(gl) <- pop.data$pop

# Calculate allele frequencies for each site
allele_freqs <- (tab(gl, freq = FALSE, NA.method = "zero"))
allele_freqs
# Keep only the first of the two alleles for each SNP (since p=1-q).
#allele_freqs <- allele_freqs[, seq(1, dim(allele_freqs)[2], 2)]

# Export allele frequencies
write.table(allele_freqs, file = "allele_freqs.tsv", row.names = TRUE, sep = "\t")

