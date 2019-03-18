# Set working directory
setwd('~/Documents/CPI/cpi_paper/4.2_Real_Data/4.2.2_Breast_Cancer')

# Load libraries
library(GEOquery)
library(limma)
library(dplyr)

# Import breast cancer data from Herschkowitz et al., 2007
# (This is the data used by Wu & Smyth)
eset <- getGEO('GSE3165')[[3]]

# Extract expression data
x <- exprs(eset)

# Phenotype data
clin <- pData(eset)

# Add basal column
clin <- clin %>%
  mutate(Basal = if_else(grepl('Basal', characteristics_ch2.11), 1, 0))

# Reduce matrix to gene symbols, remove NAs, mean-center, transpose
x <- avereps(x, ID = fData(eset)$GENE_SYMBOL)
x <- na.omit(x[rownames(x) != '', ])
x <- t(x - rowMeans(x))

# Export data
out <- list(x = x, y = clin$Basal)
saveRDS(out, 'GSE165.rds')
