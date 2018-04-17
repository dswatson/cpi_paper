### NULL SIMULATION ###
# We simulate an RNA-seq experiment comparing gene expression between two groups
# of equal size, with an expected log fold change of 0 for each gene.
# We should see t-statistics following t(n - 1) and p-values following U(0, 1)
# ...but we don't :(

# Load libraries, register cores
library(ranger)
library(DESeq2paper)
library(DESeq2)
library(doMC)
registerDoMC(4)

# Load CPI, simulation materials
source('./Scripts/cpi.R')
source('./Scripts/ranger/R/infinitesimalJackknife.R')
source('./DESeq2paper/inst/script/makeSim.R')
load('./DESeq2paper/data/meanDispPairs.RData')

# Hyper-parameters
set.seed(1)
n <- 100
p <- 5000
grp <- factor(rep(c('A', 'B'), each = n / 2))
x <- model.matrix(~ grp)
beta <- rep(0, p)  # This imposes the null hypothesis: 0 LFC between groups

# Simulate counts using joint distribution from Pickrell et al. (2010) study
cnts <- makeSim(p, n, x, beta, meanDispPairs)$mat

# Make DESeqDataSet object for normalization 
dds <- DESeqDataSetFromMatrix(cnts, colData = data.frame(grp), design = ~ grp)
dds <- estimateSizeFactors(dds)
mat <- t(counts(dds, normalized = TRUE))
dimnames(mat) <- list(paste0('s', seq_len(n)), paste0('x', seq_len(p)))

# Remove singleton genes
keep <- rowSums(cnts) > 1
mat <- mat[, keep]

# Preliminaries
p <- ncol(mat)
mtry <- floor(sqrt(p))
B <- 1e4
y <- ifelse(grp == 'A', 0, 1)  # Annoying thing where I have to code responses
                               # as numeric for rf_split(). I'll fix this later

# CPI
# Round 1: Classification forest
cpi <- rf_split(mat, y, type = 'classification', test = 't', weights = FALSE,
                mtry = mtry, B = B, replace = FALSE, 
                n.cores = 4, seed = 1)
# These t-statistics look inflated...maybe we're underestimating the standard
# error by not taking predictive variance into account? 
# We can rerun with precision weights if we set type = 'probability'

# Round 2: Probability forest
cpi_wtd <- rf_split(mat, y, type = 'probability', test = 't', weights = TRUE,
                    mtry = mtry, B = B, replace = FALSE, 
                    n.cores = 4, seed = 1)
# ...except that looks even worse. Does brute_force have the same issues?

# Round 3: Brute force
cpi_bf <- brute_force(mat, y, type = 'classification', test = 't', 
                      weights = FALSE, mtry = mtry, B = B, replace = FALSE,
                      n.cores = 4, seed = 1)
# I haven't run this because it takes a longgg time. Maybe you can try if 
# you have more cores? I suspect it will generate very similar results

# My next idea is to try standardizing the residuals, i.e. dividing each
# by its SE, to see if that improves comparisons...


