### SIMULATED DATA ###
### Evaluate error rates and nominal coverage
### Compare with and without precision weights
### Test behavior on correlated predictors
### Measure AUROC and AUPR vs. different VIMs

# Load libraries, register cores
library(mlbench)
library(ranger)
library(party)
library(tidyverse)
library(doMC)
registerDoMC(4)

# Load cpi
source('cpi.R')

# Competitor functions
mdi <- function(df, mtry, B, b) {
  f <- ranger(data = df, dependent.variable.name = 'y',
              mtry = mtry, num.trees = B, replace = FALSE, 
              importance = 'impurity', num.threads = 1, seed = b)
  importance(f)
}
mda <- function(df, mtry, B, b) {
  f <- ranger(data = df, dependent.variable.name = 'y',
              mtry = mtry, num.trees = B, replace = FALSE, 
              importance = 'permutation', num.threads = 1, seed = b)
  importance(f)
}
mda_c <- function(df, mtry, B, b, ctrl) {  
  set.seed(b)
  f <- cforest(y ~ ., data = df, controls = ctrl)
  varimp(f, conditional = TRUE)
}

# Friedman1 benchmark
n <- 100
p <- 10

# CPI loop
loop <- function(b) {
  
  # Preliminaries
  set.seed(b)
  dat <- mlbench.friedman1(n = n)
  colnames(dat$x) <- paste0('x', seq_len(p))
  mtry <- floor(p / 3)
  B <- 1000
  ctrl <- cforest_unbiased(mtry = mtry, ntree = B)
  
  # t-test
  cpi_t <- brute_force(dat$x, dat$y, test = NULL, weights = FALSE,
                       replace = FALSE, mtry = mtry, B = B, 
                       n.cores = 1, seed = b) 
  
  # Weighted t-test
  cpi_wtd <- brute_force(dat$x, dat$y, test = 't', weights = TRUE,
                         replace = FALSE, mtry = mtry, B = B, 
                         n.cores = 1, seed = b) 
  cpi_wtd <- cpi_wtd$CPI
  
  # MDI
  mdi_out <- mdi(df, mtry, B, b)
  
  # MDA
  mda_out <- mda(df, mtry, B, b)
  
  # MDA-C
  mda_c_out <- mda_c(df, mtry, B, b, ctrl)
  
  # Altogether now
  data.frame(
    Feature = colnames(dat$x),
    CPI = cpi_t,
    CPI_wtd = CPI_wtd,
    MDI = mdi_out,
    MDA = mda_out,
    MDA_C = mda_c_out,
    Run = b
  )
}
sims <- foreach(b = seq_len(100), .combine = rbind) %dopar% loop(b) 
saveRDS(sims, './Results/Simulations/Friedman1.rds')

# RNA-seq simulation (from Love et al., 2014)
library(DESeq2paper)
library(DESeq2)
source('./DESeq2paper/inst/script/makeSim.R')
load('./DESeq2paper/data/meanDispPairs.RData')

# Hyper-parameters
n <- 100
p <- 1e4
lfc <- log2(2)
grp <- factor(rep(c('A', 'B'), each = n / 2))
x <- model.matrix(~ grp)

loop <- function(b) {
  
  set.seed(b)

  # Simulate counts
  beta <- c(rep(0, p * 8/10), 
            sample(c(-lfc, lfc), p * 2/10, replace = TRUE))
  de_idx <- abs(beta) > 0
  cnts <- makeSim(p, n, x, beta, meanDispPairs)$mat
  
  # Make DESeqDataSet object for normalization 
  dds <- DESeqDataSetFromMatrix(cnts, colData = data.frame(grp), design = ~ grp)
  dds <- estimateSizeFactors(dds)
  mat <- t(counts(dds, normalized = TRUE))
  dimnames(mat) <- list(paste0('s', seq_len(n)), paste0('x', seq_len(p)))
  
  # Remove singleton genes
  keep <- rowSums(cnts) > 1
  mat <- mat[, keep]
  de_idx <- de_idx[keep]
  
  # Preliminaries
  df <- data.frame(mat, y = grp)
  p <- ncol(mat)
  mtry <- floor(sqrt(p))
  B <- 2e4
  
  # CPI
  y <- ifelse(grp == 'A', 0, 1)  # Unclear why this isn't the other way around
  cpi <- rf_split(mat, y, type = 'classification', test = NULL,
                  mtry = mtry, B = B, replace = FALSE, 
                  n.cores = 1, seed = b)
  
  # MDI
  mdi_out <- mdi(df, mtry, B, b)
  
  # MDA
  mda_out <- mda(df, mtry, B, b)
  
  # Altogether now
  data.frame(
    Feature = colnames(mat),
    CPI = cpi,
    MDI = mdi_out, 
    MDA = mda_out,
    DE = de_idx,
    Run = b
  )
  
}
sims <- foreach(b = seq_len(10), .combine = rbind) %dopar% loop(b) 
saveRDS(sims, './Results/Simulations/RNA_seq.rds')

# Strobl et al. simulation
sigma <- matrix(0, nrow = 12, ncol = 12)
sigma[1:4, 1:4] <- 0.9
diag(sigma) <- 1
x <- mvrnorm(n = 100 * 12, mu = rep(0, 12), Sigma = sigma)
colnames(x) <- paste0('x', 1:12)
beta <- c(5, 5, 2, 0, -5, -5, -2, rep(0, times = 5))
y <- drop(tcrossprod(beta, x)) + rnorm(100, mean = 0, sd = 0.5)




