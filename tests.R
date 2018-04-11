### REAL DATA ###
### 100 runs with each dataset
### Compare cpi with permutation, purity, and party::varimp(conditional = TRUE)
### Test: internal stability, mutual consistency, and processing time
### For largest differences in ranking, compare with correlation matrix
### Also worth running Altmann permutations to see p-value differences?

# Load libraries, register cores
library(mlbench)
library(ranger)
library(party)
library(matrixStats)
library(microbenchmark)
library(tidyverse)
library(doMC)
registerDoMC(4)

# Load cpi
source('cpi.R')

eval <- function(dat, x, y, type, high, runs) {
  
  # Preliminaries
  p <- ncol(x)
  if (type %in% c('regression', 'survival')) {
    df <- data.frame(x, y = y)
    mtry <- floor(p / 3)
  } else {
    df <- data.frame(x, y = as.factor(y))
    mtry <- floor(sqrt(p))
  }
  B <- ifelse(high, 1e4L, 1e3L)
  ctrl <- cforest_unbiased(mtry = mtry, ntree = B)
  
  # Brute force algorithm
  loop <- function(b) {
    imp <- brute_force(x, y, type = type, test = NULL, weights = FALSE, 
                       replace = FALSE, mtry = mtry, B = B,
                       n.cores = 1, seed = b)
    data.frame(Feature = names(imp),
                   VIM = 'CPI',
                    VI = as.numeric(imp)) %>%
      mutate(Rank = min_rank(desc(VI))) 
  }
  cpi_out <- foreach(b = seq_len(runs), .combine = rbind) %dopar% loop(b) 
  
  # MDI
  loop <- function(b) {
    f <- ranger(data = df, dependent.variable.name = 'y',
                mtry = mtry, num.trees = B, replace = FALSE, 
                importance = 'impurity', num.threads = 1, seed = b)
    imp <- importance(f)
    data.frame(Feature = names(imp), 
                   VIM = 'MDI',
                    VI = as.numeric(imp)) %>%
      mutate(Rank = min_rank(desc(VI)))
  }
  mdi_out <- foreach(b = seq_len(runs), .combine = rbind) %dopar% loop(b)
  
  # MDA
  loop <- function(b) {
    f <- ranger(data = df, dependent.variable.name = 'y',
                mtry = mtry, num.trees = B, replace = FALSE, 
                importance = 'permutation', num.threads = 1, seed = b)
    imp <- importance(f)
    data.frame(Feature = names(imp), 
                   VIM = 'MDA',
                    VI = as.numeric(imp)) %>%
      mutate(Rank = min_rank(desc(VI)))
  }
  mda_out <- foreach(b = seq_len(runs), .combine = rbind) %dopar% loop(b) 
  
  # Strobl
  loop <- function(b) {
    set.seed(b)
    f <- cforest(y ~ ., data = df, controls = ctrl)
    imp <- varimp(f, conditional = TRUE)
    data.frame(Feature = names(imp), 
                   VIM = 'MDA-C',
                    VI = as.numeric(imp)) %>%
      mutate(Rank = min_rank(desc(VI)))
  }
  strobl_out <- foreach(b = seq_len(runs), .combine = rbind) %dopar% loop(b)
  
  # Altogether now
  rbind(cpi_out, mdi_out, mda_out, strobl_out) %>%
    saveRDS(paste0('./Results/', dat, '_vim.rds'))
  
  
  
}


# Test internal consistency
vim <- readRDS('./Results/BostonHousing_vim.rds')
cpi_res <- vim %>%
  filter(VIM == 'CPI') %>%
  mutate(Run = rep(seq_len(100), each = p))
mat <- matrix(nrow = 100, ncol = 100,
              dimnames = list(paste0('r', seq_len(100)), 
                              paste0('r', seq_len(100))))
# Rank consistency
for (i in 2:100) {
  for (j in 1:(i - 1)) {
    mat[i, j] <- cor(cpi_res$Rank[cpi_res$Run == i],
                     cpi_res$Rank[cpi_res$Run == j])
  }
}
mean(mat, na.rm = TRUE)

# VIM consistency
for (i in 2:100) {
  for (j in 1:(i - 1)) {
    mat[i, j] <- cor(cpi_res$CPI[cpi_res$Run == i],
                     cpi_res$CPI[cpi_res$Run == j])
  }
}
mean(mat, na.rm = TRUE)

# External consistency
vim_summary <- vim %>%
  group_by(VIM, Feature) %>%
  mutate(Mean_VI = mean(VI),
           SD_VI = sd(VI),
       Mean_Rank = mean(Rank),
         SD_Rank = sd(Rank)) %>%
  ungroup(.) %>%
  select(Feature, VIM, Mean_VI, SD_VI, Mean_Rank, SD_Rank) %>%
  unique(.)
cor(vim_summary$Mean_Rank[vim_summary$VIM == 'CPI'],
    vim_summary$Mean_Rank[vim_summary$VIM == 'MDI'])









