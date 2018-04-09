### REAL DATA ###
### 100 runs with each dataset
### Compare cpi with permutation, purity, and party::varimp(conditional = TRUE)
### Test: internal stability, mutual consistency, and processing time
### For largest differences in ranking, compare with correlation matrix
### Also worth running Altmann permutations to see p-value differences?

# Load libraries, register cores
library(mlbech)
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
  
  # Prelimz
  p <- ncol(x)
  df <- ifelse(type == 'regression', data.frame(x, y = y), 
               data.frame(x, y = as.factor(y)))
  mtry <- ifelse(type == 'regression', floor(p/3), floor(sqrt(p)))
  B <- ifelse(high, 1e4L, 1e3L)
  ctrl <- cforest_control(teststat = 'quad', testtype = 'Univ',
                          mtry = mtry, ntree = B, replace = FALSE)
  
  # Brute force algorithm
  loop <- function(b) {
    brute_force(x, y, type = type, weights = FALSE, 
                replace = FALSE, mtry = mtry, B = B,
                n.cores = 1, seed = b) %>%
      rename(VI = Delta) %>%
      mutate(Rank = min_rank(desc(VI)), VIM = 'CPI') %>%
      select(Feature, VIM, VI, Rank)
  }
  cpi_out <- foreach(b = seq_len(runs), .combine = rbind) %dopar% loop(b) 
  
  # MDI
  loop <- function(b) {
    f <- ranger(data = df, dependent.variable.name = 'y',
                mtry = mtry, num.trees = B, replace = FALSE, 
                importance = 'impurity', num.threads = 1, seed = b)
    imp <- importance(rf)
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
    imp <- importance(rf)
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
  
  # Timing?  
  mdi <- function(df, mtry, B, par) {
    f <- ranger(data = df, dependent.variable.name = 'y', 
                mtry = mtry, num.trees = B, replace = FALSE,
                importance = 'impurity', num.threads = ifelse(par, 4, 1))
    importance(f)
  }
  mda <- function(df, mtry, B, par) {
    f <- ranger(data = df, dependent.variable.name = 'y',
                mtry = mtry, num.trees = B, replace = FALSE,
                importance = 'permutation', num.threads = ifelse(par, 4, 1))
    importance(f)
  }
  strobl <- function(df, mtry, B) {  
    f <- cforest(y ~ ., data = df, controls = ctrl)
    varimp(f, conditional = TRUE)
  }
  times <- microbenchmark(
    bf_par = brute_force(x, y, type = type, weights = FALSE,
                         replace = FALSE, mtry = mtry, B = B, 
                         n.cores = 4, seed = NULL),
    bf_ser = brute_force(x, y, type = type, weights = FALSE,
                         replace = FALSE, mtry = mtry, B = B, 
                         n.cores = 1, seed = NULL),
    mdi_par = mdi(df, mtry, B, par = TRUE),
    mdi_ser = mdi(df, mtry, B, par = FALSE),
    mda_par = mda(df, mtry, B, par = TRUE),
    mda_ser = mda(df, mtry, B, par = FALSE),
    strobl - strobl(df, mtry, B),
    times = runs, unit = 's', control = list(warmup = 5)
  )
  data.frame(VIM = times$expr,
            Time = times$time) %>%
    saveRDS(paste0('./Results/', dat, '_times.rds'))
  
}



# Test internal consistency
vim <- readRDS('./Results/BostonHousing_vim.rds')
cpi_res <- vim %>%
  filter(VIM == 'CPI') %>%
  mutatee(Run = rep(seq_len(100), each = p))
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









