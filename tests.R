### REAL DATA ###
### 100 runs with each dataset
### Compare cpi with permutation, purity, and party::varimp(conditional = TRUE)

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

eval <- function(dat, x, y, type, high) {
  
  # Prelimz
  p <- ncol(x)
  df <- data.frame(x, y = y)
  mtry <- ifelse(type == 'regression', floor(p/3), floor(sqrt(p)))
  B <- ifelse(high, 1e4L, 1e3L)
  ctrl <- cforest_control(teststat = 'quad', testtype = 'Univ',
                          mtry = mtry, ntree = B, replace = FALSE)
  
  # Brute force algorithm
  loop <- function(b) {
    brute_force(x, y, type = type, weights = FALSE, 
                replace = FALSE, mtry = mtry, B = B,
                n.cores = 1, seed = b)
  }
  foreach(b = seq_len(100), .combine = rbind) %dopar% loop(b) %>%
    mutate(Run = rep(seq_len(100), each = p)) %>%
    saveRDS(paste0('./Results/', dat, '_BruteForce.rds'))
  
  # MDI
  loop <- function(b) {
    f <- ranger(y ~ ., data = df, 
                mtry = mtry, num.trees = B, replace = FALSE, 
                importance = 'impurity', num.threads = 1, seed = b)
    importance(f)
  }
  res <- foreach(b = seq_len(100), .combine = rbind) %dopar% loop(b)
  res %>%
    as_tibble(.) %>%
    gather(key = 'Feature', value = 'MDI') %>%
    mutate(Run = rep(seq_len(100), p)) %>%
    saveRDS(paste0('./Results/', dat, '_MDI.rds'))
  
  # MDA
  loop <- function(b) {
    f <- ranger(y ~ ., data = df, 
                mtry = mtry, num.trees = B, replace = FALSE, 
                importance = 'permutation', num.threads = 1, seed = b)
    importance(f)
  }
  res <- foreach(b = seq_len(100), .combine = rbind) %dopar% loop(b)
  res %>%
    as_tibble(.) %>%
    gather(key = 'Feature', value = 'MDA') %>%
    mutate(Run = rep(seq_len(100), p)) %>%
    saveRDS(paste0('./Results/', dat, '_MDA.rds'))
  
  # Strobl
  loop <- function(b) {
    f <- cforest(medv ~ ., data = Boston, controls = ctrl)
    varimp(f, conditional = TRUE)
  }
  res <- foreach(b = seq_len(100), .combine = rbind) %dopar% loop(b)
  res %>%
    as_tibble(.) %>%
    gather(key = 'Feature', value = 'Strobl') %>%
    mutate(Run = rep(seq_len(100), p)) %>%
    saveRDS(paste0('./Results/', dat, '_Strobl.rds'))
  
  # Timing?  
  mdi <- function(df, mtry, B, par) {
    f <- ranger(y ~ ., data = df, mtry = mtry, num.trees = B, replace = FALSE,
                importance = 'impurity', num.threads = ifelse(par, 4, 1))
    importance(f)
  }
  mda <- function(df, mtry, B, par) {
    f <- ranger(y ~ ., data = df, mtry = mtry, num.trees = B, replace = FALSE,
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
    times = 100L, unit = 's', control = list(warmup = 5)
  )
  data.frame(VIM = times$expr,
            Time = times$time) %>%
    saveRDS(paste0('./Results/', dat, '_times.rds'))
  
}


# Given some model:
dat <- mlbench.friedman1(n, sd = 0)

# We can increase the SNR as follows
snr <- 3
var_res <- var(dat$y) / snr
noisy_y <- dat$y + rnorm(n, sd = sqrt(var_res))



y <- rnorm(10)
noise <- rnorm(10)
k <- sqrt(var(y) / (snr * var(noise)))
noisy_y <- y + k * noise












