### NULL SIMULATION ###

# Load libraries, register cores
library(tidyverse)
library(doMC)
registerDoMC(4)

# Load CPI
source('cpi.R')

# Simulation function
sim <- function(b, n, p) {
  set.seed(b)
  x <- matrix(runif(n * p), ncol = p,
              dimnames = list(NULL, paste0('x', seq_len(p))))
  y <- rnorm(n)
  mtry <- floor(sqrt(p))
  B <- 2 * p
  rf_split(x, y, mtry = mtry, B = B, seed = b, n.cores = 4) %>%
    mutate(Run = b)
}

# Some hyperparameters
n <- 100
p <- 5000

# Run 10 times
delta <- foreach(b = seq_len(10), .combine = rbind) %do% sim(b, n, p)

# What proportion of p-values are <= 0.05?
p_props <- sapply(seq_len(10), function(b) {
  sum(delta$p.value[delta$Run == b] <= 0.05) / p
})
p_props

# Let's see those t-statistics
ggplot(delta, aes(t)) + 
  geom_histogram(aes(y = ..density..), bins = 100, color = 'black') + 
  stat_function(fun = dt, color = 'red', args = list(df = n - 1)) + 
  theme_bw()

# Pretty obvious bias here...how about if we use cpi_full.R? In other words,
# comparing f0 to the full ensemble, rather than a randomly selected 
# subsample with equal # of trees

# Redefine simulation function
source('cpi_full.R')
sim <- function(b, n, p) {
  set.seed(b)
  x <- matrix(runif(n * p), ncol = p,
              dimnames = list(NULL, paste0('x', seq_len(p))))
  y <- rnorm(n)
  mtry <- floor(sqrt(p))
  B <- 2 * p
  rf_split(x, y, mtry = mtry, B = B, seed = b, n.cores = 4) %>%
    mutate(Run = b)
}
# Run 10 times
delta <- foreach(b = seq_len(10), .combine = rbind) %do% sim(b, n, p)

# What proportion of p-values are <= 0.05?
p_props <- sapply(seq_len(10), function(b) {
  sum(delta$p.value[delta$Run == b] <= 0.05) / p
})
p_props

# Let's see those t-statistics
ggplot(delta, aes(t)) + 
  geom_histogram(aes(y = ..density..), bins = 100, color = 'black') + 
  stat_function(fun = dt, color = 'red', args = list(df = n - 1)) + 
  theme_bw()

# Results from the latter look more stable and less biased
# ...but still biased!



