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
(p_props <- sapply(seq_len(10), function(b) {
  sum(delta$p.value[delta$Run == b] <= 0.05) / p
}))

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
(p_props <- sapply(seq_len(10), function(b) {
  sum(delta$p.value[delta$Run == b] <= 0.05) / p
}))

# Let's see those t-statistics
ggplot(delta, aes(t)) + 
  geom_histogram(aes(y = ..density..), bins = 100, color = 'black') + 
  stat_function(fun = dt, color = 'red', args = list(df = n - 1)) + 
  theme_bw()

# Results from the latter look more stable and less biased
# ...but still biased!

# Another helpful visualization is a QQ plot of p-values
data.frame(Observed = -log10(sort(delta$p.value)),
           Expected = -log10(ppoints(nrow(delta)))) %>%
  ggplot(aes(Expected, Observed)) + 
  geom_point(size = 0.5) + 
  geom_abline(intercept = 0, slope = 1, color = 'red') + 
  labs(x = expression('Expected'~-log[10](italic(p))),
       y = expression('Observed'~-log[10](italic(p)))) +
  theme_bw()

# Inflation factor?
chisq <- qchisq(p = 1 - delta$p.value, df = 1)
(lambda <- median(chisq) / qchisq(p = 0.5, df = 1))

