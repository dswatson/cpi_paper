# Load libraries
library(mlbench)
library(MASS)
library(relaimpo)
library(tidyverse)

# Consider a plot like Fig. 3 from Strobl et al. (2008):
# x-axis is feature ID, y-axis is VI, and we see different color curves 
# (with SE) for different measures of VI

# Probably worth writing a stripped down version of c_imp that just spits out
# deltas, then looping it with different seeds like 1k times to get a spread
# Will also need 1k runs of gini, perm, and strobl_cond importance

# Friedman
df <- as.data.frame(mlbench.friedman1(n = 100))
out <- c_imp(df[, 1:10], df$y)

# Strobl et al.
sigma <- matrix(0, nrow = 12, ncol = 12)
sigma[1:4, 1:4] <- 0.9
diag(sigma) <- 1
x <- mvrnorm(n = 100 * 12, mu = rep(0, 12), Sigma = sigma)
colnames(x) <- paste0('x', 1:12)
beta <- c(5, 5, 2, 0, -5, -5, -2, 0, 0, 0, 0, 0)
y <- drop(tcrossprod(beta, x)) + rnorm(100, mean = 0, sd = 0.5)
out <- c_imp(x, y)

# Compare with LMG method (ground truth for linear model?)
m <- lm(y ~ x)
lmg <- calc.relimp(m, type = 'lmg')