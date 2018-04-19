### HOLD OUT ###

# Load libraries, register cores
library(ranger)
library(matrixStats)
library(tidyverse)
library(doMC)
registerDoMC(20)

# Simulation function
sim <- function(b, n, p) {
  
  # Simulate data
  set.seed(b)
  x <- matrix(runif(n * p), ncol = p,
              dimnames = list(NULL, paste0('x', seq_len(p))))
  y <- rnorm(n)
  df <- data.frame(x, y = y)
  split <- sample.int(n, n / 2)
  df1 <- df[split, ]
  df2 <- df[-split, ]
  rm(x, y, df)
  
  # Grow forests, calculate predictions and loss
  mtry <- floor(sqrt(p))
  B <- 2 * p
  rf1 <- ranger(data = df1, dependent.variable.name = 'y', 
                mtry = mtry, num.trees = B, replace = TRUE,
                num.threads = 20, seed = b)
  rf2 <- ranger(data = df2, dependent.variable.name = 'y', 
                mtry = mtry, num.trees = B, replace = TRUE,
                num.threads = 20, seed = b)
  preds1 <- predict(rf1, data = df2, num.threads = 20, 
                    predict.all = TRUE)$predictions
  preds2 <- predict(rf2, data = df1, num.threads = 20, 
                    predict.all = TRUE)$predictions
  y_hat1 <- rowMeans2(preds1)
  y_hat2 <- rowMeans2(preds2)
  loss1 <- (df2$y - y_hat1)^2
  loss2 <- (df1$y - y_hat2)^2
  loss <- c(loss1, loss2)
  
  # Drop function
  drop <- function(j) {
    f0_idx1 <- !rf1$forest$variable.selected[, j]
    f0_idx2 <- !rf2$forest$variable.selected[, j]
    y_hat01 <- rowMeans2(preds1[, f0_idx1, drop = FALSE])
    y_hat02 <- rowMeans2(preds2[, f0_idx2, drop = FALSE])
    loss01 <- (df2$y - y_hat01)^2
    loss02 <- (df1$y - y_hat02)^2
    loss0 <- c(loss01, loss02)
    t_test <- t.test(loss0, loss, paired = TRUE, alternative = 'greater')
    out <- c(t_test$estimate, t_test$estimate / t_test$statistic,
             t_test$statistic, t_test$p.value)
    return(out)
  }
  
  # Execute in parallel, export
  delta <- foreach(j = seq_len(p), .combine = rbind) %dopar% drop(j)
  dimnames(delta) <- list(NULL, c('CPI', 'SE', 't', 'p.value'))
  delta <- data.frame(
    Feature = paste0('x', seq_len(p)), 
    delta,
    Run = b
  )
  return(delta)
  
}

# Run 10 times
n <- 100
p <- 5000
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

# Very different results across the various runs
(bias_p <- sapply(seq_len(10), function(b) {
  t.test(delta$CPI[delta$Run == b], alternative = 'greater')$p.value
}))

# But overall, the bias is significant
t.test(delta$CPI, alternative = 'greater')

# Also evident in the p-value distribution
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




