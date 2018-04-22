### CPI vs. collection of random forest subsets ###

# Load libraries, register cores
library(ranger)
library(matrixStats)
library(tidyverse)
library(doMC)
registerDoMC(20)

# Simulation function
sim <- function(b, n, p) {
  
  set.seed(b)
  x <- matrix(runif(n * p), ncol = p,
              dimnames = list(NULL, paste0('x', seq_len(p))))
  y <- rnorm(n)
  
  # Grow forest
  mtry <- floor(sqrt(p))
  B <- 2 * p
  df <- data.frame(x, y = y)
  rf <- ranger(data = df, dependent.variable.name = 'y', 
               mtry = mtry, num.trees = B, replace = FALSE,
               keep.inbag = TRUE, num.threads = 20, seed = b)
  
  # Create oob_preds object
  oob_idx <- ifelse(simplify2array(rf$inbag.counts) == 0, TRUE, NA)
  preds <- predict(rf, df, num.threads = 20, predict.all = TRUE)$predictions
  oob_preds <- oob_idx * preds
  
  # How many trees expected in sub-forest?
  splits <- foreach(b = seq_len(B), .combine = c) %dopar% 
    sum(!treeInfo(rf, b)$terminal)
  exp_B0 <- round(B * (1 - mtry/p)^mean(splits))
  
  # Random sub-forests
  loop <- function(j) {
    set.seed(j)
    sub_idx <- sample.int(B, exp_B0)
    y_hat <- rowMeans2(oob_preds[, sub_idx, drop = FALSE], na.rm = TRUE)
    loss <- (y_hat - y)^2
    return(loss)
  }
  loss <- rowMeans2(foreach(j = seq_len(100), .combine = cbind) %dopar% loop(j))
  
  # Null sub-forests
  drop <- function(j) {
    f0_idx <- !rf$forest$variable.selected[, j]
    y_hat0 <- rowMeans2(oob_preds[, f0_idx, drop = FALSE], na.rm = TRUE)
    loss0 <- (y_hat0 - y)^2
    t_test <- t.test(loss0, loss, paired = TRUE, alternative = 'greater')
    out <- c(t_test$estimate, t_test$estimate / t_test$statistic,
             t_test$statistic, t_test$p.value)
  }
  delta <- foreach(j = seq_len(p), .combine = rbind) %dopar% drop(j)
  dimnames(delta) <- list(NULL, c('CPI', 'SE', 't', 'p.value'))
  delta <- data.frame(Feature = colnames(x), delta)
  return(delta)
  
}

# Simulate data
n <- 100
p <- 5000
out <- foreach(b = seq_len(100), .combine = rbind) %do% sim(b, n, p)
saveRDS(out, 'cpi_rdm_subset.rds')

# This looks pretty darn good to me
ggplot(out, aes(t)) + 
  geom_histogram(aes(y = ..density..), bins = 100, color = 'black') + 
  stat_function(fun = dt, color = 'red', args = list(df = 99)) + 
  theme_bw()

# Although a QQ plot suggests there's a slight mismatch at the tails
n_out <- nrow(out)
data_frame(Observed = quantile(out$t, probs = ppoints(n_out)),
           Expected = qt(ppoints(n_out), df = 99)) %>%
  ggplot(aes(Expected, Observed)) + 
  geom_point(size = 0.25) + 
  geom_abline(intercept = 0, slope = 1, color = 'red') + 
  labs(x = 'Expected Quantiles', y = 'Observed Quantiles') + 
  theme_bw()

# Annoyingly, a t-test still indicates a slight upward bias in CPI
t.test(out$CPI, alternative = 'greater')

# But a KS-test fails to reject the null hypothesis that our observed t-stats
# were sampled from the expected distribution
ks.test(out$t, 'pt', df = 99, alternative = 'greater')

# What proportion of p-values are <= 0.05?
sum(out$p.value <= 0.05) / n_out

# Not bad. In fact, the p-values seem pretty uniform
ggplot(out, aes(p.value)) + 
  geom_histogram(aes(y = ..density..), bins = 100, color = 'black') + 
  stat_function(dunif, color = 'red') + 
  theme_bw()

# QQ plot is probably more informative here
chisq <- qchisq(p = 1 - dat$p.value, df = 1)
lambda_val <- median(chisq) / qchisq(p = 0.5, df = 1)
lambda_lbl <- paste('lambda ==',  round(lambda_val, 2))
df <- data_frame(Observed = -log10(sort(dat$p.value)),
                 Expected = -log10(ppoints(nrow(dat))))
ggplot(df, aes(Expected, Observed)) + 
  geom_point(size = 0.25) + 
  geom_abline(intercept = 0, slope = 1, color = 'red') + 
  labs(x = expression('Expected'~-log[10](italic(p))),
       y = expression('Observed'~-log[10](italic(p)))) +
  annotate('text', x = max(df$Expected), y = 0, size = 5,
           hjust = 1, label = lambda_lbl, parse = TRUE) +
  theme_bw()
  
# I can live with that inflation factor, personally...
