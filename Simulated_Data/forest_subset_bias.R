### IS BIAS DUE TO FOREST SUBSETTING OR *SELECTIVE* FOREST SUBSETTING? ###

# Load libraries, register cores
library(ranger)
library(matrixStats)
library(tidyverse)
library(doMC)
registerDoMC(4)

# Simulation function
sim <- function(b, n, p) {
  
  # Simulate data
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
               keep.inbag = TRUE, num.threads = 4, seed = b)
  
  # Create oob_preds object
  oob_idx <- ifelse(simplify2array(rf$inbag.counts) == 0, TRUE, NA)
  preds <- predict(rf, df, num.threads = 4, predict.all = TRUE)$predictions
  oob_preds <- oob_idx * preds
  
  # Extract predictions, loss
  y_hat <- rf$predictions
  emp_risk <- mean((y_hat - y)^2)
  
  # Proportion of trees expected in sub-forest?
  splits <- sapply(seq_len(B), function(b) {
    sum(!treeInfo(rf, b)$terminal)
  })
  pr_null <- (1 - mtry/p)^mean(splits)
  
  # Random sub-forests
  loop <- function(j) {
    set.seed(j)
    f0_idx <- sample(c(TRUE, FALSE), size = B, replace = TRUE, 
                     prob = c(pr_null, 1 - pr_null))
    y_hat0 <- rowMeans2(oob_preds[, f0_idx, drop = FALSE], na.rm = TRUE)
    emp_risk0 <- mean((y_hat0 - y)^2)
    return(emp_risk0 - emp_risk)
  }
  out <- data.frame(
    Random = foreach(j = seq_len(p), .combine = c) %dopar% loop(j)
  )
  
  # Drop-one sub-forests
  loop <- function(j) {
    f0_idx <- !rf$forest$variable.selected[, j]
    y_hat0 <- rowMeans2(oob_preds[, f0_idx, drop = FALSE], na.rm = TRUE)
    emp_risk0 <- mean((y_hat0 - y)^2)
    return(emp_risk0 - emp_risk)
  }
  out$Drop <- foreach(j = seq_len(p), .combine = c) %dopar% loop(j)
  
  # Export
  out$Run <- b
  return(out)
  
}

# Run 10 iterations
n <- 100
p <- 5000
out <- foreach(b = seq_len(10), .combine = rbind) %do% sim(b, n, p)

# Test each run separately
(pvals <- sapply(seq_len(10), function(b) {
  with(filter(out, Run == b), 
       t.test(Drop, Random, alternative = 'greater')$p.value)
}))

# Global test
with(out, t.test(Drop, Random, alternative = 'greater'))

### Though the drop method appears to have *slightly* greater empirical risk
### on average, the difference is not significant. This suggests that the 
### bias evident in null.R is due to forest subsetting in general, not 
### anything specific to the drop-one method.

# Boxplots
box <- function(b) {
  out %>%
    filter(Run == b) %>% 
    select(-Run) %>%
    gather(key = 'Method', value = 'Emp_Risk') %>%
    ggplot(aes(Method, Emp_Risk, fill = Method)) + 
    geom_boxplot() + 
    theme_bw()
}
box(1)
box(2) # etc.

# Overall boxplot
out %>%
  select(-Run) %>%
  gather(key = 'Method', value = 'Emp_Risk') %>%
  ggplot(aes(Method, Emp_Risk, fill = Method)) + 
  geom_boxplot() + 
  theme_bw()



