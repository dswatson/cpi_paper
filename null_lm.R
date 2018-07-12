# Load libraries, register cores
library(tidyverse)
library(doMC)
registerDoMC(12)

# Things that'll come up a lot
n <- 1000
p <- 10
xnames <- paste0('x', seq_len(p))

# Define loop function
loop <- function(b, k, n, p) {
  # Simulate data
  x <- matrix(runif(n * p), ncol = p, 
              dimnames = list(NULL, xnames))
  y <- rnorm(n)
  df <- data.frame(x, y)
  # Create folds
  flds <- split(sample.int(n), rep(seq_len(k), length.out = n))
  # Define cv function
  cv <- function(m) {
    idx <- flds[[m]]
    test <- df[idx, ]
    train <- df[-idx, ]
    f <- lm(y ~ ., data = train)
    y_hat <- predict(f, test)
    loss <- (test$y - y_hat)^2
    # Define drop function
    drop <- function(j) {
      test0 <- test[, -j]
      train0 <- train[, -j]
      f0 <- lm(y ~ ., data = train0)
      y_hat0 <- predict(f0, test0)
      loss0 <- (test$y - y_hat0)^2
      return(loss0)
    }
    out <- foreach(j = seq_len(p), .combine = cbind) %do% drop(j)
    out <- cbind(loss, out)
    dimnames(out) <- list(NULL, c('full', xnames))
    return(out)
  }
  # Compute over folds
  out <- foreach(m = seq_len(k), .combine = rbind) %do% cv(m) 
  as_tibble(out) %>%
    mutate(Run = b) %>%
    return(.)
}

# Iterate
out <- foreach(b = seq_len(1e4), .combine = rbind) %dopar% 
  loop(b, k = 10, n = n, p = p)

# Run t-tests
t_test <- function(b, j) {
  tmp <- out %>% filter(Run == b)
  delta <- tmp[[j]] - tmp$full
  cpi <- mean(delta)
  sigma <- sd(delta) / sqrt(n)
  t_test <- t.test(delta, alternative = 'greater')
  out <- c(cpi, sigma, t_test$statistic, t_test$p.value)
  return(out)
}

# Save to RDS
out <- foreach(b = seq_len(1e4), .combine = rbind) %:%
  foreach(j = xnames, .combine = rbind) %dopar%
  t_test(b, j)
dimnames(out) <- list(NULL, c('cpi', 'se', 't', 'p.value'))
out <- as.data.frame(out)
saveRDS(out, 'null_lm.rds')
