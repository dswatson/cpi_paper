# Load libraries, register cores
library(tidyverse)
library(doMC)
registerDoMC(20)

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
  cv <- function(k) {
    idx <- flds[[k]]
    test <- df[idx, ]
    train <- df[-idx, ]
    f <- lm(y ~ ., data = train)
    y_hat <- predict(f, test)
    rmse <- sqrt(mean((test$y - y_hat)^2))
    ll <- dnorm(test$y, mean = y_hat, sd = rmse, log = TRUE)
    # Define drop function
    drop <- function(j) {
      test0 <- test[, -j]
      train0 <- train[, -j]
      f0 <- lm(y ~ ., data = train0)
      y_hat0 <- predict(f0, test0)
      rmse0 <- sqrt(mean((test0$y - y_hat0)^2))
      ll0 <- dnorm(test0$y, mean = y_hat0, sd = rmse0, log = TRUE)
      return(ll0)
    }
    out <- foreach(j = seq_len(p), .combine = cbind) %do% drop(j)
    out <- cbind(ll, out)
    dimnames(out) <- list(NULL, c('full', xnames))
    return(out)
  }
  # Compute over folds
  out <- foreach(k = seq_len(k), .combine = rbind) %do% cv(k) 
  as_tibble(out) %>%
    mutate(Run = b) %>%
    return(.)
}
  
# Iterate
out <- foreach(b = seq_len(1e4), .combine = rbind) %dopar% 
  loop(b, k = 10, n = n, p = p)
  
# Calculate likelihood ratios
lr <- function(b, j) {
  tmp <- out %>% filter(Run == b)
  ll <- sum(tmp$full)
  ll0 <- sum(tmp[[j]])
  lambda <- ll - ll0 
  return(lambda)
}

# Save to RDS
foreach(b = seq_len(1e4), .combine = c) %:%
  foreach(j = xnames, .combine = c) %dopar% 
  lr(b, j) %>%
  saveRDS('lambdas.rds')

# And t-tests
t_test <- function(b, j) {
  tmp <- out %>% filter(Run == b)
  delta <- tmp$full - tmp[[j]]
  mu <- mean(delta)
  sigma <- sd(delta) / sqrt(100)
  t_test <- t.test(delta, alternative = 'greater')
  out <- c(mu, sigma, t_test$statistic, t_test$p.value)
  return(out)
}

# Save to RDS
out <- foreach(b = seq_len(1e4), .combine = rbind) %:%
  foreach(j = xnames, .combine = rbind) %dopar%
  t_test(b, j)
dimnames(out) <- list(NULL, c('mu', 'sigma', 't', 'p.value'))
saveRDS(as.data.frame(out), 'll_t.rds')
  
  
  
  