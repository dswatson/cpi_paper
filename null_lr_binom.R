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
  x <- matrix(rnorm(n * p), ncol = p, 
              dimnames = list(NULL, xnames))
  y <- rbinom(n, size = 1, prob = 0.5)
  df <- data.frame(x, y)
  # Create folds
  flds <- split(sample.int(n), rep(seq_len(k), length.out = n))
  # Define cv function
  cv <- function(k) {
    idx <- flds[[k]]
    test <- df[idx, ]
    train <- df[-idx, ]
    f <- glm(y ~ ., data = train, family = 'binomial')
    y_hat <- predict(f, test, type = 'response')
    ll <- dbinom(test$y, size = 1, prob = y_hat, log = TRUE)
    # Define drop function
    drop <- function(j) {
      test0 <- test[, -j]
      train0 <- train[, -j]
      f0 <- glm(y ~ ., data = train0, family = 'binomial')
      y_hat0 <- predict(f0, test0, type = 'response')
      ll0 <- dbinom(test0$y, size = 1, prob = y_hat0, log = TRUE)
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
  saveRDS('lambdas_binom.rds')




