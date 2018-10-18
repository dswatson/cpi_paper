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
saveRDS(out, 'lm_test.rds')

mu <- out %>%
  group_by(Run) %>%
  summarise(full = mean(full),
            x1 = mean(x1),
            x2 = mean(x2),
            x3 = mean(x3),
            x4 = mean(x4),
            x5 = mean(x5),
            x6 = mean(x6),
            x7 = mean(x7),
            x8 = mean(x8),
            x9 = mean(x9),
            x10 = mean(x10))
sigma <- out %>%
  group_by(Run) %>%
  summarise(full = sd(full),
            x1 = sd(x1),
            x2 = sd(x2),
            x3 = sd(x3),
            x4 = sd(x4),
            x5 = sd(x5),
            x6 = sd(x6),
            x7 = sd(x7),
            x8 = sd(x8),
            x9 = sd(x9),
            x10 = sd(x10))
rss <- out %>%
  group_by(Run) %>%
  summarise(full = sum(full),
            x1 = sum(x1),
            x2 = sum(x2),
            x3 = sum(x3),
            x4 = sum(x4),
            x5 = sum(x5),
            x6 = sum(x6),
            x7 = sum(x7),
            x8 = sum(x8),
            x9 = sum(x9),
            x10 = sum(x10))


# Run t-tests
t_test <- function(b, j) {
  tmp <- out %>% filter(Run == b)
  delta <- tmp[[j]] - tmp$full
  cpi <- mean(delta)
  se <- sd(delta) / sqrt(n)
  t_test <- t.test(delta, alternative = 'greater')
  out <- c(cpi, se, t_test$statistic, t_test$p.value)
  return(out)
}

# Save to RDS
out <- foreach(b = seq_len(1e4), .combine = rbind) %:%
  foreach(j = xnames, .combine = rbind) %dopar%
  t_test(b, j)
dimnames(out) <- list(NULL, c('cpi', 'se', 't', 'p.value'))
out <- as.data.frame(out)
saveRDS(out, 'null_lm.rds')
