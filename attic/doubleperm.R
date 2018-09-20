# Load libraries, register cores
library(mlbench)
library(data.table)
library(doMC)
registerDoMC(12)

# Loop
loop <- function(n, b) {
  
  # Simulate data
  dat <- mlbench.friedman1(3 * n)
  dat <- dat01 <- dat02 <- data.frame(dat$x, y = dat$y)
  
  # Independent permutations
  dat01[, 10] <- dat[sample.int(3 * n), 10]
  dat02[, 10] <- dat[sample.int(3 * n), 10]
  
  # Split into training and test sets
  train <- dat[seq_len(2 * n), ]
  train01 <- dat01[seq_len(2 * n), ]
  train02 <- dat02[seq_len(2 * n), ]
  test <- dat[(2 * n + 1):(3 * n), ]
  test01 <- dat01[(2 * n + 1):(3 * n), ]
  test02 <- dat02[(2 * n + 1):(3 * n), ]
  
  # Train models
  f <- lm(y ~ ., data = train)
  f01 <- lm(y ~ ., data = train01)
  f02 <- lm(y ~ ., data = train02)
  
  # Test models
  y_hat <- predict(f, test)
  y_hat01 <- predict(f01, test01)
  y_hat02 <- predict(f02, test02)
  loss <- (test$y - y_hat)^2
  loss01 <- (test01$y - y_hat01)^2
  loss02 <- (test02$y - y_hat02)^2
  
  # Define and test delta vectors
  delta1 <- loss01 - loss
  delta0 <- loss02 - loss01
  delta <- delta1 - delta0
  t_test <- t.test(delta, alternative = 'greater')
  
  # Export results
  out <- data.table(
    'CPI' = mean(delta),
     'SE' = sd(delta) / sqrt(n),
      't' = t_test$statistic,
'p.value' = t_test$p.value,
     'Run' = b
  )
  return(out)
  
}

# Execute in parallel
out <- foreach(b = seq_len(1e5), .combine = rbind) %dopar% 
  loop(n = 1000, b = b)
saveRDS(out, 'dperm.rds')





