### Experiments in inference for delta vectors ### 

# Load libraries, register cores
library(mlbench)
library(weights)
library(tidyverse)
library(ggsci)
library(doMC)
registerDoMC(12)

# PART I
# As a baseline, run un-weighted t-tests on log ratios of loss vectors

# Loop
loop <- function(n, b) {
  
  # Simulate data
  dat <- mlbench.friedman1(3 * n)
  dat <- dat0 <- data.frame(dat$x, y = dat$y)
  dat0[, 10] <- dat[sample.int(3 * n), 10]
  
  # Split into training and test sets
  idx <- sample.int(3 * n, size = 2 * n)
  train <- dat[idx, ]
  train0 <- dat0[idx, ]
  test <- dat[-idx, ]
  test0 <- dat0[-idx, ]
  
  # Train models
  f <- lm(y ~ ., data = train)
  f0 <- lm(y ~ ., data = train0)
  
  # Test models
  y_hat <- predict(f, test)
  y_hat0 <- predict(f0, test0)
  loss <- (test$y - y_hat)^2
  loss0 <- (test0$y - y_hat0)^2
  
  # Define and test delta vector
  delta <- log(loss / loss0)
  t_test <- t.test(delta, alternative = 'greater')
  
  # Calculate empirical coverage at alpha = 0.05
  se <- sd(delta) / sqrt(n)
  cvg <- 1 - (sum(loss >= loss0 + qnorm(0.975) * se) + 
              sum(loss <= loss0 - qnorm(0.975) * se)) / n
  
  # Export results
  out <- data.frame(
         'CPI' = mean(delta),
          'SE' = se,
           't' = t_test$statistic,
     'p.value' = t_test$p.value,
    'Coverage' = cvg,
    'Variance_Model' = 'Homoskedastic'
  )
  return(out)
  
}

# Execute in parallel
out1 <- foreach(b = seq_len(1e5), .combine = rbind) %dopar% 
  loop(n = 1000, b = b)
saveRDS(out1, 'unwtd.rds')


# PART II
# We know that a loss vs. loss0 scatterplot is heteroskedastic under the null
# So here we impose a particular variance model and run weighted t-tests

# Simulate data
n <- 1e4
dat <- mlbench.friedman1(3 * n)
dat <- dat0 <- data.frame(dat$x, y = dat$y)
dat0[, 10] <- dat[sample.int(3 * n), 10]

# Split into training and test sets
idx <- sample.int(3 * n, size = 2 * n)
train <- dat[idx, ]
train0 <- dat0[idx, ]
test <- dat[-idx, ]
test0 <- dat0[-idx, ]

# Train, test models
f <- lm(y ~ ., data = train)
f0 <- lm(y ~ ., data = train0)
y_hat <- predict(f, test)
y_hat0 <- predict(f0, test0)
loss <- (test$y - y_hat)^2
loss0 <- (test0$y - y_hat0)^2

# Model variance with a LOWESS curve
delta <- log(loss / loss0)
l <- lowess(x = log(loss0), y = sqrt(abs(delta)))
sigma_mod <- approxfun(l, rule = 2)

# Now we use this sigma_mod to simulate correlated loss vectors
# with predetermined heteroskedasticity

# Loop
loop <- function(n, b) {
  
  # Simulate data
  dat <- mlbench.friedman1(3 * n)
  dat <- data.frame(dat$x, y = dat$y)
  
  # Split into training and test sets
  idx <- sample.int(3 * n, size = 2 * n)
  train <- dat[idx, ]
  test <- dat[-idx, ]
  
  # Train and test models
  f <- lm(y ~ ., data = train)
  y_hat <- predict(f, test)
  loss <- log((test$y - y_hat)^2)
  
  # Expected variance 
  sigma_hat <- sigma_mod(loss)^2
  
  # Simulate loss0 via Monte Carlo
  loss0 <- loss + rnorm(n, sd = sigma_hat)
  
  # Define delta vector and precision weights
  delta <- loss - loss0
  wts <- 1 / sigma_hat^2
  
  # Run test
  t_test <- wtd.t.test(delta, weight = wts, alternative = 'greater')
  
  # Calculate empirical coverage at alpha = 0.05
  cvg <- 1 - (sum(loss >= loss0 + qnorm(0.975) * sigma_hat) + 
              sum(loss <= loss0 - qnorm(0.975) * sigma_hat)) / n
  
  # Export results
  out <- data.frame(
         'CPI' = t_test$additional[1],
          'SE' = t_test$additional[4],
           't' = t_test$coefficients[1],
     'p.value' = t_test$coefficients[3],
    'Coverage' = cvg,
    'Variance_Model' = 'Known'
  )
  return(out)
  
}

# Execute in parallel
out2 <- foreach(b = seq_len(1e5), .combine = rbind) %dopar% 
  loop(n = 1000, b = b)
saveRDS(out2, 'wtd_fixed.rds')


# PART III
# That's all well and good with *fixed* sigma_mod, but what if we need to 
# estimate variance directly from the data?

loop <- function(n, b) {
  
  # Simulate data
  dat <- mlbench.friedman1(3 * n)
  dat <- dat0 <- data.frame(dat$x, y = dat$y)
  dat0[, 10] <- dat[sample.int(3 * n), 10]
  
  # Split into training and test sets
  idx <- sample.int(3 * n, size = 2 * n)
  train <- dat[idx, ]
  train0 <- dat0[idx, ]
  test <- dat[-idx, ]
  test0 <- dat0[-idx, ]
  
  # Train models
  f <- lm(y ~ ., data = train)
  f0 <- lm(y ~ ., data = train0)
  
  # Test models
  y_hat <- predict(f, test)
  y_hat0 <- predict(f0, test0)
  loss <- log((test$y - y_hat)^2)
  loss0 <- log((test0$y - y_hat0)^2)
  
  # Fit sigma_mod, estimate weights
  delta <- loss - loss0
  l <- lowess(x = loss0, y = sqrt(abs(delta)))
  sigma_mod <- approxfun(l, rule = 2)
  sigma_hat <- sigma_mod(loss0)^2
  wts <- 1 / sigma_hat^2
  
  # Run test
  t_test <- wtd.t.test(delta, weight = wts, alternative = 'greater')
  
  # Calculate empirical coverage at alpha = 0.05
  cvg <- 1 - (sum(loss >= loss0 + qnorm(0.975) * sigma_hat) + 
              sum(loss <= loss0 - qnorm(0.975) * sigma_hat)) / n
  
  # Export results
  out <- data.frame(
         'CPI' = t_test$additional[1],
          'SE' = t_test$additional[4],
           't' = t_test$coefficients[1],
     'p.value' = t_test$coefficients[3],
    'Coverage' = cvg,
    'Variance_Model' = 'Estimated'
  )
  return(out)
  
}

# Execute in parallel
out3 <- foreach(b = seq_len(1e5), .combine = rbind) %dopar% 
  loop(n = 1000, b = b)
saveRDS(out3, 'wtd_rdm.rds')


# PART IV
# Check the t-statistics
ggplot(out1, aes(t)) +
  geom_histogram(aes(y = ..density..), bins = 100, color = 'black') +
  stat_function(fun = dt, color = 'red', args = list(df = 999)) +  
  theme_bw()
ggplot(out2, aes(t)) +
  geom_histogram(aes(y = ..density..), bins = 100, color = 'black') +
  stat_function(fun = dt, color = 'red', args = list(df = 999)) +  
  theme_bw()
ggplot(out3, aes(t)) +
  geom_histogram(aes(y = ..density..), bins = 100, color = 'black') +
  stat_function(fun = dt, color = 'red', args = list(df = 999)) +  
  theme_bw()

# Check the coverage
rbind(out1, out2, out3) %>%
  ggplot(aes(Variance_Model, Coverage, fill = Variance_Model)) + 
  geom_boxplot() + 
  geom_hline(yintercept = 0.95, color = 'red') + 
  theme_bw() + 
  scale_fill_d3()