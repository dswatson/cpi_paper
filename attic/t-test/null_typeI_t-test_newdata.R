
library(ggplot2)

n <- 1000
p <- 10

res <- replicate(10000, {
  # Training data 1
  x <- matrix(runif(n * p), ncol = p,
              dimnames = list(NULL, paste0('x', seq_len(p))))
  y <- rnorm(n)
  dat_train1 <- data.frame(y = y, x)
  
  # Training data 2
  x <- matrix(runif(n * p), ncol = p,
              dimnames = list(NULL, paste0('x', seq_len(p))))
  y <- rnorm(n)
  dat_train2 <- data.frame(y = y, x)
  
  # Test data 1
  x <- matrix(runif(n * p), ncol = p,
              dimnames = list(NULL, paste0('x', seq_len(p))))
  y <- rnorm(n)
  dat_test1 <- data.frame(y = y, x)
  
  # Test data 2
  x <- matrix(runif(n * p), ncol = p,
              dimnames = list(NULL, paste0('x', seq_len(p))))
  y <- rnorm(n)
  dat_test2 <- data.frame(y = y, x)
  
  # Fit full model (Training data 1), predict on test data 1 
  fit_full <- lm(y ~ ., dat_train1)
  pred_full <- predict(fit_full, dat_test1)
  err_full <- (pred_full - dat_test1$y)^2
  
  # Fit reduced model (Training data 2), predict on test data 2
  fit_reduced <- lm(y ~ ., dat_train2)
  pred_reduced <- predict(fit_reduced, dat_test2)
  err_reduced <- (pred_reduced - dat_test2$y)^2
  
  # Fisher test
  x <- err_reduced
  y <- err_full
  dif <- x - y
  t.test(dif)$p.value
})

# Plot histogram
res_long <- data.frame(pval = res)
ggplot(res_long, aes(pval)) + 
  geom_histogram(aes(y = ..density..), bins = 100, position="identity") 
ggsave("null_typeI_t-test_newdata.pdf")

# Type I error
mean(res < 0.05)

