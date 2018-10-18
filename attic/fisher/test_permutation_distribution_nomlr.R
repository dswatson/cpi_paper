
library(ggplot2)

n <- 100
p <- 2

res <- replicate(10000, {
  # Training data
  x <- matrix(runif(n * p), ncol = p,
              dimnames = list(NULL, paste0('x', seq_len(p))))
  y <- rnorm(n)
  dat_train <- data.frame(y = y, x)
  
  # Test data
  x <- matrix(runif(n * p), ncol = p,
              dimnames = list(NULL, paste0('x', seq_len(p))))
  y <- rnorm(n)
  dat_test <- data.frame(y = y, x)

  # Permuted data
  sample_idx <- sample(1:n)
  dat_train_reduced <- dat_train
  dat_train_reduced$x1 <- dat_train_reduced$x1[sample_idx]
  
  # Fit full model
  fit_full <- lm(y ~ ., dat_train)
  pred_full <- predict(fit_full, dat_test)
  err_full <- (pred_full - dat_test$y)^2
  
  # Fit reduced/permuted model
  fit_reduced <- lm(y ~ ., dat_train_reduced)
  pred_reduced <- predict(fit_reduced, dat_test)
  err_reduced <- (pred_reduced - dat_test$y)^2

  # Fisher test
  x <- err_reduced
  y <- err_full
  dif <- x - y
  orig_mean <- mean(dif)
  
  signs <- sample(c(-1, 1), length(dif), replace = TRUE)
  perm_mean <- mean(signs * dif)
  
  c(orig = orig_mean, perm = perm_mean)
})

res_long <- data.table::melt(res)
colnames(res_long) <- c("Source", "Repl", "CPI")

# Plot original and permutation distribution
ggplot(res_long, aes(x = Source, y = CPI)) + 
  geom_boxplot() + 
  geom_hline(yintercept = 0, col = "red")
#ggsave("perm_boxplot.pdf")

# ggplot(res_long, aes(CPI, fill = Source)) + 
#   geom_histogram(aes(y = ..density..), bins = 100, alpha = 0.8, position="identity") + 
#   geom_vline(xintercept = 0, col = "red")
# ggsave("perm_hist.pdf")
