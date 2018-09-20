
library(ggplot2)

n <- 100
p <- 2

res <- replicate(10000, {
  # Simulate loss
  err_full <- rchisq(n, df = 1)
  err_reduced <- err_full + rnorm(n, sd = .01)
  
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
#ggsave("perm_boxplot_chi2.pdf")

# ggplot(res_long, aes(CPI, fill = Source)) + 
#   geom_histogram(aes(y = ..density..), bins = 100, alpha = 0.8, position="identity") + 
#   geom_vline(xintercept = 0, col = "red")
# ggsave("perm_hist.pdf")
