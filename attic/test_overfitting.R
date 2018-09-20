
library(data.table)
library(ggplot2)
source("cpi_mlr.R")

n <- 100
p <- 10
num_replicates <- 1000

sim_data <- function(n, p, ...) {
  x <- matrix(runif(n * p), ncol = p,
              dimnames = list(NULL, paste0('x', seq_len(p))))
  y <- rnorm(n)
  dat <- data.frame(y = y, x)
  makeRegrTask(data = dat, target = "y")
}

cpi_fun <- function(task, learner, resampling, measure, p_reduce) {
  # Create resampling instance
  resample_instance <- makeResampleInstance(desc = resampling, task = task)
  
  # Fit learner and compute performance
  pred_full <- fit_learner(learner = learner, task = task, resampling = resample_instance, measure = measure, verbose = FALSE)
  aggr_full <- performance(pred_full, measure)
  
  # Fit reduced learner and compute performance
  if (p_reduce == 0) {
    reduced_task <- task
  } else {
    reduced_task <- subsetTask(task, features = getTaskFeatureNames(task)[-(1:p_reduce)])
  }
  pred_reduced <- fit_learner(learner = learner, task = reduced_task, resampling = resample_instance, measure = measure, verbose = FALSE)
  aggr_reduced <- performance(pred_reduced, measure)
  
  # Return CPI
  unname(aggr_reduced - aggr_full)
}

compare_fun <- function(n, p) {
  task <- sim_data(n, p)
  sapply(0:p, cpi_fun, task = task, learner = makeLearner("regr.lm"), 
         resampling = makeResampleDesc("CV", iters = 5), measure = mse)
}

res <- replicate(num_replicates, compare_fun(n, p))
res_long <- melt(res)
colnames(res_long) <- c("Vars_dropped", "repl", "CPI")
ggplot(res_long, aes(x = factor(Vars_dropped-1), y = CPI)) + 
  geom_boxplot() + 
  geom_hline(yintercept = 0, col = "red") + 
  xlab("Number of dropped variables") + ylab("CPI value")
ggsave("test_overfitting.pdf")
