
library(mlr)
library(ggplot2)

fit_learner <- function(learner, task, test_data) {
  # Compute error on training data
  mod <- train(learner, task)
  pred <- predict(mod, newdata = test_data)
  pred
}

compute_loss <- function(pred) {
  (pred$data$truth - pred$data$response)^2
}

n <- 100
p <- 10

res <- replicate(1000, {
  x <- matrix(runif(n * p), ncol = p,
              dimnames = list(NULL, paste0('x', seq_len(p))))
  y <- rnorm(n)
  dat <- data.frame(y = y, x)
  task_train <- makeRegrTask(data = dat, target = "y")
  
  x <- matrix(runif(n * p), ncol = p,
              dimnames = list(NULL, paste0('x', seq_len(p))))
  y <- rnorm(n)
  dat <- data.frame(y = y, x)
  task_test <- makeRegrTask(data = dat, target = "y")
  test_data <- getTaskData(task_test)
  
  learner <-  makeLearner("regr.lm")
  
  pred_full <- fit_learner(learner = learner, task = task_train, test_data = test_data)
  err_full <- compute_loss(pred_full)
  
  i <- 1
  reduced_data <- getTaskData(task_train)
  reduced_data[, getTaskFeatureNames(task_train)[i]] <- sample(reduced_data[, getTaskFeatureNames(task_train)[i]])
  reduced_task <- changeData(task_train, reduced_data)
  
  pred_reduced <- fit_learner(learner = learner, task = reduced_task,  test_data = test_data)
  err_reduced <- compute_loss(pred_reduced)
  
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

ggplot(res_long, aes(x = Source, y = CPI)) + 
  geom_boxplot() + 
  geom_hline(yintercept = 0, col = "red")
ggsave("perm_boxplot.pdf")

ggplot(res_long, aes(CPI, fill = Source)) + 
  geom_histogram(aes(y = ..density..), bins = 100, alpha = 0.8, position="identity") + 
  geom_vline(xintercept = 0, col = "red")
ggsave("perm_hist.pdf")
