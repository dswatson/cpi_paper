
library(mlr)
library(foreach)

# TODO: Add confidence intervals
# TODO: How to compute loss for classification?
brute_force_mlr <- function(task, learner, 
                            resampling = makeResampleDesc("CV", iters = 5), 
                            test = NULL,
                            verbose = FALSE, 
                            cores = 1) {
  
  # Full model resampling
  full_model <- resample(learner, task, resampling = resampling, show.info = verbose)
  err_full <- compute_loss(full_model)
  
  # For each feature, fit reduced model and return difference in error
  cpi_fun <- function(i) {
    reduced_task <- subsetTask(task, features = getTaskFeatureNames(task)[-i])
    reduced_model <- resample(learner, reduced_task, resampling = resampling, show.info = verbose)
    err_reduced <- compute_loss(reduced_model)
    
    res <- data.frame(Variable = getTaskFeatureNames(task)[i],
                      CPI = mean(err_reduced - err_full), 
                      stringsAsFactors = FALSE)

    # Statistical testing
    if (!is.null(test)) {
      if (test == "t") {
        test_result <- t.test(err_reduced, err_full, paired = TRUE, alternative = 'greater')
      } else if (test == "wilcox") {
        test_result <- wilcox.test(err_reduced, err_full, paired = TRUE, alternative = 'greater')
      } else {
        stop("Unknown test.")
      }
      res$statistic <- test_result$statistic
      res$p.value <- test_result$p.value
    }
    res
  }
  
  # Run in parallel if >1 cores
  if (cores == 1) {
    foreach(j = seq_len(getTaskNFeats(task)), .combine = rbind) %do% cpi_fun(j)
  } else {
    foreach(j = seq_len(getTaskNFeats(task)), .combine = rbind) %dopar% cpi_fun(j)
  }
}

compute_loss <- function(model) {
  if (getTaskType(model) == "regr") {
    # MSE
    loss <- (model$pred$data$truth - model$pred$data$response)^2
  } else if (getTaskType(task) == "classif") {
    if (hasLearnerProperties(model$learner, "prob")) {
      # 0/1 loss
      y <- model$pred$data$truth
      y_hat <- model$pred$data$response
      loss <- -(y * log(y_hat) + (1 - y) * log(1 - y_hat))
    } else {
      stop("For classification the learner requires probability support.")
    }
  } else {
    stop("Unknown task type.")
  }
  loss
}