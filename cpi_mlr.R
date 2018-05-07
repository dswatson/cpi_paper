
library(mlr)
library(foreach)

# TODO: Add confidence intervals
# TODO: Don't offer choice of measure but always use MSE? Classification?
brute_force_mlr <- function(task, learner, 
                            resampling = makeResampleDesc("CV", iters = 5), 
                            measures = NULL, 
                            test = NULL,
                            verbose = FALSE, 
                            cores = 1) {
  # Full model resampling
  full_model <- resample(learner, task, resampling = resampling, measures = measures, show.info = verbose)
  if (!is.null(test)) {
    err_full <- (full_model$pred$data$truth - full_model$pred$data$response)^2
  }
  
  # For each feature, fit reduced model and return difference in error
  cpi_fun <- function(i) {
    reduced_task <- subsetTask(task, features = getTaskFeatureNames(task)[-i])
    reduced_model <- resample(learner, reduced_task, resampling = resampling, measures = measures, show.info = verbose)
    
    res <- data.frame(Variable = getTaskFeatureNames(task)[i],
                      CPI = unname(reduced_model$aggr - full_model$aggr), 
                      stringsAsFactors = FALSE)

    # Statistical testing
    if (!is.null(test)) {
      err_reduced <- (reduced_model$pred$data$truth - reduced_model$pred$data$response)^2
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

