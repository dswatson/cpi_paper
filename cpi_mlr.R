
library(mlr)
library(foreach)

# TODO: Add confidence intervals
brute_force_mlr <- function(task, learner, 
                            resampling = makeResampleDesc("CV", iters = 5), 
                            measure = NULL,
                            test = NULL,
                            verbose = FALSE, 
                            cores = 1) {
  if (is.null(measure)) {
    if (getTaskType(task) == "regr") {
      measure <- mse
    } else if (getTaskType(task) == "classif") {
      measure <- logloss
    } else {
      stop("Unknown task type.")
    }
  }
  
  if (!is.null(test)) {
    if (!(measure$id %in% c("mse", "logloss"))) {
      stop("Statistical testing currently only implemented for 'mse' and 'logloss' measure.")
    }
  }
  
  # Full model resampling
  full_model <- resample(learner, task, resampling = resampling, measure = measure, show.info = verbose)
  if (!is.null(test)) {
    err_full <- compute_loss(full_model)
  }
  
  # For each feature, fit reduced model and return difference in error
  cpi_fun <- function(i) {
    reduced_task <- subsetTask(task, features = getTaskFeatureNames(task)[-i])
    reduced_model <- resample(learner, reduced_task, resampling = resampling, measure = measure, show.info = verbose)
    
    res <- data.frame(Variable = getTaskFeatureNames(task)[i],
                      CPI = unname(reduced_model$aggr - full_model$aggr), 
                      stringsAsFactors = FALSE)

    # Statistical testing
    if (!is.null(test)) {
      err_reduced <- compute_loss(reduced_model)
      if (test == "t") {
        test_result <- t.test(err_reduced, err_full, paired = TRUE, alternative = 'greater')
        res$statistic <- test_result$statistic
        res$p.value <- test_result$p.value
      } else if (test == "lrt") {
        if (measure$id == "mse") {
          ll_full <- sum(log(dnorm(full_model$pred$data$truth, mean = full_model$pred$data$response, sd = sqrt(full_model$aggr))))
          ll_reduced <- sum(log(dnorm(reduced_model$pred$data$truth, mean = reduced_model$pred$data$response, sd = sqrt(reduced_model$aggr))))
        } else if (measure$id == "logloss") {
          ll_full <- sum(-err_full)
          ll_reduced <- sum(-err_reduced)
        } else {
          stop("Likelihood ratio test not implemented for this measure.")
        }
        d <- 2 * (ll_full - ll_reduced)
        res$statistic <- d
        res$p.value <- pchisq(d, df = 1, lower = FALSE)
      } else {
        stop("Unknown test.")
      }
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
  } else if (getTaskType(model) == "classif") {
    if (hasLearnerProperties(model$learner, "prob")) {
      # logloss (taken from mlr::measureLogloss)
      probabilities <- model$pred$data[, paste("prob", model$pred$task.desc$class.levels, sep = ".")]
      
      #let's confine the predicted probabilities to [eps,1 - eps], so logLoss doesn't reach infinity under any circumstance
      eps <- 1e-15
      probabilities[probabilities > 1 - eps] <- 1 - eps
      probabilities[probabilities < eps] <- eps
      
      truth <- match(as.character(model$pred$data$truth), model$pred$task.desc$class.levels)
      p <- mlr:::getRowEls(probabilities, truth)
      loss <- -log(p)
    } else {
      stop("For classification the learner requires probability support.")
    }
  } else {
    stop("Unknown task type.")
  }
  loss
}

