
library(mlr)
library(foreach)

# TODO: Add confidence intervals
brute_force_mlr <- function(task, learner, 
                            resampling = makeResampleDesc("CV", iters = 5), 
                            measure = NULL,
                            test = NULL,
                            permute = FALSE,
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
  
  if (getTaskType(task) == "classif") {
    if (!hasLearnerProperties(learner, "prob")) {
      stop("For classification the learner requires probability support.")
    }
  }
  
  # Create resampling instance
  if (resampling %in% c("oob", "none")) {
    resample_instance <- resampling
  } else {
    resample_instance <- makeResampleInstance(desc = resampling, task = task)
  }
  
  # Fit learner and compute performance
  pred_full <- fit_learner(learner = learner, task = task, resampling = resample_instance, measure = measure, verbose = verbose)
  aggr_full <- performance(pred_full, measure)
  if (!is.null(test)) {
    err_full <- compute_loss(pred_full)
  }
  
  # For each feature, fit reduced model and return difference in error
  cpi_fun <- function(i) {
    if (permute) {
      reduced_data <- getTaskData(task)
      reduced_data[, getTaskFeatureNames(task)[i]] <- sample(reduced_data[, getTaskFeatureNames(task)[i]])
      reduced_task <- changeData(task, reduced_data)
    } else {
      reduced_task <- subsetTask(task, features = getTaskFeatureNames(task)[-i])
    }
    
    pred_reduced <- fit_learner(learner = learner, task = reduced_task, resampling = resample_instance, measure = measure, verbose = verbose)
    aggr_reduced <- performance(pred_reduced, measure)
    
    res <- data.frame(Variable = getTaskFeatureNames(task)[i],
                      CPI = unname(aggr_reduced - aggr_full), 
                      stringsAsFactors = FALSE)

    # Statistical testing
    if (!is.null(test)) {
      err_reduced <- compute_loss(pred_reduced)
      if (test == "t") {
        test_result <- t.test(err_reduced, err_full, paired = TRUE, alternative = 'greater')
        res$statistic <- test_result$statistic
        res$p.value <- test_result$p.value
      } else if (test == "U") {
        test_result <- wilcox.test(err_reduced, err_full, paired = TRUE, alternative = 'greater', exact = FALSE)
        res$statistic <- test_result$statistic
        res$p.value <- test_result$p.value
      } else if (test == "lrt") {
        if (measure$id == "mse") {
          ll_full <- sum(log(dnorm(pred_full$data$truth, mean = pred_full$data$response, sd = sqrt(aggr_full))))
          ll_reduced <- sum(log(dnorm(pred_reduced$data$truth, mean = pred_reduced$data$response, sd = sqrt(aggr_reduced))))
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

fit_learner <- function(learner, task, resampling, measure, test, verbose) {
  if (is.list(resampling)) {
    # Full model resampling
    pred <- resample(learner, task, resampling = resampling, measure = measure, show.info = verbose)$pred
  } else if (resampling == "none") {
    # Compute error on training data
    mod <- train(learner, task)
    pred <- predict(mod, task)
  } else if (resampling == "oob") {
    # Use OOB predictions if available
    if (!hasLearnerProperties(learner, "oobpreds")) {
      stop("OOB error not available for this learner.")
    }
    mod <- train(learner, task)
    pred <- getOOBPreds(mod, task)
  } else {
    stop("Unknown value for 'resampling'.")
  }
  pred
}

compute_loss <- function(pred) {
  if (getTaskType(pred) == "regr") {
    # MSE
    loss <- (pred$data$truth - pred$data$response)^2
  } else if (getTaskType(pred) == "classif") {
    # logloss (taken from mlr::measureLogloss)
    probabilities <- pred$data[, paste("prob", pred$task.desc$class.levels, sep = ".")]
    
    #let's confine the predicted probabilities to [eps,1 - eps], so logLoss doesn't reach infinity under any circumstance
    eps <- 1e-15
    probabilities[probabilities > 1 - eps] <- 1 - eps
    probabilities[probabilities < eps] <- eps
    
    truth <- match(as.character(pred$data$truth), pred$task.desc$class.levels)
    p <- mlr:::getRowEls(probabilities, truth)
    loss <- -log(p)
  } else {
    stop("Unknown task type.")
  }
  loss
}

