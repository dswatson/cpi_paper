
library(mlr)
library(foreach)

# TODO: Add confidence intervals
brute_force_mlr <- function(task, learner, 
                            resampling = NULL,
                            test_data = NULL,
                            measure = NULL,
                            test = "t",
                            permute = TRUE,
                            log = TRUE,
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
  
  if (is.character(measure)) {
    measure <- eval(parse(text = measure))
  }
  
  if (!is.null(test)) {
    if (!(measure$id %in% c("mse", "mae", "mmce", "logloss"))) {
      stop("Statistical testing currently only implemented for 'mse', 'mae', 'mmce' and 'logloss' measures.")
    }
  }
  
  if (getTaskType(task) == "classif") {
    if (!hasLearnerProperties(learner, "prob")) {
      stop("For classification the learner requires probability support.")
    }
  }
  
  # Create resampling instance
  if (is.null(resampling)) {
    if (is.null(test_data)) {
      stop("Either resampling or test_data argument required.")
    }
  } else if (is.list(resampling)) {
    resample_instance <- makeResampleInstance(desc = resampling, task = task)
  } else if (resampling %in% c("oob", "none")) {
    resample_instance <- resampling
  } else {
    stop("Unknown resampling value.")
  }
  
  # Fit learner and compute performance
  pred_full <- fit_learner(learner = learner, task = task, resampling = resample_instance, measure = measure, test_data = test_data, verbose = verbose)
  aggr_full <- performance(pred_full, measure)
  if (!is.null(test)) {
    err_full <- compute_loss(pred_full, measure)
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
    
    pred_reduced <- fit_learner(learner = learner, task = reduced_task, resampling = resample_instance, measure = measure, test_data = test_data, verbose = verbose)
    aggr_reduced <- performance(pred_reduced, measure)
    
    if (log) {
      cpi <- log(aggr_reduced) - log(aggr_full)
    } else {
      cpi <- err_reduced - err_full
    }
    
    res <- data.frame(Variable = getTaskFeatureNames(task)[i],
                      CPI = unname(cpi), 
                      stringsAsFactors = FALSE)

    # Statistical testing
    if (!is.null(test)) {
      err_reduced <- compute_loss(pred_reduced, measure)
      if (log) {
        dif <- log(err_reduced) - log(err_full)
      } else {
        dif <- err_reduced - err_full
      }
      if (test == "fisher") {
        orig_mean <- mean(dif)
        
        # B permutations
        B <- 10000
        perm_means <- replicate(B, {
          signs <- sample(c(-1, 1), length(dif), replace = TRUE)
          mean(signs * dif)
        })
        res$p.value <- sum(perm_means >= orig_mean)/B
      } else if (test == "t") {
        test_result <- t.test(dif, alternative = 'greater')
        res$statistic <- test_result$statistic
        res$p.value <- test_result$p.value
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

fit_learner <- function(learner, task, resampling = NULL, measure = NULL, test_data = NULL, verbose = FALSE) {
  if (!is.null(test_data)) {
    # Compute error on test data
    mod <- train(learner, task)
    pred <- predict(mod, newdata = test_data)
  } else if (is.list(resampling)) {
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

compute_loss <- function(pred, measure) {
  if (getTaskType(pred) == "regr") {
    if (measure$id == "mse") {
      # Squared errors
      loss <- (pred$data$truth - pred$data$response)^2
    } else if (measure$id == "mae") {
      # Absolute errors
      loss <- abs(pred$data$truth - pred$data$response)
    } else {
      stop("Unknown measure.")
    }
  } else if (getTaskType(pred) == "classif") {
    if (measure$id == "logloss") {
      # logloss (taken from mlr::measureLogloss)
      probabilities <- pred$data[, paste("prob", pred$task.desc$class.levels, sep = ".")]
      
      #let's confine the predicted probabilities to [eps,1 - eps], so logLoss doesn't reach infinity under any circumstance
      eps <- 1e-15
      probabilities[probabilities > 1 - eps] <- 1 - eps
      probabilities[probabilities < eps] <- eps
      
      truth <- match(as.character(pred$data$truth), pred$task.desc$class.levels)
      p <- mlr:::getRowEls(probabilities, truth)
      loss <- -log(p)
    } else if (measure$id == "mmce") {
      # Misclassification error
      loss <- 1*(pred$data$truth == pred$data$response)
      
      # Avoid 0 and 1
      eps <- 1e-15
      loss[loss > 1 - eps] <- 1 - eps
      loss[loss < eps] <- eps
    } else {
      stop("Unknown measure.")
    }
    
    
  } else {
    stop("Unknown task type.")
  }
  loss
}

