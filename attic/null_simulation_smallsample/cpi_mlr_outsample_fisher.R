
library(mlr)
library(foreach)

# TODO: Add confidence intervals
brute_force_mlr <- function(task, learner, 
                            test_data,
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
  
  # Fit learner and compute performance
  pred_full <- fit_learner(learner = learner, task = task, test_data = test_data)
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
    
    pred_reduced <- fit_learner(learner = learner, task = reduced_task, test_data = test_data)
    aggr_reduced <- performance(pred_reduced, measure)
    
    res <- data.frame(Variable = getTaskFeatureNames(task)[i],
                      CPI = unname(aggr_reduced - aggr_full), 
                      stringsAsFactors = FALSE)

    # Statistical testing
    if (!is.null(test)) {
      err_reduced <- compute_loss(pred_reduced)
      if (test == "fisher") {
        dif <- err_reduced - err_full
        orig_mean <- mean(dif)
        
        # B permutation
        B <- 10000
        perm_means <- replicate(B, {
          signs <- sample(c(-1, 1), length(dif), replace = TRUE)
          mean(signs * dif)
        })
        res$p.value <- sum(perm_means >= orig_mean)/B
      } else if (test == "t") {
        test_result <- t.test(err_reduced, err_full, paired = TRUE, alternative = 'greater')
        res$statistic <- test_result$statistic
        res$p.value <- test_result$p.value
      } else if (test == "U") {
        test_result <- wilcox.test(err_reduced, err_full, paired = TRUE, alternative = 'greater', exact = FALSE)
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

fit_learner <- function(learner, task, test_data) {
  # Compute error on training data
  mod <- train(learner, task)
  pred <- predict(mod, newdata = test_data)
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

