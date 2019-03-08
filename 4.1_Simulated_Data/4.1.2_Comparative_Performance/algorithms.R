
# CPI ----------------------------------------------------------------
cpi_fn <- function(data, job, instance, learner_name, ...) {
  par.vals <- switch(learner_name, 
                     regr.ranger = list(num.trees = 50), 
                     regr.nnet = list(size = 20, decay = .1, trace = FALSE), 
                     regr.svm = list(kernel = "radial"), 
                     list())
  test_data <- instance$test$env$data      
  cpi(task = instance$train, learner = makeLearner(learner_name, par.vals = par.vals), 
      test_data = test_data, ...)[, c("Variable", "p.value")]
}

# Williamson et al.'s nonparametric ANOVA ---------------------------
anova_fn <- function(data, job, instance, learner_name, ...) {
  par.vals <- switch(learner_name, 
                     regr.ranger = list(num.trees = 50), 
                     regr.nnet = list(size = 20, decay = .1, trace = FALSE), 
                     regr.svm = list(kernel = "radial"), 
                     list())
  test_data <- instance$test$env$data
  
  # Fit full model
  learner <- makeLearner(learner_name, par.vals = par.vals)
  fit <- train(learner, instance$train)
  pred <- predict(fit, newdata = test_data)
  
  # For each variable calculate ANOVA test
  pvals <- sapply(2:ncol(test_data), function(j) {
    train0 <- test_data[, -j]
    train0$y <- getPredictionResponse(pred)
    train0_task <- makeRegrTask(data = train0, target = "y")
    fit0 <- train(learner, train0_task)
    pred0 <- predict(fit0, newdata = test_data)
    vimp <- vimp_regression(test_data$y, f1 = getPredictionResponse(pred), 
                            f2 = getPredictionResponse(pred0), indx = j, 
                            run_regression = FALSE)
    pnorm(vimp$est / vimp$se, lower = FALSE)
  })
  
  # Return p-values
  data.frame(Variable = colnames(test_data)[-1], p.value = pvals)
}

# Chalupka et al.'s fast conditional independence test (FCIT) ---------------------------
fcit_fn <- function(data, job, instance, learner_name, ...) {
  par.vals <- switch(learner_name, 
                     regr.ranger = list(num.trees = 50), 
                     regr.nnet = list(size = 20, decay = .1, trace = FALSE), 
                     regr.svm = list(kernel = "radial"), 
                     list())
  test_data <- instance$test$env$data
  train_data <- instance$train$env$data
  
  # Fit full model
  learner <- makeLearner(learner_name, par.vals = par.vals)
  fit <- train(learner, instance$train)
  pred <- predict(fit, newdata = test_data)
  loss <- (test_data$y - getPredictionResponse(pred))^2
  
  # For each variable calculate FCIT
  pvals <- sapply(2:ncol(test_data), function(j) {
    train0_task <- makeRegrTask(data = train_data[, -j], target = "y")
    fit0 <- train(learner, train0_task)  
    pred0 <- predict(fit0, newdata = test_data[, -j])                         
    loss0 <- (test_data$y - getPredictionResponse(pred0))^2
    delta <- loss0 - loss                          
    t.test(delta, alternative = 'greater')$p.value
  })
  
  # Return p-values
  data.frame(Variable = colnames(train_data)[-1], p.value = pvals)
}

# Shah & Peters's generalised covariance measure (GCM) ---------------------------
gcm_fn <- function(data, job, instance, learner_name, ...) {
  par.vals <- switch(learner_name, 
                     regr.ranger = list(num.trees = 50), 
                     regr.nnet = list(size = 20, decay = .1, trace = FALSE), 
                     regr.svm = list(kernel = "radial"), 
                     list())
  test_data <- instance$test$env$data
  train_data <- instance$train$env$data
  
  # Fit full model
  learner <- makeLearner(learner_name, par.vals = par.vals)
  fit <- train(learner, instance$train)
  pred <- predict(fit, newdata = test_data)
  eps <- test_data$y - getPredictionResponse(pred)
  
  # For each variable calculate GCM
  pvals <- sapply(2:ncol(test_data), function(j) {
    train_task_f <- makeRegrTask(data = train_data[, -j], target = "y")
    fit_f <- train(learner, train_task_f)  
    pred_f <- predict(fit_f, newdata = test_data[, -j])     
    eps_f <- test_data$y - getPredictionResponse(pred_f)
    
    train_data_g <- train_data[, -1]
    colnames(train_data_g)[j - 1] <- 'y'
    test_data_g <- test_data[, -1]
    colnames(test_data_g)[j - 1] <- 'y'
                        
    train_task_g <- makeRegrTask(data = train_data_g, target = "y")
    fit_g <- train(learner, train_task_g)  
    pred_g <- predict(fit_g, newdata = test_data_g)                         
    eps_g <- test_data_g$y - getPredictionResponse(pred_g)
    r <- eps_f * eps_g
    gcm <- abs((sqrt(nrow(test_data)) * mean(r)) / 
                    sqrt(mean(r^2) - (mean(r))^2))                       
    2 * pnorm(gcm, lower = FALSE)
  })
  
  # Return p-values
  data.frame(Variable = colnames(train_data)[-1], p.value = pvals)
}

