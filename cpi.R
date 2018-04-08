### BRUTE FORCE ###
brute_force <- function(x, y, type = 'regression', weights = FALSE,
                        mtry = NULL, B = NULL, replace = TRUE,
                        n.cores = 1, seed = NULL) {

  n <- nrow(x)
  p <- ncol(x)
  df <- data.frame(x, 'y' = y)
  if (type != 'regression') {
    df$y <- as.factor(df$y)
  }
  if (is.null(mtry)) {
    if (type == 'regression') {
      mtry <- floor(p/3)
    } else {
      mtry <- floor(sqrt(p))
    }
  }
  if (is.null(B)) {
    B <- 500
  }
  
  # Define drop function
  drop <- function(j) {
    
    # Drop j
    df0 <- df[, -j]
    
    # Build null models
    if (type == 'regression') {
      if (weights) {
        # Grow forest
        rf0 <- ranger(data = df0, dependent.variable.name = 'y', 
                      mtry = mtry, num.trees = B, replace = replace,
                      keep.inbag = TRUE, num.threads = 1, seed = seed)
        # Calculate weights
        wts0 <- 1 / predict(rf0, data = df0, num.threads = 1,
                            type = 'se')$se^2
      } else {
        # Grow forest
        rf0 <- ranger(data = df0, dependent.variable.name = 'y', 
                      mtry = mtry, num.trees = B, replace = replace,
                      num.threads = 1, seed = seed)
      }
      # Calculate sample-wise loss
      loss0 <- (rf0$predictions - y)^2
    } else {
      if (type == 'probability') {
        if (weights) {
          # Grow forest
          rf0 <- ranger(data = df0, dependent.variable.name = 'y', 
                        mtry = mtry, num.trees = B, replace = replace,
                        keep.inbag = TRUE, num.threads = 1, seed = seed, 
                        probability = TRUE) 
          # Calculate weights
          wts0 <- 1 / (predict(rf0, data = df0, num.threads = 1,
                               type = 'se')$se[, 1])^2
        } else {
          # Grow forest
          rf0 <- ranger(data = df0, dependent.variable.name = 'y', 
                        mtry = mtry, num.trees = B, replace = replace,
                        num.threads = 1, seed = seed, 
                        probability = TRUE) 
        }
        # Extract probabilities
        y_hat0 <- sapply(f0$predictions, '[[', 1)
      } else if (type == 'classification') {
        # Grow forest
        rf0 <- ranger(data = df0, dependent.variable.name = 'y',
                      mtry = mtry, num.trees = B, replace = replace,
                      keep.inbag = TRUE, num.threads = 1, seed = seed,
                      classification = TRUE)
        # Extract probabilities
        oob_idx <- ifelse(simplify2array(rf0$inbag.counts) == 0, TRUE, NA)
        preds <- predict(rf0, df0, predict.all = TRUE)$predictions
        y_hat0 <- rowMeans2(oob_idx * preds, na.rm = TRUE)
      }
      # Calculate sample-wise loss
      if (type == 'regression') {
        loss0 <- (y_hat0 - y)^2
      } else {
        loss0 <- -(y * log(y_hat0) + (1 - y) * log(1 - y_hat0))
      }
    }
    
    # Perform paired one-sided t-test with optional precision weights
    if (weights) {
      df <- data.frame(Loss = c(loss, loss0),
                       Wts = c(wts, wts0),
                       Model = rep(c('Full', 'Null'), each = n),
                       Obs = rep(paste0('n', seq_len(n)), times = 2))
      s <- summary(lm(Loss ~ Model + Obs, weights = Wts, data = df))
      out <- c(coef(s)[2, 1:3], pt(coef(s)[2, 3], s$df[2], lower = FALSE))
    } else {
      t_test <- t.test(loss0, loss, paired = TRUE, alternative = 'greater')
      out <- c(t_test$estimate, t_test$estimate / t_test$statistic,
               t_test$statistic, t_test$p.value)
    }
    
    # Export
    return(out)
    
  }
  
  # Define test
  if (type == 'regression') {
    if (weights) {
      # Grow full forest
      rf <- ranger(data = df, dependent.variable.name = 'y', 
                   mtry = mtry, num.trees = B, replace = replace,
                   keep.inbag = TRUE, num.threads = n.cores, seed = seed)
      # Calculate precision
      wts <- 1 / predict(rf, data = df, num.threads = n.cores,
                         type = 'se')$se^2
    } else {
      # Grow full forest
      rf <- ranger(data = df, dependent.variable.name = 'y', 
                   mtry = mtry, num.trees = B, replace = replace,
                   num.threads = n.cores, seed = seed)
    }
    # Calculate sample-wise loss
    loss <- (rf$predictions - y)^2
  } else {
    if (type == 'probability') {
      if (weights) {
        # Grow full forest
        rf <- ranger(data = df, dependent.variable.name = 'y', 
                     mtry = mtry, num.trees = B, replace = replace,
                     keep.inbag = TRUE, num.threads = n.cores, seed = seed,
                     probability = TRUE)
        # Calculate sample-wise precision
        wts <- 1 / predict(rf, data = df, num.threads = n.cores,
                           type = 'se')$se[, 1]^2
      } else {
        # Grow full forest
        rf <- ranger(data = df, dependent.variable.name = 'y', 
                     mtry = mtry, num.trees = B, replace = replace,
                     num.threads = n.cores, seed = seed,
                     probability = TRUE)
      }
      # Extract probabilities
      y_hat <- sapply(rf$predictions, '[[', 1)
    } else if (type == 'classification') {
      if (weights) {
        stop('Weights cannot be applied with classification forests. ',
             'Rerun with type = "probability".')
      }
      # Grow full forest
      rf <- ranger(data = df, dependent.variable.name = 'y', 
                   mtry = mtry, num.trees = B, replace = replace,
                   keep.inbag = TRUE, num.threads = n.cores, seed = seed,
                   classification = TRUE)
      # Extract probabilities
      oob_idx <- ifelse(simplify2array(rf$inbag.counts) == 0, TRUE, NA)
      preds <- predict(rf, df, predict.all = TRUE)$predictions
      y_hat <- rowMeans2(oob_idx * preds, na.rm = TRUE)
    }
    # Calculate sample-wise loss
    if (loss == 'mse') {
      loss0 <- (y_hat0 - y)^2
    } else {
      loss0 <- -(y * log(y_hat0) + (1 - y) * log(1 - y_hat0))
    }
  }
  
  # Optionally execute in parallel
  if (n.cores > 1) {
    delta <- foreach(j = seq_len(p), .combine = rbind) %dopar% drop(j)
  } else {
    delta <- foreach(j = seq_len(p), .combine = rbind) %do% drop(j)
  }
  dimnames(delta) <- list(NULL, c('Delta', 'SE', 't', 'p.value'))
  return(data.frame(Feature = colnames(x), delta))
  
}











cpi <- function(rf, dat, loss = 'mse') {
  
  p <- rf$num.independent.variables
  y_hat <- f$predictions
  loss <- (f$predictions - dat$y)^2
  tree_preds <- predict(f, data = dat, predict.all = TRUE)$predictions
  oob_idx <- ifelse(simplify2array(rf$inbag.counts) == 0, TRUE, NA)
  oob_tree_preds <- oob_idx * tree_preds
  vi <- sapply(seq_len(p), function(j) {
    f0_idx <- rf$forest$variable.selected[, j]
    y_hat0 <- rowMeans(oob_tree_preds[, f0_idx], na.rm = TRUE)
    loss0 <- (y_hat0 - dat$y)^2
    loss0 - loss
  })
  return(vi)
  
}
