### BRUTE FORCE ###
brute_force <- function(x, 
                        y, 
                        type = 'regression', 
                        test = 't',
                        conf.int = FALSE,
                        conf.level = 0.95,
                        weights = FALSE,
                        p.adj = NULL, 
                        mtry = NULL, 
                        B = NULL, 
                        replace = TRUE,
                        n.cores = 1, 
                        seed = NULL) {
  
  # Preliminaries
  require(ranger)
  require(matrixStats)
  require(foreach)
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
  if (weights) {
    if (type == 'classification') {
      stop('Weights cannot be applied with classification forests. ',
           'Rerun with type = "probability".')
    }
    if (test == 'wilcox') {
      stop('Weights cannot be applied with the wilcox test. ',
           'Rerun with test = "t".')
    }
  }
  
  ### Part I: Full forest ###
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
      y_hat <- rf$predictions[, 1]
    } else if (type == 'classification') {
      # Grow full forest
      rf <- ranger(data = df, dependent.variable.name = 'y', 
                   mtry = mtry, num.trees = B, replace = replace,
                   keep.inbag = TRUE, num.threads = n.cores, seed = seed,
                   classification = TRUE)
      # Extract probabilities
      oob_idx <- ifelse(simplify2array(rf$inbag.counts) == 0, TRUE, NA)
      preds <- predict(rf, df, num.threads = n.cores,
                       predict.all = TRUE)$predictions - 1
      y_hat <- rowMeans2(oob_idx * preds, na.rm = TRUE)
    }
    # Calculate sample-wise loss
    loss <- -(y * log(y_hat) + (1 - y) * log(1 - y_hat))
  }
  
  ### Part II: Null forests ###
  drop <- function(j) {
    # Drop j
    df0 <- df[, -j]
    # Build null model
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
        y_hat0 <- f0$predictions[, 1]
      } else if (type == 'classification') {
        # Grow forest
        rf0 <- ranger(data = df0, dependent.variable.name = 'y',
                      mtry = mtry, num.trees = B, replace = replace,
                      keep.inbag = TRUE, num.threads = 1, seed = seed,
                      classification = TRUE)
        # Extract probabilities
        oob_idx <- ifelse(simplify2array(rf0$inbag.counts) == 0, TRUE, NA)
        preds <- predict(rf0, df0, num.threads = 1,
                         predict.all = TRUE)$predictions - 1
        y_hat0 <- rowMeans2(oob_idx * preds, na.rm = TRUE)
      }
      # Calculate sample-wise loss
      loss0 <- -(y * log(y_hat0) + (1 - y) * log(1 - y_hat0))
    }
    if (is.null(test)) {
      out <- mean(loss0 - loss)
    } else if (test == 't') {
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
    } else if (test == 'wilcox') {
      w_test <- wilcox.test(loss0, loss, paired = TRUE, alternative = 'greater')
      out <- c(mean(loss0 - loss), w_test$statistic, w_test$p.value)
    }
    if (conf.int) {
      q <- qnorm(1 - (1 - conf.level) / 2)
      if (weights) {
        out <- c(out[1:2], out[1] - out[2] * q, out[1] + out[2] * q, out[3:4])
      } else {
        delta <- loss0 - loss
        se <- sd(delta) / sqrt(n)
        lower <- mean(delta) - se * q
        upper <- mean(delta) + se * q
        if (is.null(test)) {
          out <- c(out, lower, upper)
        } else if (test == 't') {
          out <- c(out[1:2], lower, upper, out[3:4])
        } else if (test == 'wilcox') {
          out <- c(out[1], lower, upper, out[2:3])
        }
      }
    }
    
    # Export
    return(out)
  }
  
  # Optionally execute in parallel
  if (is.null(test) & !conf.int) {
    if (n.cores > 1) {
      res <- foreach(j = seq_len(p), .combine = c) %dopar% drop(j)
    } else {
      res <- foreach(j = seq_len(p), .combine = c) %do% drop(j)
    }
    names(res) <- colnames(x)
  } else {
    if (n.cores > 1) {
      res <- foreach(j = seq_len(p), .combine = rbind) %dopar% drop(j)
    } else {
      res <- foreach(j = seq_len(p), .combine = rbind) %do% drop(j)
    }
    if (is.null(test) & conf.int) {
      dimnames(res) <- list(NULL, c('CPI', 'CI.L', 'CI.R'))
    }
    if (test == 't') {
      if (conf.int) {
        dimnames(res) <- list(NULL, c('CPI', 'SE', 'CI.L', 'CI.R', 't', 'p.value'))
      } else {
        dimnames(res) <- list(NULL, c('CPI', 'SE', 't', 'p.value'))
      }
    } else if (test == 'wilcox') {
      if (conf.int) {
        dimnames(res) <- list(NULL, c('CPI', 'CI.L', 'CI.R', 'W', 'p.value'))
      } else {
        dimnames(res) <- list(NULL, c('CPI', 'W', 'p.value'))
      }
    }
    res <- data.frame(Feature = colnames(x), res)
  }
  
  # Adjust p-values?
  if (!is.null(p.adj) & !is.null(test)) {
    if (p.adj %in% c('fdr', 'BH')) {
      res$q.value <- p.adjust(res$p.value, method = 'fdr')
    } else {
      res$adj.p.value <- p.adjust(res$p.value, method = p.adj)
    }
  }
  return(res)
  
}

### FOREST SPLITTING ###
rf_split <- function(x, 
                     y, 
                     type = 'regression', 
                     test = 't',
                     conf.int = FALSE,
                     conf.level = 0.95,
                     weights = FALSE,
                     p.adj = NULL,
                     mtry = NULL, 
                     B = NULL, 
                     B0 = NULL, 
                     n.sub = 500,
                     replace = TRUE,
                     n.cores = 1, 
                     seed = NULL) {
  
  # Preliminaries
  require(ranger)
  require(matrixStats)
  require(foreach)
  n <- nrow(x)
  p <- ncol(x)
  if (type == 'regression') {
    df <- data.frame(x, 'y' = y)
  } else {
    df <- data.frame(x, 'y' = as.factor(y))
  }
  if (is.null(mtry)) {
    if (type == 'regression') {
      mtry <- floor(p / 3)
    } else {
      mtry <- floor(sqrt(p))
    }
  }
  if (is.null(B)) {
    B <- 500
  }
  if (is.null(B0)) {
    B0 <- round(B / 2)
  }
  if (weights) {
    if (type == 'classification') {
      stop('Weights cannot be applied with classification forests. ',
           'Rerun with type = "probability".')
    }
    if (test == 'wilcox') {
      stop('Weights cannot be applied with the wilcox test. ',
           'Rerun with test = "t".')
    }
  }
  
  ### Part I: Full forest ###
  if (type == 'regression') {
    # Grow forest
    rf <- ranger(data = df, dependent.variable.name = 'y', 
                 mtry = mtry, num.trees = B, replace = replace,
                 keep.inbag = TRUE, num.threads = n.cores, seed = seed)
    # Create n x B matrix of predictions
    preds <- predict(rf, data = df, num.threads = n.cores,
                     predict.all = TRUE)$predictions 
  } else if (type == 'probability') {
    # Grow forest
    rf <- ranger(data = df, dependent.variable.name = 'y', 
                 mtry = mtry, num.trees = B, replace = replace,
                 keep.inbag = TRUE, num.threads = n.cores, seed = seed,
                 probability = TRUE)
    # Create n x B matrix of predictions
    preds <- predict(rf, data = df, num.threads = n.cores,
                     predict.all = TRUE)$predictions[, 1, ]
  } else if (type == 'classification') {
    # Grow forest
    rf <- ranger(data = df, dependent.variable.name = 'y', 
                 mtry = mtry, num.trees = B, replace = replace,
                 keep.inbag = TRUE, num.threads = n.cores, seed = seed,
                 classification = TRUE)
    # Create n x B matrix of predictions
    preds <- predict(rf, df, num.threads = n.cores,
                     predict.all = TRUE)$predictions - 1
  }
  # Create oob_preds object
  oob_idx <- ifelse(simplify2array(rf$inbag.counts) == 0, TRUE, NA)
  oob_preds <- oob_idx * preds

  # Determine if forest splitting is appropriate
  if (n.cores > 1) {
    splits <- foreach(b = seq_len(B), .combine = c) %dopar% 
      sum(!treeInfo(rf, b)$terminal)
  } else {
    splits <- sapply(seq_len(B), function(b) {
      sum(!treeInfo(rf, b)$terminal)
    })
  }
  exp_B0 <- round(B * (1 - mtry/p)^mean(splits))
  if (exp_B0 < B0) {
    stop('Expected number of trees in each sub-forest is ', exp_B0, 
         ', less than the minimum of ', B0, '.')
  }
  
  # Calculate loss from random sub-forests
  sub_forest <- function(b) {
    set.seed(b)
    sub_idx <- sample.int(B, exp_B0)
    y_hat <- rowMeans2(oob_preds[, sub_idx, drop = FALSE], na.rm = TRUE)
    if (type == 'regression') {
      loss <- (y_hat - y)^2
    } else if (type %in% c('probability', 'classification')) {
      loss <- -(y * log(y_hat) + (1 - y) * log(1 - y_hat))
    }
    if (weights) {
      inbag_cts <- simplify2array(rf$inbag.counts)
      wts <- 1 / rInfJack(pred = preds[, sub_idx, drop = FALSE], 
                          inbag = inbag_cts[, sub_idx, drop = FALSE])$var.hat
    } else {
      wts <- rep(NA_real_, n)
    }
    return(c(loss, wts))
  }
  
  # Optionally execute in parallel
  if (n.cores > 1) {
    subs <- rowMeans2(foreach(b = seq_len(n.sub), .combine = cbind) %dopar% sub_forest(b))
  } else {
    subs <- rowMeans2(foreach(b = seq_len(n.sub), .combine = cbind) %do% sub_forest(b))
  }
  loss <- subs[seq_len(n)]
  if (weights) {
    wts <- subs[(n + 1):(2 * n)]
  }
  
  ### Part II: Null forests ###
  drop <- function(j) {
    f0_idx <- !rf$forest$variable.selected[, j]
    y_hat0 <- rowMeans2(oob_preds[, f0_idx, drop = FALSE], na.rm = TRUE)
    if (type == 'regression') {
      loss0 <- (y_hat0 - y)^2
    } else {
      loss0 <- -(y * log(y_hat0) + (1 - y) * log(1 - y_hat0))
    }
    if (weights) {
      inbag_cts <- simplify2array(rf$inbag.counts)
      wts0 <- 1 / rInfJack(pred = preds[, f0_idx, drop = FALSE], 
                           inbag = inbag_cts[, f0_idx, drop = FALSE])$var.hat
    }
    if (is.null(test)) {
      out <- mean(loss0 - loss)
    } else if (test == 't') {
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
    } else if (test == 'wilcox') {
      w_test <- wilcox.test(loss0, loss, paired = TRUE, alternative = 'greater')
      out <- c(mean(loss0 - loss), w_test$statistic, w_test$p.value)
    }
    if (conf.int) {
      q <- qnorm(1 - (1 - conf.level) / 2)
      if (weights) {
        out <- c(out[1:2], out[1] - out[2] * q, out[1] + out[2] * q, out[3:4])
      } else {
        delta <- loss0 - loss
        se <- sd(delta) / sqrt(n)
        lower <- mean(delta) - se * q
        upper <- mean(delta) + se * q
        if (is.null(test)) {
          out <- c(out, lower, upper)
        } else if (test == 't') {
          out <- c(out[1:2], lower, upper, out[3:4])
        } else if (test == 'wilcox') {
          out <- c(out[1], lower, upper, out[2:3])
        }
      }
    }
    return(out)
  }
  
  # Optionally execute in parallel
  if (is.null(test) & !conf.int) {
    if (n.cores > 1) {
      res <- foreach(j = seq_len(p), .combine = c) %dopar% drop(j)
    } else {
      res <- foreach(j = seq_len(p), .combine = c) %do% drop(j)
    }
    names(res) <- colnames(x)
  } else {
    if (n.cores > 1) {
      res <- foreach(j = seq_len(p), .combine = rbind) %dopar% drop(j)
    } else {
      res <- foreach(j = seq_len(p), .combine = rbind) %do% drop(j)
    }
    if (is.null(test) & conf.int) {
      dimnames(res) <- list(NULL, c('CPI', 'CI.L', 'CI.R'))
    }
    if (test == 't') {
      if (conf.int) {
        dimnames(res) <- list(NULL, c('CPI', 'SE', 'CI.L', 'CI.R', 't', 'p.value'))
      } else {
        dimnames(res) <- list(NULL, c('CPI', 'SE', 't', 'p.value'))
      }
    } else if (test == 'wilcox') {
      if (conf.int) {
        dimnames(res) <- list(NULL, c('CPI', 'CI.L', 'CI.R', 'W', 'p.value'))
      } else {
        dimnames(res) <- list(NULL, c('CPI', 'W', 'p.value'))
      }
    }
    res <- data.frame(Feature = colnames(x), res)
  }
  
  # Adjust p-values?
  if (!is.null(p.adj) & !is.null(test)) {
    if (p.adj %in% c('fdr', 'BH')) {
      res$q.value <- p.adjust(res$p.value, method = 'fdr')
    } else {
      res$adj.p.value <- p.adjust(res$p.value, method = p.adj)
    }
  }
  return(res)
  
}






# Problems, ideas:
# Not sure how to extend this to multi-class problems (k > 2)?
# Extend to survival forests?
# Need to revise cross entropy formula for factor inputs
