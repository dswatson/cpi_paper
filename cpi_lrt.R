### Likelihood ratio tests for CPI? ###

cpi_lrt <- function(x, 
                    y, 
                    type = 'regression',
                    mtry = NULL, 
                    B = NULL, 
                    replace = TRUE,
                    n.cores = 1, 
                    seed = NULL) {
  
  # Preliminaries
  n <- nrow(x)
  p <- ncol(x)
  df <- data.frame(x, 'y' = y)
  if (type != 'regression') {
    df$y <- as.factor(df$y)
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
  
  ### Part I: Full forest ###
  if (type == 'regression') {
    # Grow full forest
    rf <- ranger(data = df, dependent.variable.name = 'y', 
                 mtry = mtry, num.trees = B, replace = replace,
                 num.threads = n.cores, seed = seed)
    # Calculate log likelihood
    y_hat <- rf$predictions
    rmse <- sqrt(rf$prediction.error)
    ll <- sum(log(dnorm(y, mean = y_hat, sd = rmse)))
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
    # Calculate log likelihood
    ll <- sum(y * log(y_hat) + (1 - y) * log(1 - y_hat))
  }
  
  ### Part II: Null forests ###
  drop <- function(j) {
    # Drop j
    df0 <- df[, -j]
    # Build null model
    if (type == 'regression') {
      # Grow forest
      rf0 <- ranger(data = df0, dependent.variable.name = 'y', 
                    mtry = mtry, num.trees = B, replace = replace,
                    num.threads = 1, seed = seed)
      # Calculate log likelihood
      y_hat0 <- rf0$predictions
      rmse0 <- sqrt(rf0$prediction.error)
      ll0 <- sum(log(dnorm(y, mean = y_hat0, sd = rmse0)))
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
      # Calculate log likelihood
      ll0 <- sum(y * log(y_hat0) + (1 - y) * log(1 - y_hat0))
    }
    # Test
    cpi <- rf0$prediction.error - rf$prediction.error
    d <- 2 * (ll - ll0)
    p <- pchisq(d, df = 1, lower = FALSE)
    # Export
    out <- c(cpi, d, p)
    return(out)
  }
  
  # Execute in parallel
  res <- foreach(j = seq_len(p), .combine = rbind) %dopar% drop(j)
  dimnames(res) <- list(NULL, c('CPI', 'Deviance', 'p.value'))
  res <- data.frame(Feature = colnames(x), res)
  return(res)
  
}
  
  
  
  
  
  
  
  
  
  
  