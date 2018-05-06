# Calculate log-likelihood of a random forest fit with ranger

log_lik <- function(x, y, rf) {
  
  if (rf$treetype == 'Regression') {
    y_hat <- rf$predictions
    rmse <- sqrt(rf$prediction.error)
    ll <- sum(log(dnorm(y, mean = y_hat, sd = rmse)))
  } else if (rf$treetype == 'Classification') {
    oob_idx <- ifelse(simplify2array(rf$inbag.counts) == 0, TRUE, NA)
    preds <- predict(rf, x, predict.all = TRUE)$predictions - 1
    y_hat <- rowMeans2(oob_idx * preds, na.rm = TRUE)
    ll <- sum(y * log(y_hat) + (1 - y) * log(1 - y_hat))
  }
  return(ll)
  
} 