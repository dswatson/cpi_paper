### Likelihood ratio tests for CPI? ###

cpi_lrt <- function(x, y, f, f0, df, type = 'regression') {
  
  if (type == 'regression') {
    # Alternative
    y_hat <- predict(f, x)
    rmse <- sqrt(mean((y - y_hat)^2))
    ll <- sum(log(dnorm(y, mean = y_hat, sd = rmse)))
    # Null
    y_hat0 <- predict(f0, x)
    rmse0 <- sqrt(mean((y - y_hat0)^2))
    ll0 <- sum(log(dnorm(y, mean = y_hat0, sd = rmse0)))
  } else if (type == 'classification') {
    # Alternative
    y_hat <- predict(f, x)  # Make sure these are probabilites
    ll <- sum(y * log(y_hat) + (1 - y) * log(1 - y_hat))
    # Null
    y_hat0 <- predict(f0, x)
    ll0 <- sum(y * log(y_hat0) + (1 - y) * log(1 - y_hat0))
  }
  # Test
  d <- 2 * (ll - ll0)
  p <- pchisq(d, df = df)
  out <- c('deviance' = d, 'p.value' = p)
  return(out)
  
}