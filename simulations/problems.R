
# Linear data ----------------------------------------------------------------
linear_data <- function(data, job, n, p, outcome = "regr", ...) {
  beta <- rep(c(0, 0, -1, 1, -2, 2, -3, 3, -4, 4), each = p/10)
  beta0 <- 0
  
  x <- matrix(runif(n * p), ncol = p,
              dimnames = list(NULL, paste0('x', seq_len(p))))
  lp <- x %*% beta + beta0 
  
  if (outcome == "regr") {
    y <- lp + rnorm(n)
    dat <- data.frame(y = y, x)
    makeRegrTask(data = dat, target = "y")
  } else if (outcome == "classif") {
    y <- as.factor(rbinom(n, size = 1, prob = plogis(lp)))
    dat <- data.frame(y = y, x)
    makeClassifTask(data = dat, target = "y")
  }
}