
# Linear data ----------------------------------------------------------------
linear_data <- function(data, job, n, p, outcome = "regr", ...) {
  beta <- rep(c(0, 0, -.5, .5, -1, 1, -1.5, 1.5, -2, 2), each = p/10)
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

# Non-linear data ----------------------------------------------------------------
nonlinear_data <- function(data, job, n, p, outcome = "regr", ...) {
  beta <- rep(c(0, 0, -.5, .5, -1, 1, -1.5, 1.5, -2, 2), each = p/10)
  beta0 <- 0
  
  x <- matrix(runif(n * p), ncol = p,
              dimnames = list(NULL, paste0('x', seq_len(p))))
  idx <- x < .25 | x > .75
  xx <- matrix(0, nrow = n, ncol = p)
  xx[idx] <- 0
  xx[!idx] <- 1
  lp <- xx %*% beta + beta0
  
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