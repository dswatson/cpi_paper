
# Linear data ----------------------------------------------------------------
linear_data <- function(n, p, outcome = "regr", cov_base = 0, ...) {
  beta <- rep(seq(0, .9, length.out = 10), each = p/10)
  beta0 <- 0
  
  sigma <- toeplitz(cov_base^(0:(p-1)))
  
  x <- matrix(rmvnorm(n = n, sigma = sigma), ncol = p,
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
nonlinear_data <- function(n, p, outcome = "regr", cov_base = 0, ...) {
  beta <- rep(seq(0, .9, length.out = 10), each = p/10)
  beta0 <- 0
  
  sigma <- toeplitz(cov_base^(0:(p-1)))
  
  x <- matrix(rmvnorm(n = n, sigma = sigma), ncol = p,
              dimnames = list(NULL, paste0('x', seq_len(p))))
  idx <- x < -qnorm(0.75) | x > qnorm(0.75)
  xx <- matrix(0, nrow = n, ncol = p)
  xx[idx] <- -1
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

data_train_test <- function(data, job, type = "linear", ...) {
  if (type == "linear") {
    list(train = linear_data(...), 
         test = linear_data(...))      
  } else if (type == "nonlinear") {
    list(train = nonlinear_data(...), 
         test = nonlinear_data(...))     
  } else {
    stop("Unknown data type.")
  }
}

