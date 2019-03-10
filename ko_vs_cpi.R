# Set working directory
setwd('~/Documents/CPI/cpi_paper/4.1_Simulated_Data')

# Set seed
set.seed(123, kind = "L'Ecuyer-CMRG")

# Load libraries, register cores
library(data.table)
library(knockoff)
library(glmnet)
library(doMC)
registerDoMC(8)

# Outer loop
loop <- function(n = 3000, p, rho, amplitude) {
  # Simulate predictors
  x <- matrix(rnorm(n * p), ncol = p)
  if (rho == 0) {
    Sigma <- diag(p)
  } else {
    Sigma <- toeplitz(rho^(0:(p - 1)))
    x <- x %*% chol(Sigma)
  }
  dimnames(x) <- list(NULL, paste0('x', seq_len(p)))
  # Simulate signal
  k <- 60
  nonzero <- sample(p, k)
  signs <- sample(c(1, -1), size = p, replace = TRUE)
  beta <- amplitude * (seq_len(p) %in% nonzero) / sqrt(n) * signs
  signal <- x %*% beta
  # Gaussian MX knockoff parameters
  mu <- rep(0, p)
  diag_s = create.solve_asdp(Sigma)
  # Inner loop
  test <- function(b, type) {
    # Generate knockoffs
    x_tilde <- create.gaussian(x, mu, Sigma, diag_s = diag_s)
    ### Knockoff filter ###
    if (type == 'regression') {
      y <- signal + rnorm(n)
      w <- stat.glmnet_coefdiff(x, x_tilde, y, family = 'gaussian', cores = 1)
    } else if (type == 'classification') {
      y <- as.factor(rbinom(n, size = 1, prob = plogis(signal)))
      w <- stat.glmnet_coefdiff(x, x_tilde, y, family = 'binomial', cores = 1)
    }
    tau <- knockoff.threshold(w, fdr = 0.1, offset = 0)
    pos <- which(w > tau)
    neg <- which(w <= tau)
    ko_fdr <- sum(beta[pos] == 0) / max(1, length(pos))
    ko_pwr <- 1 - sum(beta[neg] != 0) / p
    ### CPI ###
    # Generate test dataset
    x_test <- matrix(rnorm(n * p), ncol = p)
    if (rho > 0) {
      x_test <- x_test %*% chol(Sigma)
    }
    dimnames(x_test) <- list(NULL, paste0('x', seq_len(p)))
    signal_test <- x_test %*% beta
    # Fit model, compute loss
    if (type == 'regression') {
      f <- cv.glmnet(x, y, family = 'gaussian', nlambda = 500, parallel = FALSE)
      y_test <- signal_test + rnorm(n)
      y_hat <- predict(f, newx = x_test, s = 'lambda.min')
      loss <- (y_test - y_hat)^2
    } else if (type == 'classification') {
      f <- cv.glmnet(x, y, family = 'binomial', nlambda = 500, parallel = FALSE)
      y_test <- rbinom(n, size = 1, prob = plogis(signal_test))
      y_hat <- predict(f, newx = x_test, s = 'lambda.min', type = 'response')
      loss <- -(y_test * log(y_hat) + (1 - y_test) * log(1 - y_hat))
    }
    p_values <- rep(1, p)
    cpi_fn <- function(j) {
      x_test[, j] <- x_tilde[, j]
      if (type == 'regression') {
        y_hat0 <- predict(f, newx = x_test, s = 'lambda.min')
        loss0 <- (y_test - y_hat0)^2
      } else if (type == 'classification') {
        y_hat0 <- predict(f, newx = x_test, s = 'lambda.min', type = 'response')
        loss0 <- -(y_test * log(y_hat0) + (1 - y_test) * log(1 - y_hat0))
      }
      delta <- loss0 - loss
      t_test <- t.test(delta, alternative = 'greater')
      return(t_test$p.value)
    }
    nonzero <- predict(f, s = 'lambda.min', type = 'nonzero')$X1
    p_values[nonzero] <- foreach(j = nonzero, .combine = c) %do% cpi_fn(j)
    q_values <- p.adjust(p_values, method = 'fdr')
    pos <- which(q_values <= 0.1)
    neg <- which(q_values > 0.1)
    cpi_fdr <- sum(beta[pos] == 0) / max(1, length(pos))
    cpi_pwr <- 1 - sum(beta[neg] != 0) / p
    out <- data.table(
         Run = b, 
        type = type,
      method = c('Knockoff', 'CPI'), 
         FDR = c(ko_fdr, cpi_fdr),
       Power = c(ko_pwr, cpi_pwr)
    )
    return(out)
  }
  # Execute inner loop in parallel
  out <- foreach(b = seq_len(200), .combine = rbind) %:%
    foreach(type = c('regression', 'classification'), .combine = rbind) %dopar% 
    test(b, type)
  out[, p := p 
    ][, rho := rho
    ][, amplitude := amplitude]
  return(out)
}

# Execute outer loop in serial, export to RDS
out <- foreach(p = c(1000, 6000), .combine = rbind) %:%
  foreach(rho = seq(from = 0, to = 0.8, by = 0.1), .combine = rbind) %do%
  loop(n = 3000, p = p, rho = rho, amplitude = 10)
saveRDS(out, 'out_rho.rds')
rm(out)
out <- foreach(p = c(1000, 6000), .combine = rbind) %:%
  foreach(amplitude = seq(from = 5, to = 13, by = 1), .combine = rbind) %do%
  loop(n = 3000, p = p, rho = 0, amplitude = amplitude)
saveRDS(out, 'out_amp.rds')


