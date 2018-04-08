




# Test on Boston housing data
library(MASS)
data(Boston)
delta <- c_imp(x = Boston[, -14], y = Boston[, 14])


# What about...
x <- matrix(rnorm(1e4), ncol = 10, 
            dimnames = list(NULL, paste0('x', seq_len(10))))
y1 <- drop(tcrossprod(1:10, x))
y2 <- ifelse(y1 > 0, 1, 0)

# (1) Divide sample-wise loss by sample-wise variance
# This method fails for simple x, y 

# (2) Linear model with precision weights
# Actually does pretty well

# (3) No precision weights
# Gets the rankings right, altho values are a little odd

# Wait also compare OLS and straightup t.test results, cuz the latter's definitely faster
# Just did and it turns out they're equivalent yayyy

# REGRESSION:
# Weighted and unweighted models both performed well (perfect ranking, anyway)
# Adding weights tends to deflate deltas and SEs, ultimately leading to 
# more conservative p-values 

# CLASSIFICATION:
# Deltas and SEs all look super tiny, presumably because the cross entropy
# loss function tends to have a smaller range than L2-loss? Might be worth
# exploring alternative loss functions (squared Pearson or deviance residuals?)
# but once again, rankings are almost perfect in both cases (one mistake in weighted version)
# Again, weights make for more conservative p-values

# Thoughts on weights: not sure what they're really meant to add? We don't typically
# evaluate RF performance with precision weights, so it's unclear why we should
# bring that in here


# At least three options for classification: 
# (1) Straightup L2 loss
# (2) "Pearson residuals", aka L2 loss divided by variance
# (3) L2 loss in linear model with precision weights
loss <- (y - y_hat)^2

loss <- (y - y_hat) / sqrt(y_hat * (1 - y_hat))

pos_loss <- function(y, y_hat) 2 * y * log(y / y_hat)
neg_loss <- function(y, y_hat) 2 * (1 - y) * log((1 - y) / 1 - y_hat)
dev_loss <- sapply(seq_len(length(y)), function(i) {
  ifelse(y[i] == 1, pos_loss(y[i], y_hat[i]), neg_loss(y[i], y_hat[i]))
})

loss <- 2 * y * log(y / y_hat) + 2 * (1 - y) * log((1 - y) / (1 - y_hat))
  
loss0 <- -(d[, 'y'] * log(y_hat0) + (1 - d[, 'y']) * log(1 - y_hat0))

# THIS JUST IN: SQUARED DEVIANCE IS EXACTLY 2 * CROSS-ENTROPY
  
  
# Register cores
library(doMC)
registerDoMC(4)

c_imp <- function(x, y, weights = TRUE,
                  mtry = floor(sqrt(ncol(x))), B = 2000, 
                  n.cores = 4, seed = 123) {
  
  # Define drop_1 function
  drop_1 <- function(j) {
    
    # Drop j
    d0 <- d[, -j]
    
    # Build models, calculate loss (and optionally wts)
      if (weights) {
        f0 <- ranger(data = d0, dependent.variable.name = 'y', 
                     mtry = mtry, num.trees = B, keep.inbag = TRUE, 
                     num.threads = 1, seed = seed, 
                     probability = TRUE) 
        wts0 <- 1 / (predict(f0, data = d0, num.threads = 1,
                             type = 'se')$se[, 1])^2
      } else {
        f0 <- ranger(data = d0, dependent.variable.name = 'y', 
                     mtry = mtry, num.trees = B, 
                     num.threads = 1, seed = seed, 
                     probability = TRUE) 
      }
      y_hat0 <- sapply(f0$predictions, '[[', 1)
      loss0 <- (y_hat0 - y)^2 / (y_hat0 * (1 - y_hat0))
    
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
  
  # Prelimz
  require(ranger)
  n <- nrow(x)
  p <- ncol(x)
  d <- cbind(x, 'y' = y)
  
  # Define test
    if (weights) {
      # Grow full forest
      f <- ranger(data = d, dependent.variable.name = 'y', 
                  mtry = mtry, num.trees = B, keep.inbag = TRUE, 
                  num.threads = n.cores, seed = seed,
                  probability = TRUE)
      # Calculate sample-wise precision
      wts <- 1 / predict(f, data = d, num.threads = n.cores,
                         type = 'se')$se[, 1]^2
    } else {
      # Grow full forest
      f <- ranger(data = d, dependent.variable.name = 'y', 
                  mtry = mtry, num.trees = B, 
                  num.threads = n.cores, seed = seed,
                  probability = TRUE)
    }
    # Calculate sample-wise loss
    y_hat <- sapply(f$predictions, '[[', 1)
    loss <- (y_hat - y)^2 / (y_hat * (1 - y_hat))
  
  # Execute in parallel
  delta <- foreach(j = seq_len(p), .combine = rbind) %dopar% drop_1(j)
  dimnames(delta) <- list(NULL, c('Delta', 'SE', 't', 'p.value'))
  return(data.frame(Feature = colnames(x), delta))
  
}












  
  
  
  
  
