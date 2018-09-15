
n <- 10000
p <- 10

res <- replicate(1000, {
  # Training data 1
  x <- matrix(runif(n * p), ncol = p,
              dimnames = list(NULL, paste0('x', seq_len(p))))
  y <- rnorm(n)
  dat_train1 <- data.frame(y = y, x)
  
  # Training data 2
  x <- matrix(runif(n * p), ncol = p,
              dimnames = list(NULL, paste0('x', seq_len(p))))
  y <- rnorm(n)
  dat_train2 <- data.frame(y = y, x)
  
  # Test data
  x <- matrix(runif(n * p), ncol = p,
              dimnames = list(NULL, paste0('x', seq_len(p))))
  y <- rnorm(n)
  dat_test <- data.frame(y = y, x)
  
  # Fit on training data 1, predict on test data
  fit1 <- lm(y ~ ., dat_train1)
  pred1 <- predict(fit1, dat_test)
  loss1 <- (pred1 - dat_test$y)^2
  
  # Fit on training data 2, predict on test data
  fit2 <- lm(y ~ ., dat_train2)
  pred2 <- predict(fit2, dat_test)
  loss2 <- (pred2 - dat_test$y)^2
  
  cbind(loss1, loss2)
})

loss1 <- res[, 1, ]
loss2 <- res[, 2, ]
dif <- loss1 - loss2
log_dif <- log(loss1) - log(loss2)

# Same results as before
plot(loss1[, 1], loss2[, 2])
plot(log(loss1[, 1]), log(loss2[, 1]))

# Heteroscedasticity in the two losses
plot(apply(loss1, 1, mean), apply(loss1, 1, var))
plot(apply(loss2, 1, mean), apply(loss2, 1, var))

# But not in the difference
plot(apply(dif, 1, mean), apply(dif, 1, var))

# Type I error inflated anyway
mean(apply(dif, 2, function(x) {t.test(x)$p.value < 0.05}))

# Type I error of log difference (ratio test) OK
mean(apply(log_dif, 2, function(x) {t.test(x)$p.value < 0.05}))



