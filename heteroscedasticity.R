
n <- 10000

res <- replicate(1000, {
  a <- rnorm(n, mean = 0, sd = .1)
  b <- rnorm(n, mean = 0, sd = .1)
  c <- rnorm(n, mean = 0, sd = 1)
  
  loss1 <- (a - c)^2
  loss2 <- (b - c)^2
  
  cbind(loss1, loss2)
})

loss1 <- res[, 1, ]
loss2 <- res[, 2, ]
dif <- loss1 - loss2

# Same results as before (but without any model)
plot(loss1[, 1], loss2[, 2])
plot(log(loss1[, 1]), log(loss2[, 1]))

# Heteroscedasticity in the two losses
plot(apply(loss1, 1, mean), apply(loss1, 1, var))
plot(apply(loss2, 1, mean), apply(loss2, 1, var))

# But not in the difference
plot(apply(dif, 1, mean), apply(dif, 1, var))

# Type I error OK
mean(apply(dif, 2, function(x) {t.test(x)$p.value < 0.05}))

# Type I error of log difference (ratio test) about the same
mean(apply(log_dif, 2, function(x) {t.test(x)$p.value < 0.05}))


