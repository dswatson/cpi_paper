### NULL SIMULATION: Why all the bias? ###

# Load libraries, register cores
library(ranger)
library(matrixStats)
library(tidyverse)
library(doMC)
registerDoMC(4)

# Load CPI
source('cpi.R')

# Simulate
set.seed(1)
n <- 100
p <- 5000
x <- matrix(runif(n * p), ncol = p)
colnames(x) <- paste0('x', seq_len(p))
y <- rnorm(n)

# Build model
mtry <- floor(sqrt(p))
B <- 1e4
df <- data.frame(x, y = y)  
rf <- ranger(data = df, dependent.variable.name = 'y', 
             mtry = mtry, num.trees = B, replace = FALSE,
             keep.inbag = TRUE, num.threads = 4, seed = 1)

# Create oob_preds object
oob_idx <- ifelse(simplify2array(rf$inbag.counts) == 0, TRUE, NA)
preds <- predict(rf, df, num.threads = 4,
                 predict.all = TRUE)$predictions
oob_preds <- oob_idx * preds

# Extract predictions, loss
y_hat <- rf$predictions
loss <- (y_hat - y)^2

# How many trees expected in sub-forest?
splits <- foreach(b = seq_len(B), .combine = c) %dopar% 
  sum(!treeInfo(rf, b)$terminal)
pr_null <- (1 - mtry/p)^mean(splits)

### Random subsampling: 
### If we make predictions using random subsamples of trees and subtract the 
### predictions of the full ensemble, we should find a 0-centered normal
### distribution as n approaches infinity

# Random sub-forests
loop <- function(j) {
  set.seed(j)
  f0_idx <- sample(c(TRUE, FALSE), size = B, replace = TRUE, 
                   prob = c(pr_null, 1 - pr_null))
  y_hat0 <- rowMeans2(oob_preds[, f0_idx, drop = FALSE], na.rm = TRUE)
  return(y_hat0 - y_hat)
}
res <- data.frame(
  Delta_Random = foreach(j = seq_len(p), .combine = c) %dopar% loop(j)
)
ggplot(res, aes(Delta_Random)) + 
  geom_histogram(aes(y = ..density..), bins = 100, color = 'black') + 
  stat_function(fun = dnorm, color = 'red', 
                args = list(mean = 0, sd = sd(res$Delta_Random))) + 
  theme_bw()

### If the forest splitting algorithm is biased, then approximating null
### forests through selective subsetting of the ensemble should produce 
### different results to the random method above

# Drop-one sub-forests
loop <- function(j) {
  f0_idx <- !rf$forest$variable.selected[, j]
  y_hat0 <- rowMeans2(oob_preds[, f0_idx, drop = FALSE], na.rm = TRUE)
  return(y_hat0 - y_hat)
}
res$Delta_Drop <- foreach(j = seq_len(p), .combine = c) %dopar% loop(j)
ggplot(res, aes(Delta_Drop)) + 
  geom_histogram(aes(y = ..density..), bins = 100, color = 'black') + 
  stat_function(fun = dnorm, color = 'red', 
                args = list(mean = 0, sd = sd(res$Delta_Drop))) + 
  theme_bw()

# Specifically, we're curious about each variable's variance
lapply(res, sd)

# And indeed, it appears the drop-one method produces substantially larger
# deltas (in absolute terms) than the random subsampling method

# But these null models aren't just *different* than their random counterparts 
# ...they're worse
risk <- data.frame(
  y_hat0_rdm = res$Delta_Random + rep(y_hat, p),
  y_hat0_drp = res$Delta_Drop + rep(y_hat, p)
) %>%
  mutate(loss0_rdm = (rep(y, p) - y_hat0_rdm)^2,
         loss0_drp = (rep(y, p) - y_hat0_drp)^2,
             batch = rep(seq_len(p), each = n)) %>% 
  group_by(batch) %>%
  summarize(emp_risk_rdm = mean(loss0_rdm),
            emp_risk_drp = mean(loss0_drp))
with(risk, t.test(emp_risk_drp, emp_risk_rdm, alternative = 'greater'))


# Random sub-forests with actual selection frequencies
loop <- function(j) {
  #set.seed(j)
  pp <- 1 - mean(rf$forest$variable.selected[, j])
  f0_idx <- sample(c(TRUE, FALSE), size = B, replace = TRUE, 
                   prob = c(pp, 1 - pp))
  y_hat0 <- rowMeans2(oob_preds[, f0_idx, drop = FALSE], na.rm = TRUE)
  return(y_hat0 - y_hat)
}
res$Delta_Random2 <- foreach(j = seq_len(p), .combine = c) %dopar% loop(j)

# Looks OK now
lapply(res, sd)

# Boxplot
myres <- rbind(data.frame(Method = "Random", value = res$Delta_Random),
               data.frame(Method = "Drop", value = res$Delta_Drop),
               data.frame(Method = "Random2", value = res$Delta_Random2))
ggplot(myres, aes(x = Method, y = value)) + 
  geom_boxplot()