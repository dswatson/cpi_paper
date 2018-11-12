# Load libraries, register cores
library(data.table)
library(mlbench)
library(ranger)
library(nnet)
library(e1071)
library(vimp)
library(tidyverse)
library(ggsci)
library(doMC)
registerDoMC(4)

# Set seed
set.seed(123, kind = "L'Ecuyer-CMRG")

# Hyperparameters
n <- 100
p <- 10

# These will come up a lot
models <- c('LM', 'RF', 'NN', 'SVM')

# Loss function
compute_loss <- function(trn, tst) {
  
  # Linear model
  lm_fit <- lm(y ~ ., data = trn)
  lm_yhat <- predict(lm_fit, newdata = tst)
  lm_eps <- tst$y - lm_yhat
  
  # Random forest
  rf_fit <- ranger(y ~ ., data = trn, num.trees = 50)
  rf_yhat <- predict(rf_fit, tst)$predictions
  rf_eps <- tst$y - rf_yhat
  
  # Neural network
  nn_fit <- nnet(y ~ ., data = trn, size = 3, decay = 0.2, 
                 linout = TRUE, trace = FALSE)
  nn_yhat <- predict(nn_fit, newdata = tst)
  nn_eps <- tst$y - nn_yhat
  
  # Support vector machine
  svm_fit <- svm(y ~ ., data = trn, kernel = 'radial', fitted = FALSE)
  svm_yhat <- predict(svm_fit, newdata = tst)
  svm_eps <- tst$y - svm_yhat
  
  # Export
  out <- data.frame(
    lm_yhat, lm_eps, rf_yhat, rf_eps, 
    nn_yhat, nn_eps, svm_yhat, svm_eps
  ) %>%
    mutate(lm_loss = lm_eps^2, rf_loss = rf_eps^2,
           nn_loss = nn_eps^2, svm_loss = svm_eps^2)
  return(out)
  
}

run <- function(sim) {
  
  # Simulate data
  if (sim == 'friedman') {
    dat <- mlbench.friedman1(3 * n)
    dat <- dat0 <- data.frame(dat$x, y = dat$y)
  } else if (sim == 'linear') {
    beta <- c(seq_len(5), rep(0, times = 5))
    x <- matrix(rnorm(3 * n * 10), ncol = 10, 
                dimnames = list(NULL, paste0('x', seq_len(10))))
    y <- x %*% beta + rnorm(3 * n)
    dat <- dat0 <- data.frame(x, y)
  }
  
  # Split into training and test sets
  idx <- seq_len(2 * n)
  train <- dat[idx, ]
  test <- dat[-idx, ]
  
  # Fit and evaluate original models
  orig <- compute_loss(train, test)
  
  ### CPI ###
  cpi_fn <- function(j) {
    dat0[, j] <- dat[sample.int(3 * n), j]
    train0 <- dat0[idx, ]
    test0 <- dat0[-idx, ]
    null <- compute_loss(train0, test0)
    delta <- data.table(
      'LM' = null$lm_loss - orig$lm_loss,
      'RF' = null$rf_loss - orig$rf_loss,
      'NN' = null$nn_loss - orig$nn_loss,
      'SVM' = null$svm_loss - orig$svm_loss
    )
    lambda <- data.table(
      'LM' = log(null$lm_loss / orig$lm_loss),
      'RF' = log(null$rf_loss / orig$rf_loss),
      'NN' = log(null$nn_loss / orig$nn_loss),
      'SVM' = log(null$svm_loss / orig$svm_loss)
    )
    out <- data.table(
      'Model' = rep(models, times = 2),
      'Feature' = j,
      'VIM' = rep(c('cpi_d', 'cpi_l'), each = length(models)),
      'VI' = c(
        mean(delta$LM), mean(delta$RF), mean(delta$NN), mean(delta$SVM),
        mean(lambda$LM), mean(lambda$RF), mean(lambda$NN), mean(lambda$SVM)
      ),
      'p.value' = c(
        t.test(delta$LM, alternative = 'greater')$p.value,
        t.test(delta$RF, alternative = 'greater')$p.value,
        t.test(delta$NN, alternative = 'greater')$p.value,
        t.test(delta$SVM, alternative = 'greater')$p.value,
        t.test(lambda$LM, alternative = 'greater')$p.value,
        t.test(lambda$RF, alternative = 'greater')$p.value,
        t.test(lambda$NN, alternative = 'greater')$p.value,
        t.test(lambda$SVM, alternative = 'greater')$p.value
      )
    )
    return(out)
  }
  cpi_df <- foreach(j = seq_len(p), .combine = rbind) %do% cpi_fn(j)
  
  ### Williamson et al.'s nonparametric ANOVA ###
  anova_fn <- function(j) {
    train0 <- test[, -c(j, 11)]
    # Linear model
    train0$y <- orig$lm_yhat
    lm_fit0 <- lm(y ~ ., data = train0)
    lm_yhat0 <- predict(lm_fit0, newdata = test)
    vimp_lm <- vimp_regression(test$y, f1 = orig$lm_yhat, f2 = lm_yhat0, indx = j, 
                               run_regression = FALSE)
    # Random forest
    train0$y <- orig$rf_yhat
    rf_fit0 <- ranger(y ~ ., data = train0, num.trees = 50)
    rf_yhat0 <- predict(rf_fit0, data = test)$predictions
    vimp_rf <- vimp_regression(test$y, f1 = orig$rf_yhat, f2 = rf_yhat0, indx = j, 
                               run_regression = FALSE)
    # Neural network
    train0$y <- orig$nn_yhat
    nn_fit0 <- nnet(y ~ ., data = train0, size = 3, decay = 0.2, 
                    linout = TRUE, trace = FALSE)
    nn_yhat0 <- predict(nn_fit0, newdata = test)
    vimp_nn <- vimp_regression(test$y, f1 = orig$nn_yhat, f2 = nn_yhat0, indx = j, 
                               run_regression = FALSE)
    # SVM
    train0$y <- orig$svm_yhat
    svm_fit0 <- svm(y ~ ., data = train0, kernel = 'radial', fitted = FALSE)
    svm_yhat0 <- predict(svm_fit0, newdata = test)
    vimp_svm <- vimp_regression(test$y, f1 = orig$svm_yhat, f2 = svm_yhat0, indx = j, 
                                run_regression = FALSE)
    # Compute nonparametric ANOVA for each model
    data.table(
      'Model' = models,
      'Feature' = j,
      'VIM' = 'anova',
      'VI' = c(
        vimp_lm$est, vimp_rf$est, vimp_nn$est, vimp_svm$est
      ),
      'p.value' = c(
        pnorm(vimp_lm$est / vimp_lm$se, lower = FALSE),
        pnorm(vimp_rf$est / vimp_rf$se, lower = FALSE),
        pnorm(vimp_nn$est / vimp_nn$se, lower = FALSE),
        pnorm(vimp_svm$est / vimp_svm$se, lower = FALSE)
      )
    )
  }
  anova_df <- foreach(j = seq_len(p), .combine = rbind) %do% anova_fn(j)
  
  ### Chalupka et al.'s nested models (delta and lambda),           ###
  ### Ramsey's conditional correlation independence (CCI) test, and ###
  ### Shah & Peters's generalised covariance measure (GCM) test     ###
  other_fn <- function(j) {
    f <- compute_loss(train[, -j], test[, -j])
    delta <- data.table(
      'LM' = f$lm_loss - orig$lm_loss,
      'RF' = f$rf_loss - orig$rf_loss,
      'NN' = f$nn_loss - orig$nn_loss,
      'SVM' = f$svm_loss - orig$svm_loss
    )
    lambda <- data.table(
      'LM' = log(f$lm_loss / orig$lm_loss),
      'RF' = log(f$rf_loss / orig$rf_loss),
      'NN' = log(f$nn_loss / orig$nn_loss),
      'SVM' = log(f$svm_loss / orig$svm_loss)
    )
    train0 <- select(train, -y)
    colnames(train0)[j] <- 'y'
    test0 <- select(test, -y)
    colnames(test0)[j] <- 'y'
    g <- compute_loss(train0, test0)
    cor_lm <- cor.test(f$lm_eps, g$lm_eps)
    cor_rf <- cor.test(f$rf_eps, g$rf_eps)
    cor_nn <- cor.test(f$nn_eps, g$nn_eps)
    cor_svm <- cor.test(f$svm_eps, g$svm_eps)
    r_lm <- f$lm_eps * g$lm_eps
    r_rf <- f$rf_eps * g$rf_eps
    r_nn <- f$nn_eps * g$nn_eps
    r_svm <- f$svm_eps * g$svm_eps
    gcm_lm <- abs((sqrt(n) * mean(r_lm)) / sqrt(mean(r_lm^2) - (mean(r_lm))^2))
    gcm_rf <- abs((sqrt(n) * mean(r_rf)) / sqrt(mean(r_rf^2) - (mean(r_rf))^2))
    gcm_nn <- abs((sqrt(n) * mean(r_nn)) / sqrt(mean(r_nn^2) - (mean(r_nn))^2))
    gcm_svm <- abs((sqrt(n) * mean(r_svm)) / sqrt(mean(r_svm^2) - (mean(r_svm))^2))
    data.table(
      'Model' = rep(models, times = 4),
      'Feature' = j,
      'VIM' = rep(c('nested_d', 'nested_l', 'cci', 'gcm'), 
                  each = length(models)),
      'VI' = c(
        mean(delta$LM), mean(delta$RF), mean(delta$NN), mean(delta$SVM),
        mean(lambda$LM), mean(lambda$RF), mean(lambda$NN), mean(lambda$SVM),
        cor_lm$estimate^2, cor_rf$estimate^2, 
        cor_nn$estimate^2, cor_svm$estimate^2, 
        gcm_lm, gcm_rf, gcm_nn, gcm_svm
      ),
      'p.value' = c(
        t.test(delta$LM, alternative = 'greater')$p.value,
        t.test(delta$RF, alternative = 'greater')$p.value,
        t.test(delta$NN, alternative = 'greater')$p.value,
        t.test(delta$SVM, alternative = 'greater')$p.value,
        t.test(lambda$LM, alternative = 'greater')$p.value,
        t.test(lambda$RF, alternative = 'greater')$p.value,
        t.test(lambda$NN, alternative = 'greater')$p.value,
        t.test(lambda$SVM, alternative = 'greater')$p.value,
        cor_lm$p.value, cor_rf$p.value, 
        cor_nn$p.value, cor_svm$p.value,
        2 * pnorm(gcm_lm, lower = FALSE), 2 * pnorm(gcm_rf, lower = FALSE),
        2 * pnorm(gcm_nn, lower = FALSE), 2 * pnorm(gcm_svm, lower = FALSE)
      )
    )
  }
  other_df <- foreach(j = seq_len(p), .combine = rbind) %do% other_fn(j)
  
  # Put it all together
  out <- rbind(cpi_df, anova_df, other_df)
  out[Model == 'LM', MSE := mean(orig$lm_loss)
    ][Model == 'RF', MSE := mean(orig$rf_loss)
    ][Model == 'NN', MSE := mean(orig$nn_loss)
    ][Model == 'SVM', MSE := mean(orig$svm_loss)
    ][, Simulation := sim]
  return(out)

}

# Analysis loop
loop <- function(b) {
  
  fri <- run('friedman')
  lin <- run('linear')
  out <- rbind(fri, lin)
  out[, Run := b]
  return(out)
  
}

# Execute in parallel
out <- foreach(b = seq_len(1e4), .combine = rbind) %dopar% loop(b)
saveRDS(out, 'comp_sim.rds')

# Empirical error rates
FPR <- out[Feature %in% 6:10, sum(Positive) / 5e4, by = .(Model, VIM)]$V1
TPR <- out[Feature %in% 1:5, sum(Positive) / 5e4, by = .(Model, VIM)]$V1
FNR <- out[Feature %in% 1:5, sum(!Positive) / 5e4, by = .(Model, VIM)]$V1
TNR <- out[Feature %in% 6:10, sum(!Positive) / 5e4, by = .(Model, VIM)]$V1
test_pos <- out[, sum(Positive), by = .(Model, VIM)]$V1
FDR <- out[Feature %in% 6:10, sum(Positive), by = .(Model, VIM)]$V1 / test_pos
test_neg <- out[, sum(!Positive), by = .(Model, VIM)]$V1
FOR <- out[Feature %in% 1:5, sum(!Positive), by = .(Model, VIM)]$V1 / test_neg
df <- data.table(
  Model = c('Linear Model', 'Random Forest', 
            'Neural Network', 'Support Vector Machine'),
  VIM = rep(c('delta', 'lambda', 'anova', 'lrt', 'cci', 'gcm'),
            each = length(models)),  # Check order
  FPR, TPR, FNR, TNR, FDR
)
ggplot(df, aes(FPR, TPR, color = VIM, shape = VIM)) +
  geom_point(size = 3) + 
  scale_color_d3() +
  geom_vline(xintercept = 0.05) + 
  labs(x = 'False Positive Rate', y = 'True Positive Rate') + 
  xlim(0, 1) + ylim(0, 1) + 
  theme_bw() + 
  facet_wrap(Simulation ~ Model)
ggplot(df, aes(FPR, 1 - FNR, color = VIM, shape = VIM)) +
  geom_point(size = 3) + 
  scale_color_d3() +
  geom_vline(xintercept = 0.05) + 
  labs(x = 'False Positive Rate', y = 'Power') + 
  xlim(0, 1) + ylim(0, 1) + 
  theme_bw() + 
  facet_wrap(Simulation ~ Model)


