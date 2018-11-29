# Load libraries, register cores
library(data.table)
library(ranger)
library(nnet)
library(e1071)
library(vimp)
library(ggsci)
library(tidyverse)
library(doMC)
registerDoMC(8)

# Set seed
set.seed(123, kind = "L'Ecuyer-CMRG")

# Hyperparameters
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
  nn_fit <- nnet(y ~ ., data = trn, size = 20, decay = 0.1, 
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

loop <- function(b, sim, test_n) {
  
  # Simulate data
  beta <- c(0, 0, -0.5, 0.5, -1, 1, -1.5, 1.5, -2, 2)
  x <- matrix(runif(3 * test_n * p), ncol = p,
              dimnames = list(NULL, paste0('x', seq_len(p))))
  if (sim == 'linear') {
    y <- x %*% beta + rnorm(3 * test_n)
    dat <- dat0 <- data.frame(x, y)
  } else if (sim == 'nonlinear') {
    idx <- x < .25 | x > .75
    xx <- matrix(0, nrow = 3 * test_n, ncol = p)
    xx[idx] <- 0
    xx[!idx] <- 1
    y <- xx %*% beta + rnorm(3 * test_n)
    dat <- dat0 <- data.frame(x, y)
  }
  
  # Split into training and test sets
  idx <- seq_len(2 * test_n)
  train <- dat[idx, ]
  test <- dat[-idx, ]
  
  # Fit and evaluate original models
  orig <- compute_loss(train, test)
  
  ### CPI ###
  cpi_fn <- function(j) {
    dat0[, j] <- dat[sample.int(3 * test_n), j]
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
      'Simulation' = sim,
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
    train0 <- test[, -j]
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
    nn_fit0 <- nnet(y ~ ., data = train0, size = 20, decay = 0.1, 
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
      'Simulation' = sim,
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
  
  ### Chalupka et al.'s fast conditional independence test (FCIT) and ###
  ### Shah & Peters's generalised covariance measure (GCM)            ###
  other_fn <- function(j) {
    f <- compute_loss(train[, -j], test[, -j])
    delta <- data.table(
      'LM' = f$lm_loss - orig$lm_loss,
      'RF' = f$rf_loss - orig$rf_loss,
      'NN' = f$nn_loss - orig$nn_loss,
      'SVM' = f$svm_loss - orig$svm_loss
    )
    train0 <- select(train, -y)
    colnames(train0)[j] <- 'y'
    test0 <- select(test, -y)
    colnames(test0)[j] <- 'y'
    g <- compute_loss(train0, test0)
    r_lm <- f$lm_eps * g$lm_eps
    r_rf <- f$rf_eps * g$rf_eps
    r_nn <- f$nn_eps * g$nn_eps
    r_svm <- f$svm_eps * g$svm_eps
    gcm_lm <- abs((sqrt(test_n) * mean(r_lm)) / 
                    sqrt(mean(r_lm^2) - (mean(r_lm))^2))
    gcm_rf <- abs((sqrt(test_n) * mean(r_rf)) / 
                    sqrt(mean(r_rf^2) - (mean(r_rf))^2))
    gcm_nn <- abs((sqrt(test_n) * mean(r_nn)) / 
                    sqrt(mean(r_nn^2) - (mean(r_nn))^2))
    gcm_svm <- abs((sqrt(test_n) * mean(r_svm)) / 
                     sqrt(mean(r_svm^2) - (mean(r_svm))^2))
    data.table(
      'Simulation' = sim,
      'Model' = rep(models, times = 2),
      'Feature' = j,
      'VIM' = rep(c('fcit', 'gcm'), each = length(models)),
      'VI' = c(
        mean(delta$LM), mean(delta$RF), mean(delta$NN), mean(delta$SVM),
        gcm_lm, gcm_rf, gcm_nn, gcm_svm
      ),
      'p.value' = c(
        t.test(delta$LM, alternative = 'greater')$p.value,
        t.test(delta$RF, alternative = 'greater')$p.value,
        t.test(delta$NN, alternative = 'greater')$p.value,
        t.test(delta$SVM, alternative = 'greater')$p.value,
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
    ][, Run := b
    ][, Size := paste0('n = ', 2 * test_n)]
  return(out)

}

# Execute in parallel
out <- foreach(b = seq_len(1e4), .combine = rbind) %:%
  foreach(sim = c('linear', 'nonlinear')) %:%
  foreach(test_n = c(50, 250, 500), .combine = rbind) %dopar% 
  loop(b, sim, test_n)
saveRDS(out, 'comp_big_sim.rds')














# Empirical error rates
out[, Positive := p.value <= 0.05]
FPR <- out[Feature %in% 6:10, sum(Positive) / 5e4, 
           by = .(Model, VIM, Size)]$V1
TPR <- out[Feature %in% 1:5, sum(Positive) / 5e4, 
           by = .(Model, VIM, Size)]$V1
test_pos <- out[, sum(Positive), 
                by = .(Model, VIM, Size)]$V1
FDR <- out[Feature %in% 6:10, sum(Positive), 
           by = .(Model, VIM, Size)]$V1 / test_pos
test_neg <- out[, sum(!Positive), 
                by = .(Model, VIM, Size)]$V1
FOR <- out[Feature %in% 1:5, sum(!Positive), 
           by = .(Model, VIM, Size)]$V1 / test_neg
df <- data.table(
  Model = c('Linear Model', 'Random Forest', 
            'Neural Network', 'Support Vector Machine'),
  VIM = rep(out[, unique(VIM)], each = length(models)),  
  Size = factor(rep(c('n = 50', 'n = 100', 'n = 200'), 
                    each = length(models) * out[, length(unique(VIM))]),
                levels = c('n = 50', 'n = 100', 'n = 200')),
  FPR, TPR
)

# Plot results
lbls <- c(bquote('CPI'[Delta]), bquote('CPI'[lambda]),
          'ANOVA', 'FCIT', 'GCM')
ggplot(df, aes(FPR, TPR, color = VIM, shape = VIM)) +
  geom_point(size = 3) + 
  scale_color_manual(name = 'CI Test',
                     breaks = unique(df$VIM),
                     labels = lbls,
                     values = pal_d3()(5)) +
  scale_shape_manual(name = 'CI Test',
                     breaks = unique(df$VIM),
                     labels = lbls, 
                     values = c(16:18, 3:5)) +
  geom_vline(xintercept = 0.05) + 
  labs(x = 'False Positive Rate', y = 'True Positive Rate') + 
  xlim(0, 1) + ylim(0, 1) + 
  theme_bw() + 
  facet_grid(Size ~ Model)



