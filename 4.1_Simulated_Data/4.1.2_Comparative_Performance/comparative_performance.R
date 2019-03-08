# Load libraries, register cores
library(data.table)
library(knockoff)
library(ranger)
library(nnet)
library(e1071)
library(vimp)
library(ggsci)
library(tidyverse)
library(doMC)
registerDoMC(20)

# Set seed
set.seed(123, kind = "L'Ecuyer-CMRG")

# Hyperparameters
p <- 10

# These will come up a lot
models <- c('LM', 'RF', 'NN', 'SVM')

loop <- function(b, rho, sim, test_n) {
  
  # Simulate data
  beta <- c(0, 0, -0.5, 0.5, -1, 1, -1.5, 1.5, -2, 2)
  Sigma <- toeplitz(rho^(0:(p - 1)))
  x <- matrix(rnorm(3 * test_n * p), ncol = p) %*% chol(Sigma)
  dimnames(x) <- list(NULL, paste0('x', seq_len(p)))
  if (sim == 'linear') {
    y <- x %*% beta + rnorm(3 * test_n)
  } else if (sim == 'nonlinear') {
    idx <- x < -qnorm(0.75) | x > qnorm(0.75)
    xx <- matrix(0, nrow = 3 * test_n, ncol = p)
    xx[idx] <- 0
    xx[!idx] <- 1
    y <- xx %*% beta + rnorm(3 * test_n)
  }
  dat <- dat0 <- data.frame(x, y)
  
  # Split into training and test sets
  idx <- seq_len(2 * test_n)
  train <- dat[idx, ]
  test <- dat[-idx, ]
  
  # Fit models
  lm_fit <- lm(y ~ ., data = train)
  rf_fit <- ranger(y ~ ., data = train, num.trees = 50)
  nn_fit <- nnet(y ~ ., data = train, size = 20, decay = 0.1, 
                 linout = TRUE, trace = FALSE)
  svm_fit <- svm(y ~ ., data = train, kernel = 'radial', fitted = FALSE)
  
  # Compute predictions, residuals, and squared errors
  lm_yhat <- predict(lm_fit, newdata = test)
  lm_eps <- test$y - lm_yhat
  lm_loss <- lm_eps^2
  rf_yhat <- predict(rf_fit, test)$predictions
  rf_eps <- test$y - rf_yhat
  rf_loss <- rf_eps^2
  nn_yhat <- predict(nn_fit, newdata = test)
  nn_eps <- test$y - nn_yhat
  nn_loss <- nn_eps^2
  svm_yhat <- predict(svm_fit, newdata = test)
  svm_eps <- test$y - svm_yhat
  svm_loss <- svm_eps^2
  
  # Generate knockoff matrix
  x_tilde <- create.second_order(x)

  ### CPI ###
  cpi_fn <- function(j) {
    dat0[, j] <- x_tilde[, j]
    test0 <- dat0[-idx, ]
    lm_yhat0 <- predict(lm_fit, newdata = test0)
    rf_yhat0 <- predict(rf_fit, test0)$predictions
    nn_yhat0 <- predict(nn_fit, newdata = test0)
    svm_yhat0 <- predict(svm_fit, newdata = test0)
    lm_loss0 <- (test$y - lm_yhat0)^2
    rf_loss0 <- (test$y - rf_yhat0)^2
    nn_loss0 <- (test$y - nn_yhat0)^2
    svm_loss0 <- (test$y - svm_yhat0)^2
    delta <- data.table(
      'LM' = lm_loss0 - lm_loss, 
      'RF' = rf_loss0 - rf_loss,
      'NN' = nn_loss0 - nn_loss, 
      'SVM' = svm_loss0 - svm_loss
    )
    lambda <- data.table(
      'LM' = log(lm_loss0 / lm_loss),
      'RF' = log(rf_loss0 / rf_loss),
      'NN' = log(nn_loss0 / nn_loss),
      'SVM' = log(svm_loss0 / svm_loss)
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
    train0$y <- lm_yhat
    lm_fit0 <- lm(y ~ ., data = train0)
    lm_yhat0 <- predict(lm_fit0, newdata = test)
    vimp_lm <- vimp_regression(test$y, f1 = lm_yhat, f2 = lm_yhat0, indx = j, 
                               run_regression = FALSE)
    # Random forest
    train0$y <- rf_yhat
    rf_fit0 <- ranger(y ~ ., data = train0, num.trees = 50)
    rf_yhat0 <- predict(rf_fit0, data = test)$predictions
    vimp_rf <- vimp_regression(test$y, f1 = rf_yhat, f2 = rf_yhat0, indx = j, 
                               run_regression = FALSE)
    # Neural network
    train0$y <- nn_yhat
    nn_fit0 <- nnet(y ~ ., data = train0, size = 20, decay = 0.1, 
                    linout = TRUE, trace = FALSE)
    nn_yhat0 <- predict(nn_fit0, newdata = test)
    vimp_nn <- vimp_regression(test$y, f1 = nn_yhat, f2 = nn_yhat0, indx = j, 
                               run_regression = FALSE)
    # SVM
    train0$y <- svm_yhat
    svm_fit0 <- svm(y ~ ., data = train0, kernel = 'radial', fitted = FALSE)
    svm_yhat0 <- predict(svm_fit0, newdata = test)
    vimp_svm <- vimp_regression(test$y, f1 = svm_yhat, f2 = svm_yhat0, indx = j, 
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
    # Fit models
    lm_fit0 <- lm(y ~ ., data = train[, -j])
    rf_fit0 <- ranger(y ~ ., data = train[, -j])
    nn_fit0 <- nnet(y ~ ., data = train[, -j], size = 20, decay = 0.1, 
                   linout = TRUE, trace = FALSE)
    svm_fit0 <- svm(y ~ ., data = train[, -j], kernel = 'radial', 
                    fitted = FALSE)
    # Compute predictions, residuals, and squared errors
    lm_yhat0 <- predict(lm_fit0, newdata = test[, -j])
    lm_eps0 <- test$y - lm_yhat0
    lm_loss0 <- lm_eps0^2
    rf_yhat0 <- predict(rf_fit0, test[, -j])$predictions
    rf_eps0 <- test$y - rf_yhat0
    rf_loss0 <- rf_eps0^2
    nn_yhat0 <- predict(nn_fit0, newdata = test[, -j])
    nn_eps0 <- test$y - nn_yhat0
    nn_loss0 <- nn_eps0^2
    svm_yhat0 <- predict(svm_fit0, newdata = test[, -j])
    svm_eps0 <- test$y - svm_yhat0
    svm_loss0 <- svm_eps0^2
    # FCIT
    delta <- data.table(
      'LM' = lm_loss0 - lm_loss,
      'RF' = rf_loss0 - rf_loss,
      'NN' = nn_loss0 - nn_loss,
      'SVM' = svm_loss0 - svm_loss
    )
    # GCM
    train0 <- select(train, -y)
    colnames(train0)[j] <- 'y'
    test0 <- select(test, -y)
    colnames(test0)[j] <- 'y'
    # Fit models
    lm_fit_g <- lm(y ~ ., data = train0)
    rf_fit_g <- ranger(y ~ ., data = train0)
    nn_fit_g <- nnet(y ~ ., data = train0, size = 20, decay = 0.1, 
                    linout = TRUE, trace = FALSE)
    svm_fit_g <- svm(y ~ ., data = train0, kernel = 'radial', fitted = FALSE)
    # Compute predictions, residuals, and squared errors
    lm_yhat_g <- predict(lm_fit_g, newdata = test0)
    lm_eps_g <- test$y - lm_yhat_g
    lm_loss_g <- lm_eps_g^2
    rf_yhat_g <- predict(rf_fit_g, test0)$predictions
    rf_eps_g <- test$y - rf_yhat_g
    rf_loss_g <- rf_eps_g^2
    nn_yhat_g <- predict(nn_fit_g, newdata = test0)
    nn_eps_g <- test$y - nn_yhat_g
    nn_loss_g <- nn_eps_g^2
    svm_yhat_g <- predict(svm_fit_g, newdata = test0)
    svm_eps_g <- test$y - svm_yhat_g
    svm_loss_g <- svm_eps_g^2
    r_lm <- lm_eps0 * lm_eps_g 
    r_rf <- rf_eps0 * rf_eps_g
    r_nn <- nn_eps0 * nn_eps_g
    r_svm <- svm_eps0 * svm_eps_g
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
    ][, rho := rho
    ][, Run := b
    ][, Size := paste0('n = ', 2 * test_n)]
  return(out)

}

# Execute in parallel
out <- foreach(b = seq_len(1e4), .combine = rbind) %:%
  foreach(rho = c(0, 0.5)) %:%
  foreach(sim = c('linear', 'nonlinear'), .combine = rbind) %:%
  foreach(test_n = c(50, 250, 500), .combine = rbind) %dopar% 
  loop(b, rho, sim, test_n)
saveRDS(out, 'comp_big_sim.rds')

# Empirical error rates
out[, Rejected := p.value <= 0.05
  ][Feature %in% 1:2, EffSize := 0
  ][Feature %in% 3:4, EffSize := 0.5
  ][Feature %in% 5:6, EffSize := 1
  ][Feature %in% 7:8, EffSize := 1.5
  ][Feature %in% 9:10, EffSize := 2]
df <- out[, .(RejRate = mean(Rejected)), 
          by = .(Simulation, Model, EffSize, VIM, Size)]
df[Simulation == 'linear', Simulation := 'Linear data'
  ][Simulation == 'nonlinear', Simulation := 'Nonlinear data'
  ][Model == 'LM', Model := 'Linear model'
  ][Model == 'RF', Model := 'Random forest'
  ][Model == 'NN', Model := 'Neural network'
  ][Model == 'SVM', Model := 'Support vector machine'
  ][, Model := factor(Model, levels = c('Linear model', 'Support vector machine',
                                        'Random forest', 'Neural network'))]

# Plot results
lbls <- c(bquote('CPI'[Delta]), bquote('CPI'[lambda]),
          'ANOVA', 'FCIT', 'GCM')
plot_fn <- function(n) {
  if (n == 100) {
    df <- df %>% filter(Size == 'n = 100')
    title <- 'n = 100'
  } else if (n == 500) {
    df <- df %>% filter(Size == 'n = 500')
    title <- 'n = 500'
  } else if (n == 1000) {
    df <- df %>% filter(Size == 'n = 1000')
    title <- 'n = 1000'
  }
  ggplot(df, aes(EffSize, RejRate, group = VIM, color = VIM, shape = VIM)) +
    geom_point() + 
    geom_line() + 
    geom_hline(yintercept = 0.05, linetype = 'dashed') +
    scale_y_continuous(breaks = c(0, .05, .25, .5, .75, 1), limits = c(0, 1)) + 
    labs(x = 'Effect size', 
         y = 'Rejection proportion') + 
    theme_bw() + 
    theme(legend.position = 'bottom') + 
    scale_color_manual(name = 'CI Test',
                       breaks = unique(df$VIM),
                       labels = lbls,
                       values = pal_npg()(5)) +
    facet_grid(Simulation ~ Model)
}


