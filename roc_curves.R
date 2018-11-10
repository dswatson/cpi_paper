# Load libraries, register cores
library(data.table)
library(mlbench)
library(ranger)
library(nnet)
library(e1071)
library(vimp)
library(precrec)
library(tidyverse)
library(gridExtra)
library(ggsci)
library(doMC)
registerDoMC(4)

# Set seed
set.seed(123, kind = "L'Ecuyer-CMRG")

# Friedman1 benchmark hyperparameters
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

# Analysis loop
loop <- function(b) {
  
  # Simulate data
  dat <- mlbench.friedman1(3 * n)
  dat <- dat0 <- data.frame(dat$x, y = dat$y)
  
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
    out <- data.table(
      'Model' = models,
      'Feature' = j,
      'VIM' = rep(c('delta', 'lambda'), each = length(models)),
      'VI' = c(
        mean(null$lm_loss - orig$lm_loss),
        mean(null$rf_loss - orig$rf_loss),
        mean(null$nn_loss - orig$nn_loss),
        mean(null$svm_loss - orig$svm_loss),
        mean(log(null$lm_loss / orig$lm_loss)),
        mean(log(null$rf_loss / orig$rf_loss)),
        mean(log(null$nn_loss / orig$nn_loss)),
        mean(log(null$svm_loss / orig$svm_loss))
      )
    )
    return(out)
  }
  cpi_df <- foreach(j = seq_len(p), .combine = rbind) %do% cpi_fn(j)
  
  ### Chalupka et al.'s LRT ###
  lrt_fn <- function(j) {
    dat0 <- dat[, -j]
    train0 <- dat0[idx, ]
    test0 <- dat0[-idx, ]
    null <- compute_loss(train0, test0)
    out <- data.table(
      'Model' = models,
      'Feature' = j,
      'VIM' = 'lrt',
      'VI' = c(
        mean(null$lm_loss - orig$lm_loss),
        mean(null$rf_loss - orig$rf_loss),
        mean(null$nn_loss - orig$nn_loss),
        mean(null$svm_loss - orig$svm_loss)
      )
    )
    return(out)
  }
  lrt_df <- foreach(j = seq_len(p), .combine = rbind) %do% lrt_fn(j)
  
  ### Williamson et al.'s nonparametric ANOVA ###
  anova_fn <- function(j) {
    train0 <- test[, -c(j, 11)]
    # Linear model
    train0$y <- orig$lm_yhat
    lm_fit0 <- lm(y ~ ., data = train0)
    lm_yhat0 <- predict(lm_fit0, newdata = test)
    # Random forest
    train0$y <- orig$rf_yhat
    rf_fit0 <- ranger(y ~ ., data = train0, num.trees = 50)
    rf_yhat0 <- predict(rf_fit0, data = test)$predictions
    # Neural network
    train0$y <- orig$nn_yhat
    nn_fit0 <- nnet(y ~ ., data = train0, size = 3, decay = 0.2, 
                    linout = TRUE, trace = FALSE)
    nn_yhat0 <- predict(nn_fit0, newdata = test)
    # SVM
    train0$y <- orig$svm_yhat
    svm_fit0 <- svm(y ~ ., data = train0, kernel = 'radial', fitted = FALSE)
    svm_yhat0 <- predict(svm_fit0, newdata = test)
    # Compute nonparametric ANOVA for each model
    data.table(
      'Model' = models,
      'Feature' = j,
      'VIM' = 'anova',
      'VI' = c(
        vimp_regression(test$y, f1 = orig$lm_yhat, f2 = lm_yhat0, indx = j, 
                        run_regression = FALSE)$est,
        vimp_regression(test$y, f1 = orig$rf_yhat, f2 = rf_yhat0, indx = j, 
                        run_regression = FALSE)$est,
        vimp_regression(test$y, f1 = orig$nn_yhat, f2 = nn_yhat0, indx = j, 
                        run_regression = FALSE)$est,
        vimp_regression(test$y, f1 = orig$svm_yhat, f2 = svm_yhat0, indx = j, 
                        run_regression = FALSE)$est
      )
    )
  }
  anova_df <- foreach(j = seq_len(p), .combine = rbind) %do% anova_fn(j)
  
  ### Ramsey's conditional correlation independence (CCI) test and ###
  ### Shah & Peter's generalised covariance measure (GCM) test     ###
  pcor_fn <- function(j) {
    f <- compute_loss(train[, -j], test[, -j])
    train0 <- select(train, -y)
    colnames(train0)[j] <- 'y'
    test0 <- select(test, -y)
    colnames(test0)[j] <- 'y'
    g <- compute_loss(train0, test0)
    r_lm <- f$lm_eps * g$lm_eps
    r_rf <- f$rf_eps * g$rf_eps
    r_nn <- f$nn_eps * g$nn_eps
    r_svm <- f$svm_eps * g$svm_eps
    data.table(
      'Model' = rep(models, times = 2),
      'Feature' = j,
      'VIM' = rep(c('cci', 'gcm'), each = 4),
      'VI' = c(
        cor(f$lm_eps, g$lm_eps)^2, cor(f$rf_eps, g$rf_eps)^2,
        cor(f$nn_eps, g$nn_eps)^2, cor(f$svm_eps, g$svm_eps)^2,
        abs((sqrt(n) * mean(r_lm)) / sqrt(mean(r_lm^2) - (mean(r_lm))^2)),
        abs((sqrt(n) * mean(r_rf)) / sqrt(mean(r_rf^2) - (mean(r_rf))^2)),
        abs((sqrt(n) * mean(r_nn)) / sqrt(mean(r_nn^2) - (mean(r_nn))^2)),
        abs((sqrt(n) * mean(r_svm)) / sqrt(mean(r_svm^2) - (mean(r_svm))^2))
      )
    )
  }
  pcor_df <- foreach(j = seq_len(p), .combine = rbind) %do% pcor_fn(j)
  
  # Put it all together
  out <- rbind(cpi_df, anova_df, lrt_df, pcor_df)
  out[, Run := b]
  return(out)
  
}

# Execute in parallel
out <- foreach(b = seq_len(100), .combine = rbind) %dopar% loop(b)

# Prepare to plot
obs <- rep(rep(c(1, 0), each = 5), times = 100)
auroc <- function(model) {
  pred <- list(CPI_d = out[VIM == 'delta' & Model == model, VI],
               CPI_l = out[VIM == 'lambda' & Model == model, VI],
               ANOVA = out[VIM == 'anova' & Model == model, VI],
               LRT = out[VIM == 'lrt' & Model == model, VI],
               CCI = out[VIM == 'cci' & Model == model, VI],
               GCM = out[VIM == 'gcm' & Model == model, VI])
  rocs <- evalmod(scores = pred, labels = obs)$rocs
  df <- seq_along(pred) %>%
    map_df(~ data_frame(FPR = rocs[[.x]]$x,
                        TPR = rocs[[.x]]$y,
                        VIM = names(pred)[.x],
                        AUC = rocs[[.x]] %>% attr('auc'))) %>%
    mutate(Model = model)
  return(df)
}
df <- foreach(m = models, .combine = rbind) %dopar% auroc(m)

# Plot away
plot_fn <- function(model) {
  auc_vals <- df %>% 
    filter(Model == model) %>%
    select(AUC) %>%
    unique(.) %>% 
    unlist(.) %>%
    as.numeric(.) %>%
    round(., 3)
  lbls <- c(bquote('CPI'[Delta]*', AUC = '~.(auc_vals[1])),
            bquote('CPI'[lambda]*', AUC = '~.(auc_vals[2])),
            paste0('ANOVA, AUC = ', auc_vals[3]),
            paste0('LRT, AUC = ', auc_vals[4]),
            paste0('CCI, AUC = ', auc_vals[5]),
            paste0('GCM, AUC = ', auc_vals[6]))
  if (model == 'LM') {
    title <- 'Linear Model'
  } else if (model == 'RF') {
    title <- 'Random Forest'
  } else if (model == 'NN') {
    title <- 'Neural Network'
  } else if (model == 'SVM') {
    title <- 'Support Vector Machine'
  }
  p <- ggplot(filter(df, Model == model),
              aes(FPR, TPR, color = VIM)) + 
    geom_line() +
    geom_abline(intercept = 0, slope = 1,
                linetype = 'dashed', color = 'black') +
    scale_color_manual(name = 'VIM',
                       breaks = unique(df$VIM),
                       labels = lbls,
                       values = pal_d3()(6)) + 
    lims(x = c(0, 1), y = c(0, 1)) +
    labs(title = title, 
         x = 'False Positive Rate', 
         y = 'True Positive Rate') +
    theme_bw() + 
    theme(plot.title = element_text(hjust = 0.5),
          legend.justification = c(0.99, 0.01),
          legend.position = c(0.99, 0.01))
  return(p)
}
p_lm <- plot_fn('LM')
p_rf <- plot_fn('RF')
p_nn <- plot_fn('NN')
p_svm <- plot_fn('SVM')
grid.arrange(p_lm, p_rf, p_nn, p_svm, ncol = 2)

# Remove AUC values from plots, create table with AUC and associated SE

