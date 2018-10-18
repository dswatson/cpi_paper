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

# Friedman1 benchmark hyperparameters
n <- 100
p <- 10

# These will come up a lot
models <- c('LM', 'RF', 'NN', 'SVM')
vims <- c('CPI', 'ANOVA', 'LRT')

# Loss function
compute_loss <- function(trn, tst) {
  
  # Linear model
  lm_fit <- lm(y ~ ., data = trn)
  lm_yhat <- predict(lm_fit, newdata = tst)
  lm_loss <- (tst$y - lm_yhat)^2
  
  # Random forest
  rf_fit <- ranger(y ~ ., data = trn, num.trees = 50)
  rf_yhat <- predict(rf_fit, tst)$predictions
  rf_loss <- (tst$y - rf_yhat)^2
  
  # Neural network
  nn_fit <- nnet(y ~ ., data = trn, size = 3, decay = 0.2, 
                 linout = TRUE, trace = FALSE)
  nn_yhat <- predict(nn_fit, newdata = tst)
  nn_loss <- (tst$y - nn_yhat)^2
  
  # Support vector machine
  svm_fit <- svm(y ~ ., data = trn, kernel = 'radial', fitted = FALSE)
  svm_yhat <- predict(svm_fit, newdata = tst)
  svm_loss <- (tst$y - svm_yhat)^2
  
  # Export
  out <- data.frame(
    lm_yhat, lm_loss, rf_yhat, rf_loss, 
    nn_yhat, nn_loss, svm_yhat, svm_loss
  )
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
      'VIM' = 'cpi',
      'VI' = c(
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

  # Put it all together
  out <- rbind(cpi_df, anova_df, lrt_df)
  out[, Run := b]
  return(out)
  
}

# Execute in parallel
out <- foreach(b = seq_len(100), .combine = rbind) %dopar% loop(b)

# Prepare to plot
obs <- rep(rep(c(1, 0), each = 5), times = 100)
auroc <- function(model) {
  pred <- list(CPI = out[VIM == 'cpi' & Model == model, VI],
               ANOVA = out[VIM == 'anova' & Model == model, VI],
               LRT = out[VIM == 'lrt' & Model == model, VI])
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
  lbls <- paste0(vims, ', AUC = ', auc_vals)
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
                       breaks = vims,
                       labels = lbls,
                       values = pal_d3()(3)) + 
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














