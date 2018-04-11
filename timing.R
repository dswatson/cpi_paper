### Processing Time Tests ###

# Load libraries
library(ranger)
library(party)
library(microbenchmark)
library(doMC)
registerDoMC(4)

# Competitor functions
mdi <- function(df, mtry, B, par) {
  f <- ranger(data = df, dependent.variable.name = 'y', 
              mtry = mtry, num.trees = B, replace = FALSE,
              importance = 'impurity', num.threads = ifelse(par, 4, 1))
  importance(f)
}
mda <- function(df, mtry, B, par) {
  f <- ranger(data = df, dependent.variable.name = 'y',
              mtry = mtry, num.trees = B, replace = FALSE,
              importance = 'permutation', num.threads = ifelse(par, 4, 1))
  importance(f)
}
strobl <- function(df, mtry, B) {  
  f <- cforest(y ~ ., data = df, controls = ctrl)
  varimp(f, conditional = TRUE)
}

eval <- function(dat, x, y, type, mtry, B) {
  df <- data.frame(x, y = y)
  times <- microbenchmark(
    bf_par = brute_force(x, y, type = type, test = NULL,
                         replace = FALSE, mtry = mtry, B = B, 
                         n.cores = 4, seed = NULL),
    bf_ser = brute_force(x, y, type = type, test = NULL, 
                         replace = FALSE, mtry = mtry, B = B, 
                         n.cores = 1, seed = NULL),
    mdi_par = mdi(df, mtry, B, par = TRUE),
    mdi_ser = mdi(df, mtry, B, par = FALSE),
    mda_par = mda(df, mtry, B, par = TRUE),
    mda_ser = mda(df, mtry, B, par = FALSE),
    strobl - strobl(df, mtry, B),
    times = runs, unit = 's', control = list(warmup = 5)
  )
  data.frame(VIM = times$expr,
             Time = times$time) %>%
    saveRDS(paste0('./Results/', dat, '_times.rds'))
}





