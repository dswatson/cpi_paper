library(ggplot2)
library(tidyr)

source("cpi_mlr.R")

n <- 100
p <- 10
num_cats <- rep(seq(2, 11), each = p/10) # 2:11 unique values each
num_replicates <- 100

## Data simulation
simulate_data <- function(n, num_cats) {
  x <- round(sapply(num_cats, runif, n = n, min = 1))
  y <- runif(n)
  data.frame(y = y, x)
}

## Run comparison
res <- replicate(num_replicates, {
  dat <- simulate_data(n, num_cats)
  task <- makeRegrTask(data = dat, target = "y")
  
  bf_lm <- brute_force_mlr(task = task, learner = makeLearner("regr.lm"), measures = list(mse))
  bf_rpart <- brute_force_mlr(task = task, learner = makeLearner("regr.rpart"), measures = list(mse))
  bf_nnet <- brute_force_mlr(task = task, learner = makeLearner("regr.nnet", size = 5, decay = 0.1, trace = FALSE), measures = list(mse))
  bf_svm <- brute_force_mlr(task = task, learner = makeLearner("regr.svm"), measures = list(mse))
  
  c(lm = bf_lm, 
    rpart = bf_rpart, 
    nnet = bf_nnet, 
    svm = bf_svm)
})

df <- data.frame(t(res))
df_long <- gather(df, value = "Importance")
df_long$Method <- factor(gsub("\\.X\\d+", "", df_long$key), 
                         levels = c("lm",
                                    "rpart", 
                                    "nnet",
                                    "svm"))
df_long$Variable <- factor(gsub("^\\w+\\.", "", df_long$key), 
                           levels = paste0("X", 1:p))
levels(df_long$Variable) <- num_cats

## Plot
ggplot(df_long, aes(x = Variable, y = Importance)) + 
  geom_boxplot(outlier.size = .5) + 
  geom_hline(yintercept = 0, col = "red") + 
  facet_wrap(~ Method, scales = "free")
ggsave("split_bias_mlr.pdf")
