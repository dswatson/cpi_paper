library(ggplot2)
library(tidyr)

source("cpi_mlr.R")

# Regression
brute_force_mlr(task = bh.task, learner = makeLearner("regr.lm"))
brute_force_mlr(task = bh.task, learner = makeLearner("regr.lm"), measure = mae)
brute_force_mlr(task = bh.task, learner = makeLearner("regr.lm"), test = "t")
brute_force_mlr(task = bh.task, learner = makeLearner("regr.lm"), test = "lrt")

# Classification (binary)
brute_force_mlr(task = pid.task, learner = makeLearner("classif.logreg", predict.type = "prob"))
brute_force_mlr(task = pid.task, learner = makeLearner("classif.logreg", predict.type = "prob"), measure = mmce)
brute_force_mlr(task = pid.task, learner = makeLearner("classif.logreg", predict.type = "prob"), test = "t")
brute_force_mlr(task = pid.task, learner = makeLearner("classif.logreg", predict.type = "prob"), test = "lrt")

# Classification (multiclass)
brute_force_mlr(task = iris.task, learner = makeLearner("classif.glmnet", predict.type = "prob"))
brute_force_mlr(task = iris.task, learner = makeLearner("classif.glmnet", predict.type = "prob"), measure = mmce)
brute_force_mlr(task = iris.task, learner = makeLearner("classif.glmnet", predict.type = "prob"), test = "t")
brute_force_mlr(task = iris.task, learner = makeLearner("classif.glmnet", predict.type = "prob"), test = "lrt")

# OOB error
brute_force_mlr(task = bh.task, learner = makeLearner("regr.ranger", num.trees = 50), resampling = "oob")
brute_force_mlr(task = bh.task, learner = makeLearner("regr.ranger", num.trees = 50), resampling = "oob", test = "t")
brute_force_mlr(task = bh.task, learner = makeLearner("regr.ranger", num.trees = 50), resampling = "oob", test = "lrt")

# Training error
brute_force_mlr(task = bh.task, learner = makeLearner("regr.lm"), resampling = "none")
brute_force_mlr(task = bh.task, learner = makeLearner("regr.lm"), resampling = "none", test = "t")
brute_force_mlr(task = bh.task, learner = makeLearner("regr.lm"), resampling = "none", test = "lrt")
