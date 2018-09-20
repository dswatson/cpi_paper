library(ggplot2)
library(tidyr)

source("cpi_mlr.R")

# Regression
brute_force_mlr(task = bh.task, learner = makeLearner("regr.lm"), 
                resampling = makeResampleDesc("CV", iters = 5))
brute_force_mlr(task = bh.task, learner = makeLearner("regr.lm"), 
                resampling = makeResampleDesc("CV", iters = 5), measure = mae)
brute_force_mlr(task = bh.task, learner = makeLearner("regr.lm"), 
                resampling = makeResampleDesc("CV", iters = 5), test = "t")
brute_force_mlr(task = bh.task, learner = makeLearner("regr.lm"), 
                resampling = makeResampleDesc("CV", iters = 5), test = "t", measure = mae)

# Classification (binary)
brute_force_mlr(task = pid.task, learner = makeLearner("classif.logreg", predict.type = "prob"), 
                resampling = makeResampleDesc("CV", iters = 5))
brute_force_mlr(task = pid.task, learner = makeLearner("classif.logreg", predict.type = "prob"), 
                resampling = makeResampleDesc("CV", iters = 5), measure = acc, test = NULL)
brute_force_mlr(task = pid.task, learner = makeLearner("classif.logreg", predict.type = "prob"), 
                resampling = makeResampleDesc("CV", iters = 5), test = "t")
brute_force_mlr(task = pid.task, learner = makeLearner("classif.logreg", predict.type = "prob"), 
                resampling = makeResampleDesc("CV", iters = 5), test = "t", measure = acc)

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
