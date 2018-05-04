
library(mlr)

brute_force_mlr <- function(task, learner, 
                            resampling = makeResampleDesc("CV", iters = 5), 
                            measures = NULL, 
                            verbose = FALSE) {
  # Full model resampling
  full_model <- resample(learner, task, resampling = resampling, measures = measures, show.info = verbose)
  
  # For each feature, fit reduced model and return difference in error
  cpi <- sapply(1:getTaskNFeats(task), function(i) {
    reduced_task <- subsetTask(task, features = getTaskFeatureNames(task)[-i])
    reduced_model <- resample(learner, reduced_task, resampling = resampling, measures = measures, show.info = verbose)
    unname(reduced_model$aggr - full_model$aggr)
  })
  names(cpi) <- getTaskFeatureNames(task)
  cpi
}

