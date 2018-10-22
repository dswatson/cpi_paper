
library(data.table)
library(batchtools)
library(ggplot2)
library(ggsci)

# Simulation parameters ----------------------------------------------------------------
num_replicates <- 10000
n <- 1000
p <- 10

# Algorithm parameters ----------------------------------------------------------------
learners <- c("classif.logreg", "classif.ranger", "classif.nnet", "classif.svm")
tests <- c("t", "fisher")
measures <- c("mmce", "logloss")

# Registry ----------------------------------------------------------------
reg_name <- "cpi_sim_holdout_classif"
reg_dir <- file.path("registries", reg_name)
dir.create("registries", showWarnings = FALSE)
unlink(reg_dir, recursive = TRUE)
makeExperimentRegistry(file.dir = reg_dir, 
                       packages = "mlr", 
                       source = c("../cpi_mlr.R", "problems.R"))

# Problems ----------------------------------------------------------------
addProblem(name = "linear", fun = linear_data)
addProblem(name = "nonlinear", fun = nonlinear_data)

# Algorithms ----------------------------------------------------------------
cpi <- function(data, job, instance, learner_name, ...) {
  par.vals <- switch(learner_name, 
                     classif.ranger = list(num.trees = 50), 
                     classif.nnet = list(size = 20, decay = .1, trace = FALSE), 
                     classif.svm = list(kernel = "radial"), 
                     list())
  as.list(brute_force_mlr(task = instance, learner = makeLearner(learner_name, par.vals = par.vals, predict.type = "prob"), 
                          resampling = makeResampleDesc("Holdout"), ...))
}
addAlgorithm(name = "cpi", fun = cpi)

# Experiments -----------------------------------------------------------
prob_design <- list(linear = expand.grid(n = n, p = p, outcome = "classif",
                                         stringsAsFactors = FALSE), 
                    nonlinear = expand.grid(n = n, p = p, outcome = "classif",
                                            stringsAsFactors = FALSE))
algo_design <- list(cpi = expand.grid(learner_name = learners,
                                      test = tests,
                                      measure = measures,
                                      permute = TRUE,
                                      log = TRUE, 
                                      stringsAsFactors = FALSE))
addExperiments(prob_design, algo_design, repls = num_replicates)
summarizeExperiments()
#testJob(1)

# Submit -----------------------------------------------------------
if (grepl("node\\d{2}|bipscluster", system("hostname", intern = TRUE))) {
  ids <- findNotStarted()
  ids[, chunk := chunk(job.id, chunk.size = 400)]
  submitJobs(ids = ids, # walltime in seconds, 10 days max, memory in MB
             resources = list(name = reg_name, chunks.as.arrayjobs = TRUE, 
                              ncpus = 1, memory = 6000, walltime = 10*24*3600, 
                              max.concurrent.jobs = 400))
} else {
  submitJobs()
}
waitForJobs()

# Get results -------------------------------------------------------------
res_wide <- flatten(flatten(ijoin(reduceResultsDataTable(), getJobPars())))
res <- melt(res_wide, measure.vars = patterns("^Variable*", "^CPI*", "^statistic*", "^p.value*"), 
            value.name = c("Variable", "CPI", "Statistic", "p.value"))
res[, Variable := factor(Variable,
                         levels = paste0("x", 1:unique(p)), 
                         labels = paste0("X", 1:unique(p)))]
res[, Learner := factor(learner_name, 
                        levels = c("classif.logreg", "classif.svm", "classif.ranger", "classif.nnet"), 
                        labels = c("Logistic regression", "Support vector machine", "Random forest", "Neural network"))]
res[, Problem := factor(problem, 
                levels = c("linear", "nonlinear"), 
                labels = c("Linear data", "Non-linear data"))]
saveRDS(res, "simulation_holdout_classif.Rds")

# Plots -------------------------------------------------------------
# Boxplots of CPI values per variable
lapply(unique(res$measure), function(m) {
  ggplot(res[measure == m, ], aes(x = Variable, y = CPI)) + 
    geom_boxplot(outlier.size = .01) + 
    facet_grid(Problem ~ Learner, scales = "free") + 
    geom_hline(yintercept = 0, col = "red") + 
    xlab("Variable") + ylab("CPI value")
  ggsave(paste0("holdout_classif_CPI_", m, ".pdf"), width = 10, height = 5)
  #ggsave(paste0("holdout_classif_CPI_", m, ".png"), width = 10, height = 5, dpi = 300)
})

# Histograms of t-test statistics (only null variables)
lapply(unique(res$measure), function(m) {
  ggplot(res[measure == m & test == "t" & Variable %in% c("X1", "X2"), ], aes(Statistic)) +
    geom_histogram(aes(y = ..density..), bins = 100) +
    facet_grid(Problem ~ Learner) +
    stat_function(fun = dt, color = 'red', args = list(df = unique(res$n) - 1)) +
    xlab("Test statistic") + ylab("Density")
  ggsave(paste0("holdout_classif_tstat_", m, ".pdf"), width = 10, height = 5)
})

# Power (mean over replications)
res[, reject := p.value <= 0.05]
res_mean <- res[, .(power = mean(reject, na.rm = TRUE)), by = .(Problem, algorithm, Learner, test, Variable, measure)]
levels(res_mean$Variable) <- rep(c(0, 0, -.5, .5, -1, 1, -1.5, 1.5, -2, 2), each = p/10)
res_mean[, Variable := abs(as.numeric(as.character(Variable)))]
res_mean[, power := mean(power), by = list(Problem, algorithm, Learner, test, Variable, measure)]
res_mean[, Test := factor(test, levels = c("fisher", "t"), labels = c("Fisher", "t-test"))]
lapply(unique(res$measure), function(m) {
  ggplot(res_mean[measure == m, ], aes(x = Variable, y = power, col = Learner, shape = Learner, linetype = Test)) +
    geom_line() + geom_point() +
    facet_wrap(~ Problem) +
    geom_hline(yintercept = 0.05, col = "black", linetype = "dashed") +
    scale_color_npg() +
    scale_y_continuous(breaks = c(0, .05, .25, .5, .75, 1), limits = c(0, 1)) + 
    xlab("Effect size") + ylab("Rejected hypotheses")
  ggsave(paste0("holdout_classif_power_", m, ".pdf"), width = 10, height = 5)
})
