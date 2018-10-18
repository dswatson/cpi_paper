
library(data.table)
library(batchtools)
library(ggplot2)

# Simulation parameters ----------------------------------------------------------------
num_replicates <- 1000
n <- 100
p <- 10

# Algorithm parameters ----------------------------------------------------------------
learners <- c("regr.lm", "regr.ranger", "regr.nnet", "regr.svm")
tests <- c("t", "fisher")

# Registry ----------------------------------------------------------------
reg_name <- "cpi_null_outsample"
reg_dir <- file.path("registries", reg_name)
dir.create("registries", showWarnings = FALSE)
unlink(reg_dir, recursive = TRUE)
makeExperimentRegistry(file.dir = reg_dir, 
                       packages = c("mlr"), 
                       source = "cpi_mlr_outsample_fisher.R")

# Problems ----------------------------------------------------------------
sim_data <- function(data, job, n, p, ...) {
 x <- matrix(runif(n * p), ncol = p,
             dimnames = list(NULL, paste0('x', seq_len(p))))
 y <- rnorm(n)
 dat <- data.frame(y = y, x)
 train <- makeRegrTask(data = dat, target = "y")

 x <- matrix(runif(n * p), ncol = p,
             dimnames = list(NULL, paste0('x', seq_len(p))))
 y <- rnorm(n)
 dat <- data.frame(y = y, x)
 test <- makeRegrTask(data = dat, target = "y")

 list(train = train, test = test)
}
addProblem(name = "sim", fun = sim_data)

# Algorithms ----------------------------------------------------------------
cpi <- function(data, job, instance, learner_name, ...) {
  par.vals <- switch(learner_name, 
                     regr.ranger = list(num.trees = 50), 
                     regr.nnet = list(size = 5, decay = 0.1, trace = FALSE), 
                     regr.svm = list(kernel = "radial"), 
                     list())
  brute_force_mlr(task = instance$train, learner =  makeLearner(learner_name, par.vals = par.vals), 
                  test_data = getTaskData(instance$test), ...)
}
addAlgorithm(name = "cpi", fun = cpi)

# Experiments -----------------------------------------------------------
prob_design <- list(sim = expand.grid(n = n, p = p, 
                                      stringsAsFactors = FALSE))
algo_design <- list(cpi = expand.grid(learner_name = learners,
                                      test = tests,
                                      permute = TRUE,
                                      stringsAsFactors = FALSE))
addExperiments(prob_design, algo_design, repls = num_replicates)
summarizeExperiments()

# Submit -----------------------------------------------------------
if (grepl("bipscluster", system("hostname", intern = TRUE))) {
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
res[, Variable := factor(Variable, levels = paste0("x", 1:unique(p)))]
saveRDS(res, "null_simulation_outsample.Rds")

# Plots -------------------------------------------------------------
# Boxplots of CPI values per variable
ggplot(res, aes(x = Variable, y = CPI)) + 
  geom_boxplot() + 
  facet_wrap(~ learner_name, scales = "free") + 
  geom_hline(yintercept = 0, col = "red") + 
  xlab("Variable") + ylab("CPI value")
ggsave("outsample_CPI.pdf")

# Type 1 error (mean over replications)
res[, reject := p.value <= 0.05]
res_mean <- res[, .(power = mean(reject, na.rm = TRUE)), by = .(problem, algorithm, learner_name, test, variable)]
ggplot(res_mean[test != "fisher_exactRankTests"], aes(x = learner_name, y = power)) + 
  geom_boxplot() + 
  #facet_wrap(~ learner_name) + 
  geom_hline(yintercept = 0.05, col = "red") + 
  xlab("Test") + ylab("Type I error") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("outsample_typeI.pdf")

# Histograms of t-test statistics (over all variables)
ggplot(res[test == "t", ], aes(Statistic)) +
  geom_histogram(aes(y = ..density..), bins = 100) +
  facet_wrap(~ learner_name) +
  stat_function(fun = dt, color = 'red', args = list(df = unique(n) - 1)) +
  xlab("Test statistic") + ylab("Density")
ggsave("outsample_t.pdf")

# Histograms of CPI values (over all variables)
ggplot(res, aes(CPI)) + 
  geom_histogram(aes(y = ..density..), bins = 100) + 
  facet_wrap(~ learner_name, scales = "free") + 
  geom_vline(xintercept = 0, col = "red") + 
  xlab("CPI value") + ylab("Density")
ggsave("outsample_CPI_hist.pdf")

# Histograms of p-values (over all variables)
ggplot(res[test == "fisher", ], aes(p.value)) + 
  geom_histogram(aes(y = ..density..), bins = 100) + 
  facet_wrap(~ learner_name) + 
  xlab("p-value") + ylab("Density")
ggsave("outsample_pvalues.pdf")
