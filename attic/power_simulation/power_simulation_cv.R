
library(data.table)
library(batchtools)
library(ggplot2)

# Simulation parameters ----------------------------------------------------------------
num_replicates <- 1000
n <- 100
p <- 10

# Algorithm parameters ----------------------------------------------------------------
learners <- c("regr.lm", "regr.ranger", "regr.nnet", "regr.svm")
tests <- c("t", "fisher", "U")

# Registry ----------------------------------------------------------------
reg_name <- "cpi_power_cv"
reg_dir <- file.path("registries", reg_name)
dir.create("registries", showWarnings = FALSE)
unlink(reg_dir, recursive = TRUE)
makeExperimentRegistry(file.dir = reg_dir, 
                       packages = c("mlr"), 
                       source = "cpi_mlr.R")

# Problems ----------------------------------------------------------------
sim_data <- function(data, job, n, p, ...) {
  beta <- rep(c(0, 0, -1, 1, -2, 2, -3, 3, -4, 4), each = p/10)
  beta0 <- 0
  
  x <- matrix(runif(n * p), ncol = p,
              dimnames = list(NULL, paste0('x', seq_len(p))))
  y <- x %*% beta + beta0 + rnorm(n)
  dat <- data.frame(y = y, x)
  makeRegrTask(data = dat, target = "y")
}
addProblem(name = "sim", fun = sim_data)

# Algorithms ----------------------------------------------------------------
cpi <- function(data, job, instance, learner_name, ...) {
  par.vals <- switch(learner_name, 
                     regr.ranger = list(num.trees = 50), 
                     regr.nnet = list(size = 3, decay = 1, trace = FALSE), 
                     regr.svm = list(kernel = "radial"), 
                     list())
  brute_force_mlr(task = instance, learner = makeLearner(learner_name, par.vals = par.vals), 
                  resampling = makeResampleDesc("CV", iters = 5), ...)
}
addAlgorithm(name = "cpi", fun = cpi)

# Experiments -----------------------------------------------------------
prob_design <- list(sim = expand.grid(n = n, p = p, 
                                      stringsAsFactors = FALSE))
algo_design <- list(cpi = expand.grid(learner_name = learners,
                                      test = tests,
                                      permute = TRUE,
                                      log = TRUE,
                                      stringsAsFactors = FALSE))
addExperiments(prob_design, algo_design, repls = num_replicates)
summarizeExperiments()

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
res[, Variable := factor(Variable, levels = paste0("x", 1:unique(p)))]
saveRDS(res, "power_simulation_cv.Rds")

# Plots -------------------------------------------------------------
# Boxplots of CPI values per variable
ggplot(res, aes(x = Variable, y = CPI)) + 
  geom_boxplot() + 
  facet_wrap(~ learner_name, scales = "free") + 
  geom_hline(yintercept = 0, col = "red") + 
  xlab("Variable") + ylab("CPI value")
ggsave("cv_CPI.pdf")

# Power (mean over replications)
res[, reject := p.value <= 0.05]
res_mean <- res[, .(power = mean(reject, na.rm = TRUE)), by = .(problem, algorithm, learner_name, test, Variable)]
levels(res_mean$Variable) <- rep(c(0, 0, -1, 1, -2, 2, -3, 3, -4, 4), each = p/10)
res_mean[, Variable := abs(as.numeric(as.character(Variable)))]
res_mean[, power := mean(power), by = list(problem, algorithm, learner_name, test, Variable)]
ggplot(res_mean, aes(x = Variable, y = power, col = test, shape = test)) + 
  geom_line() + geom_point() + 
  facet_wrap(~ learner_name) + 
  geom_hline(yintercept = 0.05, col = "black") + 
  scale_color_brewer(palette = "Set1") + 
  xlab("Effect size") + ylab("Rejected hypotheses")
ggsave("cv_power.pdf")
