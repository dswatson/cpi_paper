
library(data.table)
library(batchtools)
library(ggplot2)
library(ggsci)

set.seed(42)

# Simulation parameters ----------------------------------------------------------------
num_replicates <- 10000
n <- 1000
p <- c(20, 50, 100)
type <- c("linear", "nonlinear")
cov_base <- .5

# Algorithm parameters ----------------------------------------------------------------
learners <- c("regr.lm", "regr.ranger", "regr.nnet", "regr.svm")

# Registry ----------------------------------------------------------------
reg_name <- "comparative_performance"
reg_dir <- file.path("registries", reg_name)
dir.create("registries", showWarnings = FALSE)
unlink(reg_dir, recursive = TRUE)
makeExperimentRegistry(file.dir = reg_dir, 
                       packages = c("mlr", "cpi", "mvtnorm", "vimp"),
                       source = c("problems.R", "algorithms.R"))

# Problems ----------------------------------------------------------------
addProblem(name = "data", fun = data_train_test, seed = 43)

# Algorithms ----------------------------------------------------------------
addAlgorithm(name = "cpi", fun = cpi_fn)
addAlgorithm(name = "anova", fun = anova_fn)
addAlgorithm(name = "loco", fun = loco_fn)
addAlgorithm(name = "gcm", fun = gcm_fn)

# Experiments -----------------------------------------------------------
prob_design <- list(data = expand.grid(n = n, p = p, outcome = "regr",
                                       type = type, cov_base = cov_base,
                                       stringsAsFactors = FALSE))
algo_design <- list(cpi = expand.grid(learner_name = learners,
                                      test = "t", measure = "mse",
                                      stringsAsFactors = FALSE), 
                    anova = expand.grid(learner_name = learners,
                                        stringsAsFactors = FALSE),
                    loco = expand.grid(learner_name = learners,
                                       stringsAsFactors = FALSE),
                    gcm = expand.grid(learner_name = learners,
                                     stringsAsFactors = FALSE))
addExperiments(prob_design, algo_design, repls = num_replicates)
summarizeExperiments()
#testJob(1)

# Submit -----------------------------------------------------------
if (grepl("node\\d{2}|bipscluster", system("hostname", intern = TRUE))) {
  #ids <- findNotStarted()
  #ids <- findExperiments(findNotSubmitted(), algo.name = "cpi") 
  ids <- findExperiments(findNotSubmitted(), algo.pars = (learner_name == "regr.nnet")) 
  
  #ids <- findExperiments(findNotDone(), algo.pars = (learner_name == "regr.svm")) 
  
  #aa <- flatten(getJobTable(findNotDone()))
  #aa[, table(algorithm, learner_name)]
  
  #bb <- flatten(getJobTable(findDone()))
  #bb[, mean(time.running), by = .(algorithm, learner_name)]
  
  ids[, chunk := chunk(job.id, chunk.size = 200)]
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
res <- melt(res_wide, measure.vars = patterns("^Variable*", "^p.value*"),
            value.name = c("Variable", "p.value"))
res <- res[!is.na(Variable), ]
res[, Variable := factor(Variable,
                         levels = paste0("x", 1:max(p)),
                         labels = paste0("X", 1:max(p)))]
res[, Learner := factor(learner_name,
                        levels = c("regr.lm", "regr.svm", "regr.ranger", "regr.nnet"),
                        labels = c("Linear model", "Support vector machine", "Random forest", "Neural network"))]
res[, Type := factor(type,
                     levels = c("linear", "nonlinear"),
                     labels = c("Linear data", "Nonlinear data"))]
res[, Test := factor(algorithm,
                     levels = c("cpi", "anova", "loco", "gcm"),
                     labels = c("CPI", "ANOVA", "LOCO", "GCM"))]
saveRDS(res, paste0(reg_name, ".Rds"))

# Plots -------------------------------------------------------------
res <- readRDS(paste0(reg_name, ".Rds"))

# Merge with SVM/anova results
res[, test := NULL]
res[, measure := NULL]
res_svm <- readRDS("comparative_performance_highdim_svm.Rds")
res <- rbind(res, res_svm)

# Power (mean over replications)
lapply(unique(res$p), function(pp) {
  res_pp <- res[p == pp, ]
  res_pp[, reject := p.value <= 0.05]
  res_mean <- res_pp[, .(power = mean(reject, na.rm = TRUE)), by = .(Type, Test, Learner, Variable, cov_base, n)]
  res_mean[, Variable := factor(Variable)]
  levels(res_mean$Variable) <- rep(seq(0, .9, length.out = 10), each = pp/10)
  res_mean[, Variable := abs(as.numeric(as.character(Variable)))]
  res_mean[, power := mean(power), by = list(Type, Test, Learner, Variable, cov_base, n)]
  
  ggplot(res_mean, aes(x = Variable, y = power, col = Test, shape = Test)) +
    geom_line() + geom_point() +
    facet_grid(Type ~ Learner) +
    geom_hline(yintercept = 0.05, col = "black", linetype = "dashed") +
    scale_color_npg() +
    scale_y_continuous(breaks = c(0, .05, .25, .5, .75, 1), limits = c(0, 1)) +
    xlab("Effect size") + ylab("Rejection proportion") +
    theme_bw() + 
    theme(legend.position = 'bottom')
  ggplot2::ggsave(paste0(reg_name, "_toeplitz_highdim_", pp, ".pdf"), width = 10, height = 6)
  #ggplot2::ggsave(paste0(reg_name, "_toeplitz_highdim_", pp, ".png"), width = 10, height = 6, dpi = 300)
})
