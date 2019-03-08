
library(data.table)
<<<<<<< HEAD
library(knockoff)
library(ranger)
library(nnet)
library(e1071)
library(vimp)
=======
library(batchtools)
library(ggplot2)
>>>>>>> c9e5e406c740332f762dfd8930a243e4b8862ed9
library(ggsci)

set.seed(42)

# Simulation parameters ----------------------------------------------------------------
num_replicates <- 1000 #10000
n <- c(100, 500, 1000)
p <- 10
type <- c("linear", "nonlinear")
cov_base <- c(0, .5)

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
addAlgorithm(name = "fcit", fun = fcit_fn)
addAlgorithm(name = "gcm", fun = gcm_fn)

<<<<<<< HEAD
loop <- function(b, rho, sim, test_n) {
  
  # Simulate data
  beta <- c(0, 0, -0.5, 0.5, -1, 1, -1.5, 1.5, -2, 2)
  Sigma <- toeplitz(rho^(0:(p - 1)))
  x <- matrix(rnorm(3 * test_n * p), ncol = p) %*% chol(Sigma)
  dimnames(x) <- list(NULL, paste0('x', seq_len(p)))
  if (sim == 'linear') {
    y <- x %*% beta + rnorm(3 * test_n)
  } else if (sim == 'nonlinear') {
    idx <- x < -qnorm(0.75) | x > qnorm(0.75)
    xx <- matrix(0, nrow = 3 * test_n, ncol = p)
    xx[idx] <- 0
    xx[!idx] <- 1
    y <- xx %*% beta + rnorm(3 * test_n)
  }
  dat <- dat0 <- data.frame(x, y)
  
  # Split into training and test sets
  idx <- seq_len(2 * test_n)
  train <- dat[idx, ]
  test <- dat[-idx, ]
  
  # Fit models
  lm_fit <- lm(y ~ ., data = train)
  rf_fit <- ranger(y ~ ., data = train, num.trees = 50)
  nn_fit <- nnet(y ~ ., data = train, size = 20, decay = 0.1, 
                 linout = TRUE, trace = FALSE)
  svm_fit <- svm(y ~ ., data = train, kernel = 'radial', fitted = FALSE)
  
  # Compute predictions, residuals, and squared errors
  lm_yhat <- predict(lm_fit, newdata = test)
  lm_eps <- test$y - lm_yhat
  lm_loss <- lm_eps^2
  rf_yhat <- predict(rf_fit, test)$predictions
  rf_eps <- test$y - rf_yhat
  rf_loss <- rf_eps^2
  nn_yhat <- predict(nn_fit, newdata = test)
  nn_eps <- test$y - nn_yhat
  nn_loss <- nn_eps^2
  svm_yhat <- predict(svm_fit, newdata = test)
  svm_eps <- test$y - svm_yhat
  svm_loss <- svm_eps^2
  
  # Generate knockoff matrix
  x_tilde <- create.second_order(x)
=======
# Experiments -----------------------------------------------------------
prob_design <- list(data = expand.grid(n = n, p = p, outcome = "regr",
                                       type = type, cov_base = cov_base,
                                       stringsAsFactors = FALSE))
algo_design <- list(cpi = expand.grid(learner_name = learners,
                                      test = "t", measure = "mse",
                                      stringsAsFactors = FALSE), 
                    anova = expand.grid(learner_name = learners,
                                        stringsAsFactors = FALSE),
                    fcit = expand.grid(learner_name = learners,
                                       stringsAsFactors = FALSE),
                    gcm = expand.grid(learner_name = learners,
                                      stringsAsFactors = FALSE))
addExperiments(prob_design, algo_design, repls = num_replicates)
summarizeExperiments()
#testJob(1)
>>>>>>> c9e5e406c740332f762dfd8930a243e4b8862ed9

# aa <- testJob(1)
# bb <- testJob(113)
# cc <- testJob(227)
# dd <- testJob(339)
#testJob(6400)

# Submit -----------------------------------------------------------
if (grepl("node\\d{2}|bipscluster", system("hostname", intern = TRUE))) {
  ids <- findNotStarted()
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
res[, Variable := factor(Variable,
                         levels = paste0("x", 1:unique(p)),
                         labels = paste0("X", 1:unique(p)))]
res[, Learner := factor(learner_name,
                        levels = c("regr.lm", "regr.svm", "regr.ranger", "regr.nnet"),
                        labels = c("Linear model", "Support vector machine", "Random forest", "Neural network"))]
res[, Type := factor(type,
                     levels = c("linear", "nonlinear"),
                     labels = c("Linear data", "Nonlinear data"))]
res[, Test := factor(algorithm,
                     levels = c("cpi", "anova", "fcit", "gcm"),
                     labels = c("CPI", "ANOVA", "FCIT", "GCM"))]
saveRDS(res, paste0(reg_name, ".Rds"))

# Plots -------------------------------------------------------------
res <- readRDS(paste0(reg_name, ".Rds"))

# Power (mean over replications)
res[, reject := p.value <= 0.05]
res_mean <- res[, .(power = mean(reject, na.rm = TRUE)), by = .(Type, Test, Learner, Variable, cov_base, n)]
levels(res_mean$Variable) <- rep(seq(0, .9, length.out = 10), each = 1)
res_mean[, Variable := abs(as.numeric(as.character(Variable)))]
res_mean[, power := mean(power), by = list(Type, Test, Learner, Variable, cov_base, n)]

# Plot uncorrelated
lapply(unique(res$n), function(nn) {
  ggplot(res_mean[cov_base == 0 & n == nn, ], aes(x = Variable, y = power, col = Test, shape = Test)) +
    geom_line() + geom_point() +
    facet_grid(Type ~ Learner) +
    geom_hline(yintercept = 0.05, col = "black", linetype = "dashed") +
    scale_color_npg() +
    scale_y_continuous(breaks = c(0, .05, .25, .5, .75, 1), limits = c(0, 1)) +
    xlab("Effect size") + ylab("Rejection proportion") +
    theme_bw() + 
    theme(legend.position = 'bottom')
  ggplot2::ggsave(paste0(reg_name, "_diag_", nn, ".pdf"), width = 10, height = 5)
})

# Plot correlated
lapply(unique(res$n), function(nn) {
  ggplot(res_mean[cov_base == .5 & n == nn, ], aes(x = Variable, y = power, col = Test, shape = Test)) +
    geom_line() + geom_point() +
    facet_grid(Type ~ Learner) +
    geom_hline(yintercept = 0.05, col = "black", linetype = "dashed") +
    scale_color_npg() +
    scale_y_continuous(breaks = c(0, .05, .25, .5, .75, 1), limits = c(0, 1)) +
    xlab("Effect size") + ylab("Rejection proportion") +
    theme_bw() + 
    theme(legend.position = 'bottom')
  ggplot2::ggsave(paste0(reg_name, "_toeplitz_", nn, ".pdf"), width = 10, height = 5)
})

