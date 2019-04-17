
library(data.table)
library(batchtools)
library(ggplot2)
library(cowplot)
library(ggsci)
library(knockoff)
library(glmnet)

set.seed(42)

# Simulation parameters ----------------------------------------------------------------
num_replicates <- 10000

# Registry ----------------------------------------------------------------
reg_name <- "ko_vs_cpi"
reg_dir <- file.path("registries", reg_name)
dir.create("registries", showWarnings = FALSE)
unlink(reg_dir, recursive = TRUE)
makeExperimentRegistry(file.dir = reg_dir, 
                       packages = c("mvtnorm", "knockoff", "glmnet", "mlr", "cpi"))

# Problems ----------------------------------------------------------------
create_data <- function(data, job, n, p, rho, beta) {
  # Simulate predictors
  sigma <- toeplitz(rho^(0:(p - 1)))
  x <- matrix(rmvnorm(n = n, sigma = sigma), ncol = p,
              dimnames = list(NULL, paste0('x', seq_len(p))))
  
  # Simulate signal
  k <- 60
  nonzero <- sample(p, k)
  signs <- sample(c(1, -1), size = p, replace = TRUE)
  beta <- beta * (seq_len(p) %in% nonzero) * signs
  signal <- x %*% beta
  y <- signal + rnorm(n)
  
  # Gaussian MX knockoff parameters
  mu <- rep(0, p)
  
  # Create solver
  #diag_s <- create.solve_asdp(sigma)
  diag_s <- readRDS(paste0("solvers/", p, "_", rho, ".Rds"))
  
  # Generate knockoffs
  x_tilde <- create.gaussian(x, mu, sigma, diag_s = diag_s)
  
  list(x = x, x_tilde = x_tilde, y = y, beta = beta)
}

addProblem(name = "data", fun = create_data, seed = 43)

# Algorithms ----------------------------------------------------------------
knockoff_filter <- function(data, job, instance) {
  x <- instance$x
  x_tilde <- instance$x_tilde
  y <- instance$y
  beta <- instance$beta
  
  p <- ncol(x)
  
  w <- stat.glmnet_coefdiff(x, x_tilde, y, family = 'gaussian', cores = 1)
  tau <- knockoff.threshold(w, fdr = 0.1, offset = 0)
  pos <- which(w > tau)
  neg <- which(w <= tau)
  ko_fdr <- sum(beta[pos] == 0) / max(1, length(pos))
  
  out <- c(1 * (w[beta != 0] > tau), 
           fdr = ko_fdr)
           
  return(out)
}
addAlgorithm(name = "knockoff_filter", fun = knockoff_filter)

cpi <- function(data, job, instance) {
  x <- instance$x
  x_tilde <- instance$x_tilde
  y <- instance$y
  beta <- instance$beta
  n <- nrow(x)
  p <- ncol(x)
  
  # Prepare mlr
  task <- makeRegrTask(data = data.frame(y = y, x), target = "y")
  learner <- makeLearner("regr.cvglmnet", nlambda = 500)
  resample_instance <- makeResampleInstance(desc = makeResampleDesc("CV", iters = 10), task = task)
  measure <- mse
  
  # Error on original data
  fit_full <- cpi:::fit_learner(learner = learner, task = task, resampling = resample_instance, measure = measure)
  pred_full <- cpi:::predict_learner(fit_full, task, resampling = resample_instance)
  err_full <- cpi:::compute_loss(pred_full, measure)
  
  # Function for CPI
  cpi_fn <- function(j) {
    reduced_data <- getTaskData(task)
    reduced_data[, getTaskFeatureNames(task)[j]] <- x_tilde[, getTaskFeatureNames(task)[j]]
    reduced_task <- changeData(task, reduced_data)
    
    pred_reduced <- cpi:::predict_learner(fit_full, reduced_task, resampling = resample_instance)
    err_reduced <- cpi:::compute_loss(pred_reduced, measure)
    
    dif <- err_reduced - err_full
    t.test(dif, alternative = 'greater')$p.value
  }

  # p-values for all nonzero variables
  p_values <- rep(NA, p)
  nonzero <- which(rowSums(sapply(fit_full, function(x) {
    1:p %in% predict(x$learner.model, s = 'lambda.min', type = 'nonzero')$X1
  })) == length(fit_full))
  p_values[nonzero] <- foreach(j = nonzero, .combine = c) %do% cpi_fn(j)
  q_values <- p.adjust(p_values, method = 'fdr')
  pos <- which(q_values <= 0.1)
  neg <- which(q_values > 0.1)
  cpi_fdr <- sum(beta[pos] == 0) / max(1, length(pos))
  
  out <- c(1 * (q_values[beta != 0] <= 0.1), 
           fdr = cpi_fdr)
           
  return(out)
}
addAlgorithm(name = "cpi", fun = cpi)

# Experiments -----------------------------------------------------------
algo_design <- list(knockoff_filter = data.frame(), 
                    cpi = data.frame())
                     
# Varying rho (Fig. 1)
prob_design <- list(data = expand.grid(n = 300, p = 1000, beta = 1,
                                       rho = seq(from = 0, to = .8, by = .1), 
                                       stringsAsFactors = FALSE))
addExperiments(prob_design, algo_design, repls = num_replicates)

# Varying amplitude (Fig. 2)
prob_design <- list(data = expand.grid(n = 300, p = 1000, rho = 0, 
                                       beta = seq(from = .1, to = 1, by = .1),
                                       stringsAsFactors = FALSE))
addExperiments(prob_design, algo_design, repls = num_replicates)

summarizeExperiments()
#testJob(1) 

# Submit -----------------------------------------------------------
if (grepl("node\\d{2}|bipscluster", system("hostname", intern = TRUE))) {
  ids <- findNotDone()
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
res <- melt(res_wide, measure.vars = patterns("^result*", "^fdr*"), value.name = c("reject", "FDR"))
saveRDS(res, paste0(reg_name, ".Rds"))

# Plot results -------------------------------------------------------------
res <- readRDS(paste0(reg_name, ".Rds"))
res[, Method := factor(algorithm, levels = c("cpi", "knockoff_filter"), 
                       labels = c("CPI", "Knockoff\nFilter"))]

# Mean over replications
res[, Power := mean(reject, na.rm = TRUE), by = list(Method, n, p, rho, beta, variable)]

# Mean over variables
res[, Power := mean(Power, na.rm = TRUE), by = list(Method, n, p, rho, beta)]
res[, FDR := mean(FDR, na.rm = TRUE), by = list(Method, n, p, rho, beta)]

# Fig.1 from Candès et al. - Power
df <- res[beta == 1, mean(Power), by = list(Method, rho)]
p_rho_power <- ggplot(df, aes(x = rho, y = V1, col = Method, shape = Method)) + 
  geom_line() + geom_point() +  
  ylim(0, 1) + 
  theme_bw() + 
  scale_color_npg() + 
  xlab("Correlation coefficient") + ylab("Power")
#ggplot2::ggsave(paste0(reg_name, "_power_rho.pdf"), width = 10, height = 5)

# Fig.1 from Candès et al. - FDR
df <- res[beta == 1, mean(FDR), by = list(Method, rho)]
p_rho_fdr <- ggplot(df, aes(x = rho, y = V1, col = Method, shape = Method)) + 
  geom_hline(yintercept = 0.1, col = "black", linetype = "dashed") +
  geom_line() + geom_point() +
  ylim(0, 1) + 
  theme_bw() + 
  scale_color_npg() + 
  xlab("Correlation coefficient") + ylab("FDR")
#ggplot2::ggsave(paste0(reg_name, "_fdr_rho.pdf"), width = 10, height = 5)

# Fig.2 from Candès et al. - Power
df <- res[rho == 0, mean(Power), by = list(Method, beta)]
p_ampl_power <- ggplot(df, aes(x = beta, y = V1, col = Method, shape = Method)) + 
  geom_line() + geom_point() +
  ylim(0, 1) + 
  theme_bw() + 
  scale_color_npg() + 
  xlab("Effect size") + ylab("Power")
#ggplot2::ggsave(paste0(reg_name, "_power_ampl.pdf"), width = 10, height = 5)

# Fig.2 from Candès et al. - FDR
df <- res[rho == 0, mean(FDR), by = list(Method, beta)]
p_ampl_fdr <- ggplot(df, aes(x = beta, y = V1, col = Method, shape = Method)) + 
  geom_hline(yintercept = 0.1, col = "black", linetype = "dashed") +
  geom_line() + geom_point() +
  ylim(0, 1) + 
  theme_bw() + 
  scale_color_npg() + 
  xlab("Effect size") + ylab("FDR")
#ggplot2::ggsave(paste0(reg_name, "_fdr_ampl.pdf"), width = 10, height = 5)

# Plot all together
plot_grid(plot_grid(p_ampl_power + theme(legend.position = "none"), 
                    p_ampl_fdr + theme(legend.position = "none"), 
                    p_rho_power + theme(legend.position = "none"), 
                    p_rho_fdr + theme(legend.position = "none"), 
                    nrow = 2, ncol = 2),
          get_legend(p_ampl_power + guides(col = guide_legend(nrow = 2))), 
          nrow = 1, ncol = 2, rel_widths = c(.9, .1))
ggplot2::ggsave(paste0(reg_name, ".pdf"), width = 10, height = 8)
