
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
                       packages = c("knockoff", "glmnet"))

# Problems ----------------------------------------------------------------
create_data <- function(data, job, n, p, type, rho, amplitude) {
  # Simulate predictors
  x <- matrix(rnorm(n * p), ncol = p)
  if (rho == 0) {
    Sigma <- diag(p)
  } else {
    Sigma <- toeplitz(rho^(0:(p - 1)))
    x <- x %*% chol(Sigma)
  }
  dimnames(x) <- list(NULL, paste0('x', seq_len(p)))
  # Simulate signal
  k <- 60
  nonzero <- sample(p, k)
  signs <- sample(c(1, -1), size = p, replace = TRUE)
  beta <- amplitude * (seq_len(p) %in% nonzero) / sqrt(n) * signs
  signal <- x %*% beta
  # Gaussian MX knockoff parameters
  mu <- rep(0, p)
  
  # Create solver
  #diag_s <- create.solve_asdp(Sigma)
  diag_s <- readRDS(paste0("solvers/", p, "_", rho, ".Rds"))
  
  # Generate knockoffs
  x_tilde <- create.gaussian(x, mu, Sigma, diag_s = diag_s)
  
  # Generate test dataset
  x_test <- matrix(rnorm(n * p), ncol = p)
  if (rho > 0) {
    x_test <- x_test %*% chol(Sigma)
  }
  dimnames(x_test) <- list(NULL, paste0('x', seq_len(p)))
  signal_test <- x_test %*% beta
  
  if (type == 'regression') {
    y <- signal + rnorm(n)
    y_test <- signal_test + rnorm(n)
  } else if (type == 'classification') {
    y <- as.factor(rbinom(n, size = 1, prob = plogis(signal)))
    y_test <- rbinom(n, size = 1, prob = plogis(signal_test))
  }
  
  list(x = x, x_test = x_test, x_tilde = x_tilde, 
       y = y, y_test = y_test, beta = beta)
}

addProblem(name = "data", fun = create_data, seed = 43)

# Algorithms ----------------------------------------------------------------
knockoff_filter <- function(data, job, instance) {
  x <- instance$x
  x_tilde <- instance$x_tilde
  y <- instance$y
  beta <- instance$beta
  
  p <- ncol(x)
  
  if (is.factor(y)) {
    w <- stat.glmnet_coefdiff(x, x_tilde, y, family = 'binomial', cores = 1)
  } else {
    w <- stat.glmnet_coefdiff(x, x_tilde, y, family = 'gaussian', cores = 1)
  }
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
  x_test <- instance$x_test
  x_tilde <- instance$x_tilde
  y <- instance$y
  y_test <- instance$y_test
  beta <- instance$beta
  
  p <- ncol(x)
  
  # Fit model, compute loss
  if (is.factor(y)) {
    f <- cv.glmnet(x, y, family = 'binomial', nlambda = 500, parallel = FALSE)
    y_hat <- predict(f, newx = x_test, s = 'lambda.min', type = 'response')
    loss <- -(y_test * log(y_hat) + (1 - y_test) * log(1 - y_hat))
  } else {
    f <- cv.glmnet(x, y, family = 'gaussian', nlambda = 500, parallel = FALSE)
    y_hat <- predict(f, newx = x_test, s = 'lambda.min')
    loss <- (y_test - y_hat)^2
  }
  p_values <- rep(NA, p)
  cpi_fn <- function(j) {
    x_test[, j] <- x_tilde[, j]
    if (is.factor(y)) {
      y_hat0 <- predict(f, newx = x_test, s = 'lambda.min', type = 'response')
      loss0 <- -(y_test * log(y_hat0) + (1 - y_test) * log(1 - y_hat0))
    } else {
      y_hat0 <- predict(f, newx = x_test, s = 'lambda.min')
      loss0 <- (y_test - y_hat0)^2
    }
    delta <- loss0 - loss
    t_test <- t.test(delta, alternative = 'greater')
    return(t_test$p.value)
  }
  
  nonzero <- predict(f, s = 'lambda.min', type = 'nonzero')$X1
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
prob_design <- list(data = expand.grid(n = 300, p = 1000, type = "regression",
                                       rho = seq(from = 0, to = 0.8, by = 0.1), amplitude = 10,
                                       stringsAsFactors = FALSE))
addExperiments(prob_design, algo_design, repls = num_replicates)

# Varying amplitude (Fig. 2)
prob_design <- list(data = expand.grid(n = 300, p = 1000, type = "regression",
                                       rho = 0, amplitude = seq(from = 0, to = 18, by = 1),
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
res[, effect_size := amplitude / sqrt(300)]

# Mean over replications
res[, Power := mean(reject, na.rm = TRUE), by = list(Method, n, p, type, rho, amplitude, effect_size, variable)]

# Mean over variables
res[, Power := mean(Power, na.rm = TRUE), by = list(Method, n, p, type, rho, amplitude, effect_size)]
res[, FDR := mean(FDR, na.rm = TRUE), by = list(Method, n, p, type, rho, amplitude, effect_size)]

# Fig.1 from Candès et al. - Power
df <- res[type == "regression" & amplitude == 10 & p == 1000 & n == 300, mean(Power), by = list(Method, rho)]
p_rho_power <- ggplot(df, aes(x = rho, y = V1, col = Method, shape = Method)) + 
  geom_line() + geom_point() +  
  ylim(0, 1) + 
  theme_bw() + 
  scale_color_npg() + 
  xlab("Correlation coefficient") + ylab("Power")
#ggplot2::ggsave(paste0(reg_name, "_power_rho.pdf"), width = 10, height = 5)

# Fig.1 from Candès et al. - FDR
df <- res[type == "regression" & amplitude == 10 & p == 1000 & n == 300, mean(FDR), by = list(Method, rho)]
p_rho_fdr <- ggplot(df, aes(x = rho, y = V1, col = Method, shape = Method)) + 
  geom_hline(yintercept = 0.1, col = "black", linetype = "dashed") +
  geom_line() + geom_point() +
  ylim(0, 1) + 
  theme_bw() + 
  scale_color_npg() + 
  xlab("Correlation coefficient") + ylab("FDR")
#ggplot2::ggsave(paste0(reg_name, "_fdr_rho.pdf"), width = 10, height = 5)

# Fig.2 from Candès et al. - Power
df <- res[type == "regression" & rho == 0 & p == 1000 & n == 300, mean(Power), by = list(Method, amplitude, effect_size)]
p_ampl_power <- ggplot(df, aes(x = effect_size, y = V1, col = Method, shape = Method)) + 
  geom_line() + geom_point() +
  ylim(0, 1) + 
  theme_bw() + 
  scale_color_npg() + 
  xlab("Effect size") + ylab("Power")
#ggplot2::ggsave(paste0(reg_name, "_power_ampl.pdf"), width = 10, height = 5)

# Fig.2 from Candès et al. - FDR
df <- res[type == "regression" & rho == 0 & p == 1000 & n == 300, mean(FDR), by = list(Method, amplitude, effect_size)]
p_ampl_fdr <- ggplot(df, aes(x = effect_size, y = V1, col = Method, shape = Method)) + 
  geom_hline(yintercept = 0.1, col = "black", linetype = "dashed") +
  geom_line() + geom_point() +
  ylim(0, 1) + 
  theme_bw() + 
  scale_color_npg() + 
  xlab("Effect size") + ylab("FDR")
#ggplot2::ggsave(paste0(reg_name, "_fdr_ampl.pdf"), width = 10, height = 5)

# Plot all together
plot_grid(plot_grid(p_rho_power + theme(legend.position = "none"), 
                    p_rho_fdr + theme(legend.position = "none"), 
                    p_ampl_power + theme(legend.position = "none"), 
                    p_ampl_fdr + theme(legend.position = "none"), 
                    nrow = 2, ncol = 2),
          get_legend(p_ampl_power + guides(col = guide_legend(nrow = 2))), 
          nrow = 1, ncol = 2, rel_widths = c(.9, .1))
ggplot2::ggsave(paste0(reg_name, ".pdf"), width = 10, height = 8)
