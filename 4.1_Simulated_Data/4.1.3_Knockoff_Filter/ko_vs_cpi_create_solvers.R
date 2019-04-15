
library(data.table)
library(batchtools)
library(knockoff)

set.seed(42)

# Simulation parameters ----------------------------------------------------------------
p <- 1000
rho <- seq(from = 0, to = 0.8, by = 0.1)

# Registry ----------------------------------------------------------------
reg_name <- "ko_vs_cpi_solver"
reg_dir <- file.path("registries", reg_name)
dir.create("registries", showWarnings = FALSE)
unlink(reg_dir, recursive = TRUE)
makeRegistry(file.dir = reg_dir, 
             packages = c("knockoff"))

# Create solver ----------------------------------------------------------------
create_solver <- function(p, rho) {
  # Create covariance matrix
  if (rho == 0) {
    Sigma <- diag(p)
  } else {
    Sigma <- toeplitz(rho^(0:(p - 1)))
  }
  
  # Create and save solver
  diag_s <- create.solve_asdp(Sigma)
  
  diag_s
}

# Create jobs
batchMap(fun = create_solver, args = expand.grid(p = p, rho = rho,
                                                 stringsAsFactors = FALSE))

# Submit -----------------------------------------------------------
if (grepl("node\\d{2}|bipscluster", system("hostname", intern = TRUE))) {
  ids <- findNotSubmitted()
  #findExperiments(findNotSubmitted())
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
res <- reduceResultsList()
pars <- flatten(getJobPars(findDone()))

# Save
for (i in 1:length(res)) {
  saveRDS(res[[i]], paste0("solvers/", pars[i, p], "_", pars[i, rho], ".Rds"))
}
