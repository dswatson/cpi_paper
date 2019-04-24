# Set working directory
#setwd('~/Documents/CPI/cpi_paper/4.2_Real_Data/4.2.2_Breast_Cancer')

# Set seed
set.seed(123, kind = "L'Ecuyer-CMRG")

# Load libraries, register cores
library(data.table)
library(limma)
library(qusage)
library(corpcor)
library(knockoff)
library(ranger)
library(stringr)
library(tidyverse)
library(doMC)
registerDoMC(8)

# Import gene expression data
dat <- readRDS('GSE165.rds')

# C2 gene sets
c2 <- read.gmt('c2.all.v6.2.symbols.gmt')
tmp1 <- data.table(GeneSymbol = colnames(dat$x))
tmp2 <- seq_along(c2) %>%
  map_df(~ data.table(Pathway = names(c2)[.x],
                   GeneSymbol = unlist(c2[[.x]])) %>%
           merge(tmp1, by = 'GeneSymbol'))
c2 <- lapply(unique(tmp2$Pathway), function(p) tmp2[Pathway == p, GeneSymbol])
names(c2) <- unique(tmp2$Pathway)

# Remove sets with fewer than 25 genes
pway_size <- sapply(seq_along(c2), function(p) length(c2[[p]]))
keep <- pway_size >= 25
c2 <- c2[keep]

# Build original model
n <- nrow(dat$x)
p <- ncol(dat$x)
df <- data.frame(dat$x, y = dat$y)
rf <- ranger(data = df, dependent.variable.name = 'y', 
             num.trees = 1e4, mtry = floor(p / 3),
             keep.inbag = TRUE, classification = TRUE,
             num.threads = 8)

# Record OOB index
oob_idx <- ifelse(simplify2array(rf$inbag.counts) == 0, TRUE, NA)

# Cross entropy loss function
loss_fn <- function(mod, dat) {
  preds <- predict(mod, dat, predict.all = TRUE, num.threads = 1)$predictions
  y_hat <- rowMeans(oob_idx * preds, na.rm = TRUE)
  loss <- -(df$y * log(y_hat) + (1 - df$y) * log(1 - y_hat))
  return(loss)
}
loss <- loss_fn(rf, df)

# Create knockoff matrix
mu <- rep(0, p)
Sigma <- matrix(cov.shrink(dat$x, verbose = FALSE), nrow = p)
#solver <- create.solve_asdp(Sigma, max.size = 1000)
#x_tilde <- create.gaussian(aa, mu, Sigma, method = 'asdp', diag_s = solver)
x_tilde <- readRDS("x_tilde_8k.Rds")

# CPI function
cpi <- function(pway) {
  # Replace submatrix of interest
  genes <- c2[[pway]]
  x_0 <- dat$x
  x_0[, genes] <- x_tilde[, genes]
  df0 <- data.frame(x_0, y = dat$y)
  # Compute null loss
  loss0 <- loss_fn(rf, df0)
  # Test CPI
  delta <- loss0 - loss
  t_test <- t.test(delta, alternative = 'greater')
  # Export results
  out <- data.table(
    GeneSet = pway,
    N_Genes = length(genes),
        CPI = mean(delta),
         SE = sd(delta) / sqrt(n),
          t = t_test$statistic,
    p.value = t_test$p.value
  )
  return(out)
}

# Execute in parallel
res <- foreach(pway = names(c2), .combine = rbind) %dopar% cpi(pway) 
res <- res %>%
  arrange(p.value) %>%
  mutate(q.value = p.adjust(p.value, method = 'fdr'))
fwrite(res, 'BreastCancer_CPI_res_sdp.csv')






