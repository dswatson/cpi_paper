library(ggplot2)
library(tidyr)
library(ranger)

source("cpi.R")
#source("cpi_noseed.R")

n <- 100
p <- 1000
num_replicates <- 30 #100

num.trees <- 50
mtry <- 10

## Data simulation
simulate_data <- function(n, p) {
  x <- matrix(runif(n * p), ncol = p)
  y <- rnorm(n)
  data.frame(y = y, x)
}

## Run comparison
res <- replicate(num_replicates, {
  dat <- simulate_data(n, p)
  
  # Brute force
  #brute_force <- brute_force(y = dat[, 1], x = dat[, -1], B = num.trees, mtry = mtry, test = NULL)
  
  # Sub-forest
  rf_split1 <- rf_split(y = dat[, 1], x = dat[, -1], B = num.trees, mtry = mtry, test = NULL, n.sub = 1)
  
  # Repeated sub-forest 100
  rf_split100 <- rf_split(y = dat[, 1], x = dat[, -1], B = num.trees, mtry = mtry, test = NULL, n.sub = 100)
  
  # Repeated sub-forest 500
  rf_split500 <- rf_split(y = dat[, 1], x = dat[, -1], B = num.trees, mtry = mtry, test = NULL, n.sub = 500)
  
  c(#brute_force = brute_force, 
    rf_split1 = rf_split1, 
    rf_split100 = rf_split100,
    rf_split500 = rf_split500)
})

df <- data.frame(t(res), repl = factor(1:num_replicates))
df_long <- gather(df, key = "repl", value = "Importance")
colnames(df_long)[2] <- "key"
df_long$Method <- factor(gsub("\\.X\\d+", "", df_long$key), 
                         levels = c(#"brute_force", 
                           "rf_split1", 
                           "rf_split100", 
                           "rf_split500"))

## Plot
ggplot(df_long, aes(x = repl, y = Importance)) + 
  geom_boxplot(outlier.size = .5) + 
  geom_hline(yintercept = 0, col = "red") + 
  facet_wrap(~ Method, scales = "free")

## Values are actually the same for all except the first replication ...
df_long[df_long$key == "rf_split1.X1", "Importance"]
