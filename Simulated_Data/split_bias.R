library(ggplot2)
library(tidyr)
library(ranger)

#source("cpi.R")
source("cpi_noseed.R")

n <- 100
p <- 1000
num_cats <- rep(seq(2:11), each = p/10) # 2:11 unique values each
num_replicates <- 30 #100

num.trees <- 50
mtry <- 10

## Data simulation
simulate_data <- function(n, num_cats) {
  x <- round(sapply(num_cats, runif, n = n, min = 1))
  y <- runif(n)
  data.frame(y = y, x)
}

## Run comparison
res <- replicate(num_replicates, {
  dat <- simulate_data(n, num_cats)
  
  # Standard impurity
  #impurity <- ranger(y ~ ., dat, num.trees = num.trees, mtry = mtry, importance = "impurity")$variable.importance

  # Brute force
  #brute_force <- brute_force(y = dat[, 1], x = dat[, -1], B = num.trees, mtry = mtry, test = NULL)

  # Sub-forest
  rf_split1 <- rf_split(y = dat[, 1], x = dat[, -1], B = num.trees, mtry = mtry, test = NULL, n.sub = 1)
  
  # Repeated sub-forest 100
  rf_split100 <- rf_split(y = dat[, 1], x = dat[, -1], B = num.trees, mtry = mtry, test = NULL, n.sub = 100)
  
  # Repeated sub-forest 500
  rf_split500 <- rf_split(y = dat[, 1], x = dat[, -1], B = num.trees, mtry = mtry, test = NULL, n.sub = 500)

  c(#impurity = impurity,
    #brute_force = brute_force, 
    rf_split1 = rf_split1, 
    rf_split100 = rf_split100,
    rf_split500 = rf_split500)
})

df <- data.frame(t(res))
df_long <- gather(df, value = "Importance")
df_long$Method <- factor(gsub("\\.X\\d+", "", df_long$key), 
                         levels = c(#"impurity",
                                    #"brute_force", 
                                    "rf_split1", 
                                    "rf_split100", 
                                    "rf_split500"))
df_long$Variable <- factor(gsub("^\\w+\\.", "", df_long$key), 
                           levels = paste0("X", 1:p))
levels(df_long$Variable) <- num_cats

## Plot
ggplot(df_long, aes(x = Variable, y = Importance)) + 
  geom_boxplot(outlier.size = .5) + 
  geom_hline(yintercept = 0, col = "red") + 
  facet_wrap(~ Method, scales = "free")
#ggsave("split_bias.pdf")
