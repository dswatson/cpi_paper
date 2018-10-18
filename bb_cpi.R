### Bayesian bootstrap for posterior distribution of CPI ###

# Load libraries, register cores
library(mlbench)
library(data.table)
library(doMC)
registerDoMC(4)

# Hyperparameters
n <- 1e3
b <- 1e5

# Simulate data
dat <- mlbench.friedman1(3 * n)
dat <- dat0 <- data.frame(dat$x, y = dat$y)
dat0[, 10] <- dat[sample.int(3 * n), 10]

# Split into training and test sets
train <- dat[seq_len(2 * n), ]
train0 <- dat0[seq_len(2 * n), ]
test <- dat[(2 * n + 1):(3 * n), ]
test0 <- dat0[(2 * n + 1):(3 * n), ]

# Train models
f <- lm(y ~ ., data = train)
f0 <- lm(y ~ ., data = train0)

# Test models
y_hat <- predict(f, test)
y_hat0 <- predict(f0, test0)
loss <- log((test$y - y_hat)^2)
loss0 <- log((test0$y - y_hat0)^2)

# Define delta vector
delta <- loss0 - loss

# Draw random Dirichlet weights
wts <- matrix(rexp(n * b), ncol = n, byrow = TRUE)
wts <- wts / rowSums(wts)

# Output function
fn <- function(i) {
  # Would love to use the weighted.mean function, but no way to compute
  # weighted variance when the weights sum to 1 (equivalent to just a single
  # observation), so we're using a sampling approach instead
  delta_b <- sample(delta, size = n, replace = TRUE, prob = wts[i, ])
  data.table(mu = mean(delta_b), sigma = sd(delta_b))
}

# Execute in parallel
posterior <- foreach(i = seq_len(b), .combine = rbind) %dopar% fn(i)

# Compute a posterior for effect size as t-statistic
posterior[, t := mu / (sigma / sqrt(n))]

# Now compare to the theoretical null
library(ggplot2)
ggplot(posterior) +
  geom_histogram(aes(t, y = ..density..), bins = 60, color = 'black') +
  stat_function(fun = dt, args = list(df = n - 1), color = 'red') +
  geom_vline(xintercept = qt(0.95, df = n - 1), color = 'blue') + 
  theme_bw()

ggplot(posterior) +
  geom_histogram(aes(t, y = ..density..), bins = 60, color = 'black') +
  geom_vline(xintercept = qt(0.95, df = n - 1), 
             linetype = 'dashed', color = 'red') + 
  theme_bw()

# How much of the HDI lies beyond the critical value?
# FUN FACT: The HDI is not the empirical 2.5 and 97.5 quantiles, 
# altho that'll nearly work for symmetric, unimodal distributions like this one. 
# The true definition is the interval such that
# for all x such that p(x) > W, the integral = 0.95.
# This sets a lower bound on acceptable densities, which can lead to 
# asymmetric or even split HDIs. See Kruschke, 2015, p. 88 or ?hdi.
# So really just imagine starting W at the mode and slide it down
# until it captures 95% of the density. But you can just do:
bounds <- posterior[, as.numeric(hdi(t, credMass = 0.95))]  # Default credMass
posterior[, hdi := ifelse(t > bounds[1] & t < bounds[2], 
                          TRUE, FALSE)]
posterior[hdi == TRUE, sum(t > qt(0.95, df = n - 1)) / .N]

# What's the posterior error probability (PEP) AKA local FDR?
posterior[, sum(t <= qt(0.95, df = n - 1)) / .N]

# Reproduce this for each Xj and take the cumulative mean for q-values

# Mu looks just about normal, which is nice. 
# Sigma does too after log transform.

posterior[, eff_size := mu / sigma]

# Eff_size looking nice and normal as well


# For just CPI 
posterior <- foreach(i = seq_len(b), .combine = c) %dopar% 
  weighted.mean(delta, w = wts[i, ])

# One idea: define the ROPE using a frequentist interval, i.e.
rope <- c(-Inf, qt(0.95, df = n - 1))
hdi <- posterior[, quantile(eff_size, c(0.025, 0.975))]

# What we really want to know is what proportion of the posterior distro 
# is below the critical value, qt(0.95, df = n - 1)
posterior[, sum(eff_size < qt(0.95, df = n - 1)) / .N]

# This does assume our null distro is well-specified...


# It seems that the purported null distro, which runs 
# c(-Inf, qt(0.95, df = n - 1)), is much wider than the posterior eff_size
# distro. Will this be true over thousands of runs?

# Alternative way to get posterior distribution is with BEST
# Using default priors here (basically noninformative)
# See Kruschke, 2013; Kruschke, 2015
library(BEST)
posterior <- BESTmcmc(delta, parallel = TRUE)

# And then to put this on the t-scale
posterior <- data.frame(
     mu = posterior$mu,
  sigma = posterior$sigma,
      t = posterior$mu / (posterior$sigma / sqrt(n))
)























