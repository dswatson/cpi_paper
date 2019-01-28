### Experiments in difference vs. log-ratios ###

# Load libraries
library(tidyverse)

# Simulate data, calculate delta and lambda
df <- data_frame(
      loss = rchisq(1e4, df = 1),
     loss0 = rchisq(1e4, df = 1),
     delta = loss0 - loss,
    lambda = log(loss0 / loss),
   r_delta = rank(desc(delta)),
  r_lambda = rank(desc(lambda))
)

# Bowtie situation
ggplot(df, aes(delta, lambda)) + 
  geom_point(size = 0.25, alpha = 0.25) + 
  labs(x = bquote('CPI'[Delta]), y = bquote('CPI'[lambda])) +
  theme_bw()

# Even with the ranks
ggplot(df, aes(r_delta, r_lambda)) + 
  geom_point(size = 0.25, alpha = 0.25) + 
  theme_bw()

# Big deltas require big values of loss0
df %>% 
  arrange(r_delta) %>%
  head(.)

# Whereas big lambdas require tiny values of loss
df %>% 
  arrange(r_lambda) %>%
  head(.)

# Because loss is nonnegative, loss0 places an upper bound on delta.
# But it doesn't place an upper bound on lambda because we can drive loss
# to arbitrarily small amounts.

# I would argue the difference is more informative than the ratio here.
# There is practically no difference between a MSE of 0.001 and 0.0005,
# but the ratio test would consider these two "as different" as values of
# 200 and 100, which is ridiculous.

# The fact that differences and ratios are not rank equivalent means we 
# can't just use one measure for the stat and another for the p-value.

# Another way to think about the difference is as a matter of stability.
# The difference measure is stabler than the ratio; small perturbations
# don't have much of an impact on the difference, but could have a 
# huge impact on the ratio.

# Think about a Bayesian approach to this -- the MSE of a model is just
# a point estimate, but we should really think about the distribution of 
# credible MSEs for the model. This could ameliorate the damage with ratios;
# a density ratio would actually be as reliable as a difference, or even more so.
# But hold on -- isn't that the whole point is doing sample-wise loss?
# This gives us a loss distribution under orig and null















