library(ggplot2)
library(ggsci)
library(mlr)
library(data.table)
library(cpi)

seed <- 100

# Convert 2-level factor to binary
task <- createDummyFeatures(bh.task, method = "reference", cols = "chas")

# Linear model
set.seed(seed)
cpi_lm_log <- cpi(task = task, learner = makeLearner("regr.lm"), 
                  resampling = makeResampleDesc("Subsample", iters = 20), 
                  test = "t", measure = mse, log = FALSE)

# SVM
set.seed(seed)
cpi_svm_log <- cpi(task = task, learner = makeLearner("regr.svm", kernel = "radial"), 
                   resampling = makeResampleDesc("Subsample", iters = 20), 
                   test = "t", measure = mse, log = FALSE)

# Combine for plotting
res <- rbind(data.table(Learner = "Linear model", Log = "Multiplicative CPI", cpi_lm_log[, c("Variable", "CPI", "SE", "p.value")]),
             data.table(Learner = "Support vector machine", Log = "Multiplicative CPI", cpi_svm_log[, c("Variable", "CPI", "SE", "p.value")]))#, 
res[, p.adj := p.adjust(p.value, "holm")]
res[, signif := ifelse(p.adj <= .05, 1, 0)]
levels <- res[Learner == "Support vector machine", as.character(Variable)[order(CPI)]]
labels <- levels
labels[labels == "chas.1"] <- "chas"
res[, Variable := factor(Variable, levels = levels, labels = labels)]

# Plot
ggplot(res, aes(x = Variable, fill = Learner, y = CPI, alpha = signif)) + 
  geom_bar(stat = "identity", position = position_dodge(-.9)) + 
  geom_errorbar(aes(ymin = CPI - SE, ymax = CPI + SE), position = position_dodge(-.9)) +
  #facet_wrap(~ Log, scales = "free") + 
  scale_fill_npg() +
  scale_alpha_continuous(range = c(.4, 1)) + 
  coord_flip() + 
  theme_bw() + 
  theme(legend.position = "top") + 
  guides(alpha = FALSE) + 
  ylab(expression('CPI'[lambda]))
ggsave("boston_housing.pdf", width = 6, height = 8)
