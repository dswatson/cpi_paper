library(ggplot2)
library(ggsci)
library(mlr)
library(data.table)

source("../cpi_mlr.R")

seed <- 100

set.seed(seed)
cpi_lm_log <- brute_force_mlr(task = bh.task, learner = makeLearner("regr.lm"), 
                              resampling = makeResampleDesc("CV", iters = 10), 
                              test = "t", measure = mse, log = TRUE)

set.seed(seed)
cpi_lm_dif <- brute_force_mlr(task = bh.task, learner = makeLearner("regr.lm"), 
                              resampling = makeResampleDesc("CV", iters = 10), 
                              test = "t", measure = mse, log = FALSE)

set.seed(seed)
cpi_svm_log <- brute_force_mlr(task = bh.task, learner = makeLearner("regr.svm", kernel = "radial"), 
                               resampling = makeResampleDesc("CV", iters = 10), 
                               test = "t", measure = mse, log = TRUE)

set.seed(seed)
cpi_svm_dif <- brute_force_mlr(task = bh.task, learner = makeLearner("regr.svm", kernel = "radial"), 
                               resampling = makeResampleDesc("CV", iters = 10), 
                               test = "t", measure = mse, log = FALSE)

res <- rbind(data.table(Learner = "Linear model", Log = "Log ratio", cpi_lm_log[, c("Variable", "CPI", "SE", "p.value")]),
             data.table(Learner = "Linear model", Log = "Difference", cpi_lm_dif[, c("Variable", "CPI", "SE", "p.value")]), 
             data.table(Learner = "Support vector machine", Log = "Log ratio", cpi_svm_log[, c("Variable", "CPI", "SE", "p.value")]), 
             data.table(Learner = "Support vector machine", Log = "Difference", cpi_svm_dif[, c("Variable", "CPI", "SE", "p.value")]))
res[, signif := ifelse(p.value <= .05, "*", "")]

ggplot(res, aes(x = Variable, fill = Learner, y = CPI, label = signif)) + 
  geom_bar(stat = "identity", position = position_dodge()) + 
  #geom_errorbar(aes(ymin = CPI - SE, ymax = CPI + SE), position = position_dodge()) +
  geom_text(size = 7, position = position_dodge(width = 1), hjust = 1, vjust = .8, angle = 180) + 
  facet_wrap(~ Log, scales = "free") + 
  scale_fill_npg() +
  coord_flip() + 
  theme_bw() + 
  theme(legend.position = "top")
ggsave("boston_housing.pdf", width = 8, height = 8)
