library(Publish)
library(survival)
library(parallel)
library(doParallel)
library(foreach)
library(WeightIt)
library(survey)
library(tableone)

rm(list = ls())
###################################### hdPS with nested case-control ######################################
load("Data/proxy_ncc.RData")
load("Data/ncc4_data.RData")

# Investigator-specified covariates 
vars.measured <- c("sex", "age")

# Empirical covariates - top 200 based on LASSO
vars.proxy <- names(dat.proxy.lasso[,2:201])

# Covariates for hdPS
vars.all.lasso <- c(vars.measured, vars.proxy)

# PS formula
ps.formula <- as.formula(paste("I(anyDMD == 'Yes') ~ ", paste(vars.all.lasso, collapse = " + ")))

# Model matrix for fitting LASSO
X <- model.matrix(ps.formula, data = dat.all.lasso)[,-1]
A <- dat.all.lasso$anyDMD

# Detect the number of cores
n_cores <- parallel::detectCores()

# Create a cluster of cores
cl <- makeCluster(n_cores - 1)

# Register the cluster for parallel processing
registerDoParallel(cl)

# 5-fold cross-validation for selecting hyperparameters
set.seed(123)
fit.lasso <- cv.glmnet(x = X, y = A, nfolds = 5, parallel = T, family = "binomial", alpha = 1,
                       penalty.factor = c(rep(0,2), rep(1,200)))
stopCluster(cl)

# Predicting PS
dat.ncc4$ps.lasso <- predict(fit.lasso, type = "response", newx = X, s = fit.lasso$lambda.min)
summary(dat.ncc4$ps.lasso)

# Stabilized IPW
dat.ncc4$A <- ifelse(dat.ncc4$anyDMD == "Yes", 1, 0)
dat.ncc4$sipw.lasso <- with(dat.ncc4, ifelse(A == 1, mean(A)/ps.lasso, (1 - mean(A))/(1 - ps.lasso)))
summary(dat.ncc4$sipw.lasso)

# Truncated at 99th percentile
dat.ncc4$sipw.lasso_t <- WeightIt::trim(dat.ncc4$sipw.lasso, at = 0.99)
summary(dat.ncc4$sipw.lasso_t)

# Balance checking
design.ipw <- svydesign(ids = ~id, weights = ~sipw.lasso_t, data = dat.ncc4)
tab.ow <- svyCreateTableOne(vars = vars.measured, strata = "anyDMD", data = design.ipw, test = F)
print(tab.ow, smd = T) # Age and sex are imbalanced

# Outcome model with IPW
fit.hdps <- svycoxph(Surv(Time, Fail) ~ anyDMD + sex + age + strata(Set), design = design.ipw)
summary(fit.hdps)
publish(fit.hdps)
