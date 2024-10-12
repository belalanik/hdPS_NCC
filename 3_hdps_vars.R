library(autoCovariateSelection)
library(glmnet)
library(survival)
library(parallel)
library(doParallel)
library(foreach)

############ Datasets #################
rm(list = ls())
load("Data/simdata.RData")
load("Data/ncc4_data.RData")

head(dat.proxy)
table(dat.proxy$dim)

################### Empirical variable, i.e., proxy variable, identification ################
id <- simdat$id

step1 <- get_candidate_covariates(df = dat.proxy, domainVarname = "dim", 
                                  eventCodeVarname = "code", 
                                  patientIdVarname = "id", 
                                  patientIdVector = id, 
                                  n = 1000, 
                                  min_num_patients = 20)
out1 <- step1$covars_data

## Assessing recurrence of codes
all.equal(id, step1$patientIds)

step2 <- get_recurrence_covariates(df = out1, 
                                   eventCodeVarname = "code", 
                                   patientIdVarname = "id",
                                   patientIdVector = id)
out2 <- step2$recurrence_data

####################### Prioritizing covariates using Cox-LASSO ##########################

# Proxy variables
vars.proxy <- names(out2)[-1]

# Merge analytic data with proxies
dat.all <- merge(simdat[,c("id", "follow_up", "mortality_outcome")], out2, by = "id", all.x = T)

# Formula
formula.out <- as.formula(paste("Surv(follow_up, mortality_outcome) ~ ", paste(vars.proxy, collapse = " + ")))

# Model matrix for fitting LASSO
X <- model.matrix(formula.out, data = dat.all)[,-1]
Y <- as.matrix(data.frame(time = dat.all$follow_up, status = dat.all$mortality_outcome))

# Detect the number of cores
n_cores <- parallel::detectCores()

# Create a cluster of cores
cl <- makeCluster(n_cores - 1)

# Register the cluster for parallel processing
registerDoParallel(cl)

# 5-fold cross-validation for selecting hyperparameters
set.seed(123)
fit.lasso <- cv.glmnet(x = X, y = Y, nfolds = 5, parallel = T, alpha = 1, family = "cox")
stopCluster(cl)

#plot(fit.lasso)

# Variable ranking based on Cox-LASSO
#lasso.coef <- coef(fit.lasso, s = fit.lasso$lambda.min)
lasso.coef <- coef(fit.lasso, s = 0.0001234098)
head(lasso.coef)
lasso.coef <- data.frame(as.matrix(lasso.coef))
lasso.coef <- data.frame(vars = rownames(lasso.coef), coef = lasso.coef)
colnames(lasso.coef) <- c("vars", "coef")
rownames(lasso.coef) <- NULL
lasso.coef$coef.abs <- abs(lasso.coef$coef)
lasso.coef <- lasso.coef[order(lasso.coef$coef.abs, decreasing = T),]
head(lasso.coef)

## Variable section based on Cox-LASSO
dat.proxy.lasso <- out2[, c("id", lasso.coef$vars)]
dat.all.lasso <- merge(dat.ncc4, dat.proxy.lasso, by = "id")

## Save
save(dat.proxy, dat.all, dat.proxy.lasso, dat.all.lasso, lasso.coef, file = "Data/proxy_ncc.RData")
