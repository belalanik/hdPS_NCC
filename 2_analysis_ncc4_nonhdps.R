library(Epi)
library(Publish)
library(survival)
library(survminer)

################ Nested case-control with 4 controls per case ################
rm(list = ls())

# Load data
load("Data/simdata.RData")

# Matching
set.seed(100)
dat.ncc4 <- ccwc(
  origin = 0,
  entry = 0,
  exit = follow_up,
  fail = mortality_outcome,
  controls = 4, match = list(ses, cci, year),
  include = list(id, follow_up, mortality_outcome, anyDMD, yrs_anyDMD, sex, age),
  data = simdat,
  silent = T
)

# Drop those experienced the event before being exposed
dat.ncc4$anyDMD[dat.ncc4$yrs_anyDMD > dat.ncc4$Time] <- NA
dat.ncc4 <- dat.ncc4[complete.cases(dat.ncc4$anyDMD),]

# Save data
save(dat.ncc4, file = "Data/ncc4_data.RData")

# Number of unique ID
length(unique(dat.ncc4$id))

# Conditional logistic regression
fit.ncc1 <- clogit(Fail ~ anyDMD + sex + age + strata(Set), data = dat.ncc4, method = "efron")
publish(fit.ncc1)

# Cox PH 
fit.ncc2 <- coxph(Surv(Time, Fail) ~ anyDMD + sex + age + strata(Set), data = dat.ncc4)
publish(fit.ncc2)