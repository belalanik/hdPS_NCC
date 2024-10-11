library(Publish)
library(survival)

################ Time-dependent Cox on full cohort ################
rm(list = ls())

# Load data
load("Data/simdata.RData")

### Data set up ###
dat.tv <- simdat
dat.tv <- tmerge(dat.tv, dat.tv, id = id, event = event(follow_up, mortality_outcome))
dat.tv <- tmerge(dat.tv, dat.tv, id = id, anyDMD_tv = tdc(yrs_anyDMD, anyDMD))
table(dat.tv$anyDMD_tv, useNA = "always")
dat.tv$anyDMD_tv <- car::recode(dat.tv$anyDMD_tv, " 'Yes' = 'Yes'; else = 'No' ")
dat.tv$anyDMD_tv <- factor(dat.tv$anyDMD_tv, levels = c("No", "Yes"))
dat.tv <- dat.tv[order(dat.tv$id, dat.tv$tstart),]
head(dat.tv)

# Unadjusted 
fit.tdcox0 <- coxph(Surv(tstart, tstop, event) ~ anyDMD_tv, data = dat.tv)
publish(fit.tdcox0)

# Adjusted
fit.tdcox1 <- coxph(Surv(tstart, tstop, event) ~ anyDMD_tv + sex + age + ses + cci + year, data = dat.tv)
publish(fit.tdcox1)

cox.zph(fit.tdcox1) # Assumption is met for all variables