library(survival)
library(tableone)
rm(list = ls())

############# Generating a time-to-event outcome with a time-dependent exposure #############
#The sample size
n <- 19000
set.seed(229)
sex <- rbinom(n, 1, prob = 0.28) 
age <- round(rnorm(n, mean = 44.5, sd = 13.5)); age[age < 18] <- 18; age[age > 103] <- 103
ses <- rbinom(n, 1, prob = 0.41)
cci <- rbinom(n, 1, prob = 0.22)
year <- rbinom(n, 1, prob = 0.36)

confounders = cbind(sex, age, ses, cci, year)

lambda <- 0.03
alpha_0 <- 0.35
alpha_1 <- -0.55
alpha_2 <- 0.08
alpha_3 <- 0.45
alpha_4 <- 0.20
alpha_5 <- 0.55

ExpLin <- cbind(1, confounders) %*% c(alpha_0, alpha_1, alpha_2, alpha_3, alpha_4, alpha_5)
S <- -log(runif(n, 0, 1))/(lambda * exp(ExpLin))
beta_t <- log(0.7); beta_1 <- 0.5; beta_2 <- 0.08; beta_3 <- -0.23; beta_4 <- 0.30; beta_5 <- -0.32

U <- runif(n, 0, 1)
LinFix <- confounders %*% c(beta_1, beta_2, beta_3, beta_4, beta_5)
EventT <- ifelse(-log(1 - U) < lambda * exp(LinFix) * S, -log(1 - U)/(lambda * exp(LinFix)),
                 (-log(1 - U) - lambda * exp(LinFix) * S + lambda * exp(LinFix + beta_t) * S)/
                   (lambda * exp(LinFix + beta_t)))
CensorT <- runif(n, min = 0, max = 0.25)
# Merge into a data frame
simdat <- data.frame(id = 1:n, follow_up = pmin(EventT, CensorT), mortality_outcome = EventT < CensorT, 
                   age, sex, ses, cci, year,
                   anyDMD = -log(1 - U) >= lambda * exp(LinFix) * S & CensorT > S)
simdat$yrs_anyDMD <- pmin(simdat$follow_up, S)
simdat$anyDMD <- ifelse(simdat$anyDMD == TRUE, 1, 0)
simdat$anyDMD <- ifelse(simdat$anyDMD == TRUE, 1, 0)
simdat$yrs_anyDMD[simdat$anyDMD == 0] <- NA
head(simdat, 5)

# Exposure
round(prop.table(table(simdat$anyDMD))*100, 2)

# Exposure time
round(summary(simdat$yrs_anyDMD), 3)
simdat$yrs_anyDMD <- simdat$yrs_anyDMD * 3.3 / mean(simdat$yrs_anyDMD, na.rm = T)

# Outcome
round(prop.table(table(simdat$mortality_outcome))*100, 2)

# Follow-up time
round(summary(simdat$follow_up), 3)
simdat$follow_up <- simdat$follow_up * 11.7 / mean(simdat$follow_up, na.rm = T)

# Bivariate
table(Exposure = simdat$anyDMD, Outcome = simdat$mortality_outcome)
round(summary(simdat$yrs_anyDMD), 4)
round(summary(simdat$follow_up), 4)

# Recoding
simdat$sex <- car::recode(simdat$sex, " 0 = 'Female'; 1 = 'Male'; else = NA", as.factor = T)
simdat$ses <- car::recode(simdat$ses, " 0 = 'Low/middle'; 1 = 'High'; else = NA", as.factor = T,
                          levels = c("Low/middle", "High"))
simdat$cci <- car::recode(simdat$cci, " 0 = 'No'; 1 = 'Yes'; else = NA", as.factor = T)
simdat$year <- car::recode(simdat$year, " 0 = '1996-2005'; 1 = '2006-2017'; else = NA", as.factor = T)

simdat$anyDMD <- car::recode(simdat$anyDMD, " 0 = 'No'; 1 = 'Yes'; else = NA", as.factor = T)
simdat$mortality_outcome <- car::recode(simdat$mortality_outcome, " FALSE = 0; TRUE = 1; else = NA")

# Save data
save(simdat, file = "Data/simdata.RData")

# Variable
# id: Unique identifier
# follow_up: Follow-up time in years
# mortality_outcome: Binary outcome
# age: Age in years
# sex: Sex
# ses: Socioeconomic status 
# cci: Charlson Comorbidity Index
# year: Calendar year
# anyDMD: Time-dependent exposure
# yrs_anyDMD: Exposure time in years

# Table 1
vars <- c("anyDMD", "sex", "age", "ses", "cci", "year")
tab1 <- CreateTableOne(vars = vars, strata = "mortality_outcome", data = simdat, includeNA = T, 
                       test = F, addOverall = T)
print(tab1, showAllLevels = T)

########################### hdPS variables ###################################
id <- simdat$id

# Random proxies from difefrnt website
## Product identification numbers (PINs) - Gov.bc.ca
## https://www.cms.gov/medicare/coordination-benefits-recovery/overview/icd-code-lists
proxies <- read.csv("Data/Proxy.csv", header = T, stringsAsFactors = F, fileEncoding="latin1")
head(proxies)

# Fake proxies
set.seed(100)
n.proxy <- 200000

## 3-digit diagnostic codes
dat.diag <- data.frame(id = sample(id, size = n.proxy*0.05, replace = T), 
                       code = sample(proxies$diag, size = n.proxy*0.05, replace = T), 
                       dim = "diag")
dat.diag$code <- substr(dat.diag$code, start = 1, stop = 3)
dat.diag$code[dat.diag$code==""] <- NA
dat.diag <- na.omit(dat.diag)

## 3-digit procedure codes
dat.proc <- data.frame(id = sample(id, size = n.proxy*0.05, replace = T), 
                       code = sample(proxies$proc, size = n.proxy*0.05, replace = T), 
                       dim = "proc")
dat.proc$code <- substr(dat.proc$code, start = 1, stop = 3)
dat.proc$code[dat.proc$code==""] <- NA
dat.proc <- na.omit(dat.proc)

## 3-digit icd codes
dat.msp <- data.frame(id = sample(id, size = n.proxy*0.60, replace = T), 
                      code = sample(proxies$msp, size = n.proxy*0.60, replace = T), 
                      dim = "msp")
dat.msp$code <- substr(dat.msp$code, start = 1, stop = 3)
dat.msp$code[dat.msp$code==""] <- NA
dat.msp <- na.omit(dat.msp)

## DINPIN
dat.din <- data.frame(id = sample(id, size = n.proxy*0.30, replace = T), 
                      code = sample(proxies$din, size = n.proxy*0.30, replace = T), 
                      dim = "din")
dat.din$code[dat.din$code==""] <- NA
dat.din <- na.omit(dat.din)

# Proxy data
dat.proxy <- rbind(dat.diag, dat.proc, dat.msp, dat.din)

# Drop missing codes 
dat.proxy <- na.omit(dat.proxy)
table(dat.proxy$dim)

# Save data
save(simdat, dat.proxy, file = "Data/simdata.RData")

