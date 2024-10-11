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

# Table 1
vars <- c("anyDMD", "sex", "age", "ses", "cci", "year")
tab1 <- CreateTableOne(vars = vars, strata = "mortality_outcome", data = simdat, includeNA = T, 
                       test = F, addOverall = T)
print(tab1, showAllLevels = T)
