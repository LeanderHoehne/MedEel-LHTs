# Fit the binomial GLMM without environmental covariates
SRData$bMale <- as.factor(ifelse(SRData$fSex == "M", 1, 0))
SRMBase <- glmmTMB(bMale ~ 1 + (1|fSite), family = binomial, data = SRData)
summary(SRMBase)

# Check for overdispersion
E1 <- resid(SRMBase, type = "pearson")

# Get the sample size and the number of parameters.
N <- nrow(SRData)
p <- length(coef(SRMBase)) + 1 #' The '+1' is for the sigma

# Determine the dispersion statistic
DispersionStatistic <- sum(E1^2) / (N - p)
DispersionStatistic
# --> overdispersion is not a problem!

# Get population-level statistics
LogitEst <- fixef(SRMBase)$cond[1]
exp(LogitEst)   # population-level odds ratio
plogis(LogitEst)     # population-level male fraction
# Calculate the standard error around the global mean male fraction
LogitSE <- summary(SRMBase)$coefficients$cond[2]
LogitSE * plogis(LogitEst) * (1 - plogis(LogitEst))


# Check for sampling covariates that might need to be accounted for, e.g., gear type or meshsize
SRData$Gear <- as.factor(SRData$Gear)
# Fill the meshsize columm out as best as possible and test effect
#SRData$Meshsize_mm[which(SRData$Gear == "ELE")] <- 4    # assign the smallest observed meshsize to the electrofishing gear, given that its selectivity for smaller individuals
SRData$Gear[which(SRData$fSite == "Skadar")] <- "FYK"
SRMBase2 <- glmmTMB(bMale ~ Gear + (1|fSite), family = binomial, data = SRData)
summary(SRMBase2)
SRMBase3 <- glmmTMB(bMale ~ Meshsize_mm + (1|fSite), family = binomial, data = SRData, na.action = na.omit)
summary(SRMBase3)
# --> neither Gear nor mesh size have a significant effect on sex ratio! We do not have to take them into account!

save(SRMBase, file = "Output/Sex_ratio_base_model.Rdata")