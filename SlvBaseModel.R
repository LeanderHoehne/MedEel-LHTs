# Fit the basic mixed-effects model without environmental covariates
SlvBase <- glmmTMB(bSilver ~ Length_cm + fSex + (1|fSite), family = binomial, data = LHTStaged, na.action = na.fail)    # random intercepts model
SlvBase2 <- glmmTMB(bSilver ~ Length_cm + fSex + (1 + Length_cm|fSite), family = binomial, data = LHTStaged, na.action = na.fail)   # random intercepts and random slopes by site model
SlvBase3 <- glmmTMB(bSilver ~ Length_cm + fSex + (1 + Length_cm|fSite/fSex), family = binomial, data = LHTStaged, na.action = na.fail) 
SlvBase4 <- glmmTMB(bSilver ~ Length_cm * fSex + (1 + Length_cm|fSite/fSex), family = binomial, data = LHTStaged, na.action = na.fail) 

summary(SlvBase)
summary(SlvBase2)
summary(SlvBase3)
AIC(SlvBase, SlvBase2, SlvBase3, SlvBase4)    
# --> a model with random slopes and intercepts is best suited!

# Does the sampling gear influence the silvering curves?
LHTStaged$Gear[which(LHTStaged$Gear == "FYK/BAR")] <- "FYK"
LHTStaged$Gear[which(LHTStaged$CountryCode == "MNE")] <- "FYK"
LHTStaged2 <- LHTStaged %>% mutate(Gear = droplevels(as.factor(Gear)))
SlvBase5 <- glmmTMB(bSilver ~ Length_cm + fSex + Gear + (1 + Length_cm|fSite/fSex), family = binomial, data = LHTStaged2, na.action = na.fail) 
summary(SlvBase5)

table(LHTStaged$Gear)
# --> barrier catches should be excluded, but they are a few observations only

LHTStaged2 <- LHTStaged2 %>% filter(Gear != "BAR") %>% mutate(Gear = droplevels(Gear))

SlvBase6 <- glmmTMB(bSilver ~ Length_cm + fSex + Gear + (1 + Length_cm|fSite/fSex), family = binomial, data = LHTStaged2, na.action = na.fail) 
summary(SlvBase6)
# --> Fyke nets and electrofishing can both be used, as there are no systematic differences between the two gear types

# Refit initial model now with the reduced dataset
SlvBase <- glmmTMB(bSilver ~ Length_cm + fSex + (1 + Length_cm|fSite/fSex), family = binomial, data = LHTStaged, na.action = na.fail) 
  
  
# MODEL DIAGNOSTICS

# 1) Check for overdispersion
E1 <- resid(SlvBase, type = "pearson")

# Get the sample size and the number of parameters.
N <- nrow(LHTStaged)
p <- length(coef(SlvBase)) + 1 #' The '+1' is for the sigma

# Determine the dispersion statistic
DispersionStatistic <- sum(E1^2) / (N - p)
DispersionStatistic
# --> Tiny amount of overdispersion? Check again after including environmental covariates!


# 2) Residuals vs fitted values
LHTStaged$Resid <- resid(SlvBase, type = "pearson")
LHTStaged$Fitted <- fitted(SlvBase)

LHTStaged %>% ggplot(aes(x = Fitted, y = Resid)) +
  geom_point(col = "grey", alpha = .5) +
  geom_smooth() +
  theme_bw()    # looks ok

# Residuals vs covariates
LHTStaged %>% ggplot(aes(x = Length_cm, y = Resid)) +
  geom_point(col = "grey", alpha = .5) +
  geom_smooth() +
  theme_bw()    # no clear pattern of heterogeneity
LHTStaged %>% ggplot(aes(x = fSex, y = Resid)) +
  geom_boxplot() +
  theme_bw()

# Store model results
save(SlvBase, file = "./Output/SilveringBaseModel.rdata")


# Extract the coefficients
CoefsSlv <- coef(SlvBase)$cond$`fSex:fSite`
CoefsSlv$Sex <- substr(rownames(CoefsSlv),1,1)
CoefsSlv <- CoefsSlv %>% dplyr::select(-fSexM) %>% rename("Intercept" = `(Intercept)`, Slope = "Length_cm")
CoefsSlv$Site <- sub(".*:", "", rownames(CoefsSlv))
CoefsSlv$SiteAcronym <- LHT$SiteAcronym[match(CoefsSlv$Site, LHT$fSite)]
save(CoefsSlv, file = "Output/SlvCoefs.rdata")
