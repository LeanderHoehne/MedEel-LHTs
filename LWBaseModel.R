# 3.2a) Fit the base mixed-effects model of length-weight-relationship, not accounting for any covariates
LWMBase <- lmer(logWeight_g ~ logLength_cm + (1 + logLength_cm|fSite), data = LHT_HAB, control = lmerControl(optimizer = "bobyqa"))     # fit the model first in lme4 package, so we can use the outlier detection functions from the performance package

# Model diagnostics

# Detection of outliers / influential observations
OutliersLW <- check_outliers(LWMBase, method = "cook")
LHT_HAB$OutlierLW <- OutliersLW
# Visualize the detected outliers
LHT_HAB %>% ggplot(aes(x = logLength_cm)) +
  geom_point(aes(color = OutlierLW, y = logWeight_g))


# Refit the model without the outliers
LHT_HAB2 <- LHT_HAB %>% filter(OutlierLW == 0)
LWMBase2 <- glmmTMB(logWeight_g ~ logLength_cm + (1 + logLength_cm|fSite), data = LHT_HAB2)

# does the Gaussian model predict negative values?
LHT_HAB2$FittedLW <- exp(fitted(LWMBase2))
range(LHT_HAB2$FittedLW)
# --> no negative fitted values. We can keep the Gaussian model.


# Residuals vs. fitted values and vs. predictor variables
Diagnostics <- data.frame(Residuals = resid(LWMBase2, type = "pearson"),   # get the normalized residuals
                          Fitted = fitted(LWMBase2),    # get the fitted values
                          fSite = LHT_HAB2$fSite)

# Plot residuals vs fitted values
ggplot(Diagnostics, aes(x = Fitted, y = Residuals)) +
  geom_point() + geom_smooth() + theme_bw() +
  xlab("Fitted values") + ylab ("Residuals")    # looks roughly ok

# Try "fixing" a joint intercept across all sites (i.e. the length and weight at hatching), using only random slopes
LWMBase3 <- glmmTMB(logWeight_g ~ logLength_cm + (0 + logLength_cm|fSite), data = LHT_HAB2)
AIC(LWMBase2, LWMBase3)   # the model with random intercepts and slopes performs better, so we stick with it

# Refit the selected model with REML to get more stable parameter and SE estimates
LWMBase2 <- update(LWMBase2, REML = TRUE)


# VISUALIZE OUTPUT BY SITE
BetaLW <- summary(LWMBase2)$coef$cond
round(BetaLW, digits = 3)

# Sketch the fitted values
# 1. Create a grid of covariate values
MatrixRandom <- data.frame(logLength_cm = log(1:MaxLength)) 

# 2. Make a design matrix using the model.matrix function.
Xp.rs <- model.matrix(~logLength_cm, data = MatrixRandom)

# 3. Calculate the predicted values:
MatrixRandom$mu <- Xp.rs %*% BetaLW[,"Estimate"]

# Plot the random effect site
RandomCoefsLW <- ranef(LWMBase2)$cond$fSite; colnames(RandomCoefsLW) <- c("Intercept", "Slope")
RandomCoefsLW <- as.data.frame(RandomCoefsLW) %>% mutate(fSite = row.names(RandomCoefsLW))

SitePredsLW <- data.frame()
for (i in 1:nrow(RandomCoefsLW)) {
  FocalSite <- RandomCoefsLW$fSite[i]
  XVals_Site <- LHT_HAB2$logLength_cm[LHT_HAB2$Site == FocalSite]
  RangeX <- range(XVals_Site)
  SiteOutput <- MatrixRandom %>% filter(logLength_cm > RangeX[1] & logLength_cm < RangeX[2])
  SiteOutput$fSite <- as.factor(FocalSite)
  SiteOutput$fSiteLabel <- as.factor(LHT_HAB2$fSiteLabel[match(FocalSite, LHT_HAB2$fSite)])
  SiteOutput$Predicted <- SiteOutput$mu + RandomCoefsLW$Intercept[i] + SiteOutput$logLength_cm * RandomCoefsLW$Slope[i]
  SitePredsLW <- rbind(SitePredsLW, SiteOutput)
}
PopIntercept <- summary(LWMBase2)$coef$cond[1,1]; PopSlope <- summary(LWMBase2)$coef$cond[2,1]
PopIntSE <- summary(LWMBase2)$coef$cond[1,2]; PopSlopeSE <- summary(LWMBase2)$coef$cond[2,2]
OutputPop <- data.frame(logLength_cm = log(1:MaxLength))
OutputPop$Predicted <- PopIntercept + OutputPop$logLength_cm * PopSlope
OutputPop$Lower <- (PopIntercept - PopIntSE) + OutputPop$logLength_cm * (PopSlope - PopSlopeSE)
OutputPop$Upper <- (PopIntercept + PopIntSE) + OutputPop$logLength_cm * (PopSlope + PopSlopeSE)

# Calculate and store site-level predictions of the base model
FixefLW <- fixef(LWMBase2)$cond
SiteCoefsLW <- data.frame(fSite = RandomCoefsLW$fSite,
                          Intercept = exp(FixefLW[1] + RandomCoefsLW$Intercept),
                          Slope = FixefLW[2] + RandomCoefsLW$Slope)
# Bardawil lagoon consitutes a strong outlier with intercept and slope estimates very different from other sites. We exclude for further analyses.
SiteCoefsLW <- SiteCoefsLW %>% filter(fSite != "Bardawil")

# Add the site acronym to later match sites with the assessment table
SiteCoefsLW$SiteAcronym <- LHT$SiteAcronym[match(SiteCoefsLW$fSite, LHT$fSite)]

# Store site-specific coefficients
save(SiteCoefsLW, file = "./Output/LW_base_site_coefficients.rdata")
save(BetaLW, file = "./Output/LW_base_poplevel_coefs.rdata")
