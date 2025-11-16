# GROWTH MODEL SELECTION

# We use the original formulation of the von Bertalanffy model, in which the length at age 0 (length at glass eel arrival) can be specified

GM1 <- nlme(Length_cm ~ VBGFOrig(Age, Linf, k, L0 = 6.66),
            data = LengthAgeObs,
            fixed = list(Linf ~ fSex, 
                         k ~ fSex),
            random = Linf + k ~ fSex|fSite,
            method = "ML",
            start = list(fixed = c(70, -30, 0.3, 0.7)),   
            control = nlmeControl(maxIter = 10000, msMaxIter = 10000, pnlsMaxIter = 5000, niterEM = 10000, tolerance = 0.005, msTol = 0.005)
)
summary(GM1)

# Model diagnostics
Coefs <- coef(GM1)    # get VBGF coefficients for each individual
Coefs <- data.frame(Site = rownames(Coefs), 
                    Linf_Intercept = Coefs[,1],
                    Linf_M = Coefs[,2],
                    k_Intercept = Coefs[,3],
                    k_M = Coefs[,4])
Residuals <- resid(GM1, type = "normalized")
Diagnostics <- data.frame(Residuals = Residuals,   # get the normalized residuals
                          Fitted = fitted(GM1),    # get the fitted values
                          fSite = LengthAgeObs$fSite,
                          fSex = LengthAgeObs$fSex,
                          Age = LengthAgeObs$Age)
Diagnostics$Linf <- Coefs$Linf[match(Diagnostics$fSite, Coefs$fSite)]
Diagnostics$k <- Coefs$k[match(Diagnostics$fSite, Coefs$fSite)]

# Plot residuals vs fitted values
ggplot(Diagnostics, aes(x = Fitted, y = Residuals)) +
  geom_point() + geom_smooth() + theme_bw() +
  xlab("Fitted values") + ylab ("Residuals")    # slight trumpet effect and positive residuals at low fitted values, while negative residuals at high fitted values
# Plot residuals by the random term (= Site)
ggplot(Diagnostics, aes(x = fSite, y = Residuals)) +
  geom_boxplot() + geom_hline(yintercept = 0) + theme_bw() +
  xlab("Site") + ylab ("Residuals") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # looks ok

# Plot residuals by covariates (Age, Linf, k, t0)
ggplot(Diagnostics, aes(x = Age, y = Residuals)) +
  geom_point() + geom_smooth() + theme_bw() +
  xlab("Age") + ylab ("Residuals")    # roughly ok
ggplot(Diagnostics, aes(x = fSex, y = Residuals)) +
  geom_boxplot() + geom_hline(yintercept = 0) + theme_bw() +
  xlab("Site") + ylab ("Residuals")   # heterogeneous

# --> Heterogeneity is most problematic between the sexes - try to fit another model applying variance weights by sex

# Model 2: Adding variance weights by sex
Weights <- varIdent(form = ~ 1 | fSex)
GM2 <- update(GM1, weights = Weights)

AIC(GM1, GM2)   # weighted growth model has lower AIC

# Model diagnostics
Coefs <- coef(GM2)    # get VBGF coefficients for each individual
Coefs <- data.frame(Site = rownames(Coefs), 
                    Linf_Intercept = Coefs[,1],
                    Linf_M = Coefs[,2],
                    k_Intercept = Coefs[,3],
                    k_M = Coefs[,4])
Residuals <- resid(GM2, type = "normalized")
Diagnostics <- data.frame(Residuals = Residuals,   # get the normalized residuals
                          Fitted = fitted(GM2),    # get the fitted values
                          fSite = LengthAgeObs$fSite,
                          fSex = LengthAgeObs$fSex,
                          Age = LengthAgeObs$Age)
Diagnostics$Linf <- Coefs$Linf[match(Diagnostics$fSite, Coefs$fSite)]
Diagnostics$k <- Coefs$k[match(Diagnostics$fSite, Coefs$fSite)]

# Plot residuals vs fitted values
ggplot(Diagnostics, aes(x = Fitted, y = Residuals)) +
  geom_point() + geom_smooth() + theme_bw() +
  xlab("Fitted values") + ylab ("Residuals")    # trumpet effect has vanished and slightly better at the higher end of fitted values
# Plot residuals by the random term (= Fish ID)
ggplot(Diagnostics, aes(x = fSite, y = Residuals)) +
  geom_boxplot() + geom_hline(yintercept = 0) + theme_bw() +
  xlab("Site") + ylab ("Residuals") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))    # better - more homogeneity across sites

# Plot residuals by covariates (Age, Linf, k, t0)
ggplot(Diagnostics, aes(x = Age, y = Residuals)) +
  geom_point() + geom_smooth() + theme_bw() +
  xlab("Age") + ylab ("Residuals")    # similar
ggplot(Diagnostics, aes(x = fSex, y = Residuals)) +
  geom_boxplot() + geom_hline(yintercept = 0) + theme_bw() +
  xlab("Site") + ylab ("Residuals")   # looks good now!

# --> the variance-weighted model better meets the assumption of homoscedasticity and gets better goodness-of-fit statistic (AIC)


# Check if the sampling scheme (i.e. mixed stages vs. silver eel only) significantly influences Linf or k estimates
GM3 <- nlme(Length_cm ~ VBGFOrig(Age, Linf, k, L0 = 6.66),
            data = LengthAgeObs,
            fixed = list(Linf ~ fSex + StageSelectivity,
                         k ~ fSex + StageSelectivity),    
            random = Linf + k ~ fSex|fSite,
            method = "ML",
            weights = Weights,
            start = list(fixed = c(70, -30, 0, 0.3, 0.7, 0)), 
            control = nlmeControl(maxIter = 2000, msMaxIter = 2000, pnlsMaxIter = 500, niterEM = 1000, tolerance = 0.005, msTol = 0.005)
)
summary(GM3)
# significant effect of the stage selectivity variable only on k, not Linf

# Refit the model letting the sampling selectivity variable only affect k
GM4 <- nlme(Length_cm ~ VBGFOrig(Age, Linf, k, L0 = 6.66),
            data = LengthAgeObs,
            fixed = list(Linf ~ fSex,
                         k ~ fSex + StageSelectivity),    
            random = Linf + k ~ fSex|fSite,
            method = "ML",
            weights = Weights,
            start = list(fixed = c(70, -30, 0.3, 0.7, 0.2)), 
            control = nlmeControl(maxIter = 1000, msMaxIter = 500, pnlsMaxIter = 50, niterEM = 500, tolerance = 0.005, msTol = 0.005)
)
summary(GM4)

AIC(GM2, GM3, GM4) %>% arrange(AIC)

# --> given the significant influence of sampling scheme on k estimates and the lowest AIC, we need to account for sampling scheme as a covariate affecting k

# Model diagnostics
Coefs <- coef(GM4)    # get VBGF coefficients for each individual
Coefs <- data.frame(Site = rownames(Coefs), 
                    Linf_Intercept = Coefs[,1],
                    Linf_M = Coefs[,2],
                    k_Intercept = Coefs[,3],
                    k_M = Coefs[,4],
                    k_StageSelectivity = Coefs[,5])
Residuals <- resid(GM4, type = "normalized")
Diagnostics <- data.frame(Residuals = Residuals,   # get the normalized residuals
                          Fitted = fitted(GM4),    # get the fitted values
                          fSite = LengthAgeObs$fSite,
                          fSex = LengthAgeObs$fSex,
                          Age = LengthAgeObs$Age,
                          StageSelectivity = LengthAgeObs$StageSelectivity)
Diagnostics$Linf <- Coefs$Linf[match(Diagnostics$fSite, Coefs$fSite)]
Diagnostics$k <- Coefs$k[match(Diagnostics$fSite, Coefs$fSite)]

# Plot residuals vs fitted values
ggplot(Diagnostics, aes(x = Fitted, y = Residuals)) +
  geom_point() + geom_smooth() + theme_bw() +
  xlab("Fitted values") + ylab ("Residuals")    # similar
# Plot residuals by the random term (= Fish ID)
ggplot(Diagnostics, aes(x = fSite, y = Residuals)) +
  geom_boxplot() + geom_hline(yintercept = 0) + theme_bw() +
  xlab("Site") + ylab ("Residuals") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))    # looks good

# Plot residuals by covariates (Age, Linf, k, t0)
ggplot(Diagnostics, aes(x = Age, y = Residuals)) +
  geom_point() + geom_smooth() + theme_bw() +
  xlab("Age") + ylab ("Residuals")    # similar
ggplot(Diagnostics, aes(x = fSex, y = Residuals)) +
  geom_boxplot() + geom_hline(yintercept = 0) + theme_bw() +
  xlab("Site") + ylab ("Residuals")   # looks ok now!
ggplot(Diagnostics, aes(x = StageSelectivity, y = Residuals)) +
  geom_boxplot() + geom_hline(yintercept = 0) + theme_bw() +
  xlab("Site") + ylab ("Residuals")   # looks good

# --> Conclusion: the model using the original formulation of the von B. function with a fixed L0, sex-specific variance weights and accounting for stage-selectivity in sampling is the best model
# We will retain this as the base model!

GMBase <- GM4
save(GMBase, file = "./Output/GMBase.rdata")
