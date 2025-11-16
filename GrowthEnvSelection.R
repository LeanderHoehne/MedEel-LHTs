# Identify the environmental correlates of eel growth in Mediterranean sites

set.seed(2025)

# Fit alternative versions of the global model including all covariates to be considered
GMEnvBase <- nlme(Length_cm ~ VBGFOrig(Age, Linf, k, L0 = 6.66),
                       data = LengthAgeObs,
                       fixed = list(Linf ~ fSex + Chlorophyll + TP + DistanceGibraltar,
                                    k ~ fSex + logSurfaceArea + Temperature + Salinity + StageSelectivity),   # using original values of salinity
                       random = Linf + k ~ fSex|fSite,
                       method = "ML",
                       weights = varIdent(form = ~ 1 | fSex),
                       start = list(fixed = c(70, -30, 0.5, -0.5, 1, 0.3, 0.7, 0, 0, 0, 0.2)),
                       control = nlmeControl(maxIter = 10000, msMaxIter = 10000, pnlsMaxIter = 5000, niterEM = 10000, tolerance = 0.005, msTol = 0.005)
)
GMEnvBase2 <- nlme(Length_cm ~ VBGFOrig(Age, Linf, k, L0 = 6.66),
                        data = LengthAgeObs,
                        fixed = list(Linf ~ fSex + Chlorophyll + TP + DistanceGibraltar,
                                     k ~ fSex + logSurfaceArea + Temperature + sqSalinity + StageSelectivity),    # using squared values of salinity
                        random = Linf + k ~ fSex|fSite,
                        method = "ML",
                        weights = varIdent(form = ~ 1 | fSex),
                        start = list(fixed = c(70, -30, 0.5, -0.5, 1, 0.3, 0.7, 0, 0, 0, 0.2)),   
                        control = nlmeControl(maxIter = 5000, msMaxIter = 2000, pnlsMaxIter = 1000, niterEM = 5000, tolerance = 0.005, msTol = 0.005)
)
GMEnvBase3 <- nlme(Length_cm ~ VBGFOrig(Age, Linf, k, L0 = 6.66),
                  data = LengthAgeObs,
                  fixed = list(Linf ~ fSex + Chlorophyll + TP + logSurfaceArea + DistanceGibraltar,
                               k ~ fSex + Temperature + Salinity + StageSelectivity),
                  random = Linf + k ~ fSex|fSite,
                  method = "ML",
                  weights = varIdent(form = ~ 1 | fSex),
                  start = list(fixed = c(70, -30, 0.5, -0.5, 1, 0.3, 0.7, 0, 0, 0, 0.2)),  
                  control = nlmeControl(maxIter = 2000, msMaxIter = 1000, pnlsMaxIter = 500, niterEM = 1000, tolerance = 0.005, msTol = 0.005)
)
GMEnvBase4 <- nlme(Length_cm ~ VBGFOrig(Age, Linf, k, L0 = 6.66),
                   data = LengthAgeObs,
                   fixed = list(Linf ~ fSex + Chlorophyll + TP + logSurfaceArea + DistanceGibraltar,
                                k ~ fSex + Temperature + sqSalinity + StageSelectivity),
                   random = Linf + k ~ fSex|fSite,
                   method = "ML",
                   weights = varIdent(form = ~ 1 | fSex),
                   start = list(fixed = c(70, -30, 0.5, -0.5, 1, 0.3, 0.7, 0, 0, 0, 0.2)),  
                   control = nlmeControl(maxIter = 2000, msMaxIter = 1000, pnlsMaxIter = 500, niterEM = 1000, tolerance = 0.005, msTol = 0.005)
)
AIC(GMEnvBase, GMEnvBase2, GMEnvBase3, GMEnvBase4)
# --> using untransformed salinity values and modeling surface area on Linf improves model fit, therefore we will conduct the model selection using that configuration

# Does a model with log-transformed values of chlorophyll and TP perform better?
GMEnvBase5 <- nlme(Length_cm ~ VBGFOrig(Age, Linf, k, L0 = 6.66),
                   data = LengthAgeObs,
                   fixed = list(Linf ~ fSex + logChloro + logTP + logSurfaceArea + DistanceGibraltar,
                                k ~ fSex + Temperature + sqSalinity + StageSelectivity),
                   random = Linf + k ~ fSex|fSite,
                   method = "ML",
                   weights = varIdent(form = ~ 1 | fSex),
                   start = list(fixed = c(70, -30, 0.5, -0.5, 1, 0.3, 0.7, 0, 0, 0, 0.2)),  
                   control = nlmeControl(maxIter = 2000, msMaxIter = 1000, pnlsMaxIter = 500, niterEM = 1000, tolerance = 0.005, msTol = 0.005)
)
AIC(GMEnvBase3, GMEnvBase5)
# --> No! We will thus run the analysis with chlorophyll and TP values on the original scale.


# Backward variable selection

# Step 1
# GM1aObs: Full model minus Chlorophyll
GM1aObs <- update(GMEnvBase3, fixed = list(Linf ~ fSex + TP + logSurfaceArea + DistanceGibraltar,
                                        k ~ fSex + Temperature + Salinity + StageSelectivity),
                  start = list(fixed = c(70, -30, -0.5, 0, 1, 0.3, 0.7, 0, 0, 0.2)),
                  control = nlmeControl(maxIter = 5000, msMaxIter = 2000, pnlsMaxIter = 1000, niterEM = 5000, tolerance = 0.005, msTol = 0.005)
)

# GM1bObs: Full model minus TP
GM1bObs <- update(GMEnvBase3, fixed = list(Linf ~ fSex + Chlorophyll + logSurfaceArea + DistanceGibraltar,
                                           k ~ fSex + Temperature + Salinity + StageSelectivity),
                  start = list(fixed = c(70, -30, 0.5, 0, 1, 0.3, 0.7, 0, 0, 0.2)))

# GM1cObs: Full model minus DistanceGibraltar
GM1cObs <- update(GMEnvBase3, fixed = list(Linf ~ fSex + Chlorophyll + TP + logSurfaceArea,
                                           k ~ fSex + Temperature + Salinity + StageSelectivity),
                  start = list(fixed = c(70, -30, 0.5, -0.5, 0, 0.3, 0.7, 0, 0, 0.2)))

# GM1dObs: Full model minus SurfaceArea
GM1dObs <- update(GMEnvBase3, fixed = list(Linf ~ fSex + Chlorophyll + TP + DistanceGibraltar,
                                           k ~ fSex + Temperature + Salinity + StageSelectivity),
                  start = list(fixed = c(70, -30, 0.5, -0.5, 1, 0.3, 0.7, 0, 0, 0.2)))

# GM1eObs: Full model minus Temperature
GM1eObs <- update(GMEnvBase3, fixed = list(Linf ~ fSex + Chlorophyll + TP + logSurfaceArea + DistanceGibraltar,
                                           k ~ fSex + Salinity + StageSelectivity),
                  start = list(fixed = c(70, -30, 0.5, -0.5, 0, 1, 0.3, 0.7, 0, 0.2)))

# GM1fObs: Full model minus Salinity
GM1fObs <- update(GMEnvBase3, fixed = list(Linf ~ fSex + Chlorophyll + TP + logSurfaceArea + DistanceGibraltar,
                                           k ~ fSex + Temperature + StageSelectivity),
                  start = list(fixed = c(70, -30, 0.5, -0.5, 0, 1, 0.3, 0.7, 0, 0.2)))

AIC(GMEnvBase, GM1aObs, GM1bObs, GM1cObs, GM1dObs, GM1eObs, GM1fObs) %>% arrange(AIC)
# --> elimination of SurfaceArea leads to the largest drop in AIC, thus will be excluded in further selection


# Step 2
# GM2aObs: Full model minus SurfaceArea and Chlorophyll
GM2aObs <- update(GM1dObs, fixed = list(Linf ~ fSex + TP + DistanceGibraltar,
                                        k ~ fSex + Temperature + Salinity + StageSelectivity),
                     start = list(fixed = c(70, -30, -0.5, 1, 0.3, 0.7, 0, 0, 0.2)))

# GM2bObs: Full model minus SurfaceArea and TP
GM2bObs <- update(GM1dObs, fixed = list(Linf ~ fSex + Chlorophyll + DistanceGibraltar,
                                        k ~ fSex + Temperature + Salinity + StageSelectivity),
                  start = list(fixed = c(70, -30, 0.5, 1, 0.3, 0.7, 0, 0, 0.2)))

# GM2cObs: Full model minus SurfaceArea and DistanceGibraltar
GM2cObs <- update(GM1dObs, fixed = list(Linf ~ fSex + Chlorophyll + TP,
                                        k ~ fSex + Temperature + Salinity + StageSelectivity),
                  start = list(fixed = c(70, -30, 0.5, -0.5, 0.3, 0.7, 0, 0, 0.2)))

# GM2dObs: Full model minus SurfaceArea and Temperature
GM2dObs <- update(GM1dObs, fixed = list(Linf ~ fSex + Chlorophyll + TP + DistanceGibraltar,
                                        k ~ fSex + Salinity + StageSelectivity),
                  start = list(fixed = c(70, -30, 0.5, -0.5, 1, 0.3, 0.7, 0, 0.2)))

# GM2eObs: Full model minus SurfaceArea and Salinity
GM2eObs <- update(GM1dObs, fixed = list(Linf ~ fSex + Chlorophyll + TP + DistanceGibraltar,
                                        k ~ fSex + Temperature + StageSelectivity),
                  start = list(fixed = c(70, -30, 0.5, -0.5, 1, 0.3, 0.7, 0, 0.2)))

AIC(GM1dObs, GM2aObs, GM2bObs, GM2cObs, GM2dObs, GM2eObs) %>% arrange(AIC)
# --> elimination of Temperature leads to the largest drop in AIC, thus will be excluded in further selection


# Step 3
# GM3aObs: Full model minus SurfaceArea, Temperature and Chlorophyll
GM3aObs <- update(GM2dObs, fixed = list(Linf ~ fSex + TP + DistanceGibraltar,
                                        k ~ fSex + Salinity + StageSelectivity),
                  start = list(fixed = c(70, -30, -0.5, 1, 0.3, 0.7, 0, 0.2)))

# GM3bObs: Full model minus SurfaceArea, Temperature and TP
GM3bObs <- update(GM2dObs, fixed = list(Linf ~ fSex + Chlorophyll + DistanceGibraltar,
                                        k ~ fSex + Salinity + StageSelectivity),
                  start = list(fixed = c(70, -30, 0.5, 1, 0.3, 0.7, 0, 0.2)))

# GM3cObs: Full model minus SurfaceArea, Temperature and DistanceGibraltar
GM3cObs <- update(GM2dObs, fixed = list(Linf ~ fSex + Chlorophyll + TP,
                                        k ~ fSex + Salinity + StageSelectivity),
                  start = list(fixed = c(70, -30, 0.5, -0.5, 0.3, 0.7, 0, 0.2)))

# GM3dObs: Full model minus SurfaceArea, Temperature and Salinity
GM3dObs <- update(GM2dObs, fixed = list(Linf ~ fSex + Chlorophyll + TP + DistanceGibraltar,
                                        k ~ fSex + StageSelectivity),
                  start = list(fixed = c(70, -30, 0.5, -0.5, 1, 0.3, 0.7, 0.2)))

AIC(GM2dObs, GM3aObs, GM3bObs, GM3cObs, GM3dObs) %>% arrange(AIC)
# --> Chlorophyll should be dropped!


# Step 4
# GM4aObs: Full model minus SurfaceArea, Temperature, Chlorophyll and TP
GM4aObs <- update(GM3aObs, fixed = list(Linf ~ fSex + DistanceGibraltar,
                                        k ~ fSex + Salinity + StageSelectivity),
                  start = list(fixed = c(70, -30, 1, 0.3, 0.7, 0, 0.2)))

# GM4bObs: Full model minus SurfaceArea, Temperature, Chlorophyll and DistanceGibraltar
GM4bObs <- update(GM3aObs, fixed = list(Linf ~ fSex + TP,
                                        k ~ fSex + Salinity + StageSelectivity),
                  start = list(fixed = c(70, -30, -0.5, 0.3, 0.7, 0, 0.2)))

# GM4cObs: Full model minus SurfaceArea, Temperature, Chlorophyll and Salinity
GM4cObs <- update(GM3aObs, fixed = list(Linf ~ fSex + TP + DistanceGibraltar,
                                        k ~ fSex + StageSelectivity),
                  start = list(fixed = c(70, -30, -0.5, 1, 0.3, 0.7, 0.2)))

AIC(GM3aObs, GM4aObs, GM4bObs, GM4cObs) %>% arrange(AIC)
# --> Salinity should be dropped!


# Step 5
# GM5aObs: Full model minus SurfaceArea, Temperature, Chlorophyll, Salinity and TP
GM5aObs <- update(GM4cObs, fixed = list(Linf ~ fSex + DistanceGibraltar,
                                        k ~ fSex + StageSelectivity),
                  start = list(fixed = c(70, -30, 1, 0.3, 0.7, 0.2)),
                  control = nlmeControl(maxIter = 5000, msMaxIter = 2000, pnlsMaxIter = 1000, niterEM = 5000, tolerance = 0.005, msTol = 0.005)
)

# GM5bObs: Full model minus SurfaceArea, Temperature, Chlorophyll, Salinity and DistanceGibraltar
GM5bObs <- update(GM4cObs, fixed = list(Linf ~ fSex + TP,
                                        k ~ fSex + StageSelectivity),
                  start = list(fixed = c(70, -30, -0.5, 0.3, 0.7, 0.2)))

AIC(GM4cObs, GM5aObs, GM5bObs) %>% arrange(AIC)
# --> TP should be excluded!


# Step 6
GM6aObs <- update(GM5aObs, fixed = list(Linf ~ fSex,
                                        k ~ fSex + StageSelectivity),
                  start = list(fixed = c(70, -30, 0.4, 0.8, 0)))
AIC(GM5aObs, GM6aObs)
# --> only DistanceGibraltar is retained as an environmental covariate in the selected model!


# SENSITIVITY ANALYSIS: Repeat model selection with the extended dataset to account for eventual bias by selective subsampling of aged individuals

# Log-transform surface area values
LengthAgeExt$StageSelectivity <- as.factor(ifelse(LengthAgeExt$fSite %in% nObsLengthAge$fSite[which(nObsLengthAge$Selective == 1)], 1, 0))

GMEnvBaseExt <- nlme(Length_cm ~ VBGFOrig(Age, Linf, k, L0 = 6.66),
                  data = LengthAgeExt,
                  fixed = list(Linf ~ fSex + Chlorophyll + TP + logSurfaceArea + DistanceGibraltar,
                               k ~ fSex + Temperature + Salinity + StageSelectivity),   
                  random = Linf + k ~ fSex|fSite,
                  method = "ML",
                  weights = varIdent(form = ~ 1 | fSex),
                  start = list(fixed = c(70, -30, 0.5, -0.5, 0, 1, 0.3, 0.7, 0, 0, 0.2)),
                  control = nlmeControl(maxIter = 10000, msMaxIter = 10000, pnlsMaxIter = 5000, niterEM = 10000, tolerance = 0.005, msTol = 0.005)
)

# Backward variable selection

# Step 1
# GM1aExt: Full model minus Chlorophyll
GM1aExt <- update(GMEnvBaseExt, fixed = list(Linf ~ fSex + TP + logSurfaceArea + DistanceGibraltar,
                                           k ~ fSex + Temperature + Salinity + StageSelectivity),
                  start = list(fixed = c(70, -30, 0, 1, 0, 0.3, 0.7, 0, 0, 0.2)))

# GM1bExt: Full model minus TP
GM1bExt <- update(GMEnvBaseExt, fixed = list(Linf ~ fSex + Chlorophyll + logSurfaceArea + DistanceGibraltar,
                                            k ~ fSex + Temperature + Salinity + StageSelectivity),
                   start = list(fixed = c(70, -30, 0.5, 0, 1, 0.3, 0.7, 0, 0, 0.2))
)

# GM1cExt: Full model minus DistanceGibraltar
GM1cExt <- update(GMEnvBaseExt, fixed = list(Linf ~ fSex + Chlorophyll + TP + logSurfaceArea,
                                           k ~ fSex + Temperature + Salinity + StageSelectivity),
                  start = list(fixed = c(70, -30, 0.5, -0.5, 0, 0.3, 0.7, 0, 0, 0.2)))


# GM1dExt: Full model minus SurfaceArea
GM1dExt <- update(GMEnvBaseExt, fixed = list(Linf ~ fSex + Chlorophyll + TP + DistanceGibraltar,
                                           k ~ fSex + Temperature + Salinity + StageSelectivity),
                  start = list(fixed = c(70, -30, 0.5, -0.5, 1, 0.3, 0.7, 0, 0, 0.2)))

# GM1eExt: Full model minus Temperature
GM1eExt <- update(GMEnvBaseExt, fixed = list(Linf ~ fSex + Chlorophyll + TP + logSurfaceArea + DistanceGibraltar,
                                           k ~ fSex + Salinity + StageSelectivity),
                  start = list(fixed = c(70, -30, 0.5, -0.5, 0, 1, 0.3, 0.7, 0, 0.2)))

# GM1fExt: Full model minus Salinity
GM1fExt <- update(GMEnvBaseExt, fixed = list(Linf ~ fSex + Chlorophyll + TP + logSurfaceArea + DistanceGibraltar,
                                           k ~ fSex + Temperature + StageSelectivity),
                  start = list(fixed = c(70, -30, 0.5, -0.5, 0, 1, 0.3, 0.7, 0, 0.2)))

AIC(GMEnvBaseExt, GM1aExt, GM1bExt, GM1cExt, GM1dExt, GM1eExt, GM1fExt) %>% arrange(AIC)
# --> elimination of Salinity leads to the largest drop in AIC, thus will be excluded in further selection


# Step 2

# Minus Salinity and Chlorophyll
GM2aExt <- update(GM1fExt, fixed = list(Linf ~ fSex + TP + logSurfaceArea + DistanceGibraltar,
                                             k ~ fSex + Temperature + StageSelectivity),
                  start = list(fixed = c(70, -30, -0.5, 0, 1, 0.3, 0.7, 0, 0.2)))

# Minus Salinity and TP
GM2bExt <- update(GM1fExt, fixed = list(Linf ~ fSex + Chlorophyll + logSurfaceArea + DistanceGibraltar,
                                             k ~ fSex + Temperature + StageSelectivity),
                  start = list(fixed = c(70, -30, 0.5, 0, 1, 0.3, 0.7, 0, 0.2)))

# Minus Salinity and logSurfaceArea
GM2cExt <- update(GM1fExt, fixed = list(Linf ~ fSex + Chlorophyll + TP + DistanceGibraltar,
                                             k ~ fSex + Temperature + StageSelectivity),
                  start = list(fixed = c(70, -30, 0.5, -0.5, 1, 0.3, 0.7, 0, 0.2)))

# Minus Salinity and DistanceGibraltar
GM2dExt <- update(GM1fExt, fixed = list(Linf ~ fSex + Chlorophyll + TP + logSurfaceArea,
                                             k ~ fSex + Temperature + StageSelectivity),
                  start = list(fixed = c(70, -30, 0.5, -0.5, 0, 0.3, 0.7, 0, 0.2)))

# Minus Salinity and Temperature
GM2eExt <- update(GM1fExt, fixed = list(Linf ~ fSex + Chlorophyll + TP + logSurfaceArea + DistanceGibraltar,
                                             k ~ fSex + StageSelectivity),
                  start = list(fixed = c(70, -30, 0.5, -0.5, 0, 1, 0.3, 0.7, 0.2)))

AIC(GM1dExt, GM2aExt, GM2bExt, GM2cExt, GM2dExt, GM2eExt) %>% arrange(AIC)
# --> Surface area should be dropped!


# Step 3

# Minus Salinity, logSurfaceArea and Chlorophyll
GM3aExt <- update(GM2cExt, fixed = list(Linf ~ fSex + TP + DistanceGibraltar,
                                        k ~ fSex + Temperature + StageSelectivity),
                  start = list(fixed = c(70, -30, -0.5, 1, 0.3, 0.7, 0, 0.2)))

# Minus Salinity, logSurfaceArea and TP
GM3bExt <- update(GM2cExt, fixed = list(Linf ~ fSex + Chlorophyll + DistanceGibraltar,
                                        k ~ fSex + Temperature + StageSelectivity),
                  start = list(fixed = c(70, -30, 0.5, 1, 0.3, 0.7, 0, 0.2)))

# Minus Salinity, logSurfaceArea and DistanceGibraltar
GM3cExt <- update(GM2cExt, fixed = list(Linf ~ fSex + Chlorophyll + TP,
                                        k ~ fSex + Temperature + StageSelectivity),
                  start = list(fixed = c(70, -30, 0.5, -0.5, 0.3, 0.7, 0, 0.2)))

# Minus Salinity, logSurfaceArea and Temperature
GM3dExt <- update(GM2cExt, fixed = list(Linf ~ fSex + Chlorophyll + TP + DistanceGibraltar,
                                        k ~ fSex + StageSelectivity),
                  start = list(fixed = c(70, -30, 0.5, -0.5, 1, 0.3, 0.7, 0.2)))

AIC(GM2eExt, GM3aExt, GM3bExt, GM3cExt, GM3dExt) %>% arrange(AIC)
# --> Temperature should be dropped!


# Step 4

# Minus Temperature, Salinity, logSurfaceArea and Chlorophyll
GM4aExt <- update(GM3dExt, fixed = list(Linf ~ fSex + TP + DistanceGibraltar,
                                        k ~ fSex + StageSelectivity),
                  start = list(fixed = c(70, -30, -0.5, 1, 0.3, 0.7, 0.2)))

# Minus Temperature, Salinity, logSurfaceArea and TP
GM4bExt <- update(GM3dExt, fixed = list(Linf ~ fSex + Chlorophyll + DistanceGibraltar,
                                        k ~ fSex + StageSelectivity),
                  start = list(fixed = c(70, -30, 0.5, 1, 0.3, 0.7, 0.2)))

# Minus Temperature, Salinity, logSurfaceArea and DistanceGibraltar
GM4cExt <- update(GM3dExt, fixed = list(Linf ~ fSex + Chlorophyll + TP,
                                        k ~ fSex + StageSelectivity),
                  start = list(fixed = c(70, -30, 0.5, -0.5, 0.3, 0.7, 0.2)))

AIC(GM3dExt, GM4aExt, GM4bExt, GM4cExt) %>% arrange(AIC)
# --> Chlorophyll should be excluded!


# Step 5
# Minus Temperature, Salinity, logSurfaceArea, Chlorophyll and TP
GM5aExt <- update(GM4aExt, fixed = list(Linf ~ fSex + DistanceGibraltar,
                                        k ~ fSex + StageSelectivity),
                  start = list(fixed = c(70, -30, 1, 0.3, 0.7, 0)))

# Minus Temperature, Salinity, logSurfaceArea, Chlorophyll and DistanceGibraltar
GM5bExt <- update(GM4aExt, fixed = list(Linf ~ fSex + TP,
                                        k ~ fSex + StageSelectivity),
                  start = list(fixed = c(70, -30, -0.5, 0.3, 0.7, 0)))

AIC(GM4aExt, GM5aExt, GM5bExt) %>% arrange(AIC)
# --> TP should be excluded!


# Step 6
GM6aExt <- update(GM5aExt, fixed = list(Linf ~ fSex,
                                        k ~ fSex + StageSelectivity),
                  start = list(fixed = c(70, -30, 0.4, 0.8, 0)))
AIC(GM5aExt, GM6aExt)
# --> only DistanceGibraltar is retained in the selected model!


# --> Conclusion: using both the dataset of observed ages and the extended one applying the age-length key lead to the same result, i.e. the results of the observed ages dataset are robust!

GMEnvObs <- GM5aObs
save(GMEnvObs, file = "./Output/GMEnvObs.rdata")

save(LengthAgeExt, file = "Growth_env_data.rdata")


# SENSITIVITY ANALYSIS 2: CHECK ROBUSTNESS OF DISTANCE-FROM-GIBRALTAR EFFECT TO THE OUTLIER SITE "GUADIARO"

LengthAgeObs2 <- LengthAgeObs %>% filter(fSite != "Guadiaro")
GMObsSens <- update(GMEnvObs, data = LengthAgeObs2)
summary(GMObsSens)
# --> the effect of distance-from-Gibraltar is robust to the inclusion or exclusion of the site "Guadiaro"!


# MODEL VALIDATION

# 1) Check for spatial autocorrelation
LengthAgeObs$Resid <- residuals(GMEnvObs, type = "pearson")
LengthAgeObs$SignResid <- ifelse(LengthAgeObs$Resid >= 0, "positive", "negative")
LengthAgeObs$Size <- round(5 * sqrt(abs(LengthAgeObs$Resid) /max(LengthAgeObs$Resid))) + 1 

LengthAgeObs %>% ggplot(aes(x = Longitude, y = Latitude, size = Size, col = SignResid)) +
  geom_jitter(shape = 1, width = 1.5) +
  xlab("Longitude") + ylab("Latitude") +
  facet_wrap(~fSex, scales = "free") +
  theme_bw() +
  theme(legend.position = "bottom")

# --> no clear pattern of autocorrelation in residuals; we will confirm this with a variogram!

# Convert data frame to an sf object with geographic coordinates
LengthAge_sf <- st_as_sf(LengthAgeObs, coords = c("Longitude", "Latitude"), crs = 4326)  # Use EPSG:4326 (WGS 84)

# Convert to projected coordinates (e.g., UTM zone 32N) to ensure distance calculations are correct
LengthAge_sf <- st_transform(LengthAge_sf, crs = 32632)  # Change to your appropriate UTM zone

LengthAge_df <- as.data.frame(LengthAge_sf)

# Extract coordinates as a separate dataframe
LengthAge_df$X <- st_coordinates(LengthAge_sf)[,1]  # Easting
LengthAge_df$Y <- st_coordinates(LengthAge_sf)[,2]  # Northing

# Create a gstat variogram
vario <- variogram(Resid ~ 1, locations = ~X + Y, data = LengthAge_df)

# Plot the variogram
plot(vario)
# --> does not suggest any spatial autocorrelation!
