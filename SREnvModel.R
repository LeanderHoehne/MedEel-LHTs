# Explore distributions of explanatory variables
dotchart(SRData$Temperature)   
dotchart(SRData$DistanceGibraltar)   # outlier Guadiaro
dotchart(SRData$Salinity)
dotchart(SRData$Chlorophyll)
dotchart(SRData$TP)
dotchart(SRData$logSurfaceArea)
dotchart(SRData$logTP)
dotchart(SRData$logChloro)
# Chlorophyll and TP are candidates for a log-transformation. We will judge the decision whether to transform based on model fit.

# Decide whether to log-transform TP and chlorophyll or not
SRMEnvSens <- glmmTMB(bMale ~ DistanceGibraltar + Salinity + Temperature + Chlorophyll + TP + logSurfaceArea + (1|fSite), family = binomial, data = SRData)
SRMEnvSens2 <- glmmTMB(bMale ~ DistanceGibraltar + Salinity + Temperature + logChloro + logTP + logSurfaceArea + (1|fSite), family = binomial, data = SRData)
SRMEnvSens3 <- glmmTMB(bMale ~ DistanceGibraltar + Salinity + Temperature + Chlorophyll + logTP + logSurfaceArea + (1|fSite), family = binomial, data = SRData)
SRMEnvSens4 <- glmmTMB(bMale ~ DistanceGibraltar + Salinity + Temperature + logChloro + TP + logSurfaceArea + (1|fSite), family = binomial, data = SRData)
AIC(SRMEnvSens, SRMEnvSens2, SRMEnvSens3, SRMEnvSens4)

# --> lowest AIC in model 3, so we will log-transform TP, but not Chlorophyll!


# BACKWARD VARIABLE SELECTION
SRMEnv <- SRMEnvSens3

# Step 1
arrange(drop1(SRMEnv), AIC)   # Temperature should be dropped
SRMEnv1 <- update(SRMEnv, . ~ . -Temperature)

# Step 2
arrange(drop1(SRMEnv1), AIC)    # Chlorophyll should be dropped
SRMEnv2 <- update(SRMEnv1, . ~ . -Chlorophyll)

# Step 3
arrange(drop1(SRMEnv2), AIC)    # Salinity should be dropped
SRMEnv3 <- update(SRMEnv2, . ~ . -Salinity)

# Step 4
arrange(drop1(SRMEnv3), AIC) 
# --> SRMEnv3 is the best model!

SRMEnv <- SRMEnv3

# Get estimates, SEs and p-values for the selected model
summary(SRMEnv)
drop1(SRMEnv, test = "Chi")


# MODEL DIAGNOSTICS

# 1) Check for overdispersion
E1 <- resid(SRMEnv, type = "pearson")

# Get the sample size and the number of parameters.
N <- nrow(SRData)
p <- length(coef(SRMEnv)) + 1 # The '+1' is for the sigma

# Determine the dispersion statistic
DispersionStatistic <- sum(E1^2) / (N - p)
DispersionStatistic
# --> Overdispersion still not an issue!


# 2) Check for spatial correlation
SRData$ResidSlv <- residuals(SRMEnv, type = "pearson")

# Convert data frame to an sf object with geographic coordinates
SRData_sf <- st_as_sf(SRData, coords = c("Longitude", "Latitude"), crs = 4326)  # Use EPSG:4326 (WGS 84)

# Convert to projected coordinates (e.g., UTM zone 32N) to ensure distance calculations are correct
SRData_sf <- st_transform(SRData_sf, crs = 32632)  # Change to your appropriate UTM zone

SRData_df <- as.data.frame(SRData_sf)

# Extract coordinates as a separate dataframe
SRData_df$X <- st_coordinates(SRData_sf)[,1]  # Easting
SRData_df$Y <- st_coordinates(SRData_sf)[,2]  # Northing

# Create a gstat variogram
SRData_df <- SRData_df %>% select(-OutlierLW, -RunningID)
vario <- variogram(ResidSlv ~ 1, locations = ~X + Y, data = SRData_df)

# Plot the variogram
plot(vario)
# --> does not suggest any spatial autocorrelation


# 3) Plot residuals vs fitted values
SRData$Resid <- resid(SRMEnv, type = "pearson")
SRData$Fitted <- fitted(SRMEnv)

SRData %>% ggplot(aes(x = Fitted, y = Resid)) +
  geom_point(col = "grey", alpha = .8) +
  geom_smooth() +
  theme_bw()    # looks ok

# Residuals vs covariates
SRData %>% ggplot(aes(x = DistanceGibraltar, y = Resid)) +
  geom_point(col = "grey", alpha = .8) +
  geom_smooth() +
  theme_bw()    # no clear pattern of heterogeneity
SRData %>% ggplot(aes(x = logTP, y = Resid)) +
  geom_point(col = "grey", alpha = .8) +
  geom_smooth() +
  theme_bw()    # no clear pattern of heterogeneity
SRData %>% ggplot(aes(x = logSurfaceArea, y = Resid)) +
  geom_point(col = "grey", alpha = .8) +
  geom_smooth() +
  theme_bw()    # looks good


# SENSITIVITY ANALYSIS 1: IS THE SELECTED MODEL ROBUST TO SAMPLING VARIABLES (I.E., GEAR OR MESH SIZE)
SRMEnvGear <- glmmTMB(bMale ~ DistanceGibraltar + logTP + logSurfaceArea + Gear + Meshsize_mm + (1|fSite), family = binomial, data = SRData)
summary(SRMEnvGear)
# --> both gear and mesh size are not significant and the significances of the three environmental variables in the selected models are robust (surface area even significant in this model, so could be carefully discussed)!


# SENSITIVITY ANALYSIS 2: TEST ROBUSTNESS OF DISTANCE-FROM-GIBRALTAR EFFECT TO GUADIARO AS AN OUTLIER SITE
SRData2 <- SRData %>% filter(fSite != "Guadiaro")
SRMEnvGuad <- glmmTMB(bMale ~ DistanceGibraltar + logTP + logSurfaceArea + (1|fSite), family = binomial, data = SRData2)
drop1(SRMEnvGuad, test = "Chi")
# --> Distance-from-Gibraltar effect is robust to inclusion or exclusion of Guadiaro


save(SRMEnv, file = "Output/Sex_ratio_env_model.rdata")
