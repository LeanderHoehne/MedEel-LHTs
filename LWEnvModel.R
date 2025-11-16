LWEnvData <- LHT_HAB2 %>% filter(!is.na(fStage)) %>% mutate(fStage = droplevels(fStage))   # filter out eels with lacking stage info, because we want to account for it in the model

# Explore distribution of environmental data to see where log-transformation might be appropriate
dotchart(LWEnvData$Temperature)   # one slight lower outlier (Sigoulette)
dotchart(LWEnvData$DistanceGibraltar)   # one outlier at the lower end (Guadiaro)
dotchart(LWEnvData$Salinity)   # ok
dotchart(LWEnvData$Chlorophyll)   # slightly skewed distribution but continuous gradient
dotchart(LWEnvData$TP)    
dotchart(LWEnvData$logTP)   # better
dotchart(LWEnvData$SurfaceArea)
dotchart(LWEnvData$logSurfaceArea)    # better

# Delete cases with too few observations (<15)
LWEnvData %>% group_by(fSite) %>% tally() 
LWEnvData <- LWEnvData %>% filter(fSite != "Bardawil")
# minimum sample size in the dataset is 22 individuals

# Fit the global model including all environmental covariates
LWEnv <- glmmTMB(logWeight_g ~ logLength_cm + fStage + Chlorophyll + logTP + DistanceGibraltar + Temperature + Salinity + logSurfaceArea + (1 + logLength_cm|fSite), data = LWEnvData, na.action = na.fail)
summary(LWEnv)

# Step 1
LWEnvDrop1 <- drop1(LWEnv, direction = "backward")
rownames(LWEnvDrop1)[which.min(LWEnvDrop1[,2])]   # Temperature should be dropped
LWEnv2 <- glmmTMB(logWeight_g ~ logLength_cm + fStage + Chlorophyll + logTP + DistanceGibraltar + Salinity + logSurfaceArea + (1 + logLength_cm|fSite), data = LWEnvData, na.action = na.fail)

# Step 3
LWEnvDrop2 <- drop1(LWEnv2, direction = "backward")
rownames(LWEnvDrop2)[which.min(LWEnvDrop2[,2])]   # logTP should be dropped
LWEnv3 <- glmmTMB(logWeight_g ~ logLength_cm + fStage + Chlorophyll + DistanceGibraltar + Salinity + logSurfaceArea + (1 + logLength_cm|fSite), data = LWEnvData, na.action = na.fail)

# Step 4
LWEnvDrop3 <- drop1(LWEnv3, direction = "backward")
rownames(LWEnvDrop3)[which.min(LWEnvDrop3[,2])]   # --> LWEnv3 is the best model!


# Refit the selected model with REML to get parameter estimates
LWEnv <- update(LWEnv3, REML = TRUE)
summary(LWEnv)

# Calculate p-values based on the ML-fitted model and Likelihood Ratio Tests
drop1(LWEnv3, test = "Chi")


# Model diagnostics

# does the Gaussian model predict negative values?
LWEnvData$Fitted <- exp(fitted(LWEnv))
range(LWEnvData$Fitted)
# --> no negative fitted values. We can keep the Gaussian model.


# Plot residuals vs fitted values
LWEnvData$Resid <- resid(LWEnv)
LWEnvData$Fitted <- fitted(LWEnv)

LWEnvData %>% ggplot(aes(x = Fitted, y = Resid)) +
  geom_point(col = "grey", alpha = .8) +
  geom_smooth() +
  theme_bw()    # looks ok

# Residuals vs covariates
LWEnvData %>% ggplot(aes(x = logLength_cm, y = Resid)) +
  geom_point(col = "grey", alpha = .8) +
  geom_smooth() +
  theme_bw()    # no clear pattern of heterogeneity
LWEnvData %>% ggplot(aes(x = fStage, y = Resid)) +
  geom_boxplot() +
  theme_bw()    # quite ok
LWEnvData %>% ggplot(aes(x = DistanceGibraltar, y = Resid)) +
  geom_point(col = "grey", alpha = .8) +
  geom_smooth() +
  theme_bw()    # no clear pattern of heterogeneity
LWEnvData %>% ggplot(aes(x = Chlorophyll, y = Resid)) +
  geom_point(col = "grey", alpha = .8) +
  geom_smooth() +
  theme_bw()    # no clear pattern of heterogeneity
LWEnvData %>% ggplot(aes(x = logSurfaceArea, y = Resid)) +
  geom_point(col = "grey", alpha = .8) +
  geom_smooth() +
  theme_bw()    # slight "trumpet effect"
LWEnvData %>% ggplot(aes(x = Salinity, y = Resid)) +
  geom_point(col = "grey", alpha = .8) +
  geom_smooth() +
  theme_bw()    # no clear pattern of heterogeneity


# Refit with nlme and add variance weights on surface area to see if it makes a difference
LWEnv_nlme1 <- lme(
  fixed = logWeight_g ~ logLength_cm + fStage + Chlorophyll + DistanceGibraltar + logSurfaceArea,
  random = ~ logLength_cm | fSite,
  data = LWEnvData,
  na.action = na.fail
)
LWEnv_nlme2 <- lme(
  fixed = logWeight_g ~ logLength_cm + fStage + Chlorophyll + DistanceGibraltar + logSurfaceArea,
  random = ~ logLength_cm | fSite,
  weights = varFixed(~ logSurfaceArea),
  data = LWEnvData,
  na.action = na.fail
)

AIC(LWEnv_nlme1, LWEnv_nlme2)   # better without variance weights, so we keep the original model


# Check for spatial autocorrelation

LWEnvData$ResidLW <- residuals(LWEnv, type = "pearson")

# Convert data frame to an sf object with geographic coordinates
LWEnvData_sf <- st_as_sf(LWEnvData, coords = c("Longitude", "Latitude"), crs = 4326)  # Use EPSG:4326 (WGS 84)

# Convert to projected coordinates (e.g., UTM zone 32N) to ensure distance calculations are correct
LWEnvData_sf <- st_transform(LWEnvData_sf, crs = 32632)  # Change to the appropriate UTM zone

# Extract coordinates as a separate dataframe
LWEnvData_sf$X <- st_coordinates(LWEnvData_sf)[,1]  # Easting
LWEnvData_sf$Y <- st_coordinates(LWEnvData_sf)[,2]  # Northing

LWEnvData_sf <- LWEnvData_sf %>% dplyr::select(-OutlierLW, -RunningID)
LWEnvData_sf <- as.data.frame(LWEnvData_sf)

# Create a gstat variogram
vario <- variogram(ResidLW ~ 1, locations = ~X + Y, data = LWEnvData_sf)

# Plot the variogram
plot(vario)
# --> does not suggest any serious spatial autocorrelation

# Store the model output for visualization
save(LWEnv, file = "./Output/LW_env_model.Rdata")
save(LWEnvData, file = "LW_env_data.Rdata")
