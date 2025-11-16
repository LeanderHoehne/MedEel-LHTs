# Fit a mixed-effects model of the silvering process including environmental covariates

# Explore distributions of explanatory variables
dotchart(LHTStaged$Temperature)   # two slight outliers (one at the lower, one at the higher end), but still acceptable
dotchart(LHTStaged$DistanceGibraltar)   # one outlier at the lower end (Guadiaro)
dotchart(LHTStaged$Salinity)
dotchart(LHTStaged$Chlorophyll)
dotchart(LHTStaged$TP)    
dotchart(LHTStaged$logTP)   # similar; so we will compare AICs to judge whether TP should be log-transformed
dotchart(LHTStaged$SurfaceArea)
dotchart(LHTStaged$logSurfaceArea)    # better

# --> We will first run the analysis on the full dataset and then run sensitivity analysis without Guadiaro to check whether it is of influential nature

# Test whether TP should be log-transformed
SlvEnvTP1 <- glmmTMB(bSilver ~ Length_cm + fSex + DistanceGibraltar + Salinity + Temperature + Chlorophyll + TP + logSurfaceArea + (1 + Length_cm|fSite/fSex), 
                  family = binomial, data = LHTStaged, control = glmmTMBControl(optCtrl = list(iter.max = 10000, eval.max = 10000, rel.tol = 1e-4))) 
SlvEnvTP2 <- glmmTMB(bSilver ~ Length_cm + fSex + DistanceGibraltar + Salinity + Temperature + Chlorophyll + logTP + logSurfaceArea + (1 + Length_cm|fSite/fSex), 
                  family = binomial, data = LHTStaged, control = glmmTMBControl(optCtrl = list(iter.max = 10000, eval.max = 10000, rel.tol = 1e-4))) 
AIC(SlvEnvTP1, SlvEnvTP2)
# --> we will use log-transformed TP values!

# Explore the global model
SlvEnv <- SlvEnvTP2
summary(SlvEnv)


# BACKWARD VARIABLE SELECTION

# Step 1
drop1(SlvEnv)   # --> Temperature should be dropped
SlvEnv1 <- update(SlvEnv, . ~ . -Temperature)

# Step 2
drop1(SlvEnv1)    # --> logTP should be dropped
SlvEnv2 <- update(SlvEnv1, . ~ . -logTP)

# Step 3
drop1(SlvEnv2)    # --> Chlorophyll should be dropped
SlvEnv3 <- update(SlvEnv2, . ~ . -Chlorophyll)

# Step 4
drop1(SlvEnv3)
# --> further elimination of variables does not improve the model fit. SlvEnv3 is the best model!

summary(SlvEnv3)

SlvEnv <- SlvEnv3

save(SlvEnv, file = "Output/SlvEnvModel.rdata")

# MODEL DIAGNOSTICS

# 1) Check for overdispersion
E1 <- resid(SlvEnv, type = "pearson")

# Get the sample size and the number of parameters.
N <- nrow(LHTStaged)
p <- length(coef(SlvEnv)) + 1 # The '+1' is for the sigma

# Determine the dispersion statistic
DispersionStatistic <- sum(E1^2) / (N - p)
DispersionStatistic
# --> Overdispersion seems not to be a big issue: Dispersion statistic is reasonably close to 1, after inclusion of environmental covariates


# 2) Check for spatial correlation
LHTStaged$ResidSlv <- residuals(SlvEnv, type = "pearson")

# Convert data frame to an sf object with geographic coordinates
LHTStaged_sf <- st_as_sf(LHTStaged, coords = c("Longitude", "Latitude"), crs = 4326)  # Use EPSG:4326 (WGS 84)

# Convert to projected coordinates (e.g., UTM zone 32N) to ensure distance calculations are correct
LHTStaged_sf <- st_transform(LHTStaged_sf, crs = 32632)  # Change to your appropriate UTM zone

LHTStaged_df <- as.data.frame(LHTStaged_sf)

# Extract coordinates as a separate dataframe
LHTStaged_df$X <- st_coordinates(LHTStaged_sf)[,1]  # Easting
LHTStaged_df$Y <- st_coordinates(LHTStaged_sf)[,2]  # Northing

# Create a gstat variogram
LHTStaged_df <- LHTStaged_df %>% select(-OutlierLW, -RunningID)
vario <- variogram(ResidSlv ~ 1, locations = ~X + Y, data = LHTStaged_df)

# Plot the variogram
plot(vario)
# --> does not suggest any spatial autocorrelation


# 3) Plot residuals vs fitted values
LHTStaged$Resid <- resid(SlvEnv, type = "pearson")
LHTStaged$Fitted <- fitted(SlvEnv)

LHTStaged %>% ggplot(aes(x = Fitted, y = Resid)) +
  geom_point(col = "grey", alpha = .8) +
  geom_smooth() +
  theme_bw()    # looks ok

# Residuals vs covariates
LHTStaged %>% ggplot(aes(x = Length_cm, y = Resid)) +
  geom_point(col = "grey", alpha = .8) +
  geom_smooth() +
  theme_bw()    # no clear pattern of heterogeneity
LHTStaged %>% ggplot(aes(x = fSex, y = Resid)) +
  geom_boxplot() +
  theme_bw()    # slight heterogeneity between sexes, but as we are mostly interested in the p-values of environmental variables we will accept that
LHTStaged %>% ggplot(aes(x = DistanceGibraltar, y = Resid)) +
  geom_point(col = "grey", alpha = .8) +
  geom_smooth() +
  theme_bw()    # no clear pattern of heterogeneity
LHTStaged %>% ggplot(aes(x = Salinity, y = Resid)) +
  geom_point(col = "grey", alpha = .8) +
  geom_smooth() +
  theme_bw()    # no clear pattern of heterogeneity
LHTStaged %>% ggplot(aes(x = logSurfaceArea, y = Resid)) +
  geom_point(col = "grey", alpha = .8) +
  geom_smooth() +
  theme_bw()    # no clear pattern of heterogeneity


# SENSITIVITY ANALYSIS: TEST ROBUSTNESS OF DISTANCE-FROM-GIBRALTAR EFFECT TO GUADIARO AS AN OUTLIER SITE

LHTStagedSens <- LHTStaged %>% filter(fSite != "Guadiaro")

SlvEnvSens <- update(SlvEnv, data = LHTStagedSens)

drop1(SlvEnvSens, test = "Chi", control = glmmTMBControl(optCtrl = list(iter.max = 100000, eval.max = 100000)))   
# --> The effect of DistanceGibraltar (and Salinity and Surface area) are robust to Guadiaro as a univariate outlier. 


# Obtain estimates, SEs, and robust p-values of the retained model (SlvEnv)
summary(SlvEnv)
drop1(SlvEnv, test = "Chi")

# Store the selected model
save(SlvEnv, file = "Output/Slv_env_model.rdata")
save(LHTStaged, file = "Slv_env_data.rdata")
