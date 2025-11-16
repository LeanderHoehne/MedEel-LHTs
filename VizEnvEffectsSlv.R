# Visualize environmental effects on silvering

load("Output/Slv_env_model.rdata")
LengthSeqF <- seq(45, 75)

# 1) Effect of Distance-to-Gibraltar
LHTStagedHeads <- LHTStaged %>% group_by(fSite) %>% slice_head(n = 1)   # create a df with one observation per any site included in the silvering analysis (to give an equal weight to all sites when drawing quantiles for env. variables)
DistLevels <- quantile(LHTStagedHeads$DistanceGibraltar, probs = c(.05, .25, .5, .75, .95))
# DistLevels <- (seq(1000000, 3000000, by = 500000) - mean(HabitatDataImp$DistanceGibraltar)) / sd(HabitatDataImp$DistanceGibraltar)
# Create predictions dataframe
NewDataDist <- expand.grid(
  Length_cm = LengthSeqF,
  DistanceGibraltar = DistLevels,  # Different levels for separate lines
  fSex = "F",  # Fixing sex
  Salinity = mean(LHTStagedHeads$Salinity, na.rm = TRUE),
  logSurfaceArea = mean(LHTStagedHeads$logSurfaceArea, na.rm = TRUE),
  fSite = NA  # Random effect terms need to be set to NA for population-level predictions
)

PredsSlvDistanceF <- predict(SlvEnv, newdata = NewDataDist, type = "response", se.fit = TRUE)

# Add predicted values and confidence intervals to the new_data dataframe
NewDataDist$Predicted <- PredsSlvDistanceF$fit

# Back-transform z-scores to original scale for plotting
LabelsDistance <- as.character(round(((DistLevels * sd(HabitatDataImp$DistanceGibraltar)) + mean(HabitatDataImp$DistanceGibraltar)) / 1000))

# Prepare color palette
Blues <- brewer.pal("Blues", n = 9)[c(3,9)]
PalDist <- colorRampPalette(c(Blues[1], Blues[2]))(length(DistLevels))

# Draw the plot
SlvEnvPDist <- NewDataDist %>% ggplot(aes(x = Length_cm, y = Predicted*100, color = as.factor(DistanceGibraltar))) +
  geom_line(size = 1) +
  labs(x = "Length (cm)") +
  scale_x_continuous(name = "", breaks = seq(45, 75, by = 5)) +
  scale_y_continuous(name = "Proportion of individuals silver (%)") +
  scale_color_manual(name = "Distance from\n Gibraltar (km)", values = PalDist, labels = LabelsDistance) +
  theme_bw() +
  theme(legend.position = c(0.82, 0.2),
        legend.text = element_text(size = 9),
        legend.title = element_text(size = 10),
        axis.title = element_text(size = 12),
        axis.text.x = element_text(size = 11),
        axis.text.y = element_text(size = 11))
#ggsave(SlvEnvPDist, file = "./Output/Environmental effects/Slv_vs_DistGib.png", dpi = 600, height = 5, width = 6, units = "in")


# 2) Effect of Salinity 
SalLevels <- quantile(LHTStagedHeads$Salinity, probs = c(.05, .25, .5, .75, .95))

# Create predictions dataframe
NewDataSal <- expand.grid(
  Length_cm = LengthSeqF,
  DistanceGibraltar = mean(LHTStagedHeads$DistanceGibraltar, na.rm = TRUE),  
  fSex = "F",  
  Salinity = SalLevels,
  logSurfaceArea = mean(LHTStagedHeads$logSurfaceArea, na.rm = TRUE),
  fSite = NA
)

PredsSlvSalF <- predict(SlvEnv, newdata = NewDataSal, type = "response", se.fit = TRUE)

# Back-transform z-scores to original scale for plotting
LabelsSal <- as.character(round(SalLevels * sd(HabitatDataImp$Salinity) + mean(HabitatDataImp$Salinity), 1))

# Add predicted values and confidence intervals to the new_data dataframe
NewDataSal$Predicted <- PredsSlvSalF$fit

# Draw the plot
SlvEnvPSal <- NewDataSal %>% ggplot(aes(x = Length_cm, y = Predicted, color = as.factor(Salinity))) +
  geom_line(size = 1) +
  labs(x = "Length (cm)") +
  scale_x_continuous(breaks = seq(45, 75, by = 5)) +
  scale_y_continuous(labels = label_percent()) +
  #scale_color_manual(name = "salinity (PSU)", values = PalDist, labels = LabelsDistance) +
  scale_color_viridis_d(labels = LabelsSal, name = "Salinity (PSU)", direction = -1) +
  theme_bw() +
  theme(legend.position = c(0.82, 0.2),
        legend.text = element_text(size = 9),
        legend.title = element_text(size = 10),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 11),
        axis.text.y = element_blank())
#ggsave(SlvEnvPSal, file = "./Output/Environmental effects/Slv_vs_Salinity.png", dpi = 600, height = 5, width = 6, units = "in")


# 3) Effect of SurfaceArea 
SurfLevels <- quantile(LHTStagedHeads$logSurfaceArea, probs = c(.05, .25, .5, .75, .95))

# Create predictions dataframe
NewDataSurf <- expand.grid(
  Length_cm = LengthSeqF,
  DistanceGibraltar = mean(LHTStagedHeads$DistanceGibraltar, na.rm = TRUE),  
  fSex = "F",  
  Salinity = mean(LHTStagedHeads$Salinity, na.rm = TRUE),
  logSurfaceArea = SurfLevels,
  fSite = NA
)

PredsSlvSurfF <- predict(SlvEnv, newdata = NewDataSurf, type = "response", se.fit = TRUE)

# Back-transform z-scores to original scale for plotting
LabelsSurf <- as.character(SurfLevels * sd(HabitatDataImp$logSurfaceArea) + mean(HabitatDataImp$logSurfaceArea))

# Add predicted values and confidence intervals to the new_data dataframe
NewDataSurf$Predicted <- PredsSlvSurfF$fit

# Prepare color palette
Purples <- brewer.pal("PuRd", n = 9)[c(3,9)]
PalSurf <- colorRampPalette(c(Purples[1], Purples[2]))(length(SurfLevels))

# Draw the plot
SlvEnvPSurf <- NewDataSurf %>% ggplot(aes(x = Length_cm, y = Predicted, color = as.factor(logSurfaceArea))) +
  geom_line(size = 1) +
  labs(x = "Length (cm)") +
  scale_x_continuous(name = "", breaks = seq(45, 75, by = 5)) +
  scale_y_continuous(labels = label_percent()) +
  scale_color_manual(name = "Surface area (ha)", values = PalSurf, labels = round(exp(as.numeric(LabelsSurf)))) +
  theme_bw() +
  theme(legend.position = c(0.82, 0.2),
        legend.text = element_text(size = 9),
        legend.title = element_text(size = 10),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 11),
        axis.text.y = element_blank())
#ggsave(SlvEnvPSurf, file = "./Output/Environmental effects/Slv_vs_Surface.png", dpi = 600, height = 5, width = 6, units = "in")


SlvEnvPComb <- cowplot::plot_grid(SlvEnvPDist, SlvEnvPSal, SlvEnvPSurf, rel_widths = c(1.13, 1, 1), nrow = 1)
#ggsave(SlvEnvPComb, file = "./Output/Environmental effects/Slv_all_effects.png", dpi = 600, height = 5, width = 12, units = "in")
