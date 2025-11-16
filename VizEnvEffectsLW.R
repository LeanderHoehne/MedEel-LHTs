load("./Output/LW_env_model.Rdata")
load("./Output/LW_env_data.Rdata")

# Plot salinity effect
LevelsSalinityLW <- quantile(unique(LWEnvData$Salinity), probs = c(.05,.5,.95))
NewDataSal <- expand.grid(Salinity = LevelsSalinityLW, 
                          logLength_cm = seq(3.4, 4.6, by = .05),
                          fStage = "Y", Chlorophyll = mean(unique(LWEnvData$Chlorophyll)), logSurfaceArea = mean(unique(LWEnvData$logSurfaceArea)),
                          DistanceGibraltar = mean(unique(LWEnvData$DistanceGibraltar)), fSite = NA)
PredsLWSal <- data.frame(NewDataSal, Predicted = predict(LWEnv, newdata = NewDataSal)) %>% 
  mutate(SalLevel = as_factor(round((Salinity * sd(HabitatDataImp$Salinity)) + mean(HabitatDataImp$Salinity), 2)), 
         PredictedOrig = exp(Predicted),
         Length_cm = exp(logLength_cm))

# Create a color palette
PaletteSal <- brewer.pal("YlGnBu", n = 9)[c(3,6,9)]

EnvPlotLWSal <- PredsLWSal %>% ggplot(aes(x = Length_cm, y = PredictedOrig, color = SalLevel)) +
  geom_line() +
  guides(color = guide_legend(title = "Salinity (PSU)")) +
  scale_x_continuous(breaks = seq(30, 100, by = 10), name = "Length (cm)") +
  scale_y_continuous(name = "Weight (g)", labels = label_comma()) +
  scale_color_manual(values = PaletteSal) +
  theme_bw() +
  theme(legend.position = "bottom",
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14))
# ggsave(EnvPlotLWSal, file = "./Output/Environmental effects/LW_vs_Salinity.png", dpi = 600, height = 5.5, width = 6, units = "in")


# Plot DistanceGibraltar effect
LevelsGibraltarLW <- quantile(unique(LWEnvData$DistanceGibraltar), probs = c(.05,.5,.95))
NewDataGib <- expand.grid(Salinity = mean(unique(LWEnvData$Salinity)), 
                          logLength_cm = seq(3.4, 4.6, by = .05),
                          fStage = "Y", 
                          Chlorophyll = mean(unique(LWEnvData$Chlorophyll)), 
                          logSurfaceArea = mean(unique(LWEnvData$logSurfaceArea)),
                          DistanceGibraltar = LevelsGibraltarLW, 
                          fSite = NA)
PredsLWGib <- data.frame(NewDataGib, Predicted = predict(LWEnv, newdata = NewDataGib)) %>% 
  mutate(GibLevel = as_factor(round(((DistanceGibraltar * sd(HabitatDataImp$DistanceGibraltar)) + mean(HabitatDataImp$DistanceGibraltar)) / 1000)), 
         PredictedOrig = exp(Predicted),
         Length_cm = exp(logLength_cm))

# Create a color palette
PaletteGib <- brewer.pal("Blues", n = 9)[c(3,6,9)]

EnvPlotLWGib <- PredsLWGib %>% ggplot(aes(x = Length_cm, y = PredictedOrig, color = GibLevel)) +
  geom_line() +
  guides(color = guide_legend(title = "Distance from Gibraltar (km)")) +
  scale_x_continuous(breaks = seq(30, 100, by = 10), name = "Length (cm)") +
  scale_y_continuous(name = "Weight (g)", labels = label_comma()) +
  scale_color_manual(values = PaletteGib) +
  theme_bw() +
  theme(legend.position = "bottom",
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14))
ggsave(EnvPlotLWGib, file = "./Output/Environmental effects/LW_vs_DistGibraltar.png", dpi = 600, height = 5.5, width = 6, units = "in")


# Plot DistanceGibraltar effect
LevelsSurfaceLW <- quantile(unique(LWEnvData$logSurfaceArea), probs = c(.05,.5,.95))
NewDataSurf <- expand.grid(Salinity = mean(unique(LWEnvData$Salinity)), 
                          logLength_cm = seq(3.4, 4.6, by = .05),
                          fStage = "Y", 
                          Chlorophyll = mean(unique(LWEnvData$Chlorophyll)), 
                          logSurfaceArea = LevelsSurfaceLW,
                          DistanceGibraltar = mean(unique(LWEnvData$DistanceGibraltar)), 
                          fSite = NA)
PredsLWSurf <- data.frame(NewDataSurf, Predicted = predict(LWEnv, newdata = NewDataSurf)) %>% 
  mutate(SurfLevel = as_factor(round(exp((logSurfaceArea * sd(HabitatDataImp$logSurfaceArea)) + mean(HabitatDataImp$logSurfaceArea)))), 
         PredictedOrig = exp(Predicted),
         Length_cm = exp(logLength_cm))

# Create a color palette
PaletteSurf <- brewer.pal("PuRd", n = 9)[c(3,6,9)]

EnvPlotLWSurf <- PredsLWSurf %>% ggplot(aes(x = Length_cm, y = PredictedOrig, color = SurfLevel)) +
  geom_line() +
  guides(color = guide_legend(title = "Surface area (ha)")) +
  scale_x_continuous(breaks = seq(30, 100, by = 10), name = "Length (cm)") +
  scale_y_continuous(name = "Weight (g)", labels = label_comma()) +
  scale_color_manual(values = PaletteSurf) +
  theme_bw() +
  theme(legend.position = "bottom",
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14))
ggsave(EnvPlotLWSurf, file = "./Output/Environmental effects/LW_vs_Surface.png", dpi = 600, height = 5.5, width = 6, units = "in")
