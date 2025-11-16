# Visualize the raw data distributions of environmental variables with violin plots
PlotVars <- c("Latitude", "Longitude", "DistanceGibraltar", "SurfaceArea", 
              "Salinity", "Temperature", "Chlorophyll", "TP", "SurfaceArea")

HabitatDataPlot <- HabitatDataLong %>% filter(variable %in% PlotVars)
HabitatDataPlot$Habitat <- HabitatData$Habitat[match(HabitatDataPlot$Site, HabitatData$Site)]

HabitatDataPlot$value[which(HabitatDataPlot$variable == "DistanceGibraltar")] <- HabitatDataPlot$value[which(HabitatDataPlot$variable == "DistanceGibraltar")] / 1000

HabitatDataPlot <- HabitatDataPlot %>% mutate(variable = fct_recode(variable, "Salinity (g/l)" = "Salinity",
                                                                    "Surface area (ha)" = "SurfaceArea",
                                                                    "Temperature (°C)" = "Temperature",
                                                                    "Distance Gibraltar (km)" = "DistanceGibraltar",
                                                                    "Chlorophyll (  g/l)" = "Chlorophyll",
                                                                    "Total phosphorus (  g/l)" = "TP")) 

HabDistPlot <- HabitatDataPlot %>% ggplot(aes(x = 0, y = value)) +
  geom_jitter(aes(color = Habitat), width = .25) +
  geom_boxplot(fill = "transparent", width = .1, outliers = FALSE) +
  geom_violin(fill = "transparent") +
  scale_y_continuous(labels = label_comma(), name = element_blank()) +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "bottom",
        strip.text = element_text(size = 9)) +
  facet_wrap(~variable, scales = "free")
HabDistPlot
#ggsave(HabDistPlot, file = "Output/Habitat_data_raw_distributions.png", dpi = 600, height = 6, width = 6, units = "in")
