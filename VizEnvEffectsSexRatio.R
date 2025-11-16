# Extract predicted values for any environmental variable, adjusted for the pop.mean in other variables
PredSRDistance <- ggpredict(SRMEnv, terms = "DistanceGibraltar[all]")
PredSRSurface <- ggpredict(SRMEnv, terms = "logSurfaceArea[all]")
PredSRTP <- ggpredict(SRMEnv, terms = "logTP[all]")

# Back-transform z-scaled values of env. variables to plot them on the original scale
PredSRDistance$x <- (PredSRDistance$x * sd(HabitatDataImp$DistanceGibraltar) + mean(HabitatDataImp$DistanceGibraltar)) / 1000    # in km
PredSRSurface$x <- PredSRSurface$x * sd(HabitatDataImp$logSurfaceArea) + mean(HabitatDataImp$logSurfaceArea)
PredSRTP$x <- PredSRTP$x * sd(HabitatDataImp$logTP) + mean(HabitatDataImp$logTP)

SREnvP1 <- plot(PredSRDistance) +
  xlab("Distance from Gibraltar (km)") +
  ylab("Proportion of males (%)") +
  scale_x_continuous(labels = label_comma()) +
  scale_y_continuous(labels = label_number(scale = 100)) +
  ggtitle(element_blank()) +
  theme_bw() +
  theme(axis.title = element_text(size = 9),
        axis.text = element_text(size = 7))
#ggsave(SREnvP1, file = "./Output/Environmental effects/SR_vs_DistanceGibraltar.png", dpi = 400, height = 4.5, width = 6, units = "in")

SREnvP2 <- plot(PredSRSurface) +
  xlab("log(Surface area) in ha") +
  scale_y_continuous(labels = label_number(scale = 100)) +
  ggtitle(element_blank()) +
  theme_bw() +
  theme(axis.title.x = element_text(size = 9),
        axis.title.y = element_blank(),
        axis.text = element_text(size = 8))
#ggsave(SREnvP2, file = "./Output/Environmental effects/SR_vs_Surface.png", dpi = 400, height = 4.5, width = 6, units = "in")

SREnvP3 <- plot(PredSRTP) +
  xlab(expression("log(Total phosphorus) in " * mu * "g / l")) +
  scale_y_continuous(labels = label_number(scale = 100)) +
  ggtitle(element_blank()) +
  theme_bw() +
  theme(axis.title.x = element_text(size = 9),
        axis.title.y = element_blank(),
        axis.text = element_text(size = 8))
#ggsave(SREnvP3, file = "./Output/Environmental effects/SR_vs_TP.png", dpi = 400, height = 4.5, width = 6, units = "in")

SREnvPComb <- cowplot::plot_grid(SREnvP1, SREnvP2, SREnvP3, nrow = 1)
#ggsave(SREnvPComb, file = "./Output/Environmental effects/SR_all_effects.png", dpi = 600, height = 4, width = 9, units = "in")

save(SRData, file = "Sex_ratio_env_data.rdata")
