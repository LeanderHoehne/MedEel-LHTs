# Visualize the effects of environmental variables on growth
load("./Output/GMEnvObs.rdata")

# Test significance
if (Rerun) source('LRTests_Growth_env.R')

summary(GMEnvObs)

# Visualize the marginal effect of distance-from-Gibraltar on growth

# Get 5, 50, and 95% quantile values for distance-from-Gibralter in the growth dataset
DistHeads <- LengthAgeObs %>% group_by(fSite) %>% reframe(Dist = first(DistanceGibraltar))   # create a df with one observation per any site included in the growth analysis (to give an equal weight to all sites when drawing quantiles for env. variables)
LevelsDistance <- quantile(DistHeads$Dist, probs = c(.05, .5, .95))
# Create a new dataframe for predictions at different distance-from-Gibraltar levels
PredsGMDistance <- expand.grid(
  Age = seq(min(LengthAgeObs$Age), 20, length.out = 100),
  DistanceGibraltar = LevelsDistance,  
  fSex = levels(LengthAgeObs$fSex), 
  StageSelectivity = as.factor(0)
)
PredsGMDistance <- PredsGMDistance[-which(PredsGMDistance$fSex == "M" & PredsGMDistance$Age > 11),]
# Calculated predicted values
PredsGMDistance$Predicted <- predict(GMEnvObs, newdata = PredsGMDistance, level = 0)
PredsGMDistance$CaseID <- paste(PredsGMDistance$DistanceGibraltar, PredsGMDistance$fSex)

# Back-transform z-scores to original scale for plotting
LabelsDistance <- ((LevelsDistance * sd(HabitatDataImp$DistanceGibraltar)) + mean(HabitatDataImp$DistanceGibraltar)) / 1000

# Create a color palette
PaletteEnv <- brewer.pal("Blues", n = 9)[c(3,6,9)]

GREnvPDist <- ggplot(PredsGMDistance, aes(x = Age, y = Predicted, color = as.factor(DistanceGibraltar), group = CaseID)) +
  geom_line() +
  scale_x_continuous(name = "Continental age (years)", seq(2, 20, by = 2)) +
  scale_y_continuous(name = "Predicted length (cm)") +
  scale_color_manual(name = "Distance from \nGibraltar (km)", values = PaletteEnv, labels = as.character(round(LabelsDistance))) +
  theme_bw() +
  theme(legend.title = element_text(size = 10),
        legend.text = element_text(size = 9),
        axis.title = element_text(size = 13),
        axis.text = element_text(size = 11),
        legend.position = c(0.885, .18),
        legend.background = element_rect(fill = "white", color = "grey40", linewidth = .3))
#ggsave(GREnvPDist, file = "Output/Environmental effects/Growth_vs_Distance.png", dpi = 600, width = 6, height = 4.5, units = "in")