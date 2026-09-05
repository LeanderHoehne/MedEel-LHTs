load("Output/Sex_ratio_env_model.rdata")

# Calculate female fraction and add coordinates to results table
SRSummary <- SRSummary %>% mutate(Female = 1 - MalesFr,
                                  Latitude = HabitatData$Latitude[match(fSite, HabitatData$Site)],
                                  Longitude = HabitatData$Longitude[match(fSite, HabitatData$Site)]) %>% 
  rename(Male = MalesFr)

# Recode some coordinates to avoid overlapping pie charts
SRMResults <- SRSummary %>% mutate(plotLatitude = case_when(fSite == "Venezia" ~ 45.808024,
                                                            fSite == "Lesina" ~ 41.48607,
                                                            fSite == "Rasa" ~ 45.487066,
                                                            fSite == "Ghar El Melh" ~ 36.06593,
                                                            fSite == "Mellah" ~ 35.84902041471095,
                                                            fSite == "Oubeira" ~ 35.32131288068617,
                                                            fSite == "Kotychi" ~ 37.8268,
                                                            fSite == "Messolonghi-Aitoliko" ~ 39.04148546941868,
                                                            fSite == "Cetina" ~ 44.33561608894303,
                                                            fSite == "Bracciano" ~ 42.46570327337685,
                                                            fSite == "Fumemorte" ~ 44.3201654455297,
                                                            fSite == "Salses-Leucate" ~ 43.10584406940861,
                                                            fSite == "Ayrolle" ~ 42.32494, 
                                                            fSite == "Bages-Sigean" ~ 43.688442641588814,
                                                            fSite == "Canet" ~ 42.393469415435455,
                                                            fSite == "Complex_Vendres" ~ 44.091814207059514,
                                                            fSite == "Complex_Gruissan" ~ 42.74076,
                                                            fSite == "Complex_palavasien" ~ 44.73068,
                                                            fSite == "Complex_Petite_Camargue" ~ 44.35054,
                                                            fSite == "Fogliano" ~ 41.35038,
                                                            fSite == "Jadro" ~ 43.81544,
                                                            fSite == "Mafragh" ~ 36.49337,
                                                            fSite == "Mirna" ~ 45.87351,
                                                            fSite == "Or" ~ 44.668535873049045, 
                                                            fSite == "Sile" ~ 46.28764,
                                                            fSite == "Tevere" ~ 41.58863,
                                                            fSite == "Thau" ~ 44.39346,
                                                            fSite == "Tonga" ~ 35.16558,
                                                            fSite == "Trasimeno" ~ 43.418243790036364,
                                                            fSite == "Tunis North" ~ 35.86365,
                                                            fSite == "Varano" ~ 41.28134,
                                                            TRUE ~ as.numeric(Latitude)))
SRMResults <- SRMResults %>% mutate(plotLongitude = case_when(fSite == "Venezia" ~ 11.442372,
                                                              fSite == "Lesina" ~ 14.92442,
                                                              fSite == "Rasa" ~ 14.930952,
                                                              fSite == "Ghar El Melh" ~ 9.43104,
                                                              fSite == "Mellah" ~ 6.813005655457604,
                                                              fSite == "Oubeira" ~ 7.669278654461002,
                                                              fSite == "Kotychi" ~ 21.57407,
                                                              fSite == "Messolonghi-Aitoliko" ~ 21.695082072846407,
                                                              fSite == "Cetina" ~ 17.53475697526952,
                                                              fSite == "Bracciano" ~ 13.178310399192315,
                                                              fSite == "Fumemorte" ~ 5.218910188632545,
                                                              fSite == "Salses-Leucate" ~ 0.6095865386198144,
                                                              fSite == "Ayrolle" ~ 3.68982,
                                                              fSite == "Bages-Sigean" ~ 1.2506708270122253,
                                                              fSite == "Canet" ~ 1.9211749646957805,
                                                              fSite == "Complex_Vendres" ~ 1.912995882326056,
                                                              fSite == "Complex_Gruissan" ~ 4.53733,
                                                              fSite == "Complex_palavasien" ~ 3.51446,
                                                              fSite == "Complex_Petite_Camargue" ~ 5.18352,
                                                              fSite == "Fogliano" ~ 13.3132,
                                                              fSite == "Jadro" ~ 16.32679,
                                                              fSite == "Mafragh" ~ 6.21687,
                                                              fSite == "Mirna" ~ 14.15122,
                                                              fSite == "Or" ~ 4.456584010531679,
                                                              fSite == "Sile" ~ 12.78997,
                                                              fSite == "Tevere" ~ 12.39528,
                                                              fSite == "Thau" ~ 2.70306,
                                                              fSite == "Tonga" ~ 8.66751,
                                                              fSite == "Trasimeno" ~ 12.12547117934691,
                                                              fSite == "Tunis North" ~ 10.28872,
                                                              fSite == "Varano" ~ 15.69063,
                                                              TRUE ~ as.numeric(Longitude)))
SRMResults$Moved <- ifelse(SRMResults$Latitude == SRMResults$plotLatitude, 0, 1)

# Plot site-specific sex ratios on a map
SRMResultsLong <- SRMResults %>% pivot_longer(cols = c("Male", "Female"), 
                                              names_to = "fSex", 
                                              values_to = "Fraction") %>% 
  group_by(fSite) %>% mutate(CumFraction = cumsum(Fraction),
                             Start = lag(CumFraction, default = 0),
                             End = CumFraction)

SRMResultsLong$PieRadius <- .3

makePieSlices <- function(data, r_scale_x, r_scale_y) {
  pmap_dfr(data, function(plotLongitude, plotLatitude, Start, End, PieRadius, fSex, ...) {
    Angles <- seq(Start * 2 * pi, End * 2 * pi, length.out = 50) + pi/2
    tibble(
      x     = c(plotLongitude, plotLongitude + PieRadius * rScaleX * cos(Angles), plotLongitude),
      y     = c(plotLatitude,  plotLatitude  + PieRadius * rScaleY * sin(Angles), plotLatitude),
      fSex  = fSex,
      group = paste(plotLongitude, plotLatitude, fSex)
    )
  })
}

Width <- 7.8
Height <- 4
xRange <- 42.5
yRange <- 17

rScaleX  <- cos(45 * pi / 180) * 2
rScaleY  <- rScaleX * (Width / xRange) / (Height / yRange)

Pie0 <- makePieSlices(SRMResultsLong %>% filter(Moved == 0), r_scale_x, r_scale_y)
Pie1 <- makePieSlices(SRMResultsLong %>% filter(Moved == 1), r_scale_x, r_scale_y)

SRMap <- ggplot() +
  geom_sf(data = MedArea, fill = "grey98", color = "grey50", size = 0.5, alpha = .5) +
  geom_polygon(data = Pie0, aes(x = x, y = y, group = group, fill = fSex),
               color = "grey80", linewidth = .2) +
  geom_segment(data = SRMResults %>% filter(Moved == 1),
               aes(x = plotLongitude, xend = Longitude, y = plotLatitude, yend = Latitude),
               linewidth = .3) +
  geom_polygon(data = Pie1, aes(x = x, y = y, group = group, fill = fSex),
               color = "grey80", linewidth = .2) +
  scale_x_continuous(limits = c(-6, 36.5), name = element_blank()) +
  scale_y_continuous(limits = c(30, 47), name = element_blank()) +
  scale_fill_manual(values = Palette) +
  guides(fill = guide_legend(title = "Sex")) +
  coord_sf() +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.background = element_rect(fill = alpha("skyblue", .1), color = "black"),
        legend.position = c(0.93, 0.87),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 7),
        legend.key.size = unit(0.4, "cm"),
        legend.background = element_rect(color = "black", linewidth = .2),
        axis.text.y = element_text(angle = 90, hjust = .5))
# ggsave(SRMap, file = "./Output/Sitelevel/Sex_ratio_map.png", dpi = 600, height = Height, width = Width, units = "in")
