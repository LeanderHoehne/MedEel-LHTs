# Title: Analysis of life-history traits (growth, length-weight relation, sex ratio, silvering) across European eels from the Mediterranean area
# Description: This is the R-script underlying all analyses presented in Hoehne et al. "Distance from spawning grounds is the main driver of life-history variation in European eels across the Mediterranean"
# Author: Leander Höhne
# Last update: 21.10.2025


# 0) Data import and pre-procession ---------------------------------------

# Clear memory
remove(list=ls())
gc()
options(scipen = 999)   # suppress scientific notation

Rerun <- 1    # When this switch is activated (1), all statistical models in this script will be run. 
# It can be set to 0 to avoid re-running the (computationally intense) models or plots for which there is already an output in the "Output" folder

# Load libraries
library(tidyverse)
library(nlme)
library(RColorBrewer)
library(readxl)
library(glmmTMB)
library(car)
library(ggeffects)
library(scales)
library(sf)
library(rnaturalearth)
library(GGally)
library(ggh4x)
library(FSA)
library(missForest)
library(lme4)
library(performance)
library(ggforce)
library(gstat)
library(DescTools)
source('LibraryLHT.R')    # load customized functions for this project


# 0.1) Life-history data - Read in raw data
LHT <- read_excel("Input/Biometrics.xlsx", guess = 30000)
# Replace ?s and NDs with NAs
LHT[LHT == "ND"] <- NA
LHT[LHT == "NA"] <- NA
LHT[LHT == "?"] <- NA


# 0.2) Data procession

# 0.2.1) Data cleaning: Remove fish that cannot be used in analyses
# Remove fish that have too few information to be used in any analysis (minimum: length and weight information)
LHT <- LHT[!is.na(LHT$Length_mm) & !is.na(LHT$Weight_g),]

# Drop site-specific datasets with < 10 observations (additional sample size restrictions follow in the respective sections)
LHTObs <- LHT %>% group_by(SiteAcronym) %>% tally() %>% filter(n >= 10)   
LHT <- LHT %>% filter(SiteAcronym %in% LHTObs$SiteAcronym)


# 0.2.2) Amend data: Infer sex and stage from Durif indices (where possible), harmonize sampling gear types
LHT <- LHT %>% mutate(Sex = case_when(Stage == "MII" ~ "M", 
                         Stage %in% c("FII", "FIII", "FIV", "FV") ~ "F",
                         Sex == "Immature" ~ "IND",
                         TRUE ~ Sex),
                      Stage = case_when(Stage %in% c("I", "FII", "FIII", "YS", "PS") ~ "Y",
                           Stage %in% c("MII", "FIV", "FV") ~ "S",
                           Stage == "glass eel" ~ "G",
                           TRUE ~ Stage),
                      Gear = case_when(Gear %in% c("FYK (capechads)", "capechades (assemblage of fyke nets)", "FYK (capechades)", "double fyke nets", "Fyke nets", "capetchades")  ~ "FYK",
                                       Gear %in% c("Electrofishing", "electrofishing") ~ "ELE",
                                       TRUE ~ Gear),
                      # Convert characters to factors
                      fSite = as.factor(Site), fSex = as.factor(Sex), fStage = as.factor(Stage),
                      Length_cm = LHT$Length_mm / 10) %>% # Add a column containing length in cm
                      mutate(fSite = droplevels(fSite),
                             fStage = droplevels(fStage)) %>% 
                      relocate(fSite, .after = Site)
# Given the size dimorphism between sexes, we can assign sexes based on length (if they had been not or falsely recorded)
LHT$fStage[which(LHT$fSex == "IND")] <- "Y"   # if the sex is not differentiated yet, eels cannot be silver
LHT$fStage[which(LHT$Length_cm <= 30 & is.na(LHT$fStage))] <- "Y"    # male silver eels are seldom smaller than 30cm, so assign stage to yellow for all eels smaller than that (at least if stage is NA)
LHT$fSex[which(LHT$Length_cm > 46)] <- "F"
LHT$fSex[which(LHT$Length_cm < 45 & LHT$fStage == "S")] <- "M"

# Replace meshsizes reported as range by the largest diameter fished
LHT$Meshsize_mm[grep("-", LHT$Meshsize_mm)] <- sub(".*-", "", LHT$Meshsize_mm[grep("-", LHT$Meshsize_mm)])
LHT$Meshsize_mm[which(LHT$Meshsize_mm == "10mm), BAR (12mm)")] <- 12
LHT$Meshsize_mm[which(LHT$Meshsize_mm == "10mm x 20mm rectangular")] <- NA   # the rectangular barrier spacing cannot really be converted to a meshsize/spacing, so set NA and infer size selectivity based on distribution of captured lengths
LHT$Meshsize_mm <- as.numeric(LHT$Meshsize_mm)


# 0.3) Quality checks, inspection for univariate outliers, outlier removal
if (Rerun) dotchart(LHT$Length_cm)   # no clear outliers
if (Rerun) dotchart(LHT$Weight_g)    # no clear outliers
if (Rerun) dotchart(LHT$Age)   # Lake Garda has the oldest individuals, but no clear outlier

# initial check for outliers in length/weight
LHT %>% ggplot(aes(x = log(Length_cm), y = log(Weight_g))) + geom_point(color = "grey40", alpha = .5) + theme_bw()
# several apparent outliers that might imply protocoling errors. They will be treated in section 3, where weight data is relevant for the L/W relationship


# 0.4) Further preparations 

# Pre-define a color palette for sexes and stages
Palette <- c(M = "steelblue", Male = "steelblue", F = "darkorange3", Female = "darkorange3", Y = "darkgoldenrod2", S = "ivory3", ND = "transparent",
             Yellow = "darkgoldenrod2", Silver = "ivory3")

# Define the maximum lengths by sex (used later to visualize model predictions)
MaxLength <- max(LHT$Length_cm, na.rm = TRUE); MaxLength_M <- 47


# 1) Habitat data import and procession --------------------------------------------------

Habitat <- read_excel("Input/Habitat.xlsx", skip = 1) %>% arrange(SiteAcronym)

# Replace NDs by NAs
Habitat[Habitat == "ND"] <- NA
Habitat[Habitat == "NP"] <- NA
# Correct variable classes
Habitat[] <- lapply(Habitat, function(x) {    # characters to numbers
  if (is.character(x) && all(grepl("^-?\\d*\\.?\\d+$", x[!is.na(x)]))) {
    as.numeric(x)
  } else {
    x
  }
})
Habitat$reclutability <- as.factor(Habitat$reclutability)
Habitat$`Habitat type` <- as.factor(Habitat$`Habitat type`)
Habitat$salinity_max <- as.numeric(Habitat$salinity_max)
Habitat$annual_average_salinity <- as.numeric(Habitat$annual_average_salinity)
Habitat$annual_average_water_temperature <- as.numeric(Habitat$annual_average_water_temperature)
Habitat$trophic_status_chlorophill_a_concentration <- as.numeric(Habitat$trophic_status_chlorophill_a_concentration)

# Rename columns for simplicity
Habitat <- Habitat %>% dplyr::rename(Recruitability = "reclutability", Temperature = "annual_average_water_temperature",
                              Salinity = "annual_average_salinity", Chlorophyll = "trophic_status_chlorophill_a_concentration",
                              TP = "trophic_status_phosphorus_concentration", Nitrogen = "trophic_status_nitrogen_concentration",
                              SurfaceArea = "Current surface (ha)",
                              Habitat = "Habitat type")
# Calculate mean values in case only min and max is stated
IndSal <- which(is.na(Habitat$Salinity) & !is.na(Habitat$salinity_min) & !is.na(Habitat$salinity_max))
Habitat$Salinity[IndSal] <- apply(data.frame(Habitat$Salinity[IndSal], Habitat$Salinity[IndSal]),1,mean)
IndTemp <- which(is.na(Habitat$Temperature) & !is.na(Habitat$water_temperature_min) & !is.na(Habitat$water_temperature_max))
Habitat$Temperature[IndTemp] <- apply(data.frame(Habitat$water_temperature_min[IndTemp], Habitat$water_temperature_max[IndTemp]),1,mean)

# Summarize sites with multiple sampling sites/years to one value per site
HabitatN <- Habitat %>% group_by(Site) %>% tally()
Habitat1 <- Habitat %>% filter(Site %in% unique(LHT$Site) & Site %in% HabitatN$Site[HabitatN$n == 1])   # first, retain values from sites that exactly match the sampling site of eel biometrics
# For sites with multiple observations (i.e. samplings from different years), average them
Habitat2 <- Habitat %>% filter(Site %in% unique(LHT$Site) & Site %in% HabitatN$Site[HabitatN$n > 1]) %>% 
  group_by(Site) %>% summarise(Country = Country[1],
                               Site = Site[1],
                               SiteAcronym = SiteAcronym[1],
                               Habitat = Habitat[1],
                               Latitude = mean(na.omit(Latitude)),
                               Longitude = mean(na.omit(Longitude)),
                               SurfaceArea = sum(na.omit(SurfaceArea)),
                               Recruitability = Recruitability[1],
                               Salinity = mean(na.omit(Salinity)),
                               Temperature = mean(na.omit(Temperature)),
                               water_temperature_min = mean(na.omit(water_temperature_min)),
                               water_temperature_max = mean(na.omit(water_temperature_max)),
                               Chlorophyll = mean(na.omit(Chlorophyll)),
                               TP = mean(na.omit(TP)),
                               Nitrogen = mean(na.omit(Nitrogen)),
                               DistanceGibraltar = first(DistanceGibraltar))
# For site complexes, average environmental variables from the single sampling sites
Habitat3 <- Habitat %>% filter(Site %in% unique(LHT$Site[which(!LHT$Site %in% c(Habitat1$Site, Habitat2$Site))]) | 
                               SiteAcronym %in% unique(LHT$SiteAcronym[which(!LHT$SiteAcronym %in% c(Habitat1$SiteAcronym, Habitat2$SiteAcronym))]) |
                               SiteAcronym %in% unique(LHT$Site[which(!LHT$Site %in% Site)])) %>%   
  group_by(SiteAcronym) %>% summarise(Country = Country[1],
                                      Site = SiteAcronym[1],
                                      SiteAcronym = SiteAcronym[1],
                                      Habitat = Habitat[1],
                                      Latitude = mean(na.omit(Latitude)),
                                      Longitude = mean(na.omit(Longitude)),
                                      SurfaceArea = sum(na.omit(SurfaceArea)),
                                      Recruitability = Recruitability[1],
                                      Salinity = mean(na.omit(Salinity)),
                                      Temperature = mean(na.omit(Temperature)),
                                      water_temperature_min = mean(na.omit(water_temperature_min)),
                                      water_temperature_max = mean(na.omit(water_temperature_max)),
                                      Chlorophyll = mean(na.omit(Chlorophyll)),
                                      TP = mean(na.omit(TP)),
                                      Nitrogen = mean(na.omit(Nitrogen)),
                                      DistanceGibraltar = first(DistanceGibraltar)) %>% 
  ungroup()
HabitatData <- bind_rows(Habitat1, Habitat2, Habitat3)
HabitatData[HabitatData == "NaN"] <- NA


# 1.1) Habitat data availability and selection of variables to retain
HabitatData <- HabitatData %>% dplyr::select("Country", "Site", "SiteAcronym", "Latitude", "Longitude", "Habitat",
                                             "Temperature", "SurfaceArea", "DistanceGibraltar", "Salinity", "Longitude", 
                                             "Latitude", "TP", "Nitrogen", "Chlorophyll")
EnvNRaw <- round((colSums(is.na(HabitatData))/length(unique(HabitatData$Site))) * 100)
InfoVars <- c("Country", "SiteAcronym", "Habitat", "Latitude", "Longitude")
EnvN <- data.frame(Variable = names(EnvNRaw), Value = EnvNRaw) %>% 
  filter(!Variable %in% c(InfoVars, "Site"))

# Plot data availability
EnvPlot <- EnvN %>% ggplot(aes(x = 100 - Value, y = fct_reorder(Variable, -Value))) +
  geom_col(aes(fill = Value), show.legend = FALSE) +
  scale_x_continuous(name = "Data availability (% of sites)") +
  scale_y_discrete(name = element_blank()) +
  theme_bw()
EnvPlot
# ggsave(EnvPlot, file = "Output/EnvVarsSampleSizes.png", height = 4.5, width = 6, units = "in")

# 1.2) Checking for outliers
HabitatDataLong <- reshape2::melt(HabitatData[,-which(names(HabitatData) %in% InfoVars)], id.vars = c("Site")) %>% 
  mutate(value = as.numeric(value))
HabitatDataLong %>% ggplot(aes(x = value, y = Site)) +
  geom_point() +
  facet_wrap(~variable, scales = "free") +
  theme_bw() +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank())
# no clear outliers, any suspicious data were checked

# Plot sites on the map to check for errors in coordinates
load("./Input/PolygonsEMU.rdata")
# Prepare the country shapes of European and North African countries
sf::sf_use_s2(FALSE)
if (Rerun) {
MedArea <- rnaturalearth::ne_countries(scale = 'medium', type = 'map_units',
                        returnclass = 'sf', continent = c('Europe', 'Africa', 'Asia'))
MedArea <- st_crop(MedArea, xmin = -13, xmax = 43,
                   ymin = 28, ymax = 48)

MedArea %>% ggplot() +
  geom_sf(fill = "grey98", size = 0.5, alpha = .5) +
  geom_sf(data = PolygonsEMU, aes(geometry = Coordinates), 
          fill = "grey98", size = 0.5, alpha = .5) +
  geom_point(data = HabitatData, aes(x = Longitude, y = Latitude, fill = Habitat),
             pch = 21, color = "grey30", alpha = .6) #+
  #geom_text(data = HabitatData, aes(label = Site, x = Longitude, y = Latitude))
# looks good!
}

# 1.3) Checking for collinearity between environmental variables
if (Rerun) ggpairs(HabitatData[,-which(names(HabitatData) %in% c(InfoVars, "Site"))])
# Habitat might be problematic, due to apparent correlations with salinity 
# correlations are apparent but some outliers suggest double-checking of some sites
HabitatData %>% ggplot(aes(x = TP, y = Nitrogen)) +
  geom_label(aes(label = Site)) +
  geom_smooth()
# Tonga checked and confirmed
HabitatData %>% ggplot(aes(x = TP, y = Chlorophyll)) +
  geom_label(aes(label = Site)) +
  geom_smooth(se = FALSE, method = lm)
# Guadiaro checked and confirmed; Oubeira has multiple measurements that are averaged, so should be reliable
HabitatData %>% ggplot(aes(x = Nitrogen, y = Chlorophyll)) +
  geom_label(aes(label = Site)) +
  geom_smooth(se = FALSE, method = lm)
# Guadiaro checked and confirmed; Oubeira and Tonga have multiple measurements that are averaged, so should be reliable
HabitatData %>% ggplot(aes(x = Longitude, y = DistanceGibraltar)) +
  geom_label(aes(label = Site)) +
  geom_smooth(se = FALSE, method = lm)


# 1.4) Imputation of missing data
HabitatDataImp <- missForest(as.data.frame(HabitatData %>% dplyr::select(-Country,-SiteAcronym,-Site)))
HabitatDataImp <- as.data.frame(cbind(HabitatData[,1:3], HabitatDataImp$ximp)) %>% 
  mutate(logSurfaceArea = log(SurfaceArea),
         logTP = log(TP),
         logChloro = log(Chlorophyll)) 

# 1.5) z-transformation of environmental variables
InfoCols <- c(1:3,6)
zHabitatDataImp <- HabitatDataImp
zHabitatDataImp[,-InfoCols] <- scale(HabitatDataImp[,-InfoCols])
zHabitatDataImp <- zHabitatDataImp %>% rename(zLatitude = Latitude, zLongitude = Longitude)

# 1.6) Dealing with multicollinearity

# Calculate Variance Inflation Factors to see what other variables might cause trouble due to multicollinearity
corvif(zHabitatDataImp[,c(4:(ncol(zHabitatDataImp)-3))])
# Longitude and Distance-to-Gibraltar were (obviously) highly correlated, so we only keep Distance-to-Gibraltar and drop Longitude
zHabitatDataImp <- zHabitatDataImp %>% dplyr::select(-zLongitude)
corvif(zHabitatDataImp[,c(4:(ncol(zHabitatDataImp)-3))])
# Habitat causes trouble, probably due to its correlation with many env. variables (e.g., salinity)
zHabitatDataImp <- zHabitatDataImp %>% dplyr::select(-Habitat)
corvif(zHabitatDataImp[,c(4:(ncol(zHabitatDataImp)-3))])
# TP now has a high VIF (due to its strong correlation with nitrogen). We want to keep TP, however, because it had more directly observed values (fewer imputed values), so we drop Nitrogen
zHabitatDataImp <- zHabitatDataImp %>% dplyr::select(-Nitrogen)
corvif(zHabitatDataImp[,c(4:(ncol(zHabitatDataImp)-3))])
# --> all VIFs < 3 at this stage, as recommended by Zuur et al., 2010, MEE. But we will not use Latitude as predictor (because of its correlation > 0.7 with temperature, and VIF close to 3)
corvif(zHabitatDataImp[,c(4:(ncol(zHabitatDataImp)-3))] %>% dplyr::select(-zLatitude))
# --> without Latitude all VIFs < 2, so we're on the safe side


# 1.7) Join habitat and life-history tables
LHT_HAB <- left_join(LHT, zHabitatDataImp %>% dplyr::select(-Country, -SiteAcronym), by = c("Site")) %>% 
  mutate(fSiteLabel = as.factor(paste0(fSite, " (", CountryCode, ")")))
  

# 1.8) Visualize raw habitat data distributions by variable

source('VizHabitatDistributions.R')


# 2) Growth ---------------------------------------------------------

# 2.0) Data filtering and procession
EnvVars <- c("Latitude", "Longitude", "DistanceGibraltar", "SurfaceArea",
             "Salinity", "Temperature", "Chlorophyll", "TP", "logSurfaceArea", "logTP", "logChloro")
# Given the strong sex dimorphism, we can include only sexed (or confirmed indeterminate) eels in the analysis
LengthAge <- LHT_HAB %>% dplyr::select(RunningID, fSite, fSiteLabel, CountryCode, Habitat, Age, Length_cm, Weight_g, fStage, fSex, Gear, StagesSelected, EnvVars) %>%
  filter(fSex %in% c("M", "F")) %>%     
  mutate(fSex = droplevels(fSex),
         Age = floor(Age))    # round down age to the fully completed year
# Drop sites without or with insufficient age observations
LengthAge$CaseID <- paste(LengthAge$fSite, LengthAge$fSex)
AgeObsSex <- LengthAge %>% group_by(CaseID) %>% summarise(n = length(which(!is.na(Age))))   # Calculate sample sizes
# Delete cases with < 10 observations for any sex
LengthAge <- LengthAge %>% filter(CaseID %in% AgeObsSex$CaseID[which(AgeObsSex$n >= 10)]) %>% 
  mutate(fSite = droplevels(fSite))

# Inspect raw data length and age data for outliers
LengthAge %>% ggplot(aes(x = as.factor(Age), y = Length_cm)) +
  geom_point(aes(color = fSex), alpha = .5, show.legend = FALSE) +
  geom_boxplot(fill = "transparent", color = "grey30", outlier.color = alpha("grey70", .4)) +
  ggforce::facet_wrap_paginate(~fSiteLabel, ncol = 5, nrow = 5, page = 2) +
  theme_bw()
# exclude one outlier in Mellah lagoon (that was also excluded in the original paper and analysis by Tahri & Panfili, 2023). The rest looks ok.
LengthAge <- LengthAge[-which(LengthAge$Age == 3 & LengthAge$Length_cm > 100),]


# Inspect distribution of observed lengths and ages (see if there is a more or less continuous gradient)
LengthAge$ObsID <- as.factor(row.names(LengthAge))
LengthAge %>% ggplot(aes(x = Length_cm, y = ObsID)) +
  geom_point(aes(color = fSex)) +
  geom_vline(xintercept = 47, linetype = 2, linewidth = 1.5) +
  scale_x_continuous(breaks = c(20,40,60,80,100))
# looks good
LengthAge %>% ggplot(aes(x = Age, y = ObsID)) +
  geom_point(aes(color = fSex))
# looks good

# Create a subset df of only eels with observed ages, because the following mixed-effects models cannot handle NA observations
LengthAgeObs <- LengthAge %>% drop_na(Age)


# 2.1) Fit and select the hierarchical von Bertalanffy growth model(s), first without, then with environmental covariates

# add a binary variable indicating whether silver eels had been sampled exclusively or whether it was a mixed-stage sampling in any site
nObsLengthAge <- LengthAgeObs %>% group_by(fSite) %>% summarise(nTot = length(which(!is.na(Age))),
                                                                nSilver = length(which(!is.na(Age) & fStage == "S"))) %>% 
  mutate(Selective = ifelse(nTot == nSilver, 1, 0))
LengthAgeObs$StageSelectivity <- ifelse((LengthAgeObs$fSite %in% nObsLengthAge$fSite[which(nObsLengthAge$Selective == 1)]), 1, 0)
LengthAgeObs$StageSelectivity <- as.factor(LengthAgeObs$StageSelectivity)

if (Rerun) source('GrowthBaseSelection.R')
  
# extract and save predictions from the base model
source('GrowthBaseModelPredictions.R')

# Based on results for the "aged-only" dataset, assign ages to not aged fish using an age-length-key (ALK)
source('ALKAgeAssignment.R')


# 2.2) Model the environmental correlates of growth

# Explore distributions of the environmental variables, to check for outliers and the need to transform environmental data before running the models
dotchart(LengthAgeObs$Temperature)    # roughly ok, with a slight outlier at the lower end
dotchart(LengthAgeObs$Salinity)   # good
dotchart(LengthAgeObs$DistanceGibraltar)    # Guadiaro constitutes an outlier at the lower end
dotchart(LengthAgeObs$TP)   # skewed distribution, yet still a continuous gradient
dotchart(LengthAgeObs$Chlorophyll)    # skewed distribution, yet still a continuous gradient
dotchart(LengthAgeObs$SurfaceArea)    # must be log-transformed given unequal distribution with a strong outlier (Venezia lagoon)
dotchart(LengthAgeObs$logSurfaceArea)   # better
dotchart(LengthAgeObs$logTP)
dotchart(LengthAgeObs$logChloro)
# log-transformation of chlorophyll and TP could be justified, yet we still have a continuous gradient without outliers on the original scale. We will base the decision whether to transform on an AIC comparison in the global model!


# For salinity, we hypothesize a unimodal relationship (i.e., better growth at intermediate salinity levels and slow growth at extreme levels - low or high)
# Square the salinity values
LengthAgeObs$sqSalinity <- LengthAgeObs$Salinity^2

if (Rerun) source('GrowthEnvSelection.R')


# 2.3) Visualization of site-specific results from the base model (without environmental covariates)

# Get site-specific maximum ages
MaxAges <- LengthAgeObs %>% group_by(fSite) %>% summarise(MaxAge = max(Age))
Coefs$MaxAge <- MaxAges$MaxAge[match(Coefs$Site, MaxAges$fSite)]
Coefs$fSiteLabel <- LengthAgeObs$fSiteLabel[match(Coefs$Site, LengthAgeObs$fSite)]
SitePreds$Group <- paste(SitePreds$fSite, SitePreds$fSex)
Coefs2 <- Coefs
Coefs2$Linf_M <- Coefs2$Linf_Intercept+Coefs2$Linf_M
Coefs2$k_M <- Coefs2$k_Intercept+Coefs2$k_M
CoefsLong <- reshape2::melt(Coefs2, id.vars = c("Site", "SiteAcronym", "fSiteLabel", "MaxAge")) %>% mutate(fSex = fct_recode(variable, 
                                                                                                    "F" = "Linf_Intercept",
                                                                                                    "F" = "k_Intercept",
                                                                                                    "M" = "Linf_M",
                                                                                                    "M" = "k_M")) %>% 
  mutate(Group = paste(Site, fSex), value = as.numeric(value))
CoefsLong <- CoefsLong[which(CoefsLong$Group %in% SitePreds$Group),]
Coefs2 <- pivot_wider(CoefsLong %>% dplyr::select(-fSex, -Group), names_from = "variable", values_from = "value")

# Add an "ND" level for the factor "stage"
LengthAgeObs$fStage <- as.character(LengthAgeObs$fStage)
LengthAgeObs$fStage[which(is.na(LengthAgeObs$fStage))] <- "ND"
LengthAgeObs$fStage <- as.factor(LengthAgeObs$fStage)

# Prepare a key for habitat-specific background colors of the strip
StripColors <- LengthAgeObs %>% group_by(fSiteLabel) %>% reframe(Habitat = Habitat[1]) %>% 
  mutate(Color = case_when(Habitat == "RIV" ~ alpha("palegreen2", .4),
                           Habitat == "RIE" ~ alpha("palegreen4", .4),
                           Habitat == "LGN" ~ alpha("deepskyblue3", .4),
                           Habitat == "LAK" ~ alpha("gold", .4),
                           Habitat == "CNL" ~ alpha("peru", .4)))

Ind <- 1:30    # indices for the 1st plot
Ind2 <- 31:length(levels(LengthAgeObs$fSite))    # indices for the 2nd plot

Strip <- strip_themed(background_x = elem_list_rect(fill = StripColors$Color[Ind]))
Strip2 <- strip_themed(background_x = elem_list_rect(fill = StripColors$Color[Ind2]))

GRCoefs1 <- Coefs2 %>% mutate(Site = as.factor(Site)) %>% 
  filter(Site %in% levels(SitePreds$fSite)[Ind]) %>% mutate(Site = droplevels(Site))
GRCoefs2 <- Coefs2 %>% mutate(Site = as.factor(Site)) %>% 
  filter(Site %in% levels(SitePreds$fSite)[Ind2]) %>% mutate(Site = droplevels(Site))

GP <- LengthAgeObs %>% filter(fSite %in% levels(fSite)[Ind]) %>% 
  ggplot(aes(x = Age, color = fSex)) +
  geom_jitter(aes(y = Length_cm, shape = fStage), fill = "transparent", alpha = .4, width = .1, size = .6) +
  geom_line(data = SitePreds %>% filter(fSite %in% levels(fSite)[Ind]), aes(y = Predicted), alpha = .8, linewidth = .8) +
  geom_text(data = GRCoefs1 %>% filter(!is.na(Linf_Intercept)), aes(x = 0.22*MaxAge, y = 103, label = paste("L[\"inf\"]['_f'] == ", round(Linf_Intercept, 1))), parse = TRUE, color = "grey25", size = 1.5) +
  geom_text(data = GRCoefs1 %>% filter(!is.na(k_Intercept)), aes(x = 0.18*MaxAge, y = 96, label = paste("k[\"f\"] == ", round(k_Intercept, 2))), parse = TRUE, color = "grey25", size = 1.5) +
  geom_text(data = GRCoefs1 %>% filter(!is.na(Linf_M)), aes(x = 0.22*MaxAge, y = 89, label = paste("L[\"inf\"]['_m'] == ", round(Linf_M, 1))), parse = TRUE, color = "grey25", size = 1.5) +
  geom_text(data = GRCoefs1 %>% filter(!is.na(k_M)), aes(x = 0.18*MaxAge, y = 82, label = paste("k[\"m\"] == ", round(k_M, 2))), parse = TRUE, color = "grey25", size = 1.5) +
  scale_shape_manual(values = c(Y = 25, S = 24, ND = 21)) +
  scale_color_manual(values = Palette) +
  scale_x_continuous(name = "Continental age (years)", breaks = breaks_pretty(5)) +
  scale_y_continuous(name = "Length (cm)", breaks = seq(20, max(LengthAge$Length_cm), by = 20)) +
  guides(fill = guide_legend(title = "Sex", nrow = 1, override.aes = list(size = 2)), 
         color = guide_legend(title = "Sex", nrow = 1, override.aes = list(size = 2)), 
         shape = guide_legend(title = "Stage", nrow = 1, override.aes = list(size = 2))) +
  theme_bw() +
  theme(strip.text = element_text(size = 3.5, margin = margin(.05,0,.05,0, "cm")),
        strip.background = element_rect(fill = StripColors$Color[Ind]),
        axis.text.x = element_text(size = 5, margin = margin(1, 0, 0, 0)),
        axis.text.y = element_text(size = 5, margin = margin(0, 1, 0, 0)),
        axis.title = element_text(size = 8),
        axis.ticks = element_line(linewidth = .2),
        legend.title.position = "top",
        legend.box = "horizontal",
        panel.grid.minor = element_line(linewidth = .1),
        panel.grid.major = element_line(linewidth = .2)) +
  facet_wrap2(~fSiteLabel, scales = "free_x", strip = Strip, ncol = 5, nrow = 6)
  #ggforce::facet_wrap_paginate(~fSiteLabel, scales = "free_x", ncol = 7, nrow = 4, page = 1)
GPLegend <- cowplot::get_legend(GP)
GPMain <- GP + theme(legend.position = "none")
GPMain
#ggsave(GPMain, file = "./Output/Sitelevel/Growth_by_site_page1.png", dpi = 600, height = 6, width = 5, units = "in")

GP2 <- LengthAgeObs %>% filter(fSite %in% levels(fSite)[Ind2]) %>% 
  ggplot(aes(x = Age, color = fSex)) +
  geom_jitter(aes(y = Length_cm, shape = fStage), fill = "transparent", alpha = .4, width = .1, size = .6) +
  geom_line(data = SitePreds %>% filter(fSite %in% levels(fSite)[Ind2]), aes(y = Predicted), alpha = .8, linewidth = .8) +
  geom_text(data = GRCoefs2 %>% filter(!is.na(Linf_Intercept)), aes(x = 0.22*MaxAge, y = 105, label = paste("L[\"inf\"]['_f'] == ", round(Linf_Intercept, 1))), parse = TRUE, color = "grey25", size = 1.5) +
  geom_text(data = GRCoefs2 %>% filter(!is.na(k_Intercept)), aes(x = 0.18*MaxAge, y = 97.5, label = paste("k[\"f\"] == ", round(k_Intercept, 2))), parse = TRUE, color = "grey25", size = 1.5) +
  geom_text(data = GRCoefs2 %>% filter(!is.na(Linf_M)), aes(x = 0.22*MaxAge, y = 90, label = paste("L[\"inf\"]['_m'] == ", round(Linf_M, 1))), parse = TRUE, color = "grey25", size = 1.5) +
  geom_text(data = GRCoefs2 %>% filter(!is.na(k_M)), aes(x = 0.18*MaxAge, y = 82.5, label = paste("k[\"m\"] == ", round(k_M, 2))), parse = TRUE, color = "grey25", size = 1.5) +
  scale_shape_manual(values = c(Y = 25, S = 24, ND = 21)) +
  scale_color_manual(values = Palette) +
  scale_x_continuous(name = "Continental age (years)", breaks = breaks_pretty(5)) +
  scale_y_continuous(name = "Length (cm)", breaks = seq(20, max(LengthAge$Length_cm), by = 20)) +
  guides(fill = guide_legend(title = "Sex", nrow = 1, override.aes = list(size = 2)), 
         color = guide_legend(title = "Sex", nrow = 1, override.aes = list(size = 2)), 
         shape = guide_legend(title = "Stage", nrow = 1, override.aes = list(size = 2))) +
  theme_bw() +
  theme(strip.text = element_text(size = 4, margin = margin(.05,0,.05,0, "cm")),
        strip.background = element_rect(fill = StripColors$Color[Ind2]),
        axis.text.x = element_text(size = 5, margin = margin(1, 0, 0, 0)),
        axis.text.y = element_text(size = 5, margin = margin(0, 1, 0, 0)),
        axis.title = element_text(size = 8),
        axis.ticks = element_line(linewidth = .2),
        legend.title.position = "top",
        legend.box = "horizontal",
        panel.grid.minor = element_line(linewidth = .1),
        panel.grid.major = element_line(linewidth = .2)) +
  facet_wrap2(~fSiteLabel, scales = "free_x", strip = Strip2, ncol = 5, nrow = 5)
#ggforce::facet_wrap_paginate(~fSiteLabel, scales = "free_x", ncol = 7, nrow = 4, page = 1)
GPMain2 <- GP2 + theme(legend.position = "none")
GPMain2
#ggsave(GPMain2, file = "./Output/Sitelevel/Growth_by_site_page2.png", dpi = 600, height = 5, width = 5, units = "in")


# 2.4) Visualization of population-level results from the base model (without environmental covariates)
SitePreds$Group <- paste(SitePreds$fSite, SitePreds$fSex)
GrowthPlotPop <- SitePreds %>% ggplot(aes(x = Age)) +
  geom_line(aes(y = Predicted, group = Group, color = fSex), linewidth = 1.1, alpha = .3, show.legend = FALSE) +
  geom_line(data = PopPreds, aes(y = Predicted, color = fSex), linewidth = 1.7, show.legend = FALSE) +
  geom_ribbon(data = PopPreds, aes(ymin = Lower, ymax = Upper, fill = fSex), alpha = .25, show.legend = FALSE) +
  scale_x_continuous(name = "Continental age (years)", breaks = seq(0, 22, by = 2)) +
  scale_y_continuous(name = "Length (cm)", breaks = seq(10, max(LengthAge$Length_cm), by = 10)) +
  scale_color_manual(values = Palette) +
  scale_fill_manual(values = Palette) +
  # guides(fill = guide_legend(title = "Sex"), color = guide_legend(title = "Sex")) +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 13),
        legend.position = c(0.92, 0.15),
        legend.background = element_rect(color = "grey40"))
GrowthPlotPop
#ggsave(GrowthPlotPop, file = "Output/Poplevel/Growth_poplevel_sexes_combined.png", dpi = 400, height = 5, width = 6, units = "in")


# 2.5) Environmental effects - results extraction and visualization

source('VizEnvEffectsGrowth.R')


# 3) Length-weight relationship -------------------------------------------

# 3.1) Data procession

# We will run a linear model on log-transformed values. Thus log-transform length and weight values
LHT_HAB$logLength_cm <- log(LHT_HAB$Length_cm)
LHT_HAB$logWeight_g <- log(LHT_HAB$Weight_g)

# Raw data inspection & checking for outliers
LHT_HAB %>% ggplot(aes(x = Length_cm, y = Weight_g)) +
  geom_point(aes(color = fSex))   # some potential outliers apparent. We will detect and treat them after fitting the model. 


# 3.2) Fit the base mixed-effects model of length-weight-relationship, not accounting for any covariates

source('LWBaseModel.R')

# Plot results for mixed sexes and stages - population level
a <- round(exp(PopIntercept), 4); b <- round(PopSlope, 4)
LWPlotPop <- LHT_HAB2 %>% ggplot() +
  geom_ribbon(data = OutputPop, aes(x = exp(logLength_cm), ymin = exp(Lower) / 1000, ymax = exp(Upper) / 1000), color = "transparent", fill = "firebrick2", alpha = .1) +
  geom_line(data = SitePredsLW, aes(x = exp(logLength_cm), y = exp(Predicted) / 1000, group = fSite), color = "grey50") +
  geom_line(data = OutputPop, aes(x = exp(logLength_cm), y = exp(Predicted) / 1000), linewidth = 1.2, color = "firebrick") +
  scale_x_continuous(name = "Length (cm)", limits = c(0, 105), breaks = seq(20, 100, by = 20)) +
  scale_y_continuous(name = "Weight (kg)", limits = c(0, 3.2), breaks = seq(0, 3, by = .5), labels = label_number(drop0trailing = TRUE)) +
  guides(color = guide_legend(ncol = 2)) +
  #annotate("text", label = bquote(italic(W) == .(a) * italic(L)^.(b)), x = 15, y = 3300) +
  theme_bw() +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12))
LWPlotPop
#ggsave(LWPlotPop, file = "Output/Poplevel/Length-weight_poplevel.png", dpi = 600, height = 5, width = 6, units = "in")  


# 3.3) Fit a stage-specific mixed-effects model of length-weight-relationship, not accounting for environmental covariates

# Create a copy of the LHT df
LHT_HAB3 <- LHT_HAB2

if (Rerun) source('LWStageModel.R')


# 3.4) Visualize the results of the stage-specific models without environmental covariates

load("./Output/LW_preds_sitelevel_stage.rdata")
load("./Output/LWCoefsStage.rdata")

# Plot site-level results by life-stage
# Prepare a key for habitat-specific background colors of the strip
StripColors <- LHT_HAB3 %>% group_by(fSiteLabel) %>% reframe(Habitat = Habitat[1]) %>%
  mutate(Color = case_when(Habitat == "RIV" ~ alpha("palegreen2", .4),
                           Habitat == "RIE" ~ alpha("palegreen4", .4),
                           Habitat == "LGN" ~ alpha("deepskyblue3", .4),
                           Habitat == "LAK" ~ alpha("gold", .4),
                           Habitat == "CNL" ~ alpha("peru", .4)))

length(unique(LHT_HAB2$fSite))
Ind <- 1:35    # indices for the 1st plot
Ind2 <- 35:length(levels(LHT_HAB3$fSite))    # indices for the 2nd plot

Strip <- strip_themed(background_x = elem_list_rect(fill = StripColors$Color[Ind]))
Strip2 <- strip_themed(background_x = elem_list_rect(fill = StripColors$Color[Ind2]))


# First page of the site- and stage-specific length-weight plot
LWPlotSiteStage <- LHT_HAB3 %>% filter(fSite %in% levels(fSite)[Ind]) %>% 
  ggplot(aes(x = Length_cm)) +
  geom_point(aes(y = Weight_g/1000, fill = fStage), alpha = .35, stroke = .2, size = .8, pch = 21, color = "grey30") +
  geom_line(data = SitePredsLWStage %>% filter(fSite %in% levels(fSite)[Ind]), 
            aes(y = Weight_g/1000, color = fStage), linewidth = .5) +
  geom_text(data = CoefsLW %>% filter(fSite %in% levels(fSite)[Ind] & fStage == "Yellow"), aes(x = 28, y = 2.6, label = paste0("a[\"y\"] == ", round(a, 6)), color = fStage), size = 1.5, parse = TRUE, show.legend = FALSE) +
  geom_text(data = CoefsLW %>% filter(fSite %in% levels(fSite)[Ind] & fStage == "Yellow"), aes(x = 28, y = 2.3, label = paste0("b[\"y\"] == ", round(b, 6)), color = fStage), size = 1.5, parse = TRUE, show.legend = FALSE) +
  geom_text(data = CoefsLW %>% filter(fSite %in% levels(fSite)[Ind] & fStage == "Silver"), aes(x = 28, y = 2, label = paste0("a[\"s\"] == ", round(a, 6))), color = "ivory4", size = 1.5, parse = TRUE, show.legend = FALSE) +
  geom_text(data = CoefsLW %>% filter(fSite %in% levels(fSite)[Ind] & fStage == "Silver"), aes(x = 28, y = 1.7, label = paste0("b[\"s\"] == ", round(b, 6))), color = "ivory4", size = 1.5, parse = TRUE, show.legend = FALSE) +
  geom_text(data = CoefsLW %>% filter(fSite %in% levels(fSite)[Ind] & fStage == "ND"), aes(x = 28, y = 2.6, label = paste("a =", round(a, 6))), color = "grey30", size = 1.5, show.legend = FALSE) +
  geom_text(data = CoefsLW %>% filter(fSite %in% levels(fSite)[Ind] & fStage == "ND"), aes(x = 28, y = 2.3, label = paste("b =", round(b, 6))), color = "grey30", size = 1.5, show.legend = FALSE) +
  scale_fill_manual(values = Palette) +
  scale_color_manual(values = c(Silver = "ivory3", Yellow = "darkgoldenrod3", ND = "grey40")) +
  scale_x_continuous(name = "Length (cm)", breaks = seq(20, 100, by = 20)) +
  scale_y_continuous(name = "Weight (kg)", breaks = seq(.5, 3, by = .5), labels = label_number(drop0trailing = TRUE)) +
  guides(fill = guide_legend(title = "Stage", ncol = 1), color = guide_legend(title = "Stage", ncol = 1)) +
  theme_bw() +
  theme(strip.text = element_text(size = 3.8, margin = margin(.05,0,.05,0, "cm")),
        axis.text.x = element_text(size = 7, margin = margin(1, 0, 0, 0)),
        axis.text.y = element_text(size = 6, margin = margin(0, 1, 0, 0)),
        axis.title = element_text(size = 9),
        axis.ticks = element_line(linewidth = .2),
        legend.title.position = "top",
        legend.box = "horizontal",
        panel.grid.minor = element_line(linewidth = .1),
        panel.grid.major = element_line(linewidth = .2)) +
  facet_wrap2(~fSiteLabel, strip = Strip, ncol = 5, nrow = 7)
LWPlotSiteStage
LWLegendStage <- cowplot::get_legend(LWPlotSiteStage)
LWMainStage <- LWPlotSiteStage + theme(legend.position = "none")
LWMainStage
# ggsave(LWMainStage, file = "./Output/Sitelevel/LW_by_stage_by_site_page1.png", dpi = 600, height = 6, width = 5, units = "in")

# Second page plot
LWPlotSiteStage2 <- LHT_HAB3 %>% filter(fSite %in% levels(fSite)[Ind2]) %>% 
  ggplot(aes(x = Length_cm)) +
  geom_point(aes(y = Weight_g/1000, fill = fStage), alpha = .35, stroke = .2, size = .8, pch = 21, color = "grey30") +
  geom_line(data = SitePredsLWStage %>% filter(fSite %in% levels(fSite)[Ind2]), 
            aes(y = Weight_g/1000, color = fStage), linewidth = .5) +
  geom_text(data = CoefsLW %>% filter(fSite %in% levels(fSite)[Ind2] & fStage == "Yellow"), aes(x = 33, y = 2.8, label = paste0("a[\"y\"] == ", round(a, 6)), color = fStage), size = 1.5, parse = TRUE, show.legend = FALSE) +
  geom_text(data = CoefsLW %>% filter(fSite %in% levels(fSite)[Ind2] & fStage == "Yellow"), aes(x = 33, y = 2.5, label = paste0("b[\"y\"] == ", round(b, 6)), color = fStage), size = 1.5, parse = TRUE, show.legend = FALSE) +
  geom_text(data = CoefsLW %>% filter(fSite %in% levels(fSite)[Ind2] & fStage == "Silver"), aes(x = 33, y = 2.2, label = paste0("a[\"s\"] == ", round(a, 6))), color = "ivory4", size = 1.5, parse = TRUE, show.legend = FALSE) +
  geom_text(data = CoefsLW %>% filter(fSite %in% levels(fSite)[Ind2] & fStage == "Silver"), aes(x = 33, y = 1.9, label = paste0("b[\"s\"] == ", round(b, 6))), color = "ivory4", size = 1.5, parse = TRUE, show.legend = FALSE) +
  scale_fill_manual(values = Palette) +
  scale_color_manual(values = c(Silver = "ivory3", Yellow = "darkgoldenrod3", ND = "grey40")) +
  scale_x_continuous(name = "Length (cm)", breaks = seq(20, 100, by = 20)) +
  scale_y_continuous(name = "Weight (kg)", breaks = seq(.5, 3, by = .5), limits = c(0,3), labels = label_number(drop0trailing = TRUE)) +
  guides(fill = guide_legend(title = "Stage", ncol = 1), color = guide_legend(title = "Stage", ncol = 1)) +
  theme_bw() +
  theme(strip.text = element_text(size = 3.8, margin = margin(.05,0,.05,0, "cm")),
        axis.text.x = element_text(size = 6, margin = margin(1, 0, 0, 0)),
        axis.text.y = element_text(size = 7, margin = margin(0, 1, 0, 0)),
        axis.title = element_text(size = 9),
        axis.ticks = element_line(linewidth = .2),
        legend.position = "none",
        panel.grid.minor = element_line(linewidth = .1),
        panel.grid.major = element_line(linewidth = .2)) +
  facet_wrap2(~fSiteLabel, strip = Strip2, ncol = 5, nrow = 6)
# ggsave(LWPlotSiteStage2, file = "./Output/Sitelevel/LW_by_stage_by_site_page2.png", dpi = 600, height = 5.5, width = 4.7, units = "in")


# Plot population-level predictions of the stage-specific model

load("./Output/LW_preds_poplevel_stage.rdata")

# Add a case ID to the site-level predictions df
SitePredsLWStage$CaseID <- paste(SitePredsLWStage$fSite, SitePredsLWStage$fStage)

Palette2 <- Palette
Palette2[which(names(Palette2) == "Silver")] <- "#B8B8AF"

LWPlotPopStage <- SitePredsLWStage %>% filter(fStage != "ND") %>% mutate(fStage = droplevels(fStage)) %>%  
  ggplot(aes(x = exp(logLength_cm), fill = fStage, color = fStage)) +
  geom_ribbon(data = PopPredsLWStage, aes(ymin = exp(Lower) / 1000, ymax = exp(Upper) / 1000), color = "transparent", alpha = .4) +
  geom_line(aes(x = Length_cm, y = Weight_g / 1000, group = CaseID), alpha = .3, linewidth = .7) +
  geom_line(data = PopPredsLWStage, aes(y = exp(Predicted) / 1000), linewidth = 1.3, alpha = .8) +
  scale_fill_manual(values = Palette2) +
  scale_color_manual(values = Palette2) +
  scale_x_continuous(name = "Length (cm)", breaks = seq(20, 100, by = 20), limits = c(0, 105)) +
  scale_y_continuous(name = "Weight (kg)", breaks = seq(0, 3, by = .5), limits = c(0,3.2), labels = label_number(drop0trailing = TRUE)) +
  guides(fill = guide_legend(title = "Stage"), color = guide_legend(title = "Stage")) +
  theme_bw() +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.position = c(0.02, .98),         
        legend.justification = c(0, 1),
        legend.key.size = unit(.8, "cm"),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10),
        legend.background = element_rect(fill = "white", color = "grey40", linewidth = .3))    
LWPlotPopStage
#ggsave(LWPlotPopStage, file = "Output/Poplevel/Length_weight_staged_poplevel.png", dpi = 600, width = 6, height = 5, units = "in")


# 3.5) Fit the length-weight relationship accounting for environmental covariates

if (Rerun) source('LWEnvModel.R')


# 4) Sex ratio ---------------------------------------------------------------

# 4.1) Identify suitable datasets that can be used to infer sex ratio, drop indeterminate or not sexed individuals

# Inspect distribution of captured lengths of other monitorings with unknown mesh size, to decide whether they should be included or not
SRDataRaw <- LHT_HAB2
SRDataRaw$fRunningID <- as.factor(SRDataRaw$RunningID)
SRDataRaw$fLabel <- paste0(SRDataRaw$fSite, " - ", SRDataRaw$Gear, " - ", SRDataRaw$Meshsize_mm, "mm")

# Calculate no. of observations for each site and exclude datasets with an effective sample size < 20 eels
(CountsBySite <- SRDataRaw %>% group_by(fSite) %>% summarise(n = length(which(fSex %in% c("M", "F")))))
Exclude <- CountsBySite$fSite[which(CountsBySite$n < 20)]
SRDataRaw <- SRDataRaw %>% filter(!fSite %in% Exclude & !is.na(fSex)) %>% mutate(fSite = droplevels(fSite))
# Exclude sites / samplings for which sexes have been selected by the researchers
SRDataRaw <- SRDataRaw %>% filter(is.na(SexesSelected))

# Visual inspection of size-selectivity of the gear in any site, to aid decision about site inclusion in sex ratio analysis
SRDataRaw$fStage <- as.character(SRDataRaw$fStage)
SRDataRaw$fStage[which(is.na(SRDataRaw$fStage))] <- "ND"
SRDataRaw$fStage <- as.factor(SRDataRaw$fStage)
SRDataRaw %>% ggplot(aes(x = Length_cm, y = fRunningID, fill = fStage)) +
  geom_vline(aes(xintercept = 33.5), color = "darkgoldenrod3") +
  # geom_vline(aes(xintercept = 35), color = "ivory3") +
  geom_point(pch = 21, color = "grey30") +
  scale_x_continuous(breaks = seq(20, 100, by = 20)) +
  scale_y_discrete(name = element_blank()) +
  scale_fill_manual(values = Palette) +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        axis.text.y = element_blank(), 
        strip.text = element_text(size = 8)) +
  facet_wrap(~fLabel, scales = "free_y")
# Depict the full length spectrum and can be kept for analyses of sex ratio:
CasesToExclude <- c("Akgol - BAR - NAmm", "Akgol - FYK - NAmm", "Bafa - FYK - 12mm", "Burullus - FYK - NAmm", "Gala - FYK - 12mm", "Garda - FYK - NAmm",
                    "Karamenderes - FYK - 12mm", "Marano - FYK - 20mm", "Or - FYK - 6mm", "Or - FYK - 10mm", 
                    "Salses-Leucate - FYK - 12mm", "Vistonida - BAR - 30mm")
CasesToExcludeYellow <- c("Tonga - FYK - 20mm", "Mellah - FYK - 20mm")                    
# Subset the sites/samplings that meet the criteria for selectivity from the global dataframe
SRData <- SRDataRaw %>% filter(!fLabel %in% CasesToExclude & !(fLabel %in% CasesToExcludeYellow & fStage == "Yellow"))

SRDataRaw$Excluded <- 0
SRDataRaw$Excluded[which(!SRDataRaw$fLabel %in% unique(SRData$fLabel))] <- 1

StripColors <- SRDataRaw %>% group_by(fLabel) %>% reframe(Excluded = Excluded[1]) %>% 
  mutate(Color = case_when(Excluded == 0 ~ alpha("palegreen2", .4),
                           Excluded == 1 ~ alpha("firebrick3", .4)))
Strip <- strip_themed(background_x = elem_list_rect(fill = StripColors$Color))

SRDataRaw %>% ggplot(aes(x = Length_cm, y = fRunningID, fill = fStage)) +
  geom_vline(aes(xintercept = 33.5), color = "darkgoldenrod3", show.legend = FALSE) +
  # geom_vline(aes(xintercept = 35), color = "ivory3", show.legend = FALSE) +
  geom_point(pch = 21, color = "grey30", show.legend = FALSE) +
  scale_x_continuous(breaks = seq(20, 100, by = 20)) +
  scale_y_discrete(name = element_blank()) +
  scale_fill_manual(values = Palette) +
  guides(fill = guide_legend(title = "Stage", nrow = 1)) +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        axis.text.y = element_blank(), 
        axis.ticks = element_blank(),
        strip.text = element_text(size = 7)) +
  facet_wrap2(~fLabel, strip = Strip, scales = "free_y")

SRData <- SRData %>% filter(fSex %in% c("M", "F"))


# Recode sampling years reported as range
HyphenInd <- grep("-", SRData$SamplingYear)
SRData$SamplingYear[HyphenInd] <- sapply(SRData$SamplingYear[HyphenInd], function(x) {
  Parts <- strsplit(x, "-")[[1]]
  round(mean(as.numeric(Parts)))
})
SlashInd <- grep("/", SRData$SamplingYear)
SRData$SamplingYear[SlashInd] <- sapply(SRData$SamplingYear[SlashInd], function(x) {
  as.numeric(strsplit(x, "/")[[1]][1])
})
SRData$nYear <- as.numeric(SRData$SamplingYear) 


# Calculate observed sex ratios
SRSummary <- SRData %>% group_by(fSite) %>%
  summarise(fSite = head(fSite, 1), SiteAcronym = head(SiteAcronym, 1), Ntot = n(),
            Nmales = sum(fSex == "M"), Nfemales = sum(fSex == "F"), MalesFr = Nmales / Ntot,
            Meshsize_mm = mean(na.omit(Meshsize_mm)), Gear = first(Gear))


# 4.2) Fit the base binomial GLMM for sex ratio, without environmental covariates

source('SRBaseModel.R')



# 4.3) Mixed-effects Bernoulli model of sex ratio as a function of environmental covariates

if (Rerun) source('SREnvModel.R')


# 4.5) Visualize the results

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
PieRadius <- .3

SRMResultsLong$fSex <- factor(SRMResultsLong$fSex, levels = c("Male", "Female"))


# Draw sex ratios on the map as pie charts
SRMap <- ggplot() +
  geom_sf(data = MedArea, fill = "grey98", color = "grey50", size = 0.5, alpha = .5) +
  geom_sf(data = PolygonsEMU, aes(geometry = Coordinates), 
          fill = "grey98", size = 0.5, alpha = .5, color = "grey50") +
  geom_arc_bar(data = SRMResultsLong %>% filter(Moved == 0),
               aes(x0 = plotLongitude, y0 = plotLatitude, r0 = 0, r = PieRadius * cos(45 * pi / 180) * 2,
                   start = Start * 2 * pi, end = End * 2 * pi, fill = fSex),
               inherit.aes = FALSE, color = "grey80", linewidth = .2) +
  geom_segment(data = SRMResults %>% filter(Moved == 1), 
               aes(x = plotLongitude, xend = Longitude, y = plotLatitude, yend = Latitude),
               linewidth = .3) +
  geom_arc_bar(data = SRMResultsLong %>% filter(Moved == 1),
               aes(x0 = plotLongitude, y0 = plotLatitude, r0 = 0, r = PieRadius * cos(45 * pi / 180) * 2,
                   start = Start * 2 * pi, end = End * 2 * pi, fill = fSex),
               inherit.aes = FALSE, color = "grey80", linewidth = .2) +
  #geom_text(data = SRMResultsLong %>% filter(Moved == 0), aes(x = plotLongitude, y = plotLatitude, label = fSite), size = 2) +
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
SRMap
#ggsave(SRMap, file = "./Output/Sitelevel/Sex_ratio_map.png", dpi = 600, height = 3.5, width = 6, units = "in")


# 4.6) Visualize the environmental effects on sex ratio

source('VizEnvEffectsSexRatio.R')


# 5) Silvering function ---------------------------------------------------

# 5.0) Data procession and site selection

# We will work with the dataset in which length-weight regression outliers had already been excluded, as their length data might be unreliable
LHTStaged <- LHT_HAB2 %>% mutate(fStage = fct_recode(fStage, Y = "G")) %>% filter(!is.na(fStage), fSex %in% c("F", "M")) %>%   # subset only eels that were staged and sexed
                         mutate(fStage = droplevels(fStage), fSex = droplevels(fSex))

# Inspect data availability and delete sites with too low sample sizes / unsuitable sampling scheme (e.g., only silver-eel selective)
LHTStaged$bSilver <- as.factor(ifelse(LHTStaged$fStage == "S", 1, 0))

# Delete stage-selective samplings
LHTStaged <- LHTStaged[-which(LHTStaged$StagesSelected == "T"),]
SamplesSlv <- LHTStaged %>% group_by(fSite) %>% summarise(Ntot = n(), Nmales = sum(fSex == "M"), Nfemales = sum(fSex == "F"),
            Nyellow = sum(bSilver == 0), Nsilver = sum(bSilver == 1))
SamplesSlv$Exclude <- ifelse(SamplesSlv$Nsilver == SamplesSlv$Ntot | SamplesSlv$Nyellow == SamplesSlv$Ntot, 1, 0)
LHTStaged <- LHTStaged %>% filter(!fSite %in% SamplesSlv$fSite[which(SamplesSlv$Nsilver == SamplesSlv$Ntot)] &
                                  !fSite %in% SamplesSlv$fSite[which(SamplesSlv$Nyellow == SamplesSlv$Ntot)])
# Apply a minimum samples size of 15 individuals per any site and sex
LHTStaged$CaseID <- paste(LHTStaged$fSite, LHTStaged$fSex)
SamplesSlv2 <- LHTStaged %>% group_by(CaseID) %>% summarise(Ntot = n(), Nmales = sum(fSex == "M"), Nfemales = sum(fSex == "F"),
                                                          Nyellow = sum(bSilver == 0), Nsilver = sum(bSilver == 1))
LHTStaged <- LHTStaged %>% filter(CaseID %in% SamplesSlv2$CaseID[which(SamplesSlv2$Ntot >= 15)]) %>% 
  mutate(fSite = droplevels(fSite))


# Inspect representativeness of remaining samplings
LHTStaged$fRunningID <- as.factor(row.names(LHTStaged))
LHTStaged %>% ggplot(aes(x = Length_cm, y = fRunningID, color = fSex, shape = fStage)) +
  geom_point(alpha = .5) +
  scale_x_continuous(breaks = seq(20, 100, 10)) +
  scale_color_manual(values = c("pink", "skyblue")) +
  scale_shape_manual(values = c(17, 19)) +
  theme_classic() +
  theme(axis.title.y = element_blank(), axis.text.y = element_blank()) +
  facet_wrap(~fSite, scales = "free_y")


# 5.1) Fit a Bernoulli GLMM to the silvering data, without environmental covariates

source('SlvBaseModel.R')


# 5.2) Get coefficients and predictions from the base model

load('Output/SilveringBaseModel.rdata')

# Get coefficients of the fitted model
Fixef_Slv <- fixef(SlvBase)$cond
Ranef_Slv <- ranef(SlvBase)$cond$`fSex:fSite`
Ranef_Slv <- data.frame(fSex = sub(":.*", "", row.names(Ranef_Slv)), fSite = sub("*.:", "", row.names(Ranef_Slv)), 
                         Intercept = Ranef_Slv[,1], Slope = Ranef_Slv[,2])
MaxLength_Slv <- max(LHTStaged$Length_cm)
LengthSeq <- seq(1, MaxLength_Slv)
Preds_Slv <- data.frame()
for (site in levels(LHTStaged$fSite)) {
  for (sex in levels(LHTStaged$fSex)) {
    InterceptSex <- ifelse(sex == "F", 0, Fixef_Slv[3])
    RanCoefs <- Ranef_Slv %>% filter(fSex == sex & fSite == site)
    if (nrow(RanCoefs) > 0) { 
      for (l_i in LengthSeq) {
        eta <- Fixef_Slv[1] + (Fixef_Slv[2] + RanCoefs$Slope) * l_i + InterceptSex + RanCoefs$Intercept  
        p <- plogis(eta) * 100  # Convert to probability in %
        Preds_Slv <- rbind(Preds_Slv, data.frame(fSite = site, fSiteLabel = LHTStaged$fSiteLabel[which(LHTStaged$fSite == site)][1],
                                                 fSex = as.factor(sex), Length_cm = l_i, Probability = p))
      }
    }
  }
}
Preds_Slv <- Preds_Slv %>% mutate(fSite = droplevels(as.factor(fSite)), fSiteLabel = droplevels(fSiteLabel))
L50Df <- Preds_Slv %>% filter(Probability >= 50) %>% group_by(fSiteLabel, fSex) %>% 
  summarise(L50 = as.numeric(first(Length_cm))) %>% mutate(fSiteLabel = droplevels(fSiteLabel))


# 5.3) Plot the base model predictions on the site-level

# Prepare a key for habitat-specific background colors of the strip
StripColors <- LHTStaged %>% group_by(fSiteLabel) %>% reframe(Habitat = Habitat[1]) %>% 
  mutate(Color = case_when(Habitat == "RIV" ~ alpha("palegreen2", .4),
                           Habitat == "RIE" ~ alpha("palegreen4", .4),
                           Habitat == "LGN" ~ alpha("deepskyblue3", .4),
                           Habitat == "LAK" ~ alpha("gold", .4),
                           Habitat == "CNL" ~ alpha("peru", .4)))

LHTStaged <- LHTStaged %>% mutate(fSite = droplevels(fSite), fSiteLabel = droplevels(fSiteLabel))
length(levels(LHTStaged$fSite))
Ind <- 1:20    # indices for the 1st plot
Ind2 <- 21:length(levels(LHTStaged$fSite))    # indices for the 2nd plot

Strip <- strip_themed(background_x = elem_list_rect(fill = StripColors$Color[Ind]))
Strip2 <- strip_themed(background_x = elem_list_rect(fill = StripColors$Color[Ind2]))


# Plot preparation: reverse the order of the factor "sex" so that females are plotted on top of males in the end
LHTStaged$fSex <- factor(LHTStaged$fSex, levels = rev(levels(LHTStaged$fSex)))   
Preds_Slv$fSex <- factor(Preds_Slv$fSex, levels = levels(LHTStaged$fSex))
L50Df$fSex <- factor(L50Df$fSex, levels = levels(LHTStaged$fSex))

# Panel plot by site
SlvP1 <- LHTStaged %>% filter(fSiteLabel %in% levels(fSiteLabel)[Ind]) %>% 
  ggplot(aes(x = Length_cm, color = fSex, fill = fSex)) +
  geom_vline(data = L50Df %>% filter(fSiteLabel %in% levels(fSiteLabel)[Ind]), aes(xintercept = L50, color = fSex), linetype = 2, show.legend = FALSE) +
  geom_point(aes(y = as.numeric(as.character(bSilver))*100), pch = 21, alpha = .4, size = .8, show.legend = FALSE) +
  geom_line(data = Preds_Slv %>% filter(fSiteLabel %in% levels(fSiteLabel)[Ind]), aes(y = Probability), alpha = .8, show.legend = FALSE) +
  scale_x_continuous(name = "Length (cm)", limits = c(0, 100), breaks = seq(20, 100, by = 20)) +
  scale_y_continuous(name = "Proportion silver (%)") +
  scale_color_manual(values = Palette) +
  scale_fill_manual(values = Palette) +
  guides(color = guide_legend(title = "Sex"), fill = guide_legend(title = "Sex")) +
  theme_bw() +
  theme(legend.position = "bottom",
        strip.text = element_text(size = 4.7, margin = margin(.05,0,.05,0, "cm")),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 7),
        legend.key.size = unit(0.4, "cm"),
        panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_line(linewidth = .1),
        panel.grid.major = element_line(linewidth = .3),
        axis.text = element_text(size = 7),
        axis.title = element_text(size = 8)) +
  facet_wrap2(~fSiteLabel, strip = Strip, ncol = 5, nrow = 5)
SlvP1
#ggsave(SlvP1, file = "Output/Sitelevel/Silvering_sitelevel_page1.png", dpi = 600, height = 4.5, width = 6, units = "in")

SlvP2 <- LHTStaged %>% filter(fSiteLabel %in% levels(fSiteLabel)[Ind2]) %>% 
  ggplot(aes(x = Length_cm, color = fSex, fill = fSex)) +
  geom_vline(data = L50Df %>% filter(fSiteLabel %in% levels(fSiteLabel)[Ind2]), aes(xintercept = L50, color = fSex), linetype = 2) +
  geom_point(aes(y = as.numeric(as.character(bSilver))*100), pch = 21, alpha = .4, size = .8) +
  geom_line(data = Preds_Slv %>% filter(fSiteLabel %in% levels(fSiteLabel)[Ind2]), aes(y = Probability), alpha = .8) +
  scale_x_continuous(name = "Length (cm)", limits = c(0, 100), breaks = seq(20, 100, by = 20)) +
  scale_y_continuous(name = "Proportion silver (%)") +
  scale_color_manual(values = Palette) +
  scale_fill_manual(values = Palette) +
  guides(color = guide_legend(title = "Sex"), fill = guide_legend(title = "Sex")) +
  theme_bw() +
  theme(strip.text = element_text(size = 4.7, margin = margin(.05,0,.05,0, "cm")),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 7),
        legend.key.size = unit(0.4, "cm"),
        panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_line(linewidth = .1),
        panel.grid.major = element_line(linewidth = .3),
        axis.text = element_text(size = 7),
        axis.title = element_text(size = 8)) +
  facet_wrap2(~fSiteLabel, strip = Strip2, ncol = 5)
SlvLegend <- cowplot::get_legend(SlvP2)
SlvMain2 <- SlvP2 + theme(legend.position = "none")
#ggsave(SlvMain2, file = "Output/Sitelevel/Silvering_sitelevel_page2.png", dpi = 600, height = 4.5, width = 6, units = "in")


# Population-level plot
(MaxLength_Slv)   # insert maximum observed length in dataset and insert into the predicted length spectrum defined in the following line
SlvX <- seq(0,100, by = 0.1)
PopOutput_Slv <- ggpredict(SlvBase, terms = c("Length_cm[SlvX]", "fSex[all]"))
PopPreds_Slv <- data.frame(fSex = PopOutput_Slv$group, Length_cm = PopOutput_Slv$x, Probability = PopOutput_Slv$predicted*100,
                           Lower = PopOutput_Slv$conf.low*100, Upper = PopOutput_Slv$conf.high*100) %>% arrange(fSex)
L50Pop <- PopPreds_Slv %>% group_by(fSex) %>% summarise(L50 = first(Length_cm[which(Probability >= 50)]),
                                                        L5 = first(Length_cm[which(Probability >= 5)]),
                                                        L95 = first(Length_cm[which(Probability >= 95)]))
L50PopLong <- reshape2::melt(L50Pop, value.name = "Length_cm", variable.name = "Probability") %>% 
  mutate(Probability = as.numeric(as.character(sub("L","",Probability))))

# Plot preparations
L50Pop$fSex <- factor(L50Pop$fSex, levels = c("M", "F"))
PopPreds_Slv$fSex <- factor(PopPreds_Slv$fSex, levels = c("M", "F"))
Preds_Slv$fSex <- factor(Preds_Slv$fSex, levels = c("M", "F"))

Preds_Slv$Group <- paste(Preds_Slv$fSite, Preds_Slv$fSex)

SlvPlotPop <- ggplot() +
    geom_vline(data = L50Pop, aes(xintercept = L50, color = fSex), alpha = .6, linewidth = 1.5, linetype = 2, show.legend = FALSE) +
  #geom_rect(data = L50Pop, aes(xmin = L5, xmax = L95, ymin = 0, ymax = 100, fill = fSex), color = "transparent", alpha = .6, size = .7, show.legend = FALSE) +
  geom_line(data = Preds_Slv, aes(group = Group, x = Length_cm, y = Probability, color = fSex), alpha = .35) +
  geom_ribbon(data = PopPreds_Slv, aes(x = Length_cm, ymin = Lower, ymax = Upper, fill = fSex), alpha = .5) +
  geom_line(data = PopPreds_Slv, aes(x = Length_cm, y = Probability, color = fSex), linewidth = 1.3) +
  scale_x_continuous(name = "Length (cm)", limits = c(0, 100), breaks = seq(20, 100, by = 20)) +
  scale_y_continuous(name = "Proportion silver (%)") +
  scale_color_manual(values = Palette) +
  scale_fill_manual(values = Palette) +
  guides(color = guide_legend(title = "Sex"), fill = guide_legend(title = "Sex")) +
  theme_bw() +
  theme(panel.grid.minor.y = element_blank(),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.position = c(0.9, 0.14),
        legend.background = element_rect(color = "black", linewidth = .2))
SlvPlotPop
#ggsave(SlvPlotPop, file = "Output/Poplevel/Silvering_poplevel.png", dpi = 600, height = 5, width = 6, units = "in")


# Population-level mean length of a silver eel
SlvML <- LHT_HAB2 %>% filter(fStage == "S") %>% group_by(fSite, fSex) %>% 
  summarise(MeanLength = mean(Length_cm))
mean(SlvML$MeanLength[which(SlvML$fSex == "M")])
mean(SlvML$MeanLength[which(SlvML$fSex == "F")])

# 5.4) Identify environmental correlates of silvering

if (Rerun) source('SlvEnvModel.R')


# 5.5) Plot environmental effect on silvering

source('VizEnvEffectsSlv.R')



# 6) Bayesian Multivariate Multilevel model  ------------------------------

source('MultivariateMultilevelModel.R')
