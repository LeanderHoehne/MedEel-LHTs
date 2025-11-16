# Combine glass and yellow eels into one stage class
LHT_HAB3$fStage <- as.character(LHT_HAB3$fStage)
LHT_HAB3$fStage[which(is.na(LHT_HAB3$fStage))] <- "ND"
LHT_HAB3$fStage <- as.factor(LHT_HAB3$fStage)
LHT_HAB3 <- LHT_HAB3 %>% mutate(fStage = fct_recode(fStage, Yellow = "Y", Silver = "S"))

# Drop sites with insufficient length-weight observations
LHT_HAB3$CaseID <- paste(LHT_HAB3$fSite, LHT_HAB3$fStage)
LWObsStage <- LHT_HAB3 %>% group_by(CaseID) %>% tally()   # Calculate sample sizes
# Delete cases with < 15 observations for any stage
LHT_HAB3 <- LHT_HAB3 %>% filter(CaseID %in% LWObsStage$CaseID[which(LWObsStage$n >= 15)]) %>% 
  mutate(fSite = droplevels(fSite)) %>% 
  mutate(fSiteLabel = as.factor(paste0(fSite, " (", CountryCode, ")")))

# For sites that have sufficient data by stage, and unstaged individuals in addition, leave out the unstaged eels from analysis
CasesBySite <- LHT_HAB3 %>% group_by(fSite) %>% summarise(N = length(unique(CaseID)))
LHT_HAB3 <- LHT_HAB3[-which(LHT_HAB3$fSite %in% CasesBySite$fSite[which(CasesBySite$N == 3)] & LHT_HAB3$fStage == "ND"),]

# Test whether lifestage should be accounted for
LWMStage <- glmmTMB(logWeight_g ~ logLength_cm * fStage + (1 + logLength_cm|fSite), data = LHT_HAB3)
summary(LWMStage)
# --> significant effect of stage on the LW-relationship. Interaction term is significant, so random intercepts AND slopes are justified.
# In this model the effect of length on the weight of an eel is allowed to vary by site, but no so the effect of stage. We assume that this is fixed across sites.

# Refit the stage-specific model with REML before coefficients are obtained
LWMStage <- update(LWMStage, REML = TRUE)

# Extract the random and fixed effects' coefficients
RanefsLWStage <- ranef(LWMStage)$cond$fSite; colnames(RanefsLWStage) <- c("Intercept", "Slope")
RanefsLWStage <- as.data.frame(RanefsLWStage) %>% mutate(fSite = row.names(RanefsLWStage))

CoefsLW <- expand.grid(fSite = row.names(RanefsLWStage), fStage = levels(LHT_HAB3$fStage)) %>% arrange(fSite)
CoefsLW$a <- NA
CoefsLW$b <- NA
FixInt <- fixef(LWMStage)$cond["(Intercept)"]
FixLength <- fixef(LWMStage)$cond["logLength_cm"]
FixStages <- c(0, fixef(LWMStage)$cond["fStageSilver"], fixef(LWMStage)$cond["fStageYellow"])
FixInteraction <- c(0, fixef(LWMStage)$cond["logLength_cm:fStageSilver"], fixef(LWMStage)$cond["logLength_cm:fStageYellow"])

# Calculate the a and b coefficients of the lenght-weight relationship for any site-stage combination
for (i in RanefsLWStage$fSite) {
  for (s in 1:length(levels(LHT_HAB3$fStage))) {
    RanInt <- RanefsLWStage$Intercept[which(RanefsLWStage$fSite == i)]
    CoefsLW$a[which(CoefsLW$fSite == i & CoefsLW$fStage == levels(LHT_HAB3$fStage)[s])] <- exp(FixInt + FixStages[s] + RanInt)
    RanSlope <- RanefsLWStage$Slope[which(RanefsLWStage$fSite == i)]
    CoefsLW$b[which(CoefsLW$fSite == i & CoefsLW$fStage == levels(LHT_HAB3$fStage)[s])] <- FixLength + FixInteraction[s] + RanSlope
  }
}

# Subset only site-stage combinations actually occurring in the data
Cases <- LHT_HAB3 %>% group_by(fSite) %>% reframe(Stages = unique(fStage)) %>% mutate(ID = paste(fSite, Stages))
CoefsLW$ID <- paste(CoefsLW$fSite, CoefsLW$fStage)
CoefsLW <- CoefsLW %>% filter(ID %in% Cases$ID)
CoefsLW$fSiteLabel <- LHT_HAB3$fSiteLabel[match(CoefsLW$fSite, LHT_HAB3$fSite)]

# Get predicted weight values using the previously calulcated a and b coefficients of any site and stage
SitePredsLWStage <- data.frame(fSite = as.factor(NA), fSiteLabel = as.factor(NA), fStage = as.factor(NA), ID = NA, Length_cm = NA, Weight_g = NA)
for (i in 1:nrow(CoefsLW)) {
  CaseData <- CoefsLW[i,]
  LengthRange <- range(LHT_HAB3$Length_cm[which(LHT_HAB3$fSite == CaseData$fSite & LHT_HAB3$fStage == CaseData$fStage)])
  LengthSeq <- seq(LengthRange[1], LengthRange[2], by = 1)
  Preds <- CaseData$a * LengthSeq^CaseData$b
  Appendix <- data.frame(fSite = CaseData$fSite, fSiteLabel = LHT_HAB3$fSiteLabel[match(CaseData$fSite, LHT_HAB3$fSite)], fStage = CaseData$fStage, ID = CaseData$ID, Length_cm = LengthSeq, Weight_g = Preds)
  SitePredsLWStage <- rbind(SitePredsLWStage, Appendix)
}
SitePredsLWStage <- SitePredsLWStage[-1,] %>% filter(ID %in% SitePredsLWStage$ID)
SitePredsLWStage$fSite <- as.factor(SitePredsLWStage$fSite)
SitePredsLWStage$fSiteLabel <- as.factor(SitePredsLWStage$fSiteLabel)


# Alternatively, use this script to obtain site- and stage-specific predictions only for the range of observed values in any site
# source('LWStageModelPredsPopPlot.R')


# Get predicted values on the population level (extrapolating beyond observed lengths)
XLims_Yellow <- range(LHT_HAB3$logLength_cm[which(LHT_HAB3$fStage == "Yellow")])
XLims_Silver <- c(min(LHT_HAB3$logLength_cm[which(LHT_HAB3$fStage == "Silver")]), log(105)) 
PopPredsLWStage <- expand.grid(fStage = levels(LHT_HAB3$fStage), logLength_cm = log(1:MaxLength)) %>% 
  arrange(fStage) %>% 
  filter(!(fStage == "Silver" & (logLength_cm < XLims_Silver[1] | logLength_cm > XLims_Silver[2]))) %>% 
  filter(!(fStage == "Yellow" & (logLength_cm < XLims_Yellow[1] | logLength_cm > XLims_Yellow[2])))
XpLWPop <- model.matrix(~ logLength_cm * fStage, data = PopPredsLWStage)
PopPredsLWStage$Predicted <- XpLWPop %*% fixef(LWMStage)$cond    # calculate predicted values
PopPredsLWStage$SE <- sqrt(diag(XpLWPop %*% vcov(LWMStage)$cond %*% t(XpLWPop)))    # calculate standard errors
PopPredsLWStage$Upper <- PopPredsLWStage$Predicted + 1.96 * PopPredsLWStage$SE
PopPredsLWStage$Lower <- PopPredsLWStage$Predicted - 1.96 * PopPredsLWStage$SE
PopPredsLWStage <- PopPredsLWStage %>% filter(fStage != "ND") %>% mutate(fStage = droplevels(fStage))

save(CoefsLW, file = "./Output/LWCoefsStage.rdata")
save(PopPredsLWStage, file = "./Output/LW_preds_poplevel_stage.rdata")
save(SitePredsLWStage, file = "./Output/LW_preds_sitelevel_stage.rdata")


# Get population level results for the L-W curve coefficients
FixefLW <- fixef(LWMStage)$cond
aSilver <- exp(FixefLW["(Intercept)"] + FixefLW["fStageSilver"])
aYellow <- exp(FixefLW["(Intercept)"] + FixefLW["fStageYellow"])
bSilver <- FixefLW["logLength_cm"] + FixefLW["logLength_cm:fStageSilver"]
bYellow <- FixefLW["logLength_cm"] + FixefLW["logLength_cm:fStageYellow"]
