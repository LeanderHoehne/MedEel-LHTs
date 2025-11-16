# Load the output from the selected model
load("./Output/GMBase.rdata")
Coefs <- coef(GMBase)
Coefs <- data.frame(Site = rownames(Coefs), 
                    Linf_Intercept = Coefs[,1],
                    Linf_M = Coefs[,2],
                    k_Intercept = Coefs[,3],
                    k_M = Coefs[,4])

# Population-level parameter estimates and confidence intervals
MaxAge <- max(LengthAgeObs$Age); MaxAge_M <- max(LengthAgeObs$Age[which(LengthAgeObs$fSex == "M")])
SeqF <- seq(0, MaxAge, by = 0.2); SeqM <- seq(0, MaxAge_M, by = 0.2)
(Intervals <- intervals(GMBase, which = "fixed")$fixed)
# Population-level predictions
PopPreds <- data.frame(fSex = as.factor(c(rep("F", length(SeqF)), rep("M", length(SeqM)))), Age = c(SeqF, SeqM))
PopPreds <- data.frame(PopPreds, Predicted = NA, Lower = NA, Upper = NA)
PopPreds$Predicted[which(PopPreds$fSex == "F")] <- VBGFOrig(SeqF, Intervals[1,2], Intervals[3,2], L0)
PopPreds$Predicted[which(PopPreds$fSex == "M")] <- VBGFOrig(SeqM, Intervals[1,2] + Intervals[2,2], Intervals[3,2] + Intervals[4,2], L0)
PopPreds$Lower[which(PopPreds$fSex == "F")] <- VBGFOrig(SeqF, Intervals[1,1], Intervals[3,1], L0)
PopPreds$Lower[which(PopPreds$fSex == "M")] <- VBGFOrig(SeqM, Intervals[1,2] + Intervals[2,1], Intervals[3,2] + Intervals[4,1], L0)
PopPreds$Upper[which(PopPreds$fSex == "F")] <- VBGFOrig(SeqF, Intervals[1,3], Intervals[3,3], L0)
PopPreds$Upper[which(PopPreds$fSex == "M")] <- VBGFOrig(SeqM, Intervals[1,2] + Intervals[2,3], Intervals[3,2] + Intervals[4,3], L0)

# Site-specific predictions
SitePreds <- data.frame(fSite = NA, fSiteLabel = NA, Habitat = NA, fSex = NA, Age = NA, Predicted = NA)
for (i in unique(Coefs$Site)) {
  SiteCoefs <- Coefs[Coefs$Site == i,]
  SiteObs <- LengthAgeObs %>% filter(fSite == i)
  for (s in unique(SiteObs$fSex)) {
    Seq <- seq(0, max(SiteObs$Age[SiteObs$fSex == s]), by = 0.2)
    Linf <- ifelse(s == "F", SiteCoefs$Linf_Intercept, SiteCoefs$Linf_Intercept + SiteCoefs$Linf_M)
    k <- ifelse(s == "F", SiteCoefs$k_Intercept, SiteCoefs$k_Intercept + SiteCoefs$k_M)
    Predicted <- VBGFOrig(Seq, Linf, k, L0)
    Appendix <- data.frame(fSite = i, fSiteLabel = SiteObs$fSiteLabel[1], Habitat = SiteObs$Habitat[1], fSex = s, Age = Seq, Predicted)
    SitePreds <- rbind(SitePreds, Appendix)
  }
}
SitePreds <- SitePreds[-1,]
SitePreds <- SitePreds %>% mutate(fSite = as.factor(fSite), 
                                  fSex = as.factor(fSex))

Coefs$SiteAcronym <- LHT$SiteAcronym[match(Coefs$Site, LHT$fSite)]

save(Coefs, file = "./Output/BaseModelCoefs_site.rdata")
save(PopPreds, file = "./Output/BaseModelPreds_pop.rdata")
save(SitePreds, file = "./Output/BaseModelPreds_site.rdata")