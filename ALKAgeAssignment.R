# Use an age-length key to cope with potential bias through selective subsampling (following Lusk et al., 2021)
LengthAgeExt <- LengthAge[0,]
for (i in unique(LengthAge$fSite)) {
  SiteData <- LengthAge %>% filter(fSite == i)    # select site-specific observations
  for (s in unique(SiteData$fSex)) {
    # Create subsets for both the aged and not aged individuals
    AgeData <- SiteData %>% filter(!is.na(Age) & fSex == s)
    LengthData <- SiteData %>% filter(is.na(Age) & fSex == s)
    LengthAgeImpSex <- SiteData %>% filter(fSex == s)   # keep the sex-specific observations for any site, which gets overwritten in the following loop if there are observations of non-aged eels
    if (nrow(LengthData) > 0) {   # check if age has been determined only on a subsample of individuals
      rm(LengthAgeImpSex)
      StartCat <- min(AgeData$Length_cm)  # lowest length bin with observed ages (gets rounded down to the nearest 5- or 10-cm step below)
      Width <- ifelse(s == "M", 5, 10)    # "width" of the length bins in cm (5cm bins for males; 10cm bins for females)
      StartCat <- RoundTo(StartCat, multiple = Width, FUN = floor)
      AgeData2 <- lencat(~Length_cm, data = AgeData, startcat = StartCat, w = Width)   # add length-bin column to both datasets
      AgeFreqTable <- with(AgeData2,table(LCat, Age))    # age-frequencies by length-bin
      AgePropTable <- prop.table(AgeFreqTable, margin = 1)    # convert to proportions
      MaxLengthCatAged <- as.numeric(row.names(AgePropTable)[nrow(AgePropTable)]) + (Width-.01)   # define the "upper" end of the maximum length-bin with age observations
      LengthData <- LengthData %>% filter(Length_cm >= StartCat & Length_cm <= MaxLengthCatAged)    # length observations falling outside the length-bins with age observations cannot be reliably assigned an age and are therefore omitted
      LengthDataImp <- alkIndivAge(AgePropTable, Age ~ Length_cm, data = LengthData)    # apply the Age Length Key to impute ages to non-aged fish 
      LengthAgeImpSex <- rbind(AgeData, LengthDataImp)
    }
    LengthAgeExt <- rbind(LengthAgeExt, LengthAgeImpSex)
  }
  #cat(paste(i, "; Original df =", nrow(SiteData), "obs. ; Assigned df =", nrow(LengthAgeExt[which(LengthAgeExt$fSite == i),]), "obs.", "\n"))
}
