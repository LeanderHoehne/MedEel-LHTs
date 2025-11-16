# Likelihood ratio tests to obtain p-values for environmental covariates affecting growth

# Get the parameter and SE estimates
summary(GMEnvObs)

# Calculate p-values from Likelihood Ratio Tests
GMEnvObs2 <- update(GMEnvObs, fixed = list(Linf ~ fSex,  # Removed DistanceGibraltar
                                           k ~ fSex + StageSelectivity),
                    start = list(fixed = c(70, -30, 0.4, 0.8, 0)))
anova(GMEnvObs, GMEnvObs2)  # p DistanceGibraltar = 0.0009

GMEnvObs3 <- update(GMEnvObs, fixed = list(Linf ~ fSex + DistanceGibraltar,  # Removed StageSelectivity
                                           k ~ fSex),
                    start = list(fixed = c(70, -30, 0, 0.4, 0.8)))
anova(GMEnvObs, GMEnvObs3)  # p StageSelectivity = 0.0082

GMEnvObs4 <- update(GMEnvObs, fixed = list(Linf ~ DistanceGibraltar,  # Removed fSex-Linf
                                           k ~ fSex + StageSelectivity),
                    start = list(fixed = c(50, 0, 0.4, 0.8, 0)))
anova(GMEnvObs, GMEnvObs4)  # p fSex-Linf = <0.0001

GMEnvObs5 <- update(GMEnvObs, fixed = list(Linf ~ fSex + DistanceGibraltar,  # Removed fSex-k
                                           k ~ StageSelectivity),
                    start = list(fixed = c(70, -30, 0, 0.7, 0)))
anova(GMEnvObs, GMEnvObs5)  # p fSex-k = <0.0001
