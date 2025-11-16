# Fit a Multivariate Multilevel model in a Bayesian framework to the LHT data using the brms package

library(brms)   # requires Stan to be installed prior. See: https://learnb4ss.github.io/learnB4SS/articles/install-brms.html
library(loo)

# Load the datasets used for the isolated trait-specific analyses
load("LW_env_data.Rdata")
load("Sex_ratio_env_data.rdata")
load("Slv_env_data.rdata")
load("Growth_env_data.Rdata")

# Subset only eels that were included in all of the trait-specific analyses (i.e., those that have complete observations for all biometrics and only sites with sufficient sample sizes)
IDCounts <- table(c(LengthAgeExt$RunningID, LWEnvData$RunningID, SRData$RunningID, LHTStaged$RunningID))
IDsMulti <- names(IDCounts)[which(IDCounts == 4)]
LHTMulti <- LengthAgeExt %>% filter(RunningID %in% IDsMulti) %>% drop_na(Age, fStage, fSex) %>% 
  mutate(bMale = as.factor(ifelse(fSex == "M", 1, 0)), bSilver = as.factor(ifelse(fStage == "S", 1, 0)),
         L0 = 6.66, fSex = droplevels(fSex), logLength_cm = log(Length_cm), logWeight_g = log(Weight_g))   


# Define the model formulae for the multivariate model, analogous to the previous trait-specific frequentist models
silvering_formula <- bf(bSilver ~ Length_cm + fSex + DistanceGibraltar + Salinity + Temperature + Chlorophyll + TP + logSurfaceArea + (1 + Length_cm|fSite/fSex),
                        family = bernoulli())
length_weight_formula <- bf(logWeight_g ~ logLength_cm + fStage + Chlorophyll + TP + DistanceGibraltar + Temperature + Salinity + logSurfaceArea + (1 + logLength_cm|fSite),
                            family = gaussian())
sexratio_formula <- bf(bMale ~ DistanceGibraltar + Salinity + Temperature + Chlorophyll + TP + logSurfaceArea + (1|fSite),
                       family = bernoulli())
growth_formula <- bf(
  Length_cm ~ 6.66 + (Linf - 6.66) * (1 - exp(-k * Age)),
  Linf ~ fSex + logSurfaceArea + Chlorophyll + TP + DistanceGibraltar + (1 + fSex | fSite), 
  k ~ fSex + Temperature + Salinity + (1 + fSex | fSite), 
  nl = TRUE
)

# Specify priors
priors <- c(
  # GROWTH 
  prior(gamma(67, 1), class = b, nlpar = Linf, coef = Intercept, resp = Lengthcm),
  prior(normal(-27, 5), class = b, nlpar = Linf, coef = fSexM, resp = Lengthcm), 
  prior(gamma(0.5, 1.5), class = b, nlpar = k, coef = Intercept, resp = Lengthcm),   # k: baseline for females
  prior(normal(0.7, 0.3), class = b, nlpar = k, coef = fSexM, resp = Lengthcm),    # k: difference for males (male - female)
  prior(normal(0, 1), class = b, nlpar = Linf, coef = Chlorophyll, resp = Lengthcm),
  prior(normal(0, 1), class = b, nlpar = Linf, coef = TP, resp = Lengthcm),
  prior(normal(0, 1), class = b, nlpar = Linf, coef = DistanceGibraltar, resp = Lengthcm),
  prior(normal(0, 1), class = b, nlpar = Linf, coef = logSurfaceArea, resp = Lengthcm),
  prior(normal(0, 1), class = b, nlpar = k, coef = Temperature, resp = Lengthcm),
  prior(normal(0, 1), class = b, nlpar = k, coef = Salinity, resp = Lengthcm),
  prior(normal(10, 5), class = sigma, resp = Lengthcm),
  # SILVERING
  prior(normal(-14, 5), class = Intercept, resp = bSilver),
  prior(normal(6, 3), class = b, coef = fSexM, resp = bSilver),
  prior(normal(0, 1), class = b, coef = DistanceGibraltar, resp = bSilver),
  prior(normal(0, 1), class = b, coef = Length_cm, resp = bSilver),
  prior(normal(0, 1), class = b, coef = Chlorophyll, resp = bSilver),
  prior(normal(0, 1), class = b, coef = logSurfaceArea, resp = bSilver),
  prior(normal(0, 1), class = b, coef = TP, resp = bSilver),
  prior(normal(0, 1), class = b, coef = Salinity, resp = bSilver),
  prior(normal(0, 1), class = b, coef = Temperature, resp = bSilver),
  # LENGTH-WEIGHT
  prior(normal(-7, 3), class = Intercept, resp = logWeightg),
  prior(normal(0, 1), class = b, coef = DistanceGibraltar, resp = logWeightg),
  prior(normal(3.17, 1), class = b, coef = logLength_cm, resp = logWeightg),
  prior(normal(0, 1), class = b, coef = Chlorophyll, resp = logWeightg),
  prior(normal(0, 1), class = b, coef = logSurfaceArea, resp = logWeightg),
  prior(normal(0, 1), class = b, coef = TP, resp = logWeightg),
  prior(normal(0, 1), class = b, coef = Salinity, resp = logWeightg),
  prior(normal(0, 1), class = b, coef = Temperature, resp = logWeightg),
  prior(normal(0, 1), class = b, coef = fStageY, resp = logWeightg),
  # SEX RATIO
  prior(normal(-.7, .2), class = Intercept, resp = bMale),
  prior(normal(0, 1), class = b, coef = DistanceGibraltar, resp = bMale),
  prior(normal(0, 1), class = b, coef = Chlorophyll, resp = bMale),
  prior(normal(0, 1), class = b, coef = logSurfaceArea, resp = bMale),
  prior(normal(0, 1), class = b, coef = TP, resp = bMale),
  prior(normal(0, 1), class = b, coef = Salinity, resp = bMale),
  prior(normal(0, 1), class = b, coef = Temperature, resp = bMale)
)


# Initial values for some parameters, to aid convergence
inits <- list(
    b_Lengthcm_Linf = c(rnorm(1,70,10),rnorm(1,-30,5),rnorm(4,0,1)),
    b_Lengthcm_k = c(runif(1,0.2,0.6),runif(1,0.4,1),rnorm(2,0,1)),             
    sigma_Lengthcm = runif(1,5,15)
)
inits2 <- list(
  b_Lengthcm_Linf = c(rnorm(1,70,10), rnorm(1,-30,5),rnorm(4,0,1)),
  b_Lengthcm_k = c(runif(1,0.2,0.6),runif(1,0.4,1),rnorm(2,0,1)),
  sigma_Lengthcm = runif(1,5,15)
)
inits3 <- list(
  b_Lengthcm_Linf = c(rnorm(1,70,10), rnorm(1,-30,5),rnorm(4,0,1)),
  b_Lengthcm_k = c(runif(1,0.2,0.6),runif(1,0.4,1),rnorm(2,0,1)),
  sigma_Lengthcm = runif(1,5,15)
)

init_list <- list(inits, inits2, inits3)

# Configuration for the No-U-Turn-Sampler
control <- list(
  adapt_engaged = TRUE,
  adapt_delta = 0.95, 
  stepsize = 0.05, 
  max_treedepth = 15
)


# GET MCMC SAMPLES
# The following call runs the full model including all environmental covariates. 
# This step would need to be repeated excluding one env. covariate at a time from the full model, saving the results of each run to the output folder.

multi_model_full <- brm(growth_formula + silvering_formula + length_weight_formula + sexratio_formula,
                   data = LHTMulti, chains = 3, cores = 3, iter = 40000, prior = priors, init = init_list, 
                   control = control)


# MODEL VALIDATION
# For each model, inspect the traceplots for all parameters to ensure convergence
plot(multi_model_full)
# In addition, check the summary to ensure all R.hat values (a formal measure of convergence) are close to 1 (<1.1)
summary(multi_model_full)
# Only proceed to model comparison, if all models have converged well! Otherwise, update or re-run the respective model with the same no. of iterations.


# MODEL COMPARISON
# Once the full model and all leave-one-predictor-out model runs have been completed and saved, we can compare the models using PSIS-LOO

# Calculate the PSIS-LOO for any model
LOO_full <- loo(multi_model_full)   # exemplary for the full model. Do that also for all leave-one-predictor-out model's posterior samples
# then you can compare all models using loo_compare(LOO_full, LOO_minusSalinity, LOO_minusTemperature, ...)
