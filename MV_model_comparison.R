library(loo)
library(brms)
library(tidyverse)

load("C:/Users/Utente/OneDrive/PostDoc/Analysis/Life-history traits/MV model runs/MV_v2_Output/Output/MV_v2_full.rdata")
load("C:/Users/Utente/OneDrive/PostDoc/Analysis/Life-history traits/MV model runs/MV_v2_Output/Output/MV_v2_minusDist.rdata")
load("C:/Users/Utente/OneDrive/PostDoc/Analysis/Life-history traits/MV model runs/MV_v2_Output/Output/MV_v2_minusSal.rdata")
load("C:/Users/Utente/OneDrive/PostDoc/Analysis/Life-history traits/MV model runs/MV_v2_Output/Output/MV_v2_minusSurf.rdata")
load("C:/Users/Utente/OneDrive/PostDoc/Analysis/Life-history traits/MV model runs/MV_v2_Output/Output/MV_v2_minusTemp.rdata")
load("C:/Users/Utente/OneDrive/PostDoc/Analysis/Life-history traits/MV model runs/MV_v2_Output/Output/MV_v2_minusChloro.rdata")
load("C:/Users/Utente/OneDrive/PostDoc/Analysis/Life-history traits/MV model runs/MV_v2_Output/Output/MV_v2_minusTP.rdata")

# Model validation
plot(multi_model_full)
summary(multi_model_full)
summary(multi_model_dist)
summary(multi_model_sal)
summary(multi_model_surf)
summary(multi_model_temp)
plot(multi_model_temp)
summary(multi_model_chloro)
summary(multi_model_tp)

# PSIS-LOO
LOO_full <- loo(multi_model_full)
LOO_Dist <- loo(multi_model_dist)
LOO_Sal <- loo(multi_model_sal)
LOO_Surf <- loo(multi_model_surf)
LOO_Temp <- loo(multi_model_temp)
LOO_Chloro <- loo(multi_model_chloro)
LOO_TP <- loo(multi_model_tp)

loo_compare(LOO_full, LOO_Dist, LOO_Sal, LOO_Surf, LOO_Temp, LOO_Chloro, LOO_TP)

# WAIC
MVM_full2 <- add_criterion(multi_model_full, "waic")
MVM_Dist2 <- add_criterion(multi_model_dist, "waic")
MVM_Sal2 <- add_criterion(multi_model_sal, "waic")
MVM_Surf2 <- add_criterion(multi_model_surf, "waic")
MVM_Temp2 <- add_criterion(multi_model_temp, "waic")
MVM_Chloro2 <- add_criterion(multi_model_chloro, "waic")
MVM_TP2 <- add_criterion(multi_model_tp, "waic")

loo_compare(MVM_full2, MVM_Dist2, MVM_Sal2, MVM_Surf2, MVM_Temp2, MVM_Chloro2, MVM_TP2, criterion = "waic")


R2_full <- bayes_R2(MVM_full)
R2_Dist <- bayes_R2(MVM_Dist)
R2_Sal <- bayes_R2(MVM_Sal)
R2_Surf <- bayes_R2(MVM_Surf)
R2_Temp <- bayes_R2(MVM_Temp)

R2_Df <- data.frame(Model = rep(c("Full", "Dist", "Sal", "Surf", "Temp"), each = 4), 
                    Response = rep(rownames(R2_full), 5), 
                    R2 = c(R2_full[,1], R2_Dist[,1], R2_Sal[,1], R2_Surf[,1], R2_Temp[,1]))

(R2_avg_Df <- R2_Df %>% group_by(Model) %>% summarise(Mean_R2 = format(mean(R2),digits = 5)) %>% arrange(Mean_R2))

