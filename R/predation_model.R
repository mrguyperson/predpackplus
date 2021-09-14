library("tidyverse")
library("openxlsx")


##### Data on predation preventation from temperature (mortFishAqPredT) #####
# load the data
# This uses physological measures of predators to get the potential predation effect of T
predTData = read.xlsx(xlsxFile = "./inSALMO Fish Parameters.xlsx",
                      sheet = "mortAqByPredMet",
                      na.strings = "NA") %>% 
  group_by(author, year, journal, species) %>% 
  mutate(unitlessValue = 1-value/max(value)) %>% 
  ungroup() 

# do a logistic fit
predTModel = glm(predTData$unitlessValue ~ predTData$x,
                 family=quasibinomial(logit),
                 data=predTData)

# add in predictions for plotting
predTWFitData = predTData %>% 
  mutate(predict = predict.glm(predTModel, type = "response")) %>% 
  arrange(X)

# solve for inSALMO Parameters 
# convert form m to cm
mortFishAqPredT1 = -(log(1/0.1-1)+predTModel[[1]][1])/predTModel[[1]][2]
mortFishAqPredT9 = -(log(1/0.9-1)+predTModel[[1]][1])/predTModel[[1]][2]

# Plot

ggplot(predTWFitData, aes(x = X)) +
  theme_classic(base_size = 30) +
  labs(y = "Fraction present", x = "T (C)") +
  geom_point(aes(y = unitlessValue)) +
  geom_path(aes(y = predict), color = "red")

B <- (log(0.1/0.9) * 2) / (mortFishAqPredT1 - mortFishAqPredT9)

A <- log(0.1/0.9) - (B * mortFishAqPredT1)


temperature <- 3
S_Temp <- exp(A + (B * temperature)) / (1 + exp(A + (B * temperature)))

length <- 0.1
cover <- 0.1
depth <- 0.1
light <- 0.1
density <- 0.1


survival_increase_list<- c(length, 
                           cover,
                           depth, 
                           light, 
                           density, 
                           S_Temp)

reach_pred_min <- 0.9

s_t <- reach_pred_min + ((1 - reach_pred_min) * (1 - prod(1 - survival_increase_list)))

s_t