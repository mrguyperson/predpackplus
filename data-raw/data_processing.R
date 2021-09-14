library("tidyverse")
library("openxlsx")
library("rstudioapi")
library("fs")

#setwd(dirname(rstudioapi::getSourceEditorContext()$path))

data_folder  <- 'data'
file_xlsx <- 'inSALMO Fish Parameters.xlsx'
file_csv <- 'lmb_stb_combined.csv'

##### Data on predation prevention from temperature (mortFishAqPredT) #####
# load the data
# This uses physiological measures of predators to get the potential predation effect of T
pred_t_data  <- function(){
  read.xlsx(xlsxFile = path(data_folder, file_xlsx),
                      sheet = "mortAqByPredMet",
                      na.strings = "NA") %>% 
  group_by(author, year, journal, species) %>% 
  mutate(unitless_value = 1-value/max(value)) %>% 
  ungroup()  %>%
  dplyr::select(magnitude, variable, unitless_value)
}

##### Data on predation prevention from length (mortFishAqPredL) #####
# load the data and convert all metrics to daily survival

# function to calculate gape-limited prey size for a given predator size (based on species)

prey_conv <- function(a, B, pred_L){
  exp(a + B * log(pred_L))
}

# an angle value used to calculate the reaction distance of predators

angle_calc <- function(length){
  0.0167 * exp(9.14 - 2.4 * log(length) + 0.229 * log(length)^2)
}

# calculates the survival relationship with prey size
# based on the max size of prey that a given predator can consume
# and the size distribution of predators in the Sacramento River
pred_len_data <- function(){
  read.csv(path(data_folder, file_csv)) %>%
    mutate(max_prey_length = prey_conv(0.443, 0.774, length_mm),
           safety = cumulative_proportion - proportion_of_total,
           angle = angle_calc(length_mm),
           reaction_distance = max_prey_length / (2 * tan(angle / 2)),
           max_prey_length = max_prey_length / 10,
           variable = "length") %>%
    dplyr::rename(magnitude = max_prey_length, unitless_value = safety) %>%
    dplyr::select(variable, magnitude, unitless_value)
}


##### Data on predation prevention from depth (mortFishAqPredD) #####
# load the data
# this uses the fractional occurrence of small fish as a proxy for survival
# what is the max possible survival

predation_depth_pre_data <- function(){ 
  read.xlsx(xlsxFile = path(data_folder, file_xlsx),
            sheet = "mortFishByOccurence",
            na.strings = "NA") 
}

pred_depth_data  <- function(maxSurvival=0.9){ 
  predation_depth_pre_data() %>% 
    filter(fishSize_mm =="< 50",
           variable == "Depth") %>% 
    mutate(fraction = cumlitaveFraction - lag(cumlitaveFraction),
           unitless_value = fraction/max(fraction, na.rm = T)*maxSurvival) %>%
    dplyr::rename(magnitude = value) %>%
    dplyr::select(magnitude, variable, unitless_value)
}

##### Data on predation prevention from cover (mortFishAqPredH) #####
# load the data
# this uses the fractional occurrence of small fish as a proxy for survival
# what is the max possible survival
predation_cover_pre_data <- function(){
  read.xlsx(xlsxFile = path(data_folder, file_xlsx),
            sheet = "mortFishByOccurence",
            na.strings = "NA")
}

pred_cover_data <- function(maxSurvival=0.9){
  predation_cover_pre_data() %>% 
    filter(fishSize_mm =="< 50",
           variable == "Dis to Cover") %>% 
    mutate(fraction = cumlitaveFraction - lag(cumlitaveFraction),
           fraction = ifelse(is.na(fraction), cumlitaveFraction, fraction),
           unitless_value = fraction/max(fraction, na.rm = T)*maxSurvival) %>%
    dplyr::rename(magnitude = value) %>%
    dplyr::select(variable, magnitude, unitless_value)
}


# runs all the data functions and combines them into a single dataframe
full_raw_data  <- function(){
  df <- bind_rows(pred_t_data(), pred_len_data(), pred_depth_data(), pred_cover_data())
  df %>% mutate(variable = tolower(variable),
                unitless_value = as.numeric(unitless_value))
}

