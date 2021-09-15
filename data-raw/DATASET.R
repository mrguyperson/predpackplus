## code to prepare `DATASET` dataset goes here
pred_t_data  <- function(){
  readr::read_csv('data-raw/mortAqByPredMet.csv') %>%
    dplyr::group_by(author, year, journal, species) %>%
    dplyr::mutate(unitless_value = 1 - value/max(value)) %>%
    dplyr::ungroup()  %>%
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
  readr::read_csv('data-raw/lmb_stb_combined.csv') %>%
    dplyr::mutate(max_prey_length = prey_conv(0.443, 0.774, length_mm),
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

pred_depth_data  <- function(maxSurvival=0.9){
  readr::read_csv('data-raw/mortFishByOccurence.csv') %>%
    dplyr::filter(fishSize_mm =="< 50",
                  variable == "Depth") %>%
    dplyr:: mutate(fraction = cumlitaveFraction - lag(cumlitaveFraction),
                   unitless_value = fraction/max(fraction, na.rm = T)*maxSurvival) %>%
    dplyr::rename(magnitude = value) %>%
    dplyr::select(magnitude, variable, unitless_value)
}

##### Data on predation prevention from cover (mortFishAqPredH) #####
# load the data
# this uses the fractional occurrence of small fish as a proxy for survival
# what is the max possible survival

pred_cover_data <- function(maxSurvival=0.9){
  readr::read_csv('data-raw/mortFishByOccurence.csv') %>%
    dplyr::filter(fishSize_mm =="< 50",
           variable == "Dis to Cover") %>%
    dplyr::mutate(fraction = cumlitaveFraction - lag(cumlitaveFraction),
           fraction = ifelse(is.na(fraction), cumlitaveFraction, fraction),
           unitless_value = fraction/max(fraction, na.rm = T)*maxSurvival) %>%
    dplyr::rename(magnitude = value) %>%
    dplyr::select(variable, magnitude, unitless_value)
}


# runs all the data functions and combines them into a single dataframe
full_raw_data  <- function(){
  df <- dplyr::bind_rows(pred_t_data(),
                         pred_len_data(),
                         pred_depth_data(),
                         pred_cover_data())
  df %>% dplyr::mutate(variable = tolower(variable),
                       unitless_value = as.numeric(unitless_value))
}

table_of_logistic_models <- function(df){
  df %>%
    # nest the table per variable
    tidyr::nest(data = c(magnitude, unitless_value)) %>%

    # add a new column of fitted glm models for each variable
    dplyr::mutate(fit = purrr::map(data, ~ glm(unitless_value ~ magnitude,
                                               family = quasibinomial(logit),
                                               data = .)))
}

survival_prediction_table <- function(df, model_table){
  df %>%
    # nest the table per variable
    tidyr::nest(magnitude = -variable) %>%

    # join the x values with the table that has the glm models
    dplyr::inner_join(model_table, by = 'variable') %>%

    # preditct survival values based on the x values and glms
    dplyr::mutate(survival = purrr::map2(fit,
                                         magnitude,
                                         predict.glm,
                                         type='response')) %>%

    # unnest the x values and survival values
    tidyr::unnest(c(magnitude, survival)) %>%

    # save the variable name, x value, and survival prediction
    dplyr::select(variable, magnitude, survival)
}

# create a table of the variables affecting predation; the min and max values you expect to appear in the environment; and
# the increment by which the values should change

range_of_params <- function(){
  data.frame(variable = c('temp', 'length', 'dis to cover', 'depth'),
             min_x = c(0.0, 0.0, 0.0, 0.0),
             max_x = c(30, 100, 3, 2),
             interval = c(0.1, 0.01, 0.01, 0.01)
  )
}

# create a table of values for each variable over the range and increment specified above
# if someone knows how to do this via an apply function, let me know

x_value_df <- function(df){

  # create an empty list to store results from the loop
  prediction_list = list()

  # loop over each row to create a unique entry for each incremented value of a given variable; e.g., 5.3 cm, 5.4 cm, etc.
  for (i in 1:nrow(df)){
    new_df <- data.frame(variable = df[i,1], magnitude = seq(df[i,2], df[i,3], by = df[i,4]))

    # store each new data frame in the list
    prediction_list[[i]] <- new_df
  }

  # combine the list of dataframes into one large frame
  dplyr::bind_rows(prediction_list)
}

# outputs a table of predicted survival values for each variable affecting predation
# df is intended to work with functions from prediction_tables.R
# model_table is intended to work with model_table.R

survival_prediction_table <- function(df, model_table){
  df %>%
    # nest the table per variable
    tidyr::nest(magnitude = -variable) %>%

    # join the x values with the table that has the glm models
    dplyr::inner_join(model_table, by = 'variable') %>%

    # preditct survival values based on the x values and glms
    dplyr::mutate(survival = purrr::map2(fit,
                                         magnitude,
                                         predict.glm,
                                         type='response')) %>%

    # unnest the x values and survival values
    tidyr::unnest(c(magnitude, survival)) %>%

    # save the variable name, x value, and survival prediction
    dplyr::select(variable, magnitude, survival)
}

predation_survival_driver_func <- function(){
  # create a dataframe of data from the literature
  raw_data <- full_raw_data()

  # create a table of glm's fitted to the literature data
  model_table <- table_of_logistic_models(raw_data)

  # create a table of incremented x values for each variable
  x_values <- x_value_df(range_of_params())

  # make predictions for the table of x values
  survival_table <- survival_prediction_table(x_values, model_table)
}

DATASET <- predation_survival_driver_func()

usethis::use_data(DATASET, internal = TRUE)
