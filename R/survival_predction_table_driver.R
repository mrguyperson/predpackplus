library("tidyverse")
library("openxlsx")
library("rstudioapi")
library("broom")
library('fs')
#setwd(dirname(rstudioapi::getSourceEditorContext()$path))

source(path('script','functions','data_processing.R'))
source(path('script','functions','x_values_table.R'))
source(path('script','functions','model_table.R'))
source(path('script','functions','survival_prediction_table.R'))

# driver function that outputs a table of variables affecting predation survival,
# as well as the range of those values and prections

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
