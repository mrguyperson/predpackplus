library("tidyverse")

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
  bind_rows(prediction_list)
}
