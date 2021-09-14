library("rstudioapi")
library('tidyverse')
library('future.apply')
library('fs')
source(path('script','functions','survival_predction_table_driver.R'))
source(path('script','functions','generate_preds.R'))

plan(multisession)

# set-up values
# pred_pos <- get_pred_postitions()
# stream_grid_frame <- create_stream_raster_frame(pred_pos)
# enc_probs <- calc_enc_probs(stream_grid_frame)
default_base_pred_success <- 0.8
default_reach_pred_min <- 1 - default_base_pred_success
# survival_table <- predation_survival_driver_func()
# num_trans <- calculate_num_trans_traversed()

calculate_num_trans_traversed <- function(transect_length=default_transect_length, 
                                          n_transects=default_n_transects, 
                                          grid_size=default_grid_size){
  
  round(transect_length * n_transects / grid_size)

}

get_enc_prob_vector <- function(enc_prob_dataframe){
  enc_prob_dataframe %>% pull(enc_prob)
}

randomly_sample_encounter_probs <- function(enc_probs, num_trans_traversed){
  sample(enc_probs, num_trans_traversed)
}

calcule_risk_modifier <- function(survival_table, selected_variable, variable_magnitude){
  survival_table %>%
    filter(variable == as.character(selected_variable) & magnitude == as.character(variable_magnitude)) %>%
    pull(survival)
}

set_up_dataframe_of_prey_fish_lengths <- function(number_of_fish, mean = 12, sd = (1.7/14) * mean, precision = 1){
  vals <- round(rnorm(2*number_of_fish, mean, sd), precision)
  vals <- vals[vals>0][1:number_of_fish]
  data.frame(fish = (seq(number_of_fish)), length = vals)
}

calculate_encounter_prob_based_on_length <- function(fish_length_boost, enc_prob){
  encounter_check_values <- (1- fish_length_boost) * enc_prob
}

check_if_encounter_occurs <- function(encounter_prob_based_on_length){
  encounter <- sample(c(0, 1), size = 1, prob = c(1 - encounter_prob_based_on_length, encounter_prob_based_on_length))
}

evaluate_encounter <- function(pred_success = default_base_pred_success, pred_min = default_reach_pred_min){
  survival <- sample(c(0,1), 1, prob = c(pred_success, pred_min))
  if(survival == 1.0){
    outcome <- 1.0
  } else {
    outcome <- 0
  }
  
}

encounter_simulator <- function(encounter_occurrence){
  if(encounter_occurrence == 1){
    outcome <- evaluate_encounter()
  } else {
    outcome <- 1
  }
  return(outcome)
}

encounter_frame <- function(fish_number, enc_prob_vector, num_trans_traversed){
  data.frame(fish = fish_number, enc_prob = randomly_sample_encounter_probs(enc_prob_vector, num_trans_traversed))
}

add_length_survival_boost <- function(fish_frame, surv_table){
  fish_frame %>%
    mutate(survival_boost = lapply(length, FUN = calcule_risk_modifier, survival_table = surv_table, selected_variable = 'length'))
}

create_frame_of_all_cells_traversed_per_fish <- function(fish_vector, enc_probs, num_trans){
  df <- bind_rows(lapply(fish_vector, 
                         FUN = encounter_frame, 
                         enc_prob_vector = enc_probs, 
                         num_trans_traversed = num_trans))
}

combine_encounter_frame_and_fish_frame <- function(fish_frame, enc_frame){
  inner_join(fish_frame, enc_frame, by = 'fish')
}

add_encounter_simulation_columns <- function(combined_frame){
  combined_frame%>%
    mutate(survival_boost = as.numeric(survival_boost),
           enc_prob = as.numeric(enc_prob),
           modified_enc = calculate_encounter_prob_based_on_length(survival_boost, enc_prob),
           encounter = lapply(modified_enc, check_if_encounter_occurs),
           alive = as.numeric(lapply(encounter, encounter_simulator))) 
}

evaluate_final_status_of_fish <- function(outcome_total, num_trans_traversed){
  if(outcome_total < num_trans_traversed){
    final_status <- 0
  } else {
    final_status <- 1
  }
  return(as.numeric(final_status))
}

sum_encounter_outcomes <- function(df){
  df %>%
    group_by(fish) %>%
    summarize(total_alive_outcomes = sum(as.numeric(alive)))
}

calculate_number_of_survivors <- function(fish_summary_frame, num_trans_traversed){
  fish_summary_frame %>%
    mutate(final_status = lapply(total_alive_outcomes, FUN = evaluate_final_status_of_fish, num_trans_traversed)) %>%
    unnest(final_status) %>%
    summarize(num_surviving = sum(final_status))
}

calculat_proportion_of_survivors <- function(number_of_survivors, number_of_fish){
  proportion_of_survivors <- number_of_survivors[[1]] / number_of_fish
}

survival_simulation_driver <- function(number_of_fish, mean_fish_length, total_transects){
  num_trans <- calculate_num_trans_traversed(n_transects = total_transects)
  pred_pos <- get_pred_postitions(n_transects = total_transects)
  stream_grid_frame <- create_stream_raster_frame(pred_pos, n_transects = total_transects)
  enc_prob_frame <- calc_enc_probs(stream_grid_frame)
  enc_probs <- get_enc_prob_vector(enc_prob_frame)
  surv_table <- predation_survival_driver_func()
  fish_frame <- set_up_dataframe_of_prey_fish_lengths(number_of_fish, mean = mean_fish_length)
  fish_frame_with_survival_from_length <- add_length_survival_boost(fish_frame,surv_table)
  cells_traversed <- create_frame_of_all_cells_traversed_per_fish(fish_frame_with_survival_from_length$fish, enc_probs, num_trans)
  joined_fish_cell_traversed <- combine_encounter_frame_and_fish_frame(fish_frame_with_survival_from_length, cells_traversed)
  joined_frame_with_sim_cols <- add_encounter_simulation_columns(joined_fish_cell_traversed)
  fish_summary_frame <- sum_encounter_outcomes(joined_frame_with_sim_cols)
  number_of_survivors <- calculate_number_of_survivors(fish_summary_frame, num_trans)
  proprtion_of_survivors <- calculat_proportion_of_survivors(number_of_survivors, number_of_fish)
}

#### using future apply
future_apply_calculate_number_of_survivors <- function(fish_summary_frame, num_trans_traversed){
  fish_summary_frame %>%
    mutate(final_status = future_lapply(total_alive_outcomes, FUN = evaluate_final_status_of_fish, num_trans_traversed)) %>%
    unnest(final_status) %>%
    summarize(num_surviving = sum(final_status))
}


future_apply_add_encounter_simulation_columns <- function(combined_frame){
  combined_frame%>%
    mutate(survival_boost = as.numeric(survival_boost),
           enc_prob = as.numeric(enc_prob),
           modified_enc = calculate_encounter_prob_based_on_length(survival_boost, enc_prob),
           encounter = future_lapply(modified_enc, check_if_encounter_occurs, future.seed=TRUE),
           alive = as.numeric(lapply(encounter, encounter_simulator))) 
}


future_apply_add_length_survival_boost <- function(fish_frame, surv_table){
  fish_frame %>%
    mutate(survival_boost = future_lapply(length, FUN = calcule_risk_modifier, survival_table = surv_table, selected_variable = 'length',future.seed=TRUE))
}

future_apply_create_frame_of_all_cells_traversed_per_fish <- function(fish_vector, enc_probs, num_trans){
  df <- bind_rows(future_lapply(fish_vector, 
                         FUN = encounter_frame, 
                         enc_prob_vector = enc_probs, 
                         num_trans_traversed = num_trans,
                         future.seed=TRUE))
}



future_apply_survival_simulation_driver <- function(number_of_fish, mean_fish_length, total_transects){
  num_trans <- calculate_num_trans_traversed(n_transects = total_transects)
  pred_pos <- get_pred_postitions(n_transects = total_transects)
  stream_grid_frame <- create_stream_raster_frame(pred_pos, n_transects = total_transects)
  enc_prob_frame <- calc_enc_probs(stream_grid_frame)
  enc_probs <- get_enc_prob_vector(enc_prob_frame)
  surv_table <- predation_survival_driver_func()
  fish_frame <- set_up_dataframe_of_prey_fish_lengths(number_of_fish, mean = mean_fish_length)
  fish_frame_with_survival_from_length <- future_apply_add_length_survival_boost(fish_frame,surv_table)
  cells_traversed <- future_apply_create_frame_of_all_cells_traversed_per_fish(fish_frame_with_survival_from_length$fish, enc_probs, num_trans)
  joined_fish_cell_traversed <- combine_encounter_frame_and_fish_frame(fish_frame_with_survival_from_length, cells_traversed)
  joined_frame_with_sim_cols <- future_apply_add_encounter_simulation_columns(joined_fish_cell_traversed)
  fish_summary_frame <- sum_encounter_outcomes(joined_frame_with_sim_cols)
  number_of_survivors <- future_apply_calculate_number_of_survivors(fish_summary_frame, num_trans)
  proprtion_of_survivors <- calculat_proportion_of_survivors(number_of_survivors, number_of_fish)
}


#### using future map
future_map_calculate_number_of_survivors <- function(fish_summary_frame, num_trans_traversed){
  fish_summary_frame %>%
    mutate(final_status = future_map(total_alive_outcomes, evaluate_final_status_of_fish, num_trans_traversed)) %>%
    unnest(final_status) %>%
    summarize(num_surviving = sum(final_status))
}


future_map_add_encounter_simulation_columns <- function(combined_frame){
  combined_frame%>%
    mutate(survival_boost = as.numeric(survival_boost),
           enc_prob = as.numeric(enc_prob),
           modified_enc = calculate_encounter_prob_based_on_length(survival_boost, enc_prob),
           encounter = future_map(modified_enc, check_if_encounter_occurs, .options=furrr_options(seed=TRUE)),
           alive = as.numeric(lapply(encounter, encounter_simulator))) 
}


future_map_add_length_survival_boost <- function(fish_frame, surv_table){
  fish_frame %>%
    mutate(survival_boost = future_lapply(length, calcule_risk_modifier, survival_table = surv_table, selected_variable = 'length',future.seed=TRUE))
}

future_map_create_frame_of_all_cells_traversed_per_fish <- function(fish_vector, enc_probs, num_trans){
  df <- bind_rows(future_map(fish_vector, 
                                encounter_frame, 
                                enc_prob_vector = enc_probs, 
                                num_trans_traversed = num_trans,
                                .options=furrr_options(seed=TRUE)))
}



future_map_survival_simulation_driver <- function(number_of_fish, mean_fish_length, total_transects){
  num_trans <- calculate_num_trans_traversed(n_transects = total_transects)
  pred_pos <- get_pred_postitions(n_transects = total_transects)
  stream_grid_frame <- create_stream_raster_frame(pred_pos, n_transects = total_transects)
  enc_prob_frame <- calc_enc_probs(stream_grid_frame)
  enc_probs <- get_enc_prob_vector(enc_prob_frame)
  surv_table <- predation_survival_driver_func()
  fish_frame <- set_up_dataframe_of_prey_fish_lengths(number_of_fish, mean = mean_fish_length)
  fish_frame_with_survival_from_length <- future_map_add_length_survival_boost(fish_frame,surv_table)
  cells_traversed <- future_map_create_frame_of_all_cells_traversed_per_fish(fish_frame_with_survival_from_length$fish, enc_probs, num_trans)
  joined_fish_cell_traversed <- combine_encounter_frame_and_fish_frame(fish_frame_with_survival_from_length, cells_traversed)
  joined_frame_with_sim_cols <- future_map_add_encounter_simulation_columns(joined_fish_cell_traversed)
  fish_summary_frame <- sum_encounter_outcomes(joined_frame_with_sim_cols)
  number_of_survivors <- future_map_calculate_number_of_survivors(fish_summary_frame, num_trans)
  proprtion_of_survivors <- calculat_proportion_of_survivors(number_of_survivors, number_of_fish)
}


library(tictoc)

tic()

s <- survival_simulation_driver(number_of_fish = 1000, mean_fish_length = 10, total_transects = 10)
s
toc()


tic()

fs <- future_apply_survival_simulation_driver(number_of_fish = 1000, mean_fish_length = 10, total_transects = 50)
fs
toc()

tic()

fm <- future_map_survival_simulation_driver(number_of_fish = 1000, mean_fish_length = 10, total_transects = 50)
fm
toc()





library(furrr)
lengths = seq(7, 15, by = 1)
reach_length_km = seq(1, 20, by=1)
tic()

df <- expand.grid(mean_fish_length = lengths, distance_traveled = reach_length_km) %>%
  arrange(mean_fish_length) %>%
  mutate(survival = future_map2(mean_fish_length,distance_traveled,
                                future_map_survival_simulation_driver, 
                                number_of_fish = 50,
                                .options=furrr_options(seed=TRUE)))
toc()

df2 <- df %>%
  mutate(mean_fish_length = as.factor(mean_fish_length),
         distance_traveled = as.numeric(distance_traveled),
         log_survival = log(as.numeric(survival) + min(as.numeric(survival)[as.numeric(survival) > 0]) / 2))



labels = c(
  "7" = "Mean length = 7 cm",
  "8" = "Mean length = 8 cm",
  "9" = "Mean length = 9 cm",
  "10" = "Mean length = 10 cm",
  "11" = "Mean length = 11 cm",
  "12" = "Mean length = 12 cm",
  "13" = "Mean length = 13 cm",
  "14" = "Mean length = 14 cm",
  "15" = "Mean length = 15 cm"
  
)



library(ggpubr)

ggplot(df2, aes(x = distance_traveled, y = log_survival)) +
  geom_point() +
  geom_smooth(method = "lm") +
  stat_regline_equation(label.x = 1, label.y = -3.35, aes(label = ..eq.label..)) +
  stat_regline_equation(label.x = 1, label.y = -3.85, aes(label = ..rr.label..)) +
  facet_wrap(~ mean_fish_length, labeller = as_labeller(labels)) +

  labs(title = "Log survival vs. distance traveled (km) per size class (cm)",
       x = "Distance traveled (km)",
       y = "Survival (log)") +
  theme(panel.spacing=unit(.05, "lines"),
        panel.border = element_rect(color = "black", fill = NA, size = 1), 
        strip.background = element_rect(color = "black", size = 1))





# lapply(joined_frame_with_sim_cols$modified_enc, check_if_encounter_occurs)
# 
# joined_frame_with_sim_cols %>%
#   filter(is.na(modified_enc))
# 
# test <- predation_survival_driver_func() %>%
#   filter(variable == 'length')
# 
# View(test)
# 
# num_trans <- calculate_num_trans_traversed()
# 
# pred_pos <- get_pred_postitions()
# stream_grid_frame <- create_stream_raster_frame(pred_pos)
# enc_prob<- calc_enc_probs(stream_grid_frame)
# 
# enc_probs <- get_enc_prob_vector(enc_prob)
# 
# enc_prob_sample <- randomly_sample_encounter_probs(enc_probs, num_trans)
# 
# survival_table <- predation_survival_driver_func()
# 
# 
# fish_frame <- set_up_dataframe_of_fish_lengths(100, precision=0) %>%
#   mutate(survival_boost = lapply(length, FUN = calcule_risk_modifier, survival_table = survival_table, selected_variable = 'length'))
# 
# 
# 
# 
# 
# df <- bind_rows(lapply(fish_frame$fish, FUN = encounter_frame, enc_prob_vector = enc_probs, num_trans_traversed = num_trans))
# 
# 
# total_fish_frame <- inner_join(fish_frame, df, by = 'fish') %>%
#   mutate(survival_boost = as.numeric(survival_boost),
#          enc_prob = as.numeric(enc_prob),
#          modified_enc = calculate_encounter_prob_based_on_length(survival_boost, enc_prob),
#          encounter = lapply(modified_enc, check_if_encounter_occurs),
#          alive = as.numeric(lapply(encounter, encounter_simulator))) 
# fish_summary <- sum_encounter_outcomes(total_fish_frame)
# 
# fish_summary %>%
#   mutate(final_status = lapply(total_alive_outcomes, FUN = evaluate_final_status_of_fish, num_trans_traversed = num_trans)) %>%
#   unnest(final_status) %>%
#   summarize(num_surviving = sum(final_status))
# 
# 
# 
# 
# # encounter_func <- function(encounter_check, pred_success = base_pred_success, pred_min = reach_pred_min){
# #   encounter <- sample(c(0, 1), size = 1, prob = c(1 - encounter_check, encounter_check))
# #   if(encounter == 1){
# #     print("Attack!!!")
# #     print(i)
# #     survival <- sample(c(0,1), 1, prob = c(pred_success, pred_min))
# #     if(survival == 1.0){
# #       outcome <- 1.0
# #     } else {
# #       outcome <- 0
# #     }
# #   } else {
# #     outcome <- 1.0
# #   }
# # }
# 
# 
# 
# number_of_fish <- 20
# live_or_die_list <- c()
# encounters <- calc_enc_probs(stream_grid_frame) %>%
#   pull(enc_prob)
# tic()
# for(i in seq(number_of_fish)){
#   fish_length <- round(rnorm(1, mean = 7.5, sd = 4), 1)
#   
#   length_prob <- survival_table %>%
#     filter(variable == 'length' & x == fish_length) %>%
#     pull(survival)
#   new_encs <- sample(encounters, num_trans_traversed)
#   encounter_check = new_encs * length_prob
#   outcome = as.numeric(lapply(encounter_check, encounter_func))
#   outcome_sum <- sum(outcome)
#   if(outcome_sum < num_trans_traversed){
#     survival <- 0
#   } else {
#     survival <- 1
#   }
#   live_or_die_list[i] <- survival
# }
# 
# toc()
# sum(live_or_die_list)/number_of_fish
# 
# thing = "length"
# 
# survival_table %>%
#   filter(variable == thing)
# 
# calcule_risk_modifier(selected_variable = 'length', variable_magnitude = 10)
