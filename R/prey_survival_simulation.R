#' Generate Prey Fish
#'
#' creates a dataframe of prey fish at a given mean size
#'
#' @param number_of_fish number (pos integer) of prey fish desired
#' @param mean_length mean length of fish in cm; default is 14
#' @param sd_length std dev of fish length in cm; default is 1.7 and scales with mean
#' @param precision integer for level of precision for length (e.g., to the 10's or 100's place)
#'
#' @return dataframe of fish, each with a unique number and length
#' @export
#'
#' @examples prey_fish_lengths(number_of_fish = 100)
#' @importFrom rlang .data
#' @source defaults based on Steel et al. 2020. "Applying the mean free-path length model to juvenile Chinook salmon migrating in the Sacramento River, California"
#'
prey_fish_lengths <- function(number_of_fish, mean_length = 14, sd_length = (1.7/14) * mean_length, precision = 1){
  vals <- round(stats::rnorm(2*number_of_fish, mean_length, sd_length), precision)
  vals <- vals[vals>0][1:number_of_fish]
  data.frame(fish = (seq(number_of_fish)), length = vals)
}

#' Survival Boost from Length
#'
#' Checks the effect a fish's length has on the probability of encountering a predator
#'
#' @param fish_frame dataframe of fish and their lengths; use output from the prey_fish_lengths function
#'
#' @return a dataframe with fish lengths and the effect of those lengths on predator encounters
#' @export
#'
#' @examples fish_frame <- prey_fish_lengths(number_of_fish = 10)
#' fish_with_length_survival_boost(fish_frame)
#' @note this function can be parallelized; e.g., by setting plan(multisession)

fish_with_length_survival_boost <- function(fish_frame){
  fish_frame %>%
    dplyr::mutate(survival_boost = furrr::future_map(length, calculate_risk_modifier, selected_variable = 'length'))
}

#' Number of Cells Traversed
#'
#' Calculates the number of grid cells each fish traverses in the environment.
#' Assumes fish take the shortest path possible through the environment (i.e., in a straight line from left to right)
#'
#' @param transect_length length of each transect in meters; default is 1000
#' @param n_transects integer of transects in the model; default is 20
#' @param grid_size length of side of raster grid in meters; default is 15
#'
#' @return integer value of cells traversed
#' @export
#'
#' @examples calculate_path_length(n_transects = 5)
#' calculate_path_length(n_transects = 20)
calculate_path_length <- function(transect_length=1000,
                                          n_transects=20,
                                          grid_size=15){

  round(transect_length * n_transects / grid_size)

}

#' Dataframe of Cells Traversed for Each Fish
#'
#' Creates a dataframe of "paths" for each fish to "travel". A number of cells equal to path_length
#' is randomly selected for each fish and becomes that fish's travel path.
#'
#' @param number_of_fish number (pos integer) of prey fish desired; should be the same as that for the prey_fish_lengths() function
#' @param enc_prob_vector vector of all encounter probabilities calculated for each cell; see calc_enc_probs()
#' @param path_length number of cells traversed; determined by calculate_path_length()
#'
#' @return a dataframe of encounter probabilities for each fish
#' @export
#' @note this function can be parallelized; e.g., by setting plan(multisession)
#' @examples pred_pos <- get_pred_positions()
#' raster <- create_stream_raster_frame(pred_pos)
#' enc_prob_frame <- calc_enc_probs(raster)
#' enc_prob_vector <- enc_prob_frame %>% dplyr::pull(enc_prob)
#' path_length = calculate_path_length()
#' frame_of_all_cells_traversed_per_fish(number_of_fish = 10,
#'                                       enc_prob_vector = enc_prob_vector,
#'                                       path_length = path_length)
#'
#'
frame_of_all_cells_traversed_per_fish <- function(number_of_fish, enc_prob_vector, path_length){
  dplyr::bind_rows(furrr::future_map(seq(number_of_fish),
                                     encounter_frame,
                                     enc_prob_vector = enc_prob_vector,
                                     path_length = path_length,
                                     .options = furrr::furrr_options(seed=TRUE)))
}


#' Joined Table of Prey Fish, Survival Length Modifiers, and Path
#'
#' joins the outputs from frame_of_all_cells_traversed_per_fish() and fish_with_length_survival_boost()
#'
#' @param fish_frame a dataframe of fish, their lengths, and length survival boosts; see fish_with_length_survival_boost()
#' @param enc_frame a dataframe of fish and the encounter probabilities for predator in each cell on their path; see frame_of_all_cells_traversed_per_fish()
#'
#' @return a joined dataframe
#' @export

combine_encounter_frame_and_fish_frame <- function(fish_frame, enc_frame){
  dplyr::inner_join(fish_frame, enc_frame, by = 'fish')
}

#' Encounter Simulations
#'
#' Simulates encounters for each fish in each cell. First checks whether an encounter occurs
#' based on fish length. If no encounter occurs, the fish survives. If an encounter occurs,
#' them a survival check occurs based on the abilities of the predator.
#'
#' @param combined_frame dataframe of prey fish, their lengths, survival boosts, and paths; see ()
#'
#' @return dataframe with all prior fish data plus columns indicating whether an encounter occurred and whether the fish survived
#' @export
#' @importFrom rlang .data
#' @note this function can be parallelized; e.g., by setting plan(multisession)

simulate_encounters <- function(combined_frame){
  combined_frame %>%
    dplyr::mutate(survival_boost = as.numeric(.data$survival_boost),
                  enc_prob = as.numeric(.data$enc_prob),
                  modified_enc = calculate_encounter_prob_based_on_length(.data$survival_boost, .data$enc_prob),
                  encounter = furrr::future_map(.data$modified_enc, check_if_encounter_occurs, .options = furrr::furrr_options(seed=TRUE)),
                  alive = as.numeric(furrr::future_map(.data$encounter, encounter_simulator, .options = furrr::furrr_options(seed=TRUE))))
}

#' Final Status of Prey Fish
#'
#' Determines whether a fish survived the environment. All the outcomes (1 or 0) of encounters are summed. If the total is less than
#' the number of transects (i.e., a death (that is, a 0) occurred in a cell), then the fish is considered to have died.
#'
#' @param fish_frame dataframe with prey fish, lengths, survival modifiers, and simulation results; see simulate_encounters()
#' @param path_length number of cells traversed; determined by calculate_path_length()
#'
#' @return dataframe with a "final_status" column to indicate overall survival for each fish
#' @export
#' @importFrom rlang .data

calculate_final_status <- function(fish_frame, path_length){
  sum_encounter_outcomes(fish_frame) %>%
    dplyr::mutate(final_status = furrr::future_map(.data$total_alive_outcomes, evaluate_final_status_of_fish, path_length)) %>%
    tidyr::unnest(.data$final_status)
}

#' Calculate Total Survivors
#'
#' Simply calculates the total number of surviving fish based on the output of calculate_final_status()
#'
#' @param final_status_frame dataframe with final status of each fish; see calculate_final_status()
#'
#' @return an integer of the number of survivng prey fish
#' @export
#'

calculate_total_survivors <- function(final_status_frame){
  surv <- final_status_frame %>%
    dplyr::summarize(num_surviving = sum(.data$final_status))
  return(surv$num_surviving[[1]])
}

#' Calculate Proportion of Surviving fish
#'
#'  Simply calculates the proportion of survivors relative to the number of fish at the start of the model
#'
#' @param number_of_survivors an integer value of surviving fish
#' @param number_of_fish the number of fish at the start of the model
#'
#' @return a proportion
#' @export
#'
#' @examples calculate_proportion_of_survivors(number_of_survivors = 10, number_of_fish = 50)
#'
calculate_proportion_of_survivors <- function(number_of_survivors, number_of_fish){
  number_of_survivors / number_of_fish
}


#' Run a Full Simulation
#'
#' Runs a full simulation with a user-specified number of fish. Users can also adjust fish
#' mean length and sd, environment size, grid size, and predator reaction distance.
#' The model runs through the following:
#' calculates predators and their positions,
#' calculates grid cells and their encounter probabilities,
#' calculates a unique path for each fish,
#' simulates and resolves encounters for each fish in each cell,
#' determines survival after all fish have gone through each cell.
#'
#' The return value is the proportion of survivors.
#'
#'
#' @param number_of_fish number (pos integer) of prey fish desired
#' @param mean_length mean length of fish in cm; default is 14
#' @param sd_length std dev of fish length in cm; default is 1.7 and scales with mean
#' @param transect_length length of each transect in meters; default is 1000
#' @param channel_width width of the channel in meters; default is 100
#' @param lit_zone_size the size of the littoral zone (i.e., nearshore area) in meters; default is 5
#' @param grid_size length of side of raster grid in meters; default is 15
#' @param n_transects integer of transects in the model; default is 20
#' @param grid_size length of side of raster grid in meters; default is 15
#' @param reaction_dis maximum distance (in m) away from a predator that can trigger an encounter; default is 0.50
#'
#' @return the proportion of surviving fish
#' @export
#'
#' @examples survival_simulation_driver (number_of_fish = 20)
#' @source defaults based on Steel et al. 2020. "Applying the mean free-path length model to juvenile Chinook salmon migrating in the Sacramento River, California"
#' and Michel et al. 2018. "Non-native fish predator density and molecular-based diet estimates suggest differing effects of predator species on Juvenile Salmon in the San Joaquin River, California"
#' @note this function can be parallelized; e.g., by setting plan(multisession)

survival_simulation_driver <- function(number_of_fish,
                                       mean_length = 14,
                                       sd_length = (1.7/14) * mean_length,
                                       transect_length = 1000,
                                       n_transects = 20,
                                       lit_zone_size = 5,
                                       channel_width = 100,
                                       grid_size = 15,
                                       reaction_dis = 0.50){
  # determine the number of gird cells each fish must traverse
  path_length <- calculate_path_length(n_transects = n_transects)

  # create a field of predators and their positions in the stream
  pred_pos <- get_pred_positions(transect_length = transect_length,
                                 n_transects = n_transects,
                                 lit_zone_size = lit_zone_size,
                                 channel_width = channel_width)
  # create a raster to layer over a grid onto the stream and get predator counts in each grid cell
  stream_grid_frame <- create_stream_raster_frame(pred_pos,
                                                  transect_length = transect_length,
                                                  channel_width = channel_width,
                                                  grid_size = grid_size,
                                                  n_transects = n_transects)
  # calculate the probability of initiating a predator encounter in each cell
  enc_prob_frame <- calc_enc_probs(stream_grid_frame,
                                   grid_size = grid_size,
                                   reaction_dis = reaction_dis)
  # pull the vector of all encounter probabilities
  enc_probs <- get_enc_prob_vector(enc_prob_frame)

  # create a dataframe of prey fish with their own lengths
  fish_frame <- prey_fish_lengths(number_of_fish,
                                  mean_length = mean_length,
                                  sd_length = sd_length)

  # calculate the boost to survival for each fish based on its length
  fish_frame_with_survival_from_length <- fish_with_length_survival_boost(fish_frame)

  # select a random assortment of cells for each fish to traverse equal to the path_length
  cells_traversed <- frame_of_all_cells_traversed_per_fish(number_of_fish = number_of_fish,
                                                           enc_prob_vector = enc_probs,
                                                           path_length = path_length)

  # join the dataframes with fish paths and fish lengths/survival boosts
  joined_fish_cell_traversed <- combine_encounter_frame_and_fish_frame(fish_frame = fish_frame_with_survival_from_length,
                                                                       enc_frame = cells_traversed)

  # add columns with simulation results for each grid cell and fish
  joined_frame_with_sim_cols <- simulate_encounters(joined_fish_cell_traversed)

  # determine whether each fish survived the simulation
  fish_summary_frame <- calculate_final_status(fish_frame = joined_frame_with_sim_cols,
                                               path_length = path_length)

  # calculate the total number of survivors
  number_of_survivors <- calculate_total_survivors(final_status_frame = fish_summary_frame)

  # calculate the proportion of survivors
  proportion_of_survivors <- calculate_proportion_of_survivors(number_of_survivors = number_of_survivors,
                                                              number_of_fish = number_of_fish)
  return(proportion_of_survivors)
}




# Helpers ---------------------------------------------------------------------

calculate_risk_modifier <- function(selected_variable, variable_magnitude){
 survival_modifiers %>%
    dplyr::filter(.data$variable == as.character(selected_variable) & .data$magnitude == as.character(variable_magnitude)) %>%
    dplyr::pull(.data$survival)
}

encounter_frame <- function(fish_number, enc_prob_vector, path_length){
  data.frame(fish = fish_number, enc_prob = randomly_sample_encounter_probs(enc_prob_vector, path_length))
}

randomly_sample_encounter_probs <- function(enc_probs, path_length){
  sample(enc_probs, path_length)
}

#' @importFrom rlang .data
get_enc_prob_vector <- function(enc_prob_dataframe){
  enc_prob_dataframe %>% dplyr::pull(.data$enc_prob)
}

calculate_encounter_prob_based_on_length <- function(fish_length_boost, enc_prob){
  (1- fish_length_boost) * enc_prob
}

check_if_encounter_occurs <- function(encounter_prob_based_on_length){
  sample(c(0, 1), size = 1, prob = c(1 - encounter_prob_based_on_length, encounter_prob_based_on_length))
}

encounter_simulator <- function(encounter_occurrence){
  if(encounter_occurrence == 1){
    outcome <- evaluate_encounter()
  } else {
    outcome <- 1
  }
  return(outcome)
}

evaluate_encounter <- function(pred_success = 0.8, pred_min = 1 - pred_success){
  survival <- sample(c(0,1), 1, prob = c(pred_success, pred_min))
  if(survival == 1.0){
    outcome <- 1.0
  } else {
    outcome <- 0
  }

}

evaluate_final_status_of_fish <- function(outcome_total, path_length){
  if(outcome_total < path_length){
    final_status <- 0
  } else {
    final_status <- 1
  }
  return(as.numeric(final_status))
}

#' @importFrom rlang .data
sum_encounter_outcomes <- function(df){
  df %>%
    dplyr::group_by(.data$fish) %>%
    dplyr::summarize(total_alive_outcomes = sum(as.numeric(.data$alive)))
}


