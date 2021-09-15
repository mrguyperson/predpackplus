#' Generate Prey Fish
#'
#' creates a dataframe of prey fish at a given mean size
#'
#' @param number_of_fish number (pos integer) of prey fish desired
#' @param mean mean length of fish in cm; default is 14
#' @param sd std dev of fish length in cm; default is 1.7 and scales with mean
#' @param precision integer for level of precision for length (e.g., to the 10's or 100's place)
#'
#' @return dataframe of fish, each with a unique number and length
#' @export
#'
#' @examples prey_fish_lengths(number_of_fish = 100)
#' @importFrom rlang .data
#' @source defaults based on Steel et al. 2020. "Applying the mean free-path length model to juvenile Chinook salmon migrating in the Sacramento River, California"
#'
prey_fish_lengths <- function(number_of_fish, mean = 14, sd = (1.7/14) * mean, precision = 1){
  vals <- round(stats::rnorm(2*number_of_fish, mean, sd), precision)
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
#' @examples calculate_num_cells_traversed(n_transects = 5)
#' calculate_num_cells_traversed(n_transects = 20)
calculate_num_cells_traversed <- function(transect_length=1000,
                                          n_transects=20,
                                          grid_size=15){

  round(transect_length * n_transects / grid_size)

}

#' Dataframe of Cells Traversed for Each Fish
#'
#' Creates a dataframe of "paths" for each fish to "travel". A number of cells equal to num_trans_traversed
#' is randomly selected for each fish and becomes that fish's travel path.
#'
#' @param number_of_fish number (pos integer) of prey fish desired; should be the same as that for the prey_fish_lengths() function
#' @param enc_prob_vector vector of all encounter probabilities calculated for each cell; see calc_enc_probs()
#' @param num_trans_traversed number of cells traversed; determined by calculate_num_cells_traversed()
#'
#' @return a dataframe of encounter probabilities for each fish
#' @export
#' @note this function can be parallelized; e.g., by setting plan(multisession)
#' @examples pred_pos <- get_pred_positions()
#' raster <- create_stream_raster_frame(pred_pos)
#' enc_prob_frame <- calc_enc_probs(raster)
#' enc_prob_vector <- enc_prob_frame %>% dplyr::pull(enc_prob)
#' num_trans = calculate_num_cells_traversed()
#' frame_of_all_cells_traversed_per_fish(number_of_fish = 10,
#'                                       enc_prob_vector = enc_prob_vector,
#'                                       num_trans_traversed = num_trans)
#'
#'
frame_of_all_cells_traversed_per_fish <- function(number_of_fish, enc_prob_vector, num_trans_traversed){
  dplyr::bind_rows(furrr::future_map(seq(number_of_fish),
                                     encounter_frame,
                                     enc_prob_vector = enc_prob_vector,
                                     num_trans_traversed = num_trans_traversed,
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
  combined_frame%>%
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
#' @param num_trans_traversed number of cells traversed; determined by calculate_num_cells_traversed()
#'
#' @return dataframe with a "final_status" column to indicate overall survival for each fish
#' @export
#' @importFrom rlang .data

calculate_final_status <- function(fish_frame, num_trans_traversed){
  sum_encounter_outcomes(.data$sim) %>%
    dplyr::mutate(final_status = furrr::future_map(.data$total_alive_outcomes, evaluate_final_status_of_fish, num_trans_traversed)) %>%
    tidyr::unnest(.data$final_status)
}

# Helpers ---------------------------------------------------------------------

calculate_risk_modifier <- function(selected_variable, variable_magnitude){
 survival_modifiers %>%
    dplyr::filter(.data$variable == as.character(selected_variable) & .data$magnitude == as.character(variable_magnitude)) %>%
    dplyr::pull(.data$survival)
}

encounter_frame <- function(fish_number, enc_prob_vector, num_trans_traversed){
  data.frame(fish = fish_number, enc_prob = randomly_sample_encounter_probs(enc_prob_vector, num_trans_traversed))
}

randomly_sample_encounter_probs <- function(enc_probs, num_trans_traversed){
  sample(enc_probs, num_trans_traversed)
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

evaluate_final_status_of_fish <- function(outcome_total, num_trans_traversed){
  if(outcome_total < num_trans_traversed){
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


