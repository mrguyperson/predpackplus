#' Numbers of largemouth bass in a reach
#' 
#' Calculate the numbers of largemouth bass in a reach based on field data
#'
#' @param n_transects number of transects in the model 
#' @param mean mean number of fish per transect; default is 333
#' @param sd std dev of number of fish per transect; default is 195
#'
#' @return a vector of numbers of fish per transect
#' @export 
#'
#' @examples lmb_calc(n_transects = 10)
#' @source defaults based on Michel et al. 2018. "Non-native fish predator density and molecular-based diet estimates suggest differing effects of predator species on Juvenile Salmon in the San Joaquin River, California"

lmb_calc <- function(n_transects, 
                     mean = 333, 
                     sd = 195){
  lmb <- stats::rnorm(2*n_transects, mean, sd)
  lmb <- round(lmb[lmb > 0][1:n_transects])
  return(lmb)
}

#' Numbers of striped bass in a reach
#' 
#' Calculate the numbers of striped bass in a reach based on field data
#'
#' @param n_transects number of transects in the model
#' @param agg_ratio the ratio of large aggregations to low numbers of fish in a reach; default is 1/3 
#' @param mean_low mean number of fish per transect when no aggregation is present; default is 60
#' @param sd_low std dev of number of fish per transect when no aggregation is present; default is 43
#' @param mean_high mean number of fish per transect when an aggregation is present; default is 870
#' @param sd_high std dev of number of fish per transect when an aggregation is present; default is 491.5
#'
#' @return a vector of numbers of fish per transect
#' @export 
#'
#' @examples striper_calc(n_transects = 10)
#' @source defaults based on Michel et al. 2018. "Non-native fish predator density and molecular-based diet estimates suggest differing effects of predator species on Juvenile Salmon in the San Joaquin River, California"

striper_calc <- function(n_transects, 
                         agg_ratio = 1/3, 
                         mean_low = 60, 
                         sd_low = 43, 
                         mean_high = 870, 
                         sd_high = 491.5){
  aggs <- stats::rnorm(n_transects, mean_high, sd_high)
  not_aggs <- stats::rnorm(2*n_transects, mean_low, sd_low)
  
  aggs <- aggs[aggs>0][1:round(n_transects * agg_ratio)]
  not_aggs <- not_aggs[not_aggs > 0][1:(n_transects-length(aggs))]
  
  stripers <- round(sample(c(aggs, not_aggs), n_transects))
  return(stripers)
  
}

#' Predator positions across all transects
#' 
#' Calculates the position of all predators in terms of distance from shore and distance downstream for all transects
#'
#' @param transect_length length of each transect in meters; default is 1000
#' @param n_transects number of transects in the model; default is 20
#' @param lit_zone_size the size of the littoral zone (i.e., nearshore area) in meters; default is 5
#' @param channel_width width of the channel in meters; default is 100
#'
#' @return a dataframe of x, y positions for all predators
#' @export
#'
#' @examples get_pred_positions(n_transects = 5)
#' @source defaults based on Michel et al. 2018. "Non-native fish predator density and molecular-based diet estimates suggest differing effects of predator species on Juvenile Salmon in the San Joaquin River, California"

get_pred_positions <- function(transect_length = 1000, 
                                n_transects = 20, 
                                lit_zone_size = 5, 
                                channel_width = 100){
  left_bank_dis_ds <- vector(mode="list",length=n_transects)
  left_bank_dis_fr_s <- vector(mode="list",length=n_transects)
  right_bank_dis_ds <- vector(mode="list",length=n_transects)
  right_bank_dis_fr_s <- vector(mode="list",length=n_transects)
  channel_dis_ds <- vector(mode="list",length=n_transects)
  channel_dis_fr_s <- vector(mode="list",length=n_transects)
  
  lmb <- lmb_calc(n_transects)
  stripers <- striper_calc(n_transects)
  
  for(i in seq(n_transects)){
    
    lb_ds <- distance_downstream(lmb[i] / 2 + stripers[i] / 4, transect_length, i)
    lb_fs <- distance_from_shore(lmb[i] / 2 + stripers[i] / 4, 0, lit_zone_size)
    
    rb_ds <- distance_downstream(lmb[i] / 2 + stripers[i] / 4, transect_length, i)
    rb_fs <- distance_from_shore(lmb[i] / 2 + stripers[i] / 4, channel_width - lit_zone_size, channel_width)
    
    ch_ds <- distance_downstream(stripers[i] / 2, transect_length, i)
    ch_fs <- distance_from_shore(stripers[i] / 2, lit_zone_size, channel_width - lit_zone_size)
    
    left_bank_dis_ds[[i]] <- lb_ds
    left_bank_dis_fr_s[[i]] <- lb_fs
    right_bank_dis_ds[[i]] <- rb_ds
    right_bank_dis_fr_s[[i]] <- rb_fs
    channel_dis_ds[[i]] <- ch_ds
    channel_dis_fr_s[[i]] <- ch_fs
    
    
    
  }
  df <- data.frame(distance_downstream = unlist(c(left_bank_dis_ds,
                                                  right_bank_dis_ds,
                                                  channel_dis_ds)), 
                   distance_from_shore = unlist(c(left_bank_dis_fr_s,
                                                  right_bank_dis_fr_s,
                                                  channel_dis_fr_s))) %>%
    dplyr::arrange(distance_downstream)
  return(df)
  
}

#' Creates and overlays a raster onto the predator positions
#' 
#' raster size mimics the cells in netlogo
#' returns a data.frame with coordinates for the cells and counts of fish in each cell
#' 
#' @param df dataframe of predator position data; use get_pred_positions function
#' @param transect_length length of each transect in meters; default is 1000
#' @param channel_width width of the channel in meters; default is 100
#' @param grid_size length of side of raster grid in meters; default is 15
#' @param n_transects number of transects in the model; default is 20
#'
#' @return a dataframe of raster coordinates
#' @export
#' @examples pred_pos <- get_pred_positions()
#' create_stream_raster_frame(pred_pos)
#'

create_stream_raster_frame <- function(df, 
                                       transect_length = 1000, 
                                       channel_width = 100, 
                                       grid_size = 15, 
                                       n_transects = 20){
  r <- raster::raster(xmn = 0,
                      ymn = 0, 
                      xmx = transect_length * n_transects, 
                      ymx = channel_width, 
                      res = grid_size)
  r[] <- 0
  tab <- table(raster::cellFromXY(r, df))
  r[as.numeric(names(tab))] <- tab
  d <- data.frame(sp::coordinates(r), count=r[])
  return(d)
}

#' Probability that prey encounters a predator in each cell
#' 
#' based on the idea that each predator will engage prey within a certain radius
#' total area occupied by predators divided by cell area is the encounter probability
#'
#' @param df dataframe of predator counts per cell; use output from create_stream_raster_frame function
#' @param grid_size length of side of raster grid in meters; default is 15
#' @param reaction_dis maximum distance (in m) away from a predator that can trigger an encounter; default is 0.50
#'
#' @return a dataframe of encounter probabilities per grid cell for the entire reach
#' @export
#' @examples pred_pos <- get_pred_positions()
#' raster <- create_stream_raster_frame(pred_pos)
#' calc_enc_probs(raster)
#' @importFrom rlang .data
#'

calc_enc_probs <- function(df, 
                           grid_size = 15, 
                           reaction_dis = 0.50){
  data.frame(sp::coordinates(df), count = df[]) %>%
    dplyr::mutate(pred_area = .data$count * reaction_dis^2 * pi,
                  enc_prob = .data$pred_area / grid_size^2)
}


#' Positions of all predators in the stream
#'
#' graphs all predator positions
#'
#' @param df dataframe of predator positions; use output of get_pred_positions function
#'
#' @return a graph of all predator positions
#' @export
#' @examples pred_pos <- get_pred_positions(n_transects = 5)
#' graph_pred_positions(pred_pos)
#'

graph_pred_positions <- function(df){
  ggplot2::ggplot(df, ggplot2::aes(x = distance_downstream)) +
    ggplot2::theme_classic(base_size = 30) +
    ggplot2::labs(y = "Distance from left bank (m)", x = "Distance downstream (m)") +
    ggplot2::geom_point(ggplot2::aes(y = distance_from_shore))
}


#' Heatmap of encounter probabilities 
#' 
#' Graphs a heatmap of encounter probabilities 
#'
#' @param df dataframe of encounter probabilities; use output from calc_enc_probs function
#'
#' @return a heatmap of encounter probabilities 
#' @export
#' @importFrom rlang .data
#'
#' @examples pred_pos <- get_pred_positions()
#' raster <- create_stream_raster_frame(pred_pos)
#' enc_probs <- calc_enc_probs(raster)
#' graph_enc_probs(enc_probs)

graph_enc_probs <- function(df){
  pa  <- wesanderson::wes_palettes %>% 
    names()
  pal <-  wesanderson::wes_palette(name = pa[7], n = 10, type = "continuous")
  ggplot2::ggplot(df, ggplot2::aes(.data$x, .data$y, fill = .data$enc_prob)) + 
    ggplot2::labs(y = "Distance from left bank (m)", x = "Distance downstream (m)") +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_gradientn(colors=pal)
  
}

# Helpers --------------------------------------------------------------------

distance_downstream <- function(number_of_bass, 
                                transect_length, 
                                current_transect){
  stats::runif(number_of_bass, 
        min = transect_length * current_transect - transect_length, 
        max = transect_length * current_transect)
}

distance_from_shore <- function(number_of_bass, min, max){
  distance_from_shore = stats::runif(number_of_bass, 
                              min = min, 
                              max = max)
}


