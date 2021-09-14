library('raster')
library('dplyr')
library('tidyverse')

library('wesanderson')
# library('sn')


default_lit_zone_size <- 5
default_channel_width <- 100
default_transect_length  <- 1000
default_n_transects  <- 20
default_grid_size <- 15
default_reaction_dis <- 0.5

# number_of_stripers <- function(mean = 334, sd = 447, skewness = 0.9){
#   stripers <- 0
#   params <- cp2dp(c(mean, sd, skewness), "SN")
#   while(stripers < 21){
#     stripers <- rsn(1, dp = params)
#     if (stripers >= 21){
#       return(round(stripers))
#     }
#   }
# }


# calculates the number of striped bass in a reach;
# assumes that striped bass either occur in small numbers or in huge aggregations
# means and SDs are from values in Michel et al. 2018 

# number_of_stripers <- function(mean_low = 60, sd_low = 43, mean_high = 870, sd_high = 491.5, agg_prob = 1/3){
#   stripers <- 0
#   aggregation <- sample(c(TRUE,FALSE), 1, prob = c(agg_prob, 1 - agg_prob))
#   
#   if(aggregation){
#     while(stripers < 1){
#       stripers <- rnorm(1, mean = mean_high, sd = sd_high)
#       if (stripers >= 1){
#         return(round(stripers))
#         }
#       }
#     }
#   else{
#     while(stripers < 1){
#       stripers <- rnorm(1, mean = mean_low, sd = sd_low)
#       if (stripers >= 1){
#         return(round(stripers))
#       }
#     }
#   }  
# }
# 
# # calculates number of striped bass in a reach;
# # assumes a normal distribution but only picks values with 1 SD of the mean
# 
# number_of_lmb <- function(mean = 333, sd = 195){
#   lmb <- 0
#   while(lmb <= mean - sd | lmb >= mean + sd){
#     lmb <- rnorm(1, mean = mean, sd = mean)
#     if (lmb >= mean - sd & lmb <= mean + sd){
#       return(round(lmb))
#     }
#   }
# }
# 
# # creates a data.frame of predator x, y positions in the littoral zone;
# # i.e., near shore
# # calculates the position from shore and the position downstream
# # can be used for either the left or right bank
# 
# bank_preds <- function(number_of_bass, transect_length, lit_zone_size, zone_start, zone_end, current_transect){
#   data.frame(distance_downstream = runif(number_of_bass, 
#                                          min = transect_length * current_transect - transect_length, 
#                                          max = transect_length * current_transect),
#              distance_from_shore = runif(number_of_bass, 
#                                          min = zone_start, 
#                                          max = zone_end)) %>%
#     arrange(distance_downstream)
# }
# 
# # creates a data.frame of predator x, y positions in the channel area
# 
# channel_preds <- function(number_of_bass, transect_length, lit_zone_size, channel_width, current_transect){
#   data.frame(distance_downstream = runif(number_of_bass, 
#                                          min = transect_length * current_transect - transect_length, 
#                                          max = transect_length * current_transect),
#              distance_from_shore = runif(number_of_bass, 
#                                          min = lit_zone_size, 
#                                          max = channel_width - lit_zone_size)) %>%
#     arrange(distance_downstream)
# }
# 
# # calculates the number of lmb and stb per reach and their positions for each zone
# # compiles all points into one data.frame
# 
# get_pred_postitions <- function(transect_length, n_transects, lit_zone_size, channel_width){
#   left_bank_preds_list <- vector(mode="list",length=n_transects)
#   right_bank_preds_list <- vector(mode="list",length=n_transects)
#   channel_preds_list <- vector(mode="list",length=n_transects)
#   
#   for(i in seq(n_transects)){
#     lmb <-  number_of_lmb()
#     stripers <- number_of_stripers()
#     
#     left_bank_preds <- bank_preds(lmb / 2 + stripers / 4, 
#                                   transect_length, 
#                                   lit_zone_size, 
#                                   0, 
#                                   lit_zone_size, 
#                                   i)
#     right_bank_preds <- bank_preds(lmb / 2 + stripers / 4, 
#                                    transect_length, 
#                                    lit_zone_size, 
#                                    channel_width - lit_zone_size, 
#                                    channel_width, 
#                                    i)
#     channel_preds <- channel_preds(stripers / 2,
#                                    transect_length, 
#                                    lit_zone_size, 
#                                    channel_width, 
#                                    i)
#     left_bank_preds_list[[i]] <- left_bank_preds
#     right_bank_preds_list[[i]] <- right_bank_preds
#     channel_preds_list[[i]] <- channel_preds
#   }
#   df <- bind_rows(left_bank_preds_list, right_bank_preds_list, channel_preds_list) %>%
#     arrange(distance_downstream)
# 
# }

# creates and overlays a raster onto the predator positons
# raster size mimics the cells in netlogo
# returns a data.frame with coordinates for the cells and counts of fish in each cell

create_stream_raster_frame <- function(df, 
                                       transect_length = default_transect_length, 
                                       channel_width = default_channel_width, 
                                       grid_size = default_grid_size, 
                                       n_transects = default_n_transects){
  r <- raster(xmn = 0, 
              ymn = 0, 
              xmx = transect_length * n_transects, 
              ymx = channel_width, 
              res = grid_size)
  r[] <- 0
  tab <- table(cellFromXY(r, df))
  r[as.numeric(names(tab))] <- tab
  d <- data.frame(coordinates(r), count=r[])
}

# calculates encounter probability for each cell
# based on the idea that each predator will engage prey with a certain radius
# total area occupied by predators divided by cell area is the encounter probability

calc_enc_probs <- function(df, 
                           grid_size = default_grid_size, 
                           reaction_dis = default_reaction_dis){
  data.frame(coordinates(df), count=df[]) %>%
    mutate(pred_area = count * reaction_dis^2 * pi,
           enc_prob = pred_area / grid_size^2)
}

# graphs the positions of all predators in the stream

graph_pred_positions <- function(df){
  ggplot(df, aes(x = distance_downstream)) +
    theme_classic(base_size = 30) +
    labs(y = "Distance from left bank (m)", x = "Distance downstream (m)") +
    geom_point(aes(y = distance_from_shore))
}

# graphs a heatmap of encounter probabilities 

graph_enc_probs <- function(df){
  pa  <- wes_palettes %>% 
    names()
  pal <-  wes_palette(name = pa[7], n = 10, type = "continuous")
  ggplot(df, aes(x, y, fill = enc_prob)) + 
    labs(y = "Distance from left bank (m)", x = "Distance downstream (m)") +
    geom_tile() +
    scale_fill_gradientn(colors=pal)
  
}



lmb_calc <- function(n, 
                     mean = 333, 
                     sd = 195){
  draws <- rnorm(2*n, mean, sd)
  draws <- draws[draws > 0][1:n]
}

striper_calc <- function(n, 
                         agg_ratio = 1/3, 
                         mean_low = 60, 
                         sd_low = 43, 
                         mean_high = 870, 
                         sd_high = 491.5){
  aggs <- rnorm(n, mean_high, sd_high)
  not_aggs <- rnorm(2*n, mean_low, sd_low)
  
  aggs <- aggs[aggs>0][1:round(n * agg_ratio)]
  not_aggs <- not_aggs[not_aggs > 0][1:(n-length(aggs))]

  stripers <- sample(c(aggs, not_aggs), n)
  
}


distance_downstream <- function(number_of_bass, 
                                transect_length = default_transect_length, 
                                current_transect){
  runif(number_of_bass, 
        min = transect_length * current_transect - transect_length, 
        max = transect_length * current_transect)
}

distance_from_shore <- function(number_of_bass, min, max){
  distance_from_shore = runif(number_of_bass, 
                              min = min, 
                              max = max)
}



get_pred_postitions <- function(transect_length = default_transect_length, 
                                n_transects = default_n_transects, 
                                lit_zone_size = default_lit_zone_size, 
                                channel_width = default_channel_width){
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
    arrange(distance_downstream)
  
}

# 
# 
# 
# 
# tic()
# pred_pos <- get_pred_postitions()
# 
# 
# 
# stream_grid_frame <- create_stream_raster_frame(pred_pos)
# enc_probs <- calc_enc_probs(stream_grid_frame)
# 
# 
# graph_pred_positions(pred_pos)
# graph_enc_probs(enc_probs)

