% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generate_preds_and_encs.R
\name{create_stream_raster_frame}
\alias{create_stream_raster_frame}
\title{Environmental Raster}
\usage{
create_stream_raster_frame(
  df,
  transect_length = 1000,
  channel_width = 100,
  grid_size = 15,
  n_transects = 20
)
}
\arguments{
\item{df}{dataframe of predator position data; use get_pred_positions function}

\item{transect_length}{length of each transect in meters; default is 1000}

\item{channel_width}{width of the channel in meters; default is 100}

\item{grid_size}{length of side of raster grid in meters; default is 15}

\item{n_transects}{integer of transects in the model; default is 20}
}
\value{
a dataframe of raster coordinates
}
\description{
Creates and overlays a raster onto the predator positions
raster size mimics the cells in netlogo
returns a data.frame with coordinates for the cells and counts of fish in each cell
}
\examples{
pred_pos <- get_pred_positions()
create_stream_raster_frame(pred_pos)

}
