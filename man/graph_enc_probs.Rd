% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generate_preds_and_encs.R
\name{graph_enc_probs}
\alias{graph_enc_probs}
\title{Heatmap of Encounter Probabilities}
\usage{
graph_enc_probs(df)
}
\arguments{
\item{df}{dataframe of encounter probabilities; use output from calc_enc_probs function}
}
\value{
a heatmap of encounter probabilities
}
\description{
Graphs a heatmap of encounter probabilities
}
\examples{
pred_pos <- get_pred_positions()
raster <- create_stream_raster_frame(pred_pos)
enc_probs <- calc_enc_probs(raster)
graph_enc_probs(enc_probs)
}
