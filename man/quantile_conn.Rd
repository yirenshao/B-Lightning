% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/quantile_conn.R
\name{quantile_conn}
\alias{quantile_conn}
\title{Filter out candidates with low connectivity}
\usage{
quantile_conn(
  index_marker_known,
  index_marker_candidate,
  so,
  cutoff,
  num.variablefeatures
)
}
\arguments{
\item{index_marker_known}{input/known gene markers}

\item{index_marker_candidate}{candidate gene markers to be tested for their
connectivity with input/known gene markers}

\item{so}{Seurat Object}

\item{cutoff}{threshold of connectivity}
}
\value{
a char vector of genes with high connectivity
}
\description{
This function calculates each candidate's connectivity with known gene markers
and only keep the ones with high connectivity based on set threshold.
}
