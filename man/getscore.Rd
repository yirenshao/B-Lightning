% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getscore.R
\name{getscore}
\alias{getscore}
\title{getscore}
\usage{
getscore(score = "CFS", index_gene_up, index_gene_down, subcounts)
}
\arguments{
\item{index_gene_up}{names of upregulated markers}

\item{index_gene_down}{names of downregulated markers}

\item{subcounts}{count matrix (gene markers by cells)
#' @return score vector}
}
\description{
Calculate scores of cells
}
\details{
This function calculates scores of each cell based on gene markers' expression levels.
}
