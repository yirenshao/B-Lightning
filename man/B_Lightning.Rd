% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BLightning.R
\name{B_Lightning}
\alias{B_Lightning}
\title{B_Lightning}
\usage{
B_Lightning(
  so,
  upregulated.markers,
  downregulated.markers,
  score = "CFS",
  estimated.nonfeatured.proportion = 0.9,
  connectivity.cutoff = 4,
  num.variablefeatures = 2000,
  alpha.genes = 0.05,
  max.iter = 10000
)
}
\arguments{
\item{so}{Preprocessed Seurat Object}

\item{upregulated.markers}{a char/index vector of established upregulated gene markers}

\item{downregulated.markers}{a char/index vector of established downregulated gene markers}

\item{score}{the type of cellular score used to group cells: CFS/GSVA}

\item{estimated.nonfeatured.proportion}{estimated proportion of nonfeatured cells
in the whole cell population}

\item{connectivity.cutoff}{number of established genes that candidates have to be correlated with}

\item{num.variablefeatures}{number of highly variable genes used to calculate the null distribution of connectivity}

\item{alpha.genes}{TypeI error of identifying differentally expressed genes}
}
\value{
A list that contains the number of iterations, new upregulated gene markers,
new downregulated gene markers, pvalues of new upregulated gene markers,
pvaluds of new downregualted markers
}
\description{
An iterative fishing method for gene markers
}
