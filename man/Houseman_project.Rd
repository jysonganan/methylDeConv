% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/deconvAlgorithms.R
\name{Houseman_project}
\alias{Houseman_project}
\title{Houseman algorithm}
\usage{
Houseman_project(betamatrix, compTable, probes_select)
}
\arguments{
\item{betamatrix}{Beta value matrix of methylation array for mixture samples.}

\item{compTable}{Average reference profiles over each cell type.}

\item{probes_select}{Probes selected.}
}
\value{
The matrix of deconvolution results with rows as mixture samples and columns as cell types.
}
\description{
Houseman deconvolution algorithm
}
