% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/buildReference.R
\name{build_reference_450k}
\alias{build_reference_450k}
\title{Build the reference library for 450k arrays}
\usage{
build_reference_450k(extend = TRUE)
}
\arguments{
\item{extend}{If TRUE, generate the extended reference library; otherwise, generate the reference library of
only 6 immune cell types. Default value is TRUE.}
}
\value{
A list of beta value reference matrix (ref_betamatrix) and reference cell types (ref_phenotype).
}
\description{
Build the reference library for 450k to include "CD8T", "CD4T", "NK", "Bcell", "Mono", "Gran",
And extend the reference library for 450k to include "Epithelial", "Fibroblast".
}
