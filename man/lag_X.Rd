% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lag_X.R
\name{lag_X}
\alias{lag_X}
\title{Generate lag covariate from design matrix
The resulting design matrix W correspond to column-wise expansion of Phi}
\usage{
lag_X(X, ncolW = NULL)
}
\arguments{
\item{X, }{design matrix for mean vector}
}
\value{
W, design matrix for autocoefficient coefficient
}
\description{
Generate lag covariate from design matrix
The resulting design matrix W correspond to column-wise expansion of Phi
}
