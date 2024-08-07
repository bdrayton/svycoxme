% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data_generation.R
\name{one_dataset}
\alias{one_dataset}
\title{create a dataset for testing cox model stuff.}
\usage{
one_dataset(
  formula,
  dists,
  dist_args,
  error,
  censoring_time,
  coefficients = c(),
  random_effect_variance = list(),
  seed = NULL,
  random_effect_seed = NULL
)
}
\arguments{
\item{formula}{Model equation for the model you plan to fit to the data. Only the right hand side is used at the moment.}

\item{dists}{expressions for generating the variables in the formula}

\item{dist_args}{any values that are needed for the expressions in dists}

\item{coefficients}{vector of fixed effect coefficients. If there are more coefficients than fixed effects, only those up to the number of fixed effects will be used and the rest ignored.}

\item{random_effect_variance}{named list of thetas. Names to match formula terms.}

\item{seed}{passed to set.seed. Use if you need the random processes to be reproducible.}
}
\description{
create a dataset for testing cox model stuff.
}
