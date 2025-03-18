#' Synthetic population data
#'
#' This is made-up time-to-event data with properties that make it useful for testing and
#' demonstrating `svycoxme` functions. There is a single level of clustering,
#' identified with group_id, and the *X* covariates depend on *Z* covariates.
#'
#' @format
#' A data frame with 20,000 rows and 10 columns:
#' \describe{
#'    \item{X1}{Observation-level \eqn{N(\mu_{X1}, 1) } distributed covariate where \eqn{\mu_{X1}=0.5*(Z1 + 1) } }
#'    \item{X2}{Cluster-level \eqn{N(\mu_{X2}, 1) } distributed covariate where \eqn{\mu_{X2}=0.5*(Z2 + Z3)}}
#'    \item{X3}{Cluster-level binary covariate where \eqn{Pr(X3=1) = Z3}}
#'    \item{Z1}{Stratum membership. Takes the values 1 to 5}
#'    \item{Z2}{cluster-level N(0,1) distributed covariate}
#'    \item{Z3}{cluster-level Uniform(0,1) distributed covariate}
#'    \item{stat_time}{Event or Censoring time}
#'    \item{stat}{Event/Censoring indicator. Event=1; Censoring=0}
#'    \item{group_id}{Unique cluster ID}
#'    \item{obs_id}{Unique observation ID}
#'    \item{sampled}{Sampling indicator. Is this observation in \link[svycoxme]{samp_srcs}?}
#' }
"pop"

#' Simple random cluster sample of 100 clusters from synthetic population data, \link[svycoxme]{pop}.
#'
#' This is made-up time-to-event data with properties that make it useful for testing and
#' demonstrating `svycoxme` functions. There is a single level of clustering,
#' identified with group_id, and the *X* covariates depend on *Z* covariates.
#'
#' @format
#' A data frame with 20,000 rows and 10 columns:
#' \describe{
#'    \item{X1}{Observation-level \eqn{N(\mu_{X1}, 1) } distributed covariate where \eqn{\mu_{X1}=0.5*(Z1 + 1)}.}
#'    \item{X2}{Cluster-level \eqn{N(\mu_{X2}, 1) } distributed covariate where \eqn{\mu_{X2}=0.5*(Z2 + Z3)}.}
#'    \item{X3}{Cluster-level binary covariate where \eqn{Pr(X3=1) = Z3}.}
#'    \item{Z1}{Stratum membership. Takes the values 1 to 5.}
#'    \item{Z2}{cluster-level N(0,1) distributed covariate.}
#'    \item{Z3}{cluster-level Uniform(0,1) distributed covariate.}
#'    \item{stat_time}{Event or Censoring time.}
#'    \item{stat}{Event/Censoring indicator. Event=1; Censoring=0.}
#'    \item{group_id}{Unique cluster ID.}
#'    \item{obs_id}{Unique observation ID.}
#'    \item{fpc}{Total number of clusters in the population.}
#'    \item{weight}{Observation-level inverse probability of selection weight.}
#' }
"samp_srcs"
