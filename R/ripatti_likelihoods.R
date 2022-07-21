#' sort and index data
#'
#' Sorts data by one or more variables, and adds an index to the sorted tibble
#'
#' @param data A tibble to sort
#' @param sortVars tidyselected variables to sort by
#' @param index name to give the index variable. Default is "index".
#'
#' @export
#'
#' @importFrom rlang :=
#' @importFrom magrittr %>%
#' @import Matrix

sortAndIndex <- function(data, sort_vars, index = "index") {
  data %>%
    dplyr::arrange(dplyr::across({{ sort_vars }})) %>%
    dplyr::mutate({{ index }} := dplyr::row_number(), .before = everything())
}



## Generic methods

#' paritial likelhoods
#'
#' The specific method called depends on the class of the parsed data. At time of writing,
#' I have only implemented the penalised partial likelihood in Ripatti 2000. Next will be
#' the 'full' partial likelihood, that is not dropping the part that they drop.
#'
#' params are incorrect, but I will update later, when documentation is a higher priority.
#'
#' @param parms numeric vector of \eqn{\beta} and b values i.e \code{c(beta, b)}. The condition \code{length(parms) == length(c(X, Z))} must be satisfied.
#' @param X character vector containing names of X columns in data.
#' @param t column in data with time of failure/censoring data.
#' @param i column in data that identifies clusters.
#' @param dij column in data indicating if observation was censored or a failure.
#' @param D diagonal matrix with given \eqn{\theta}. Must be a square matrix of order \code{length(Z)}.
#' @param data tibble with all these columns
#'
#' @rdname lp
#' @export

lp <- function(params, formula, parsed_data, other_args){
  UseMethod("lp", parsed_data)
}

#' lp.ppl
#' @export
#' @rdname lp

lp.ppl <- function(params, formula, parsed_data, other_args) {

  # These steps are moved to improve performance of optim.
  # # This sort var could be pulled out of the formula.
  # ds_sorted <- sortAndIndex(data, sort_vars = stat_time)
  #
  # parsed_data <- lme4::lFormula(formula, data = ds_sorted)

  beta_index <- seq_len(ncol(parsed_data$X) - 1)

  beta <- params[beta_index]
  b <- params[-beta_index]

  # drop the intercept column from the X model.matrix
  risk_score <- parsed_data$X[, -1] %*% beta + Matrix::crossprod(parsed_data$reTrms$Zt, b)

  exp_risk_score <- exp(risk_score)

  # calculate risk sets
  # I do it this way to preserve the class and other slots. rev(risk_score) coerces to a numeric vector.
  # probably no point, as cumsum convert it anyway, and I need to remake the matrix with Matrix.
  rev_exp_risk_score <- exp_risk_score
  rev_exp_risk_score@x <- rev(exp_risk_score@x)

  at_risk <- Matrix::Matrix(rev(cumsum(rev_exp_risk_score)), ncol = 1)

  # would be good to peel stat out of the Surv() object in the parsed formula.
  # This is awkward because Surv() doesn't preserve the variable names, so I may
  # need to parse the formula myself, and including consideration of the "type"
  # attribute from the evaluated call. This has sort of been resolved.

  # this is passed in via other args now
  # stat <- Matrix::Matrix(ds_sorted$stat, ncol = 1)

  # penalty

  parsed_data$reTrms$Lambdat@x <- other_args$theta[parsed_data$reTrms$Lind]

  penalty <- 0.5 * t(b) %*% Matrix::solve(parsed_data$reTrms$Lambdat) %*% b

  penalised_likelihood <- sum(other_args$stat * (risk_score - log(at_risk))) - penalty

  penalised_likelihood@x

}

#' lp_gr_beta
#' @rdname lp_gr
#' @export

lp_gr_beta <- function(params, formula, parsed_data, other_args){
  UseMethod("lp_gr_beta", parsed_data)
}

#' lp_gr_beta.ppl
#' @rdname lp_gr
#' @export

lp_gr_beta.ppl <- function(params, formula, parsed_data, other_args) {

  # ds_sorted <- sortAndIndex(data, sort_vars = stat_time)
  #
  # parsed_data <- lme4::lFormula(formula, data = ds_sorted)

  beta_index <- seq_len(ncol(parsed_data$X) - 1)

  beta <- params[beta_index]
  b <- params[-beta_index]

  risk_score <- parsed_data$X[, -1] %*% beta + Matrix::crossprod(parsed_data$reTrms$Zt, b)

  exp_risk_score <- exp(risk_score)

  rev_exp_risk_score <- exp_risk_score
  rev_exp_risk_score@x <- rev(exp_risk_score@x)

  at_risk <- Matrix(rev(cumsum(rev_exp_risk_score)), ncol = 1)

  exp_risk_score_X <- exp_risk_score * parsed_data$X[,-1]

  ## I replaced this with fast_risk_sets, which is faster when there are lots of columns. Else it's about the same.
  #
  # at_risk_list <- apply(exp_risk_score_X, 2, function(column){
  #   Matrix(rev(cumsum(rev(column))), ncol = 1)
  # })
  #
  # at_risk_X <- Reduce(cbind, at_risk_list)

  at_risk_X <- fast_risk_sets(exp_risk_score_X)

  # stat <- Matrix(ds_sorted$stat, ncol = 1)

  likelihood_gradients <- colSums(other_args$stat * (parsed_data$X[, -1] - at_risk_X/at_risk))

  likelihood_gradients

}

#' lp_gr_b
#' @rdname lp_gr
#' @export

lp_gr_b <- function(params, formula, parsed_data, other_args){
  UseMethod("lp_gr_b", parsed_data)
}

#' lp_gr_b.ppl
#' @rdname lp_gr
#' @export

lp_gr_b.ppl <- function(params, formula, parsed_data, other_args) {
  #
  # ds_sorted <- sortAndIndex(data, sort_vars = stat_time)
  #
  # parsed_data <- lme4::lFormula(formula, data = ds_sorted)

  beta_index <- seq_len(ncol(parsed_data$X) - 1)

  beta <- params[beta_index]
  b <- params[-beta_index]

  risk_score <- parsed_data$X[, -1] %*% beta + crossprod(parsed_data$reTrms$Zt, b)

  exp_risk_score <- exp(risk_score)

  rev_exp_risk_score <- exp_risk_score
  rev_exp_risk_score@x <- rev(exp_risk_score@x)

  at_risk <- Matrix(rev(cumsum(rev_exp_risk_score)), ncol = 1)

  bl <- length(b)

  exp_risk_score_Z <- Matrix(rep(exp_risk_score, bl), ncol = bl) * t(parsed_data$reTrms$Zt)

  # at_risk_list <- apply(exp_risk_score_Z, 2, function(column){
  #   Matrix(rev(cumsum(rev(column))), ncol = 1)
  # })
  #
  # at_risk_Z <- Reduce(cbind, at_risk_list)

  at_risk_Z <- fast_risk_sets(exp_risk_score_Z)

  # stat <- Matrix(ds_sorted$stat, ncol = 1)

  likelihood_gradients_unpenalised <- colSums(other_args$stat * (t(parsed_data$reTrms$Zt) - at_risk_Z/at_risk))

  parsed_data$reTrms$Lambdat@x <- other_args$theta[parsed_data$reTrms$Lind]

  penalty <- t(b) %*% solve(parsed_data$reTrms$Lambdat)

  # this accessing @x is probably poor practice?
  likelihood_gradients_unpenalised - penalty@x

}

#' lp_gr
#'
#' caluculate jacobian
#'
#' @rdname lp_gr
#' @export

lp_gr <- function(params, formula, parsed_data, other_args){
  UseMethod("lp_gr", parsed_data)
}


#' lp_gr.ppl
#' @rdname lp_gr
#' @export

lp_gr.ppl <- function(params, formula, parsed_data, other_args) {

  c(lp_gr_beta.ppl(params = params, formula = formula, parsed_data = parsed_data, other_args = other_args),
    lp_gr_b.ppl(   params = params, formula = formula, parsed_data = parsed_data, other_args = other_args))

}

#' add ppl class
#'
#' for method dispatch
#'
#' @export

make_ppl <- function(data){

  # stopifnot(is.data.frame(data))

  class(data) <- c("ppl", oldClass(data))

  data

}

#' fast risk sets
#'
#' not that fast. a target for speeding up.
#'
#' @export

fast_risk_sets <- function(a_matrix){

  # flip it
  rev_index <- rev(seq_len(nrow(a_matrix)))
  m <- a_matrix[rev_index, ]

  # cumsums
  m <- matrixStats::colCumsums(as.matrix(m))

  # flip it back
  m <- m[rev_index,]

  # restore class
  Matrix(m)
}


#' likelihood for theta
#'
#' @export
#'

theta_ipl <- function(theta, formula, parsed_data, other_args){

  # set up D
  parsed_data$reTrms$Lambdat@x <- theta[parsed_data$reTrms$Lind]

  D <- parsed_data$reTrms$Lambdat

  ## This gets passed in now, to avoid recalculation
  # n_fixed <- length(attr(terms(coxme:::formula1(formula)$fixed), "order"))
  #
  # fixed_formula <- lme4:::getFixedFormula(formula)
  #
  # fit0 <- survival::coxph(fixed_formula, data = data)
  #
  # start_params <- c(coef(fit0), rep(0, length(parsed_data$reTrms$Lind)))

  ests <- optim(par = other_args$start_params,
                fn = lp,
                gr = lp_gr,
                formula = formula,
                parsed_data = make_ppl(parsed_data),
                other_args = list(theta = theta,
                                  stat = other_args$stat),
                method = "BFGS",
                control = list(fnscale = -1),
                hessian = TRUE)

  # take the part of the hessian associated with the random effects only
  Kbb <- ests$hessian[-(1:other_args$n_fixed), -(1:other_args$n_fixed)]

  c(-0.5 * determinant(D)$modulus -0.5 * determinant(Kbb)$modulus + ests$value)

}

#' this will be a proper function that tries a few sensible point on the
#' theta sample space and returns the best one, based on theta_ipl.
#' Currently just returns a vector of ones of length n_theta

get_start_theta <- function(n_theta){

  rep(1, n_theta)

}

#' estimate theta beta and b
#'
#' @export

est_parameters <- function(formula, data) {

  # this stuff is needed for each call to lp and lp_grd, but doesn't change, so
  # I calculate it once here, and pass it in.

  ds_sorted <- sortAndIndex(data, sort_vars = stat_time)

  parsed_data <- lme4::lFormula(formula, data = ds_sorted)

  theta_start <- get_start_theta(length(parsed_data$reTrms$flist))


  n_fixed <- length(attr(terms(coxme:::formula1(formula)$fixed), "order"))

  fixed_formula <- lme4:::getFixedFormula(formula)

  fit0 <- survival::coxph(fixed_formula, data = data)

  start_params <- c(coef(fit0), rep(0, length(parsed_data$reTrms$Lind)))

  # assumes that the response is of the form Surv(time, stat). Behaviour for other Surv formats is undefined.
  stat <- Matrix(unclass(parsed_data$fr[,1])[, "status"], ncol = 1)

  theta_est <- optim(par = theta_start,
                     fn = theta_ipl,
                     gr = NULL,
                     formula = formula,
                     parsed_data = parsed_data,
                     other_args = list(n_fixed = n_fixed,
                                       start_params = start_params,
                                       stat = stat),
                     method = "L-BFGS-B",
                     control = list(fnscale = -1),
                     lower = 0.00001, upper = Inf)



  beta_b_est <- optim(par = start_params,
                      fn = lp,
                      gr = lp_gr,
                      formula = formula,
                      parsed_data = make_ppl(parsed_data),
                      other_args = list(theta = theta_est$par,
                                        stat = stat),
                      method = "BFGS",
                      control = list(fnscale = -1),
                      hessian = TRUE)

  list(theta = theta_est$par,
       beta = beta_b_est$par[seq_len(n_fixed)],
       b = beta_b_est$par[-seq_len(n_fixed)])

}

