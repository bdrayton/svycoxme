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

  to_return <- penalised_likelihood@x

  attr(to_return, "penalty") <- penalty@x

  to_return

}

#' lp.pplextra
#' @export
#' @rdname lp

lp.pplextra <- function(params, formula, parsed_data, other_args) {

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
  # attribute from the evaluated call.

  # pass in via other_args
  # stat <- Matrix::Matrix(ds_sorted$stat, ncol = 1)

  # penalty

  parsed_data$reTrms$Lambdat@x <- other_args$theta[parsed_data$reTrms$Lind]

  D <- parsed_data$reTrms$Lambdat

  D_inverse <- Matrix::solve(D)

  penalty <- 0.5 * t(b) %*% D_inverse %*% b

  penalised_likelihood <- sum(other_args$stat * (risk_score - log(at_risk))) - penalty

  # terms 1 and 2 from ripatti and palmgren, which they drop.

  # determinant returns the logarithm by default. Matrix just calls the base
  # method by default. IDK if there speed advantage from sparseness.

  term1 <- -0.5 * log(det(D))

  # second term is more complicated. Need cumulative hazard (breslow), and Zt(Z)

  cumulative_hazard <- Matrix::Matrix(cumsum(other_args$stat/at_risk), ncol = 1)


  Zt <- parsed_data$reTrms$Zt

  Zt_ncol <- ncol(Zt)

  ZtZ_exp_risk_score <- vector(mode = "list", length = Zt_ncol)

  for (i in seq_len(Zt_ncol)) {

    ZtZ_exp_risk_score[[i]] <-  cumulative_hazard[i] * exp_risk_score[i] * Matrix::tcrossprod(Zt[,i, drop = FALSE])

  }

  # term2 can be NaN, as this determinant can be negative, making the log undefined.

  term2_determinant <- determinant(Reduce("+", ZtZ_exp_risk_score) - D_inverse)

  if(term2_determinant$sign == -1) {
    warning("Term 2 determinant is negative, so the likelihood is undefined")
    return(NaN)
  }


  # strip attributes with c()
  c(term1 + 0.5 * term2_determinant$modulus + penalised_likelihood@x)

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

#' lp_gr_beta.pplextra
#' @rdname lp_gr
#' @export

lp_gr.pplextra <- function(params, formula, parsed_data, other_args) {

  beta_index <- seq_len(ncol(parsed_data$X) - 1)

  beta <- params[beta_index]
  b <- params[-beta_index]

  # drop the intercept column from the X model.matrix
  X <- parsed_data$X[, -1, drop = FALSE]

  risk_score <- X %*% beta + Matrix::crossprod(parsed_data$reTrms$Zt, b)

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
  # attribute from the evaluated call.

  # pass in via other_args
  # stat <- Matrix::Matrix(ds_sorted$stat, ncol = 1)

  # penalty

  parsed_data$reTrms$Lambdat@x <- other_args$theta[parsed_data$reTrms$Lind]

  D <- parsed_data$reTrms$Lambdat

  D_inverse <- Matrix::solve(D)

  penalty <- 0.5 * t(b) %*% D_inverse %*% b

  penalised_likelihood <- sum(other_args$stat * (risk_score - log(at_risk))) - penalty

  # terms 1 and 2 from ripatti and palmgren, which they drop.

  # determinant returns the logarithm by default. Matrix just calls the base
  # method by default. IDK if there speed advantage from sparseness.

  term1 <- -0.5 * log(det(D))

  # second term is more complicated. Need cumulative hazard (breslow), and Zt(Z)

  cumulative_hazard <- Matrix::Matrix(cumsum(other_args$stat/at_risk), ncol = 1)


  Zt <- parsed_data$reTrms$Zt

  Zt_ncol <- ncol(Zt)
  Zt_nrow <- nrow(Zt)

  ZtZ_exp_risk_score <- vector(mode = "list", length = Zt_ncol)

  # each X will have a sub list of length i
  X_ZtZ_exp_risk_score <- lapply(seq_len(other_args$n_fixed), function(i){vector(mode = "list", length = Zt_ncol)})

  # each Z will also have a sublist of length i
  Z_ZtZ_exp_risk_score <- lapply(seq_len(Zt_nrow), function(i){vector(mode = "list", length = Zt_ncol)})


  for (i in seq_len(Zt_ncol)) {

    ZtZ_exp_risk_score[[i]] <-  cumulative_hazard[i] * exp_risk_score[i] * Matrix::tcrossprod(Zt[,i, drop = FALSE])

    for (j in seq_len(other_args$n_fixed)) {

      X_ZtZ_exp_risk_score[[j]][[i]] <-  X[i, j, drop = TRUE] * ZtZ_exp_risk_score[[i]]

    }

    for (k in seq_len(Zt_nrow)) {

      Z_ZtZ_exp_risk_score[[k]][[i]] <-  parsed_data$reTrms$Zt[k, i, drop = TRUE] * ZtZ_exp_risk_score[[i]]

    }

  }

  reduced_X_ZtZ <- lapply(X_ZtZ_exp_risk_score, Reduce, f = `+`)
  reduced_Z_ZtZ <- lapply(Z_ZtZ_exp_risk_score, Reduce, f = `+`)

  solved_reduced_ZtZ <- solve(Reduce("+", ZtZ_exp_risk_score) - D_inverse)

  gr_beta <- unlist(lapply(reduced_X_ZtZ, function(x){
    product <- solved_reduced_ZtZ %*% x

    -0.5 * sum(diag(product))

  }))

  gr_b <- unlist(lapply(reduced_Z_ZtZ, function(x){
    product <- solved_reduced_ZtZ %*% x

    -0.5 * sum(diag(product))

  }))

  base_gr <- lp_gr.ppl(params = params, formula = formula, parsed_data = parsed_data, other_args = other_args)

  base_gr - c(gr_beta, gr_b)

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

#' add pplextra class
#'
#' for method dispatch
#'
#' @export

make_pplextra <- function(data){

  # stopifnot(is.data.frame(data))

  class(data) <- c("pplextra", oldClass(data))

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
                parsed_data = parsed_data,
                other_args = c(list(theta = theta), other_args),
                method = "BFGS",
                control = list(fnscale = -1),
                hessian = TRUE)

  # In Ripatti Palmgren they only use this part of the Hessian matrix.
  if(other_args$re_only) {
    # take the part of the hessian associated with the random effects only
    Kbb <- ests$hessian[-(1:other_args$n_fixed), -(1:other_args$n_fixed)]
  } else {
    Kbb <- ests$hessian
  }

  c(-0.5 * determinant(D)$modulus -0.5 * determinant(Kbb)$modulus + ests$value)

}

#' likelihood for theta
#'
#' without nested optimisation of beta and b.
#'
#' @export
#'

theta_ipl2 <- function(theta, formula, parsed_data, Kbb, value){

  # set up D
  parsed_data$reTrms$Lambdat@x <- theta[parsed_data$reTrms$Lind]

  D <- parsed_data$reTrms$Lambdat

  c(-0.5 * determinant(D)$modulus -0.5 * determinant(Kbb)$modulus + value)

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

est_parameters <- function(formula, data, method, start_params = NULL, theta_start = NULL,
                           control = list(re_only = TRUE, convergence_threshold = 0.00001, max_iter = 10)) {

  # Allow start parameters to be passed in, so I can use ppl, and then refin with pplextra.

  # this stuff is needed for each call to lp and lp_grd, but doesn't change, so
  # I calculate it once here, and pass it in.

  ds_sorted <- sortAndIndex(data, sort_vars = stat_time)

  parsed_data <- lme4::lFormula(formula, data = ds_sorted)

  parsed_data <- switch(method,
         "ppl" = make_ppl(parsed_data),
         "pplextra" = make_pplextra(parsed_data))

  if (is.null(theta_start))
    theta_start <- get_start_theta(length(parsed_data$reTrms$flist))

  n_fixed <- length(attr(terms(coxme:::formula1(formula)$fixed), "order"))

  if (is.null(start_params)) {

    fixed_formula <- lme4:::getFixedFormula(formula)

    fit0 <- survival::coxph(fixed_formula, data = data)

    # start_params <- c(coef(fit0), rep(0, length(parsed_data$reTrms$Lind)))
    # the full likelihood in ripatti often evaluates to NaN when the random effects are all 0,
    # so I'll set them using rnorm given theta_start.

    start_params <- c(coef(fit0),
                      rnorm(n = length(parsed_data$reTrms$Lind),
                            mean = 0,
                            sd = sqrt(theta_start[parsed_data$reTrms$Lind])))

  }

  # assumes that the response is of the form Surv(time, stat). Behaviour for other Surv formats is undefined.
  stat <- Matrix(unclass(parsed_data$fr[,1])[, "status"], ncol = 1)

  theta_est <- optim(par = theta_start,
                     fn = theta_ipl,
                     gr = NULL,
                     formula = formula,
                     parsed_data = parsed_data,
                     other_args = list(n_fixed = n_fixed,
                                       start_params = start_params,
                                       stat = stat,
                                       re_only = control$re_only),
                     method = "L-BFGS-B",
                     control = list(fnscale = -1),
                     lower = 0.00001, upper = Inf)

  beta_b_est <- optim(par = start_params,
                      fn = lp,
                      gr = lp_gr,
                      formula = formula,
                      parsed_data = parsed_data,
                      other_args = list(theta = theta_est$par,
                                        stat = stat),
                      method = "BFGS",
                      control = list(fnscale = -1))

  # iterate between these steps until convergence.

  for (i in seq_len(control$max_iter)){

    current_ests <- c(beta_b_est$par, theta_est$par)

    theta_est <- optim(par = theta_est$par,
                       fn = theta_ipl,
                       gr = NULL,
                       formula = formula,
                       parsed_data = parsed_data,
                       other_args = list(n_fixed = n_fixed,
                                         start_params = beta_b_est$par,
                                         stat = stat,
                                         re_only = control$re_only),
                       method = "L-BFGS-B",
                       control = list(fnscale = -1),
                       lower = 0.00001, upper = Inf)

    beta_b_est <- optim(par = beta_b_est$par,
                        fn = lp,
                        gr = lp_gr,
                        formula = formula,
                        parsed_data = parsed_data,
                        other_args = list(theta = theta_est$par,
                                          stat = stat),
                        method = "BFGS",
                        control = list(fnscale = -1))

    # test convergence
    converged <- max(abs(current_ests - c(beta_b_est$par, theta_est$par))) <= control$convergence_threshold

    if(converged) break

  }

  beta_b_est$hessian <- optimHess(beta_b_est$par, fn = lp, gr = lp_gr, parsed_data = parsed_data,
                                 other_args = list(theta = theta_est$par,
                                                   stat = stat),
                                 control = list(fnscale = -1))

  theta_est$hessian <- optimHess(theta_est$par, fn = theta_ipl, gr = NULL, parsed_data = parsed_data,
                                 other_args = list(n_fixed = n_fixed,
                                                   start_params = beta_b_est$par,
                                                   stat = stat,
                                                   re_only = control$re_only),
                                 control = list(fnscale = -1))


  # # get the hessian for theta # this is about as accurate as the call for hessian4
  # theta_est$hessian2 <- optimHess(theta_est$par, fn = theta_ipl, gr = NULL, parsed_data = parsed_data,
  #                                other_args = list(n_fixed = n_fixed,
  #                                                  start_params = beta_b_est$par,
  #                                                  stat = stat,
  #                                                  re_only = control$re_only),
  #                                control = list(fnscale = -1))
  #
  # # a better hessian?
  # theta_est$hessian3 <- numDeriv::hessian(theta_ipl, x = theta_est$par,
  #                                         parsed_data = parsed_data,
  #                                         other_args = list(n_fixed = n_fixed,
  #                                                           start_params = beta_b_est$par,
  #                                                           stat = stat,
  #                                                           re_only = control$re_only))

  # the final hessian to try.
  # this works best
  # theta2 <- optim(par = theta_est$par,
  #       fn = theta_ipl,
  #       gr = NULL,
  #       formula = formula,
  #       parsed_data = parsed_data,
  #       other_args = list(n_fixed = n_fixed,
  #                         start_params = beta_b_est$par,
  #                         stat = stat,
  #                         re_only = control$re_only),
  #       method = "L-BFGS-B",
  #       control = list(fnscale = -1),
  #       lower = 0.00001, upper = Inf,
  #       hessian = TRUE)

  # theta_est$hessian4 <- theta2$hessian

  ests <- list(theta = theta_est$par,
       beta = beta_b_est$par[seq_len(n_fixed)],
       b = beta_b_est$par[-seq_len(n_fixed)])

  attr(ests, "theta_ests") <- theta_est
  attr(ests, "beta_b_est") <- beta_b_est
  attr(ests, "iterations") <- i

  ests

}


#' estimate theta beta and b
#'
#' Instead of nested optim calls, uses alternating estimation of beta and theta until convergence.
#'
#' @export

est_parameters_loop <- function(formula, data, method, start_params = NULL, theta_start = NULL, control = list(re_only = TRUE)) {

  # Allow start parameters to be passed in, so I can use ppl, and then refin with pplextra.

  # this stuff is needed for each call to lp and lp_grd, but doesn't change, so
  # I calculate it once here, and pass it in.

  ds_sorted <- sortAndIndex(data, sort_vars = stat_time)

  parsed_data <- lme4::lFormula(formula, data = ds_sorted)

  parsed_data <- switch(method,
                        "ppl" = make_ppl(parsed_data),
                        "pplextra" = make_pplextra(parsed_data))

  if (is.null(theta_start))
    theta_start <- get_start_theta(length(parsed_data$reTrms$flist))

  n_fixed <- length(attr(terms(coxme:::formula1(formula)$fixed), "order"))

  if (is.null(start_params)) {

    fixed_formula <- lme4:::getFixedFormula(formula)

    fit0 <- survival::coxph(fixed_formula, data = data)

    # start_params <- c(coef(fit0), rep(0, length(parsed_data$reTrms$Lind)))
    # the full likelihood in ripatti often evaluates to NaN when the random effects are all 0,
    # so I'll set them using rnorm given theta_start.

    start_params <- c(coef(fit0),
                      rnorm(n = length(parsed_data$reTrms$Lind),
                            mean = 0,
                            sd = sqrt(theta_start[parsed_data$reTrms$Lind])))

  }

  # assumes that the response is of the form Surv(time, stat). Behaviour for other Surv formats is undefined.
  stat <- Matrix(unclass(parsed_data$fr[,1])[, "status"], ncol = 1)

  # initial estimation of beta and b

  beta_b_est <- optim(par = start_params,
                      fn = lp,
                      gr = lp_gr,
                      formula = formula,
                      parsed_data = parsed_data,
                      other_args = list(theta = theta_est$par,
                                        stat = stat),
                      method = "BFGS",
                      control = list(fnscale = -1),
                      hessian = TRUE)

  # initial estimation of theta, given beta and b.






  theta_est <- optim(par = theta_start,
                     fn = theta_ipl,
                     gr = NULL,
                     formula = formula,
                     parsed_data = parsed_data,
                     other_args = list(n_fixed = n_fixed,
                                       start_params = start_params,
                                       stat = stat,
                                       re_only = control$re_only),
                     method = "L-BFGS-B",
                     control = list(fnscale = -1),
                     lower = 0.00001, upper = Inf,
                     hessian = TRUE)

  beta_b_est <- optim(par = start_params,
                      fn = lp,
                      gr = lp_gr,
                      formula = formula,
                      parsed_data = parsed_data,
                      other_args = list(theta = theta_est$par,
                                        stat = stat),
                      method = "BFGS",
                      control = list(fnscale = -1),
                      hessian = TRUE)

  # get the hessian for theta
  theta_est$hessian2 <- optimHess(theta_est$par, fn = theta_ipl, gr = NULL, parsed_data = parsed_data,
                                  other_args = list(n_fixed = n_fixed,
                                                    start_params = beta_b_est$par,
                                                    stat = stat,
                                                    re_only = control$re_only),
                                  control = list(fnscale = -1))

  # a better hessian?
  theta_est$hessian3 <- numDeriv::hessian(theta_ipl, x = theta_est$par,
                                          parsed_data = parsed_data,
                                          other_args = list(n_fixed = n_fixed,
                                                            start_params = beta_b_est$par,
                                                            stat = stat,
                                                            re_only = control$re_only))

  # the final hessian to try.
  theta2 <- optim(par = theta_est$par,
                  fn = theta_ipl,
                  gr = NULL,
                  formula = formula,
                  parsed_data = parsed_data,
                  other_args = list(n_fixed = n_fixed,
                                    start_params = beta_b_est$par,
                                    stat = stat,
                                    re_only = control$re_only),
                  method = "L-BFGS-B",
                  control = list(fnscale = -1),
                  lower = 0.00001, upper = Inf,
                  hessian = TRUE)

  theta_est$hessian4 <- theta2$hessian


  ests <- list(theta = theta_est$par,
               beta = beta_b_est$par[seq_len(n_fixed)],
               b = beta_b_est$par[-seq_len(n_fixed)])

  attr(ests, "theta_ests") <- theta_est
  attr(ests, "beta_b_est") <- beta_b_est

  ests

}























