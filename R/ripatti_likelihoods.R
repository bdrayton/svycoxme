


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

sortAndIndex <- function(data, sort_vars, index = "index") {
  data %>%
    dplyr::arrange(dplyr::across({{ sort_vars }})) %>%
    dplyr::mutate({{ index }} := dplyr::row_number(), .before = everything())
}


#' calculate (XB + Zb) and exp(XB + Zb)
#'
#' @param data a tibble containing X, Z. It can have other columns too.
#' @param X character vector of X columns, i.e. fixed-effect covariates
#' @param Z character vector of Z columns, i.e. random-effect covariates. Currently this is a set of indicator variables describing cluster membership
#' @param parms vector of fixed and random effects of length ncol(data[X]) + ncol(data[Z])
#'
#' @returns The tibble given to \code{data} is returned with two additional columns. lp contains (XB + Zb) and A contains exp(XB + Zb)
#'
#' @importFrom magrittr %>%
#'


calcLinearPredictor <- function(data, X, Z, parms) {
  tempCols <- glue::glue("p{c(X, Z)}")

  long_data <- data %>%
    dplyr::select(index, all_of(X), all_of(Z)) %>%
    tidyr::pivot_longer(cols = c(all_of(X), all_of(Z)))

  wide_data <- long_data %>%
    dplyr::mutate(
      effects = rep(dplyr::all_of(parms), nrow(data)),
      partial_linear_pred = value * effects
    ) %>%
    dplyr::select(index, name, partial_linear_pred) %>%
    tidyr::pivot_wider(
      values_from = partial_linear_pred, names_from = name,
      names_prefix = "p"
    )

  dplyr::inner_join(data, wide_data, by = "index") %>%
    dplyr::mutate(
      lp = rowSums(dplyr::select(., dplyr::all_of(tempCols))),
      A = exp(lp)
    ) %>%
    dplyr::select(-dplyr::all_of(tempCols))
}





#' create a dataset for testing cox model stuff.
#'
#' takes three fixed effects and one parameter, theta, that represents
#' the variance of random effects drawn from a normal distribution.
#'


one_dataset <- function(control) {

  n = control$k*control$nk

  M = rep(1:control$k, each = control$nk)

  X1 = rnorm(n, 0, 1)
  X2 = rep(rnorm(control$k, 0, 1), each = control$nk)

  # cluster level binary treatment allocation
  X3 = rep(rep(c(1, 0), ceiling(control$k/2))[1:control$k], each = control$nk)

  X = cbind(X1, X2, X3)

  b = rnorm(control$k, 0, sqrt(control$theta))

  b_rep = rep(b, each = control$nk)

  error = rexp(n, 10)

  t = exp(-X%*%control$beta - b_rep) * error

  stat = sample(rep(c(0, 1), round(n*c(0.2, 0.8))), n)

  dataset = data.frame(X, t, stat, M)

  attr(dataset, "random_effects") <- b

  dataset

}


#' calculate sum exp(XB + Zb) for each failure time. The are the sums of linear predictors for
#' the risk set at each failure time. Also multiplies these sets by X or Z or whatever when this is needed.
#'
#' @param data tibble containing all the data.
#' @param vars names of variable to multiply the risk set by. These variables will be pivoted longer. defaults to NULL.
#' @param varCol name to the product of vars and A.
#' @param A exponentiated linear predictor. i.e. exp(XB + Zb)
#' @param index numerical index for event times.
#'
#' @returns Returns the tibble provided in data. If vars = NULL, the data will have the same number of rows, and an
#' additional column cumsum_A which is the cumulative sum of A, taken from the bottom of the data up. This is related
#' to the idea of risk sets.
#'

calcRiskSets <- function(data, vars = NULL, varCol = NULL, A = A, index = index) {
  if (is.null(vars)) {
    data %>%
      dplyr::arrange(dplyr::desc({{ index }})) %>%
      dplyr::mutate(cumsum_A = cumsum(A)) %>%
      dplyr::arrange(index)
  } else {
    varColName <- rlang::sym(varCol)

    long_data <- data %>%
      tidyr::pivot_longer(cols = dplyr::all_of(vars), names_to = {{ varCol }}) %>%
      dplyr::mutate(
        A =  {{ A }},
        "{{varColName}}_A" := A * value,
        "{{varColName}}" := forcats::as_factor(!!varColName)
      )

    long_data %>%
      dplyr::arrange(dplyr::desc({{ index }})) %>%
      dplyr::group_by(!!varColName) %>%
      dplyr::mutate(
        cumsum_A = cumsum(A),
        "cumsum_{{varColName}}_A" := cumsum(!!rlang::sym(glue::glue("{varCol}_A"))),
        cumsum_A_squared = cumsum_A^2
      ) %>%
      dplyr::arrange(index)
  }
}


#' penalised partial likelihood in Ripatti 2000
#'
#' returns the PPL proposed by Ripatti and Palmgren
#'
#' @param parms numeric vector of \eqn{\beta} and b values i.e \code{c(beta, b)}. The condition \code{length(parms) == length(c(X, Z))} must be satisfied.
#' @param X character vector containing names of X columns in data.
#' @param t column in data with time of failure/censoring data.
#' @param i column in data that identifies clusters.
#' @param dij column in data indicating if observation was censored or a failure.
#' @param D diagonal matrix with given \eqn{\theta}. Must be a square matrix of order \code{length(Z)}.
#' @param data tibble with all these columns
#'
#' @export

lp <- function(parms, X, t, dij, theta, cluster, data, ...) {

  data_with_Z <- add_Z(data, cluster)
  sortedIndexedData <- sortAndIndex(data = data_with_Z, sort_vars = {{ t }})

  b <- parms[-seq_len(length.out = length(X))]

  D = theta * diag(length(b))

  penalty = c(0.5 * t(b)%*%solve(D)%*%b)

  terms1 <- calcLinearPredictor(sortedIndexedData, X = X, Z = attr(data_with_Z, "Z_names"), parms = parms)

  terms2 <- calcRiskSets(terms1)

  ll <- terms2 %>%
    dplyr::mutate(li =  {{ dij }} * (lp - log(cumsum_A))) %>%
    dplyr::summarise(ll = sum(li)) %>%
    dplyr::pull(ll)

  ppl <- ll - penalty

  attr(ppl, "penalty") <- penalty

  ppl

}

minuslp <- function(parms, X, t, dij, D, data, ...) {
  -1 * lp(
    parms = parms, X = X, t = {{ t }},
    dij = {{ dij }}, D = D, data = data, ... = ...
  )
}



#' helper for calcCrossProducts
#'
#' takes a tibble and selects everything in the character vector f, as well
#' an index variable from sortAndIndex. pivots everything in f to long format,
#' and converts the 'name' column to a factor. Levels are determined by order they appear in the data.
#'
#' @param data A tibble containing an index column called index, and all columns named in f
#' @param f character vector of columns to make longer.
#'

makeLong <- function(data, f) {
  data %>%
    dplyr::select(index, all_of(f)) %>%
    tidyr::pivot_longer(cols = -1) %>%
    dplyr::mutate(name = forcats::as_factor(name))
}

#' Get the cross products
#'
#' Calculate cross products for a set of numeric variables.
#'
#' @param data a tibble containing the numeric variables to cross multiply
#' @param f1 first set of variables.
#' @param f2 second set of variables (can be the same variables as in f1)
#' @param n1 name for the first set of variables
#' @param n2 name for the second set of variables.
#'

calcCrossProducts <- function(data, f1, f2, n1, n2) {

  if(missing(f2)) f2 <- f1

  l1 <- makeLong(data, f1)

  l2 <- makeLong(data, f2)

  dplyr::full_join(l1, l2, by = "index") %>%
    dplyr::mutate("{n1}{n2}" := value.x * value.y) %>%
    dplyr::select(
      index, "{n1}" := name.x, "{n2}" := name.y,
      !!rlang::sym(glue::glue("{n1}{n2}"))
    )
}





#' partial derivative with respect to the fixed effects, \eqn{\beta}.
#'
#' @param parms numeric vector of \eqn{\beta} and b values i.e \code{c(beta, b)}. The condition \code{length(parms) == length(c(X, Z))} must be satisfied.
#' @param X character vector containing names of X columns in data.
#' @param cluster charter vector with column name containing cluster membership
#' @param t column in data with time of failure/censoring data.
#' @param dij column in data indicating if observation was censored or a failure.
#' @param D diagonal matrix with given \eqn{\theta}. Must be a square matrix of order \code{length(Z)}.
#' @param data tibble with all these columns


dlp_beta <- function(parms, X, cluster, t, dij, data, ...) {


#   Z_formula <- formula(glue::glue(" ~ as.factor( {cluster} ) - 1"))
#
#   Z_matrix <- model.matrix(Z_formula, data = data)
#
#   colnames(Z_matrix) <- paste0("Z", seq(ncol(Z_matrix)))

  data_with_Z <- add_Z(data, cluster)

  Z_names <- attr(data_with_Z, "Z_names")

  b <- parms[-seq_len(length.out = length(X))]

  B <- parms[seq_len(length.out = length(X))]

  sortedIndexedData <- sortAndIndex(data_with_Z, sort_vars = {{ t }})

  addedLP <- calcLinearPredictor(sortedIndexedData, X = X, Z = Z_names, parms = parms)

  addedCumsums <- calcRiskSets(addedLP, X, "Xr")

  ll <- addedCumsums %>%
    dplyr::mutate(li =  {{ dij }} * (value - cumsum_Xr_A / cumsum_A)) %>%
    dplyr::summarise(ll = sum(li), .groups = "drop")

  ll

}

#' partial derivative with respect to the random effects, b.
#'
#' @param parms numeric vector of \eqn{\beta} and b values i.e \code{c(beta, b)}. The condition \code{length(parms) == length(c(X, Z))} must be satisfied.
#' @param X character vector containing names of X columns in data.
#' @param cluster charter vector with column name containing cluster membership
#' @param t column in data with time of failure/censoring data.
#' @param dij column in data indicating if observation was censored or a failure.
#' @param D diagonal matrix with given \eqn{\theta}. Must be a square matrix of order \code{length(Z)}.
#' @param data tibble with all these columns


dlp_b <- function(parms, X, cluster, t, dij, D, data, ...) {

  # Z_formula <- formula(glue::glue(" ~ as.factor( {cluster} ) - 1"))
  #
  # Z_matrix <- model.matrix(Z_formula, data = data)
  #
  # Z_colnames <- paste0("Z", seq(ncol(Z_matrix)))
  #
  # colnames(Z_matrix) <- Z_colnames
  #
  # data_with_Z <- dplyr::bind_cols(data, data.frame(Z_matrix))

  data_with_Z <- add_Z(data, cluster)

  Z_names <- attr(data_with_Z, "Z_names")

  sortedIndexedData <- sortAndIndex(data_with_Z, sort_vars = {{ t }})

  b <- parms[-seq_len(length.out = length(X))]

  B <- parms[seq_len(length.out = length(X))]

  # The penalties need names so they can be matched later. The order of the data
  # isn't guaranteed to be what you expect after this, so best to just match on
  # name.

  penalty <- data.frame(Zi = Z_names, penalty = (solve(D) %*% b[Z_names]))

  addedLP <- calcLinearPredictor(sortedIndexedData, X = X, Z = Z_names, parms = parms)

  addedCumSums <- calcRiskSets(addedLP, vars = Z_names, varCol = "Zi")

  ll_unpenalized <- addedCumSums %>%
    dplyr::mutate(li =  {{ dij }} * (value - cumsum_Zi_A / cumsum_A)) %>%
    dplyr::summarise(ll = sum(li), .groups = "drop")

  dplyr::inner_join(ll_unpenalized, penalty, by = "Zi") %>%
    dplyr::mutate(ll = ll - penalty) %>%
    dplyr::select(Zi, ll)

}

#' lp_grd
#'
#' combine the gradient function for \eqn{\beta} and b to make one vector of
#' gradients.
#'
#' @export

lp_grd <- function(parms, X, cluster, t, dij, theta, data, ...) {

  b <- parms[-seq_len(length.out = length(X))]

  D <- theta * diag(length(b))

  l4 <- dlp_beta(parms = parms, X = X, cluster = cluster,
                 t = {{ t }}, dij = {{ dij }}, data = data) %>%
    dplyr::pull(ll)

  l5 <- dlp_b(parms  = parms, X = X, cluster = cluster,
              t = {{ t }}, dij = {{ dij }}, D = D, data = data) %>%
    dplyr::pull(ll)

  c(l4, l5)

}



#' inner product
#'
#' fast version of a %*% t(b) when a and b are diagonal matrices of the same size is:
#' sum(a*a).
#'
#'

inner <- function(a, b = a) {
  sum(a * b)
}


#' Shared frailty theta
#'
#' In a shared frailty model with i.i.d. frailty terms, the solution to the
#' estimating equation for \theta has a closed form, implemented here.
#'
#'
#'
#' @export

est_theta <- function(b, K_ppl){

  c((inner(b) + sum(diag(solve(K_ppl)))) / length(b))

}

#' Shared frailty theta
#'
#' In a shared frailty model with i.i.d. frailty terms, the solution to the
#' estimating equation for \theta has a closed form, implemented here.
#'
#' In the code, the identity matrix is the derivative of D, but for more
#' complex frailty models, D and it's derivative will be more complex (and this code won't work correctly).
#'
#'
#'
#' @export

est_theta2 <- function(par, b, K_ppl){

  theta = par

  I = diag(length(b))

  D = theta * I

  D_inv = solve(D)

  # term1 <- -0.5 * sum(diag(D_inv %*% I))
  #
  # term2 <- -0.5 * sum(diag(solve(K_ppl) %*% D_inv %*% I %*% D_inv))
  #
  # term3 <- 0.5 * t(b) %*% D_inv %*% I %*% D_inv %*% b
  #
  # term1 + term2 + term3

  c(-0.5 * (
    sum(diag(D_inv))
    + sum(diag(-solve(K_ppl) %*% D_inv %*% D_inv))
    - t(b) %*% D_inv %*% D_inv %*% b
  ))


}


#' @export



est_theta3 <- function(par, parms, X, t, cluster, dij, data, b){

  theta = par

  I = diag(length(b))

  D = theta * I

  D_inv = solve(D)

  q = length(b)

  K_ppl <- bb(parms = my_params, X = X, t = t,
              cluster = "M", dij = stat, data = sample_data,
              theta = theta, return_matrix = TRUE)


  c((sum(diag(solve(K_ppl))) + t(b)%*%b)/q - theta)


}






#' estimate all parameters in a shared frailty model
#'
#' Algorithm iterates between estimating fixed and random effects for a given
#' theta, and estimating theta given the current estimates fixed and random effects,
#' until convergence is reached.
#'
#'
#' @export

estimate_parameters <- function(start_parms, theta, X, t, cluster, dij, data){

  # estimate beta and b for given theta

  nb <- dplyr::n_distinct(data[, cluster, drop = TRUE])

  fit_optim <- optim(par = start_parms,
                     fn = lp,
                     gr = lp_grd,
                     X = X,
                     t = {{ t }},
                     cluster = cluster,
                     dij = {{ dij }},
                     D = theta * diag(nb),
                     theta = theta,
                     data = data,
                     method = "BFGS",
                     control = list(fnscale = -1))

  # estimate theta for given beta and b

  b_hat <- fit_optim$par[-seq_along(X)]

  K_ppl <- bb(parms = fit_optim$par, X = X, t = {{ t }}, cluster = cluster,
              dij = {{ dij }}, data = data, theta = theta, return_matrix = TRUE)

  # new_theta = est_theta2(par = theta, b = b_hat, K_ppl = K_ppl)

  # new_theta = est_theta(b_hat, K_ppl)
  # new_theta <- var(b_hat)

  new_theta = optim(par = theta, fn = pl_theta, b = b_hat, K_ppl = K_ppl, method = "Brent", lower = 0, upper = 100,
                    control = list(fnscale = -1))$par


  list(new_theta = new_theta,
       new_parms = fit_optim$par)


}


#' loops estimate_parameters and est_theta until convergence or max iterations
#' is reached.
#'
#' @export

estimate_parameters_loop <- function(beta,
                                     theta,
                                     k,
                                     nk,
                                     cluster = "M",
                                     max_iter = 100,
                                     convergence_threshold = 0.00001,
                                     start_theta = 4,
                                     X = c("X1", "X2", "X3")){

  sample_data <- one_dataset(control = list(k = k, nk = nk, beta = beta, theta = theta))

  fit <- survival::coxph(survival::Surv(t, stat) ~ X1 + X2 + X3, data = sample_data)

  nb <- dplyr::n_distinct(sample_data[, cluster, drop = TRUE])

  start_parameters = c(coef(fit), rep(0, nb))

  names(start_parameters) <- c(X, paste0("Z", seq_len(nb)))

  current_estimates <- estimate_parameters(start_parms = start_parameters, theta = start_theta, X = X, t = "t",
                                           cluster = cluster, dij = stat, data = sample_data)

  estimate_history <- list(current_estimates)

  for (i in 1:max_iter) {

    current_estimates <- try(estimate_parameters(start_parms = current_estimates$new_parms,
                                                 theta = current_estimates$new_theta,
                                                 X = X,
                                                 t = "t",
                                                 cluster = cluster,
                                                 dij = stat,
                                                 data = sample_data))

    if("try-error" %in% class(current_estimates)) {
      converged = FALSE
      cat("fit failed on iteration", i, "\n")
      break
    }

    estimate_history[[i+1]] <- current_estimates

    biggest_diff = max(abs(c(estimate_history[[i]]$new_theta - estimate_history[[i+1]]$new_theta,
                             estimate_history[[i]]$new_parms - estimate_history[[i+1]]$new_parms)), na.rm = TRUE)

    if(biggest_diff <= convergence_threshold){
      converged = TRUE
      cat("converged in", i, "iterations", "\n")
      break
    }

  }

  if(!exists("converged")){converged = FALSE}

  list(
    estimate_history = estimate_history,
    converged = converged,
    sample_data = sample_data)

}



#' profile likelihood for \theta
#'
#' Only for shared frailty models.
#'
#' @export
#'

pl_theta <- function(theta, b, K_ppl, lp){

  # This is the IPL, removing the bit that doesn't depend on theta, which should be
  # fine, but... maybe not? Below I try adding it back in.

  # D <- theta * diag(length(b))
  #
  # c((-1/2) * ( log(det(D)) + log(det(K_ppl)) + t(b) %*% solve(D) %*% b ))

  D <- theta * diag(length(b))

  c((-1/2) * ( log(det(D)) + log(det(K_ppl)) + 2 * lp))


}


#' likelihood function for optim. Takes a value of \eqn{\theta}, plus everything
#' needed to maximise \eqn{\beta} and \eqn{b} for that \eqn{\theta}. In future,
#' need to extend this to work for a vector of thetas.
#'
#' @export

optim_ipl <- function(theta, start_parms, X, t, cluster, dij, data) {

    # estimate beta and b for given theta

    nb <- dplyr::n_distinct(data[, cluster, drop = TRUE])

    fit_optim <- optim(par = start_parms,
                       fn = lp,
                       gr = lp_grd,
                       X = X,
                       t = {{ t }},
                       cluster = cluster,
                       dij = {{ dij }},
                       D = theta * diag(nb),
                       theta = theta,
                       data = data,
                       method = "BFGS",
                       control = list(fnscale = -1))

    # return likelihood for given theta, and estimated beta and b

    b_hat <- fit_optim$par[-seq_along(X)]

    K_ppl_hat <- bb(parms = fit_optim$par, X = X, t = {{ t }}, cluster = cluster,
                    dij = {{ dij }}, data = data, theta = theta, return_matrix = TRUE)

    ipl <- pl_theta(theta = theta, b = b_hat, K_ppl = K_ppl_hat, lp = fit_optim$value)

    ipl

}









