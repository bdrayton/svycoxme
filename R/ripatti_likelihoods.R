


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

sortAndIndex <- function(data, sort_vars, index = "index") {
  data %>%
    dplyr::arrange(dplyr::across({{ sort_vars }})) %>%
    dplyr::mutate({{ index }} := dplyr::row_number(), .before = everything())
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
#' @param Z character vector containing names of Z columns in data.
#' @param t column in data with time of failure/censoring data.
#' @param i column in data that identifies clusters.
#' @param dij column in data indicating if observation was censored or a failure.
#' @param D diagonal matrix with given \eqn{\theta}. Must be a square matrix of order \code{length(Z)}.
#' @param data tibble with all these columns
#'
#' @export

lp <- function(parms, X, t, dij, D, data, ...) {

  Z_matrix <- model.matrix( ~ as.factor(M) - 1, data = data)

  colnames(Z_matrix) <- paste0("Z", seq(ncol(Z_matrix)))

  data_with_Z <- dplyr::bind_cols(data, data.frame(Z_matrix))

  b <- parms[-seq_len(length.out = length(X))]

  penalty = 0.5 * sum(diag(solve(D)*b))
  #
  # penalty <- data %>%
  #   dplyr::distinct({{ cluster }}) %>%
  #   dplyr::mutate(
  #     theta_inverse = diag(solve(D)),
  #     penalty = dplyr::all_of(b)^2 * theta_inverse
  #   ) %>%
  #   dplyr::summarise(penalty = 0.5 * sum(penalty)) %>%
  #   dplyr::pull(penalty)

  sortedIndexedData <- sortAndIndex(data_with_Z, {{ t }})

  terms1 <- calcLinearPredictor(sortedIndexedData, X = X, Z = colnames(Z_matrix), parms = parms)

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

  Z_formula <- formula(glue::glue(" ~ as.factor( {cluster} ) - 1"))

  Z_matrix <- model.matrix(Z_formula, data = data)

  colnames(Z_matrix) <- paste0("Z", seq(ncol(Z_matrix)))

  data_with_Z <- dplyr::bind_cols(data, data.frame(Z_matrix))

  b <- parms[-seq_len(length.out = length(X))]

  B <- parms[seq_len(length.out = length(X))]

  sortedIndexedData <- sortAndIndex(data_with_Z, sort_vars = {{ t }})

  addedLP <- calcLinearPredictor(sortedIndexedData, X = X, Z = colnames(Z_matrix), parms = parms)

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

  Z_formula <- formula(glue::glue(" ~ as.factor( {cluster} ) - 1"))

  Z_matrix <- model.matrix(Z_formula, data = data)

  Z_colnames <- paste0("Z", seq(ncol(Z_matrix)))

  colnames(Z_matrix) <- Z_colnames

  data_with_Z <- dplyr::bind_cols(data, data.frame(Z_matrix))

  b <- parms[-seq_len(length.out = length(X))]

  B <- parms[seq_len(length.out = length(X))]

  # The penalties need names so they can be matched later. The order of the data
  # isn't guaranteed to be what you expect after this, so best to just match on
  # name.

  penalty <- data.frame(Zi = Z_colnames, penalty = 0.5 * (solve(D) %*% b[Z_colnames]))

  sortedIndexedData <- sortAndIndex(data_with_Z, sort_vars = {{ t }})

  addedLP <- calcLinearPredictor(sortedIndexedData, X = X, Z = Z_colnames, parms = parms)

  addedCumSums <- calcRiskSets(addedLP, vars = Z_colnames, varCol = "Zi")

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

lp_grd <- function(parms, X, cluster, t, dij, D, data) {

  l4 <- dlp_beta(parms = parms, X = X, cluster = cluster,
                 t = {{ t }}, dij = {{ dij }}, data = data) %>%
    dplyr::pull(ll)

  l5 <- dlp_b(parms  = parms, X = X, cluster = cluster,
              t = {{ t }}, dij = {{ dij }}, D = D, data = data) %>%
    dplyr::pull(ll)

  c(l4, l5)

}

#' @export
#'
minus_lp_grd <- function(parms, X, Z, t, wi, i, D, wji, dij, data){

  -lp_grd(parms = parms, X = X, Z = Z,
          t = {{ t }}, wi = {{ wi }}, i = {{ i }}, D = D, wji = {{ wji }},
          dij = {{ dij }}, data = data)

}

























