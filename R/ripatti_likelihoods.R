


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
    data |>
      dplyr::arrange(dplyr::desc({{ index }})) |>
      dplyr::mutate(cumsum_A = cumsum(A)) |>
      dplyr::arrange(index)
  } else {
    varColName <- rlang::sym(varCol)

    long_data <- data |>
      tidyr::pivot_longer(cols = dplyr::all_of(vars), names_to = {{ varCol }}) |>
      dplyr::mutate(
        wi_wji_A = {{ wi }} * {{ wji }} * {{ A }},
        "{{varColName}}_A" := A * value,
        "{{varColName}}" := forcats::as_factor(!!varColName)
      )

    long_data |>
      dplyr::arrange(dplyr::desc({{ index }})) |>
      dplyr::group_by(!!varColName) |>
      dplyr::mutate(
        cumsum_A = cumsum(A),
        "cumsum_{{varColName}}_A" := cumsum(!!rlang::sym(glue::glue("{varCol}_A"))),
        cumsum_A_squared = cumsum_A^2
      ) |>
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

lp <- function(parms, X, t, cluster, dij, D, data) {

  Z_matrix <- model.matrix( ~ as.factor(M) - 1, data = data)

  colnames(Z_matrix) <- paste0("Z", seq(ncol(Z_matrix)))

  data_with_Z <- dplyr::bind_cols(data, data.frame(Z_matrix))

  b <- parms[-seq_len(length.out = length(X))]

  penalty <- data |>
    dplyr::distinct({{ cluster }}) |>
    dplyr::mutate(
      theta_inverse = diag(solve(D)),
      penalty = dplyr::all_of(b)^2 * theta_inverse
    ) |>
    dplyr::summarise(penalty = 0.5 * sum(penalty)) |>
    dplyr::pull(penalty)

  sortedIndexedData <- sortAndIndex(data_with_Z, {{ t }})

  terms1 <- calcLinearPredictor(sortedIndexedData, X = X, Z = colnames(Z_matrix), parms = parms)

  terms2 <- calcRiskSets(terms1)

  ll <- terms2 |>
    dplyr::mutate(li =  {{ dij }} * (lp - log(cumsum_A))) |>
    dplyr::summarise(ll = sum(li)) |>
    dplyr::pull(ll)

  ll - penalty

}

minuslp <- function(parms, X, t, cluster, dij, D, data) {
  -1 * lp(
    parms = parms, X = X, t = {{ t }}, cluster = {{ cluster }},
    dij = {{ dij }}, D = D, data = data
  )
}



