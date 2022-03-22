
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
  data |>
    dplyr::arrange(dplyr::across({{ sort_vars }})) |>
    dplyr::mutate({{ index }} := dplyr::row_number(), .before = everything())
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
  data |>
    dplyr::select(index, all_of(f)) |>
    tidyr::pivot_longer(cols = -1) |>
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
  l1 <- makeLong(data, f1)

  l2 <- makeLong(data, f2)

  dplyr::full_join(l1, l2, by = "index") |>
    dplyr::mutate("{n1}{n2}" := value.x * value.y) |>
    dplyr::select(
      index, "{n1}" := name.x, "{n2}" := name.y,
      !!rlang::sym(glue::glue("{n1}{n2}"))
    )
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

  long_data <- data |>
    dplyr::select(index, all_of(X), all_of(Z)) |>
    tidyr::pivot_longer(cols = c(all_of(X), all_of(Z)))

  wide_data <- long_data |>
    dplyr::mutate(
      effects = rep(dplyr::all_of(parms), nrow(data)),
      partial_linear_pred = value * effects
    ) |>
    dplyr::select(index, name, partial_linear_pred) |>
    tidyr::pivot_wider(
      values_from = partial_linear_pred, names_from = name,
      names_prefix = "p"
    )

  dplyr::inner_join(data, wide_data, by = "index") %>%
    dplyr::mutate(
      lp = rowSums(dplyr::select(., dplyr::all_of(tempCols))),
      A = exp(lp)
    ) |>
    dplyr::select(-dplyr::all_of(tempCols))
}


#' #' calculate sum wi sum wji exp(XB + Zb) for each failure time. these are the weighted sums of linear predictors for
#' #' the risk set at each failure time. Also multiplies these sets by X or Z or whatever when this is needed.
#' #'
#' #' @param data tibble containing all the data.
#' #' @param vars names of variable to multiply the risk set by. These variables will be pivoted longer. defaults to NULL.
#' #' @param varCol name to the product of vars and wi \* wji \* A.
#' #' @param wi cluster level weights
#' #' @param wji individual level weights within clusters
#' #' @param A exponentiated linear predictor. i.e. exp(XB + Zb)
#' #' @param index numerical index for event times.
#' #'
#' #' @returns Returns the tibble provided in data. If vars = NULL, the data will have the same number of rows, and two
#' #' additional columns - wi_wji_A is the product of the exponetiated linear predictor and the weighting variables.
#' #' cumsum_wi_wji_A is the cumulative sum, taken from the bottom of the data up. This is related to the idea of
#' #' risk sets.
#' #'
#'
#' calcRiskSets <- function(data, vars = NULL, varCol = NULL, wi = wi, wji = wji, A = A, index = index) {
#'   if (is.null(vars)) {
#'     data |>
#'       dplyr::mutate(wi_wji_A = {{ wi }} * {{ wji }} * {{ A }}) |>
#'       dplyr::arrange(dplyr::desc({{ index }})) |>
#'       dplyr::mutate(cumsum_wi_wji_A = cumsum(wi_wji_A)) |>
#'       dplyr::arrange(index)
#'   } else {
#'     varColName <- rlang::sym(varCol)
#'
#'     long_data <- data |>
#'       tidyr::pivot_longer(cols = dplyr::all_of(vars), names_to = {{ varCol }}) |>
#'       dplyr::mutate(
#'         wi_wji_A = {{ wi }} * {{ wji }} * {{ A }},
#'         "{{varColName}}_wi_wji_A" := wi_wji_A * value,
#'         "{{varColName}}" := forcats::as_factor(!!varColName)
#'       )
#'
#'     long_data |>
#'       dplyr::arrange(dplyr::desc({{ index }})) |>
#'       dplyr::group_by(!!varColName) |>
#'       dplyr::mutate(
#'         cumsum_wi_wji_A = cumsum(wi_wji_A),
#'         "cumsum_{{varColName}}_wi_wji_A" := cumsum(!!rlang::sym(glue::glue("{varCol}_wi_wji_A"))),
#'         cumsum_wi_wji_A_squared = cumsum_wi_wji_A^2
#'       ) |>
#'       dplyr::arrange(index)
#'   }
#' }

#'
#' #' calculate the penalty terms for dlp_b_b
#' #'
#' #' Calcualtes the penalty term for the second partial derivative of the pseudo-partial likelihood
#' #' with respect to the random effects, b.
#' #'
#' #' @param data tibble of data.
#' #' @param D_theta diagonal matrix for given theta. Row and column names must
#' #' match Z columns in the data
#' #' @param Z character vector of columns in data indicating cluster membership
#' #'
#'
#' calcPenalty <- function(data, D, Z, wi, i) {
#'
#'   # cluster weights
#'   cluster_weights <- data |>
#'     dplyr::distinct({{ i }}, {{ wi }}) |>
#'     dplyr::pull({{ wi }})
#'
#'   # diagonal matrix with nr = length(Z)
#'   penalty <- -0.5 * diag(cluster_weights) %*% solve(D)
#'
#'   rownames(penalty) <- colnames(penalty)
#'
#'   penalty |> tibble::as_tibble(rownames = "Zi") |>
#'     tidyr::pivot_longer(cols = dplyr::all_of(Z), names_to = "Zm", values_to = "penalty")
#'
#'   # DTibble <- tibble::tibble(
#'   #   Zi = colnames(D),
#'   #   tibble::as_tibble(solve(D))
#'   # ) |>
#'   #   tidyr::pivot_longer(-1, names_to = "Zm")
#'   #
#'   #
#'   # clusterWeights <- data |>
#'   #   dplyr::select(dplyr::all_of(Z), {{ wi }}) |>
#'   #   dplyr::distinct() |>
#'   #   dplyr::arrange(dplyr::across(dplyr::all_of(Z), dplyr::desc)) |>
#'   #   dplyr::mutate(Zi = dplyr::all_of(Z), .before = dplyr::all_of(Z)) |>
#'   #   tidyr::pivot_longer(cols = dplyr::all_of(Z), names_to = "Zm") |>
#'   #   dplyr::mutate(wi = wi * value) |>
#'   #   dplyr::select(Zi, Zm, wi)
#'   #
#'   # dplyr::full_join(clusterWeights, DTibble, by = c("Zi", "Zm")) |>
#'   #   dplyr::mutate(penalty = wi * value) |>
#'   #   dplyr::select(Zi, Zm, penalty) |>
#'   #   dplyr::mutate(dplyr::across(c(Zi, Zm), forcats::as_factor))
#' }
#'
#' # calcPenalty(myData, D_theta, zColumns)
#' #
#'
#'
#'
#' #' pseudo likelihood in Wang 2019
#' #'
#' #' returns the likelihood proposed in Wang 2019
#' #'
#' #' @param parms numeric vector of \eqn{\beta} and b values i.e \code{c(beta, b)}. The condition \code{length(parms) == length(c(X, Z))} must be satisfied.
#' #' @param X character vector containing names of X columns in data.
#' #' @param Z character vector containing names of Z columns in data.
#' #' @param t column in data with time of failure/censoring data.
#' #' @param i column in data that identifies clusters.
#' #' @param wi column in data containing cluster-level weights.
#' #' @param wji column in data containing within-cluster weights.
#' #' @param dij column in data indicating if observation was censored or a failure.
#' #' @param D diagonal matrix with given \eqn{\theta}. Must be a square matrix of order \code{length(Z)}.
#' #' @param data tibble with all these columns
#' #'
#' #' @export
#'
#'
#'
#' lp <- function(parms, X, Z, t, i, wi, wji, dij, D, data) {
#'   b <- parms[-seq_len(length.out = length(X))]
#'
#'   penalty <- data |>
#'     dplyr::distinct({{ i }}, {{ wi }}) |>
#'     dplyr::mutate(
#'       theta_inverse = diag(solve(D)),
#'       penalty = {{ wi }} * dplyr::all_of(b)^2 * theta_inverse
#'     ) |>
#'     dplyr::summarise(penalty = 0.5 * sum(penalty)) |>
#'     dplyr::pull(penalty)
#'
#'   sortedIndexedData <- sortAndIndex(data, {{ t }})
#'
#'   terms1 <- calcLinearPredictor(sortedIndexedData, X = X, Z = Z, parms = parms)
#'
#'   terms2 <- calcRiskSets(terms1)
#'
#'   ll <- terms2 |>
#'     dplyr::mutate(li = {{ wi }} * {{ wji }} * {{ dij }} * (lp - log(cumsum_wi_wji_A))) |>
#'     dplyr::summarise(ll = sum(li)) |>
#'     dplyr::pull(ll)
#'
#'   ll - penalty
#' }

# lp(parms = c(B, b), X = xColumns, Z = zColumns, t = t, wi = wi, wji = wji, dij = dij, D = D_theta, i = cluster, data = myData)
#'
#' #'
#' #'
#' #' pseudo likelihood in Wang 2019, multiplied by -1
#' #'
#' #' returns the likelihood proposed in Wang 2019, multiplied by -1, which is required by many optimisers because they find the min not the max.
#' #'
#' #' @param parms numeric vector of \eqn{\beta} and b values i.e \code{c(beta, b)}. The condition \code{length(parms) == length(c(X, Z))} must be satisfied.
#' #' @param X character vector containing names of X columns in data.
#' #' @param Z character vector containing names of Z columns in data.
#' #' @param t column in data with time of failure/censoring data.
#' #' @param i column in data that identifies clusters.
#' #' @param wi column in data containing cluster-level weights.
#' #' @param wji column in data containing within-cluster weights.
#' #' @param dij column in data indicating if observation was censored or a failure.
#' #' @param D diagonal matrix with given \eqn{\theta}. Must be a square matrix of order \code{length(Z)}.
#' #' @param data tibble with all these columns
#' #'
#' #'
#' #'
#' #'
#'
#' minuslp <- function(parms, X, Z, t, i, wi, wji, dij, D, data) {
#'   -1 * lp(
#'     parms = parms, X = X, Z = Z, t = {{ t }}, i = {{ i }}, wi = {{ wi }}, wji = {{ wji }},
#'     dij = {{ dij }}, D = D, data = data
#'   )
#' }


#' partial derivative with respect to the fixed effects, \eqn{\beta}.
#'
#' @param parms numeric vector of \eqn{\beta} and b values i.e \code{c(beta, b)}. The condition \code{length(parms) == length(c(X, Z))} must be satisfied.
#' @param X character vector containing names of X columns in data.
#' @param Z character vector containing names of Z columns in data.
#' @param t column in data with time of failure/censoring data.
#' @param wi column in data containing cluster-level weights.
#' @param wji column in data containing within-cluster weights.
#' @param dij column in data indicating if observation was censored or a failure.
#' @param D diagonal matrix with given \eqn{\theta}. Must be a square matrix of order \code{length(Z)}.
#' @param data tibble with all these columns


dlp_beta <- function(parms, X, Z, t, wi, wji, dij, data) {

  B <- parms[1:length(X)]

  b <- parms[-(1:length(X))]

  sortedIndexedData <- sortAndIndex(data, sort_vars = {{ t }})

  addedLP <- calcLinearPredictor(sortedIndexedData, X = X, Z = Z, parms = parms)

  addedCumsums <- calcRiskSets(addedLP, X, "Xr")

  ll <- addedCumsums |>
    dplyr::mutate(li = {{ wi }} * {{ wji }} * {{ dij }} * (value - cumsum_Xr_wi_wji_A / cumsum_wi_wji_A)) |>
    dplyr::summarise(ll = sum(li), .groups = "drop")

  ll
}


# dlp_beta(parms = c(B, b), X = xColumns, Z = zColumns, t = t, wi = wi, wji = wji, dij = dij, data = myData)

#' partial derivative with respect to the random effects, b.
#'
#' @param parms numeric vector of \eqn{\beta} and b values i.e \code{c(beta, b)}. The condition \code{length(parms) == length(c(X, Z))} must be satisfied.
#' @param X character vector containing names of X columns in data.
#' @param Z character vector containing names of Z columns in data.
#' @param t column in data with time of failure/censoring data.
#' @param wi column in data containing cluster-level weights.
#' @param wji column in data containing within-cluster weights.
#' @param dij column in data indicating if observation was censored or a failure.
#' @param D diagonal matrix with given \eqn{\theta}. Must be a square matrix of order \code{length(Z)}.
#' @param data tibble with all these columns


dlp_b <- function(parms, X, Z, t, wi, wji, dij, i, D, data) {

  B <- parms[1:length(X)]

  b <- parms[-(1:length(X))]

  penalty <- data |>
    dplyr::distinct({{ i }}, {{ wi }}) |>
    dplyr::mutate(
      Zi = forcats::as_factor(dplyr::all_of(Z)),
      theta_inverse = diag(solve(D)),
      penalty = 0.5 * {{ wi }} * dplyr::all_of(b) * theta_inverse
    ) |>
    dplyr::select(Zi, penalty)

  sortedIndexedData <- sortAndIndex(data, sort_vars = {{ t }})

  addedLP <- calcLinearPredictor(sortedIndexedData, X = X, Z = Z, parms = parms)

  addedCumSums <- calcRiskSets(addedLP, Z, "Zi")

  ll_unpenalized <- addedCumSums |>
    dplyr::mutate(li = {{ wi }} * {{ wji }} * {{ dij }} * (value - cumsum_Zi_wi_wji_A / cumsum_wi_wji_A)) |>
    dplyr::summarise(ll = sum(li), .groups = "drop")

  dplyr::inner_join(ll_unpenalized, penalty, by = "Zi") |>
    dplyr::mutate(ll = ll - penalty) |>
    dplyr::select(Zi, ll)

}

#' lp_grd
#'
#' combine the gradient function for \eqn{\beta} and b to make one vector of
#' gradients.
#'
#' @export

lp_grd <- function(parms, X, Z, t, wi, i, D, wji, dij, data) {

  l4 <- dlp_beta(parms = parms, X = X, Z = Z, t = {{t}}, wi = {{wi}},
                 wji = {{wji}}, dij = {{dij}}, data = data) %>%
    dplyr::pull(ll)

  l5 <- dlp_b(parms = parms, X = X, Z = Z, t = {{t}}, wi = {{wi}},
              i = {{i}}, D = D, wji = {{wji}}, dij = {{dij}}, data = data) %>%
    dplyr::pull(ll)

  c(l4, l5)

}

#' @export
minus_lp_grd <- function(parms, X, Z, t, wi, i, D, wji, dij, data){

  -lp_grd(parms = parms, X = X, Z = Z,
          t = {{ t }}, wi = {{ wi }}, i = {{ i }}, D = D, wji = {{ wji }},
          dij = {{ dij }}, data = data)

}


# dlp_b(parms = c(B, b), X = xColumns, Z = zColumns, t = t, wi = wi,
# i = cluster, D = D_theta, wji = wji, dij = dij, data = myData)


#' Wang 2019 likelihood differentiated twice with respect to the fixed effects, \eqn{\beta}.
#'
#' @param parms numeric vector of \eqn{\beta} and b values i.e \code{c(beta, b)}. The condition \code{length(parms) == length(c(X, Z))} must be satisfied.
#' @param X character vector containing names of X columns in data.
#' @param Z character vector containing names of Z columns in data.
#' @param t column in data with time of failure/censoring data.
#' @param wi column in data containing cluster-level weights.
#' @param wji column in data containing within-cluster weights.
#' @param dij column in data indicating if observation was censored or a failure.
#' @param D diagonal matrix with given \eqn{\theta}. Must be a square matrix of order \code{length(Z)}.
#' @param data tibble with all these columns
#'
#' @export


dlp_beta_beta <- function(parms, X, Z, t, wi, wji, dij, data) {
  require(magrittr)

  sortedIndexedData <- sortAndIndex(data, {{ t }})

  addedLP <- calcLinearPredictor(sortedIndexedData, X = X, Z = Z, parms = parms)

  termsX <- calcRiskSets(addedLP, X, "Xr")

  XCrossProducts <- calcCrossProducts(sortedIndexedData, X, X, "Xr", "Xs")

  termsX2 <- dplyr::left_join(

    termsX %>%
      dplyr::select(
        index, A, {{ wi }}, {{ wji }}, wi_wji_A, Xr, cumsum_wi_wji_A,
        cumsum_wi_wji_A_squared, cumsum_Xr_wi_wji_A
      ),
    XCrossProducts,
    by = c("index", "Xr")
  ) %>%
    dplyr::left_join(

      termsX %>% dplyr::select(index, Xr, cumsum_Xs_wi_wji_A = cumsum_Xr_wi_wji_A),
      by = c("index", "Xs" = "Xr")
    ) %>%
    dplyr::arrange(dplyr::desc(index)) %>%
    dplyr::group_by(Xr, Xs) %>%
    dplyr::mutate(
      XrXs_wi_wji_A = XrXs * wi_wji_A,
      cumsum_XrXs_wi_wji_A = cumsum(XrXs_wi_wji_A)
    ) %>%
    dplyr::arrange(index) %>%
    dplyr::group_by(Xr, Xs)

  ll <- termsX2 %>%
    dplyr::mutate(
      li = {{ wi }} * {{ wji }} * ((-1 * (cumsum_XrXs_wi_wji_A / cumsum_wi_wji_A_squared) *
        (cumsum_wi_wji_A / cumsum_wi_wji_A_squared)) +
        ((cumsum_Xr_wi_wji_A / cumsum_wi_wji_A_squared) *
          (cumsum_Xs_wi_wji_A / cumsum_wi_wji_A_squared)))
    ) %>%
    dplyr::summarise(ll = sum(li), .groups = "drop")

  ll
}

#
# dlp_beta_beta(parms = c(B, b), X = xColumns, Z = zColumns, t = t,
#               wi = wi, wji = wji, dij = dij, data = myData)



#' Wang 2019 likelihood differentiated twice with respect to the random effects, b.
#'
#' @param parms numeric vector of \eqn{\beta} and b values i.e \code{c(beta, b)}. The condition \code{length(parms) == length(c(X, Z))} must be satisfied.
#' @param X character vector containing names of X columns in data.
#' @param Z character vector containing names of Z columns in data.
#' @param t column in data with time of failure/censoring data.
#' @param wi column in data containing cluster-level weights.
#' @param wji column in data containing within-cluster weights.
#' @param dij column in data indicating if observation was censored or a failure.
#' @param D diagonal matrix with given \eqn{\theta}. Must be a square matrix of order \code{length(Z)}. row and column names must match \code{Z}
#' @param data tibble with all these columns
#'
#' @export
#'
#' @importFrom magrittr %>%
#'


dlp_b_b <- function(parms, X, Z, t, wi, wji, dij, i, D, data) {

  sortedIndexedData <- sortAndIndex(data, {{ t }})

  terms1 <- calcLinearPredictor(sortedIndexedData, X = X, Z = Z, parms = parms)

  termsZi <- calcRiskSets(terms1, Z, "Zi")

  ZCrossProducts <- calcCrossProducts(sortedIndexedData, Z, Z, "Zi", "Zm")

  key_parts_of_termsZi <- termsZi |> dplyr::select(
      index,
      A,
      {{ dij }},
      {{ wi }},
      {{ wji }},
      wi_wji_A, Zi,
      cumsum_wi_wji_A,
      cumsum_wi_wji_A_squared,
      cumsum_Zi_wi_wji_A)

  termsZi2 <-
    dplyr::left_join(
      key_parts_of_termsZi,
      ZCrossProducts,
      by = c("index", "Zi"))

  termsZm <- termsZi |> dplyr::select(index, Zi, cumsum_Zm_wi_wji_A = cumsum_Zi_wi_wji_A)

  termsZiZm <- dplyr::left_join(
      termsZi2, termsZm,
      by = c("index", "Zm" = "Zi")) |>
    dplyr::arrange(dplyr::desc(index)) %>%
    dplyr::group_by(Zi, Zm) |>
    dplyr::mutate(
      ZiZm_wi_wji_A = ZiZm * wi_wji_A,
      cumsum_ZiZm_wi_wji_A = cumsum(ZiZm_wi_wji_A)) |>
    dplyr::arrange(index) %>%
    dplyr::group_by(Zi, Zm)

  penaltyMatrix <- calcPenalty(data, D = D, Z = Z, i = {{ i }}, wi = {{ wi }})

  ll_unpenalised <- termsZiZm |>
    dplyr::mutate(
      li = {{ wi }} * {{ wji }} * {{ dij }} * (((cumsum_Zi_wi_wji_A / cumsum_wi_wji_A_squared) *
        (cumsum_Zm_wi_wji_A / cumsum_wi_wji_A_squared)) -
        ((cumsum_ZiZm_wi_wji_A / cumsum_wi_wji_A_squared) *
          (cumsum_wi_wji_A / cumsum_wi_wji_A_squared))),
    ) |>
    dplyr::summarise(unpenalisedll = sum(li), .groups = "drop")

  ll_unpenalised |>
    dplyr::inner_join(penaltyMatrix, by = c("Zi", "Zm")) |>
    dplyr::mutate(penalisedll = unpenalisedll + penalty) |>
    dplyr::select(Zi, Zm, ll = penalisedll)

}



#' Wang 2019 likelihood differentiated twice with respect to the fixed effects, \eqn{\beta}, and then the random effects \b (and vice versa, as it is equivalent.)
#'
#' @param parms numeric vector of \eqn{\beta} and b values i.e \code{c(beta, b)}. The condition \code{length(parms) == length(c(X, Z))} must be satisfied.
#' @param X character vector containing names of X columns in data.
#' @param Z character vector containing names of Z columns in data.
#' @param t column in data with time of failure/censoring data.
#' @param wi column in data containing cluster-level weights.
#' @param wji column in data containing within-cluster weights.
#' @param dij column in data indicating if observation was censored or a failure.
#' @param D diagonal matrix with given \eqn{\theta}. Must be a square matrix of order \code{length(Z)}.
#' @param data tibble with all these columns



dlp_b_beta <- function(parms, X, Z, t, wi, wji, dij, data) {
  require(magrittr)

  sortedIndexedData <- sortAndIndex(data, {{ t }})

  terms1 <- calcLinearPredictor(sortedIndexedData, X = X, Z = Z, parms = parms)

  termsX <- calcRiskSets(terms1, X, "Xr")

  termsZ <- calcRiskSets(terms1, Z, "Zi")

  XZCrossProducts <- calcCrossProducts(sortedIndexedData, Z, X, "Zi", "Xr")

  termsZX <- dplyr::inner_join(
    termsZ %>%
      dplyr::select(
        index, {{ wi }}, {{ wji }}, wi_wji_A, d, Zi, cumsum_Zi_wi_wji_A,
        cumsum_wi_wji_A, cumsum_wi_wji_A_squared
      ),
    termsX %>%
      dplyr::select(index, Xr, cumsum_Xr_wi_wji_A),
    by = "index"
  )

  addCrossProducts <- dplyr::inner_join(
    termsZX,
    XZCrossProducts,
    by = c("index", "Xr", "Zi")
  )

  calcFinalTerms <- addCrossProducts %>%
    dplyr::arrange(dplyr::desc(index)) %>%
    dplyr::group_by(Xr, Zi) %>%
    dplyr::mutate(
      ZiXr_wi_wji_A = ZiXr * wi_wji_A,
      cumsum_ZiXr_wi_wji_A = cumsum(ZiXr_wi_wji_A)
    ) %>%
    dplyr::arrange(index) %>%
    dplyr::group_by(Xr, Zi)

  ll <- calcFinalTerms %>%
    dplyr::mutate(
      li = {{ wi }} * {{ wji }} * {{ dij }} *
        (((cumsum_Zi_wi_wji_A / cumsum_wi_wji_A_squared) * (cumsum_Xr_wi_wji_A / cumsum_wi_wji_A_squared))
        - ((cumsum_ZiXr_wi_wji_A / cumsum_wi_wji_A_squared) * (cumsum_wi_wji_A / cumsum_wi_wji_A_squared))),
    ) %>%
    dplyr::summarise(ll = sum(li), .groups = "drop")

  ll
}

# dlp_b_beta(parms = c(B, b), X = xColumns, Z = zColumns, t = t, wi = wi,
#            wji = wji, dij = dij, data = myData)
#



# calculate theta hat given beta and b.
