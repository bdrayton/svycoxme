

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



BB <- function(parms, X, t, cluster, dij, data){

  Z_formula <- formula(glue::glue(" ~ as.factor( {cluster} ) - 1"))

  Z_matrix <- model.matrix(Z_formula, data = data)

  Z_names <- paste0("Z", seq(ncol(Z_matrix)))

  colnames(Z_matrix) <- Z_names

  data_with_Z <- dplyr::bind_cols(data, data.frame(Z_matrix))

  b <- parms[-seq_len(length.out = length(X))]

  B <- parms[seq_len(length.out = length(X))]

  sortedIndexedData <- sortAndIndex(data_with_Z, {{ t }})

  addedLP <- calcLinearPredictor(sortedIndexedData, X = X, Z = Z_names, parms = parms)

  termsX <- calcRiskSets(addedLP, X, "Xr")

  XCrossProducts <- calcCrossProducts(sortedIndexedData, X, X, "Xr", "Xs")

  d6 <- dplyr::left_join(termsX, XCrossProducts, by = c("index", "Xr"))

  d7 <- d6 %>%
    dplyr::arrange(desc(index)) %>%
    dplyr::group_by(Xr, Xs) %>%
    dplyr::mutate(XrXs_A = XrXs * A,
                  cumsum_XrXs_A = cumsum(XrXs_A)) %>%
    dplyr::arrange(index) %>%
    dplyr::group_by(Xr, Xs)

  ll <- d7 %>%
    dplyr::mutate(
      ll_parts = stat * ((cumsum_Xr_A^2)/(cumsum_A)^2 - cumsum_XrXs_A/cumsum_A)
    ) %>%
    dplyr::summarise(ll = sum(ll_parts), .groups = "drop")

  ll

}

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

bb <- function(parms, X, t, cluster, dij, data, theta){
  Z_formula <- formula(glue::glue(" ~ as.factor( {cluster} ) - 1"))

  Z_matrix <- model.matrix(Z_formula, data = data)

  Z_names <- paste0("Z", seq(ncol(Z_matrix)))

  colnames(Z_matrix) <- Z_names

  data_with_Z <- dplyr::bind_cols(data, data.frame(Z_matrix))

  b <- parms[-seq_len(length.out = length(X))]

  D = theta * diag(length(b))

  B <- parms[seq_len(length.out = length(X))]

  sortedIndexedData <- sortAndIndex(data_with_Z, {{ t }})

  addedLP <- calcLinearPredictor(sortedIndexedData, X = X, Z = Z_names, parms = parms)

  d4 <- calcRiskSets(addedLP, vars = Z_names, varCol = "Zr")

  d5 <- calcCrossProducts(sortedIndexedData, Z_names, Z_names, "Zr", "Zs")

  d6 <- dplyr::left_join(d4, d5, by = c("index", "Zr"))

  d7 <- d6 %>%
    dplyr::arrange(desc(index)) %>%
    dplyr::group_by(Zr, Zs) %>%
    dplyr::mutate(ZrZs_A = ZrZs * A,
                  cumsum_ZrZs_A = cumsum(ZrZs_A)) %>%
    dplyr::arrange(index) %>%
    dplyr::group_by(Zr, Zs)

  unpenalised <- d7 %>%
    dplyr::mutate(
      ll_parts = stat * ((cumsum_Zr_A^2)/(cumsum_A)^2 - cumsum_ZrZs_A/cumsum_A)
    ) %>%
    dplyr::summarise(ll_unpenalised = sum(ll_parts), .groups = "drop")

  penalty = as.data.frame(solve(D))

  names(penalty) <- Z_names

  penalty_long_df <- penalty %>%
    dplyr::mutate(Zr = tidyselect::all_of(Z_names), .before = everything()) %>%
    tidyr::pivot_longer(cols = Z_names, names_to = "Zs", values_to = "penalty")

  ll <- dplyr::left_join(unpenalised, penalty_long_df, by = c("Zr", "Zs")) %>%
    dplyr::mutate(ll = ll_unpenalised - penalty) %>%
    dplyr::select(Zr, Zs, ll)

  ll

}

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

Bb <- function(parms, X, t, cluster, dij, data){

  Z_formula <- formula(glue::glue(" ~ as.factor( {cluster} ) - 1"))

  Z_matrix <- model.matrix(Z_formula, data = data)

  Z_names <- paste0("Z", seq(ncol(Z_matrix)))

  colnames(Z_matrix) <- Z_names

  data_with_Z <- dplyr::bind_cols(data, data.frame(Z_matrix))

  b <- parms[-seq_len(length.out = length(X))]

  B <- parms[seq_len(length.out = length(X))]

  sortedIndexedData <- sortAndIndex(data_with_Z, {{ t }})

  addedLP <- calcLinearPredictor(sortedIndexedData, X = X, Z = Z_names, parms = parms)

  d4x <- calcRiskSets(addedLP, vars = X, varCol = "Xr")

  d4z <- calcRiskSets(addedLP, vars = Z_names, varCol = "Zr")

  d4 <- dplyr::left_join(
    d4x %>% dplyr::select(index, stat, A, cumsum_A, Xr, cumsum_Xr_A),
    d4z %>% dplyr::select(index, Zr, cumsum_Zr_A), by = "index")

  d5 <- calcCrossProducts(sortedIndexedData, X, Z_names, "Xr", "Zr")

  d6 <- dplyr::left_join(d4, d5, by = c("index", "Xr", "Zr"))

  d7 <- d6 %>%
    dplyr::arrange(desc(index)) %>%
    dplyr::group_by(Xr, Zr) %>%
    dplyr::mutate(XrZr_A = XrZr * A,
                  cumsum_XrZr_A = cumsum(XrZr_A)) %>%
    dplyr::arrange(index) %>%
    dplyr::group_by(Xr, Zr)

  ll <- d7 %>%
    dplyr::mutate(
      ll_parts = stat * ((cumsum_Xr_A * cumsum_Zr_A)/(cumsum_A)^2 - cumsum_XrZr_A/cumsum_A)) %>%
    dplyr::summarise(ll = sum(ll_parts), .groups = "drop")

  ll

}

#' Calculate the PPL hessian
#'
#' @export
#'

ppl_hessian <- function(parms, X, t, cluster, dij, data, theta){


  H11_long <- BB(parms = parms, X = X, t = {{ t }}, cluster = "M",
                 dij = {{ dij }}, data = data) %>%
              dplyr::rename(rows = Xr, cols = Xs)

  H22_long <- bb(parms = parms, X = X, t = {{ t }}, cluster = "M",
                 dij = {{ dij }}, data = data, theta = theta) %>%
              dplyr::rename(rows = Zr, cols = Zs)

  H21_long <- Bb(parms = parms, X = X, t = {{ t }}, cluster = "M",
                 dij = {{ dij }}, data = data) %>%
              dplyr::rename(rows = Xr, cols = Zr)

  H12_long <- dplyr::rename(H21_long, newrows = cols, newcols = rows) %>%
              dplyr::rename(rows = newrows, cols = newcols)

  H_long <- dplyr::bind_rows(H11_long, H12_long, H21_long, H22_long) %>%
    tidyr::pivot_wider(names_from = cols, values_from = ll) %>%
    tibble::column_to_rownames(var = "rows") %>%
    as.matrix()

  H_long

}


















