


lp_gr_beta <- function(params, formula, data, theta) {

  ds_sorted <- sortAndIndex(data, sort_vars = stat_time)

  parsed_data <- lme4::lFormula(formula, data = ds_sorted)

  beta_index <- seq_len(ncol(parsed_data$X) - 1)

  beta <- params[beta_index]
  b <- params[-beta_index]

  risk_score <- parsed_data$X[, -1] %*% beta + crossprod(parsed_data$reTrms$Zt, b)

  exp_risk_score <- exp(risk_score)

  rev_exp_risk_score <- exp_risk_score
  rev_exp_risk_score@x <- rev(exp_risk_score@x)

  at_risk <- Matrix(rev(cumsum(rev_exp_risk_score)), ncol = 1)

  exp_risk_score_X <- exp_risk_score * parsed_data$X[,-1]

  at_risk_list <- apply(exp_risk_score_X, 2, function(column){
    Matrix(rev(cumsum(rev(column))), ncol = 1)
  })

  at_risk_X <- Reduce(cbind, at_risk_list)

  stat <- Matrix(ds_sorted$stat, ncol = 1)

  likelihood_gradients <- colSums(stat * (parsed_data$X[, -1] - at_risk_X/at_risk))

  likelihood_gradients

}

ds <- one_dataset(~X1 + X2 + X3 + (1 | M1),
                  dists = list(X1 = ~rnorm(n),
                               X2 = ~rnorm(k * nk),
                               X3 = ~rbinom(n, 1, 0.5),
                               M1 = ~rep(1:k, each = nk),
                               error = ~rexp(n, 10),
                               stat = ~sample(rep(c(0, 1), round(n * c(0.2, 0.8))), n)),
                  dist_args = list(k = 10, nk = 4, n = 200),
                  coefficients = c(1, 1, 1),
                  random_effect_variance = list(M1 = 0.5)
)

beta = c(1, 1, 1)
b <- c(rnorm(10, mean = 0, sd = 1))
theta <- 1

dlp_beta(parms = c(beta, b),
         X = c("X1", "X2", "X3"),
         stat_time = stat_time,
         dij = stat,
         theta = theta,
         cluster = "M1", data = ds)

lp_gr_beta(params = c(beta, b), survival::Surv(stat_time, stat) ~ X1 + X2 + X3 + (1 | M1), data = ds, theta = theta)


microbenchmark::microbenchmark(
  lp(parms = c(beta, b),
     X = c("X1", "X2", "X3"), stat_time = stat_time, dij = stat, theta = theta, cluster = "M1", data = ds),
  new_lp(params = c(beta, b), survival::Surv(stat_time, stat) ~ X1 + X2 + X3 + (1 | M1), data = ds, theta = theta)
)


dlp_beta <- function(parms, X, cluster, stat_time, dij, data, ...) {

  data_with_Z <- add_Z(data, cluster)

  Z_names <- attr(data_with_Z, "Z_names")

  b <- parms[-seq_len(length.out = length(X))]

  B <- parms[seq_len(length.out = length(X))]

  sortedIndexedData <- sortAndIndex(data_with_Z, sort_vars = {{ stat_time }})

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


dlp_b <- function(parms, X, cluster, stat_time, dij, D, data, ...) {

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

  sortedIndexedData <- sortAndIndex(data_with_Z, sort_vars = {{ stat_time }})

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
