


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

microbenchmark::microbenchmark(
dlp_beta(parms = c(beta, b),
         X = c("X1", "X2", "X3"),
         stat_time = stat_time,
         dij = stat,
         theta = theta,
         cluster = "M1", data = ds),
lp_gr_beta(params = c(beta, b), survival::Surv(stat_time, stat) ~ X1 + X2 + X3 + (1 | M1), data = ds, theta = theta)
)


solve(parsed_data$reTrms$Lambdat) %*% b


lp_gr_b <- function(params, formula, data, theta) {

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

  bl <- length(b)

  exp_risk_score_Z <- Matrix(rep(exp_risk_score, bl), ncol = bl) * t(parsed_data$reTrms$Zt)

  at_risk_list <- apply(exp_risk_score_Z, 2, function(column){
    Matrix(rev(cumsum(rev(column))), ncol = 1)
  })

  at_risk_Z <- Reduce(cbind, at_risk_list)

  stat <- Matrix(ds_sorted$stat, ncol = 1)

  likelihood_gradients_unpenalised <- colSums(stat * (t(parsed_data$reTrms$Zt) - at_risk_Z/at_risk))

  penalty <- t(b) %*% solve(parsed_data$reTrms$Lambdat)

  # this accessing @x is probably poor practice?
  likelihood_gradients_unpenalised - penalty@x

}


my_params <- c(beta, b)
names(my_params) <- c(c("X1", "X2", "X3"), paste0("Z", seq_len(10)))

dlp_b(parms = my_params,
         X = c("X1", "X2", "X3"),
         stat_time = stat_time,
         dij = stat,
         theta = theta,
         D = theta * diag(10),
         cluster = "M1", data = ds)

lp_gr_b(params = c(beta, b),
        survival::Surv(stat_time, stat) ~ X1 + X2 + X3 + (1 | M1),
        data = ds, theta = theta)

microbenchmark::microbenchmark(
dlp_b(parms = my_params,
      X = c("X1", "X2", "X3"),
      stat_time = stat_time,
      dij = stat,
      theta = theta,
      D = theta * diag(10),
      cluster = "M1", data = ds)
,
lp_gr_b(params = c(beta, b),
        survival::Surv(stat_time, stat) ~ X1 + X2 + X3 + (1 | M1),
        data = ds, theta = theta)
)


lp_gr <- function(params, formula, data, theta) {

  c(lp_gr_beta(params = params, formula = formula, data = data, theta = theta),
    lp_gr_b(params = params, formula = formula, data = data, theta = theta))

}

lp_gr(params = c(beta, b),
        survival::Surv(stat_time, stat) ~ X1 + X2 + X3 + (1 | M1),
        data = ds, theta = theta)






