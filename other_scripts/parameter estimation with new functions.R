
#  rewrite theta ipl to use bb.
# write a theta ipl with nested optimisation of beta and b.

new_lp <- function(params, formula, data, theta) {

  # This sort var could be pulled out of the formula.
  ds_sorted <- sortAndIndex(data, sort_vars = stat_time)

  parsed_data <- lme4::lFormula(formula, data = ds_sorted)

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

  # would be good to peel stat out of the Surv() obeject in the parsed formula...
  stat <- Matrix::Matrix(ds_sorted$stat, ncol = 1)

  # penalty

  parsed_data$reTrms$Lambdat@x <- theta[parsed_data$reTrms$Lind]

  penalty <- 0.5 * t(b) %*% Matrix::solve(parsed_data$reTrms$Lambdat) %*% b

  penalised_likelihood <- sum(stat * (risk_score - log(at_risk))) - penalty

  penalised_likelihood@x

}




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

  ## I replaced this with fast_risk_sets, which is faster when there are lots of columns. Else it's about the same.
  #
  # at_risk_list <- apply(exp_risk_score_X, 2, function(column){
  #   Matrix(rev(cumsum(rev(column))), ncol = 1)
  # })
  #
  # at_risk_X <- Reduce(cbind, at_risk_list)

  at_risk_X <- fast_risk_sets(exp_risk_score_X)

  stat <- Matrix(ds_sorted$stat, ncol = 1)

  likelihood_gradients <- colSums(stat * (parsed_data$X[, -1] - at_risk_X/at_risk))

  likelihood_gradients

}


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

  # at_risk_list <- apply(exp_risk_score_Z, 2, function(column){
  #   Matrix(rev(cumsum(rev(column))), ncol = 1)
  # })
  #
  # at_risk_Z <- Reduce(cbind, at_risk_list)

  at_risk_Z <- fast_risk_sets(exp_risk_score_Z)

  stat <- Matrix(ds_sorted$stat, ncol = 1)

  likelihood_gradients_unpenalised <- colSums(stat * (t(parsed_data$reTrms$Zt) - at_risk_Z/at_risk))

  parsed_data$reTrms$Lambdat@x <- theta[parsed_data$reTrms$Lind]

  penalty <- t(b) %*% solve(parsed_data$reTrms$Lambdat)

  # this accessing @x is probably poor practice?
  likelihood_gradients_unpenalised - penalty@x

}

lp_gr <- function(params, formula, data, theta) {

  c(lp_gr_beta(params = params, formula = formula, data = data, theta = theta),
    lp_gr_b(params = params, formula = formula, data = data, theta = theta))

}


new_theta_ipl <- function(theta, formula, data){

  parsed_data <- lme4::lFormula(formula, data = data)

  # set up D
  parsed_data$reTrms$Lambdat@x <- theta[parsed_data$reTrms$Lind]

  D <- parsed_data$reTrms$Lambdat

  n_fixed <- length(attr(terms(coxme:::formula1(formula)$fixed), "order"))

  fixed_formula <- lme4:::getFixedFormula(formula)

  fit0 <- survival::coxph(fixed_formula, data = data)

  start_params <- c(coef(fit0), rep(0, length(parsed_data$reTrms$Lind)))

  ests <- optim(par = start_params,
                fn = new_lp,
                gr = lp_gr,
                formula = formula,
                data = data,
                theta = theta,
                method = "BFGS",
                control = list(fnscale = -1),
                hessian = TRUE)

  # take the part of the hessian associated with the random effects only
  Kbb <- ests$hessian[-(1:n_fixed), -(1:n_fixed)]

  -0.5 * log(det(D)) -0.5 * log(det(Kbb)) + ests$value

}



### end of functions ###

my_k = 10
my_nk = 10
my_n = my_k * my_nk
my_beta = c(1, -0.7, 0.5)
my_theta = 2

fit0 <- survival::coxph(survival::Surv(stat_time, stat) ~ X1 + X2 + X3, data = ds)
my_start_parameters <- c(coef(fit0), rep(0, my_k))

ds <- one_dataset(~X1 + X2 + X3 + (1 | M1),
                  dists = list(X1 = ~rnorm(n),
                               X2 = ~rnorm(k * nk),
                               X3 = ~rbinom(n, 1, 0.5),
                               M1 = ~rep(1:k, each = nk),
                               error = ~rexp(n, 10),
                               stat = ~sample(rep(c(0, 1), round(n * c(0.2, 0.8))), n)),
                  dist_args = list(k = my_k, nk = my_nk, n = my_n),
                  coefficients = my_beta,
                  random_effect_variance = list(M1 = my_theta[1])
)



new_theta_ipl(theta = my_theta,
              formula = survival::Surv(stat_time, stat) ~ X1 + X2 + X3 + (1 | M1),
              start_params = my_start_parameters,
              data = ds)

theta_est <- optim(par = 0.5,
                   fn = new_theta_ipl,
                   gr = NULL,
                   formula = survival::Surv(stat_time, stat) ~ X1 + X2 + X3 + (1 | M1),
                   start_params = my_start_parameters,
                   data = ds,
                   method = "L-BFGS-B",
                   control = list(fnscale = -1),
                   lower = 0.00001,
                   upper = 3)

coxfit <- coxme::coxme(survival::Surv(stat_time, stat) ~ X1 + X2 + X3 + (1 | M1), data = ds)

coxme::VarCorr(coxfit)

theta_est$par

# try two random effects.

my_k = 10
my_nk = 10
my_n = my_k * my_nk
my_beta = c(1, -0.7, 0.5)
my_theta = c(2, 1)

fit0 <- survival::coxph(survival::Surv(stat_time, stat) ~ X1 + X2 + X3, data = ds)
my_start_parameters <- c(coef(fit0), rep(0, my_k + 2))

ds <- one_dataset(~X1 + X2 + X3 + (1 | M1) + (1 | M2),
                  dists = list(X1 = ~rnorm(n),
                               X2 = ~rnorm(k * nk),
                               X3 = ~rbinom(n, 1, 0.5),
                               M1 = ~rep(1:k, each = nk),
                               M2 = ~rep(c("l", "r"), ceiling(n/2))[seq_len(n)],
                               error = ~rexp(n, 10),
                               stat = ~sample(rep(c(0, 1), round(n * c(0.2, 0.8))), n)),
                  dist_args = list(k = my_k, nk = my_nk, n = my_n),
                  coefficients = my_beta,
                  random_effect_variance = list(M1 = my_theta[1],
                                                M2 = my_theta[2])
)

new_lp(theta = my_theta,
       formula = survival::Surv(stat_time, stat) ~ X1 + X2 + X3 + (1 | M1) + (1 | M2),
       params = my_start_parameters,
       data = ds)


new_theta_ipl(theta = my_theta,
              formula = survival::Surv(stat_time, stat) ~ X1 + X2 + X3 + (1 | M1) + (1 | M2),
              start_params = my_start_parameters,
              data = ds)

theta_est <- optim(par = c(1, 1),
                   fn = new_theta_ipl,
                   gr = NULL,
                   formula = survival::Surv(stat_time, stat) ~ X1 + X2 + X3 + (1 | M1) + (1 | M2),
                   start_params = my_start_parameters,
                   data = ds,
                   method = "L-BFGS-B",
                   control = list(fnscale = -1),
                   lower = 0.00001,
                   upper = 3)

coxfit <- coxme::coxme(survival::Surv(stat_time, stat) ~ X1 + X2 + X3 + (1 | M1) + (1 | M2), data = ds)

theta_est$par

coxme::VarCorr(coxfit)

# try two random effects, with one nested.

my_k = 10
# should revise the my_nk to account for more than one random effect variable, or intercepts and slopes
my_nk = 10
my_n = my_k * my_nk
my_beta = c(1, -0.7, 0.5)
my_theta = c(2, 1)

ds <- one_dataset(~X1 + X2 + X3 + (1 | M1) + (1 | M2),
                  dists = list(X1 = ~rnorm(n),
                               X2 = ~rnorm(k * nk),
                               X3 = ~rbinom(n, 1, 0.5),
                               M1 = ~rep(1:k, each = nk),
                               M2 = ~rep(c("l", "r"), ceiling(n/2))[seq_len(n)],
                               error = ~rexp(n, 10),
                               stat = ~sample(rep(c(0, 1), round(n * c(0.2, 0.8))), n)),
                  dist_args = list(k = my_k, nk = my_nk, n = my_n),
                  coefficients = my_beta,
                  random_effect_variance = list(M1 = my_theta[1],
                                                M2 = my_theta[2])
)

# need 30 random effect terms
# lme4::lFormula(survival::Surv(stat_time, stat) ~ X1 + X2 + X3 + (1 | M1) + (1 | M1:M2), data = ds)

fit0 <- survival::coxph(survival::Surv(stat_time, stat) ~ X1 + X2 + X3, data = ds)
my_start_parameters <- c(coef(fit0), rep(0, 30))

new_lp(theta = my_theta,
       formula = survival::Surv(stat_time, stat) ~ X1 + X2 + X3 + (1 | M1) + (1 | M1:M2),
       params = my_start_parameters,
       data = ds)


new_theta_ipl(theta = my_theta,
              formula = survival::Surv(stat_time, stat) ~ X1 + X2 + X3 + (1 | M1) + (1 | M1:M2),
              start_params = my_start_parameters,
              data = ds)

theta_est <- optim(par = c(1, 1),
                   fn = new_theta_ipl,
                   gr = NULL,
                   formula = survival::Surv(stat_time, stat) ~ X1 + X2 + X3 + (1 | M1) + (1 | M1:M2),
                   start_params = my_start_parameters,
                   data = ds,
                   method = "L-BFGS-B",
                   control = list(fnscale = -1),
                   lower = 0.00001,
                   upper = 3)

# coxme doesn't work for this. Mine works, but looks poor. a proper simulation is needed to evaluate it.
# coxfit <- coxme::coxme(survival::Surv(stat_time, stat) ~ X1 + X2 + X3 + (1 | M1) + (1 | M1:M2), data = ds)

theta_est$par

# coxme::VarCorr(coxfit)
