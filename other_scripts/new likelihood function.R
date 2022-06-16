

# Calculate the penalised partial likelihoods using sparse matrix operations.
# use lme4 formula parsing

library(lme4)

ds <- one_dataset(~X1 + X2 + X3 + (1 | M1) + (1| M2),
              dists = list(X1 = ~rnorm(n),
                           X2 = ~rnorm(k * nk),
                           X3 = ~rbinom(n, 1, 0.5),
                           M1 = ~rep(1:k, each = nk),
                           M2 = ~rep(c("l", "r"), ceiling(n/2))[seq_len(n)],
                           error = ~rexp(n, 10),
                           stat = ~sample(rep(c(0, 1), round(n * c(0.2, 0.8))), n)),
              dist_args = list(k = 10, nk = 4, n = 200),
              coefficients = c(1, 1, 1),
              random_effect_variance = list(M1 = 1, M2 = 2)
  )

# need to sort the data first, for cumsum and riskset calcs.

ds_sorted <- sortAndIndex(ds, sort_vars = stat_time)

# pass this through the parser

parsed_data <- lme4::lFormula(survival::Surv(stat_time, stat) ~ X1 + X2 + X3 + (1 | M1) + (1 | M2), data = ds_sorted)

parsed_data

# stuff to calculate:
# linear predictor, and exponentiated form. call risk_score and exp_risk_score
# X times these, Z times these. Later, for hessian, also need exp_risk_score %*% Z %*% t(Z), or something along those lines.
# random effects and penalty, which is where I start.

# in this case there are 12 random effects. 10 drawn from one dist and 2 from another.

b <- c(rnorm(10, mean = 0, sd = 1), rnorm(2, mean = 0, sd = sqrt(2)))

# This is D. There are some weird rules about copying this I think.

t(b) %*% solve(parsed_data$reTrms$Lambdat) %*% b


### what about nested random intercepts?
parsed_data <- lme4::lFormula(survival::Surv(stat_time, stat) ~ X1 + X2 + X3 + (1 | M1) + (1 | M1:M2), data = ds_sorted)

# in this case there are 30 random effects. 10 drawn from one dist and 20 from another.

b <- c(rnorm(10, mean = 0, sd = 1), rnorm(20, mean = 0, sd = sqrt(2)))
t(b) %*% solve(parsed_data$reTrms$Lambdat) %*% b


# risk score and exp(riskscore.)
# need coef for intercept, and each covariate (per levels - 1 for each factor covariate).

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

ds_sorted <- sortAndIndex(ds, sort_vars = stat_time)

beta = c(1, 1, 1)
b <- c(rnorm(10, mean = 0, sd = 1))

parsed_data <- lme4::lFormula(survival::Surv(stat_time, stat) ~ X1 + X2 + X3 + (1 | M1), data = ds_sorted)

# drop the intercept column from the X model.matrix
risk_score <- parsed_data$X[, -1] %*% beta + crossprod(parsed_data$reTrms$Zt, b)

exp_risk_score <- exp(risk_score)

data_with_Z <- add_Z(ds_sorted, "M1")

# drop M2

r1 <- calcLinearPredictor(data_with_Z, X = c("X1", "X2", "X3"), Z = attr(data_with_Z, "Z_names"), c(beta, b))

table(r1 == as.vector(risk_score))

cbind(r1$lp, risk_score@x, r1$lp - risk_score@x)

max(r1$lp - risk_score@x) < .Machine$double.eps*10

# calculate risk sets
# I do it this way to preserve the class and other slots. rev(risk_score) coerces to a numeric vector.
# probably no point, as cumsum convert it anyway, and I need to remake the matrix with Matrix.
rev_exp_risk_score <- exp_risk_score
rev_exp_risk_score@x <- rev(exp_risk_score@x)

at_risk <- Matrix(rev(cumsum(rev_exp_risk_score)), ncol = 1)

risk_sets <- calcRiskSets(data = r1)

max(risk_sets$cumsum_A - at_risk@x)

risk_sets_X <- calcRiskSets(r1, c("X1", "X2", "X3"), "Xr")

exp_risk_score_X <- exp_risk_score * parsed_data$X[,-1]

at_risk_list <- apply(exp_risk_score_X, 2, function(column){
  Matrix(rev(cumsum(column)), ncol = 1)
})

at_risk_X <- Reduce(cbind, at_risk_list)

# would be good to peel stat out of the Surv() obeject in the parsed formula...
stat <- Matrix(ds_sorted$stat, ncol = 1)

colSums(stat * (parsed_data$X[, -1] - at_risk_X/at_risk))


ll <- addedCumsums %>%
  dplyr::mutate(li =  {{ dij }} * (value - cumsum_Xr_A / cumsum_A)) %>%
  dplyr::summarise(ll = sum(li), .groups = "drop")

ll

# unpenalised likelihood
sum(stat * (risk_score - log(at_risk)))

# penalty

# need to update lambdat to current theta estimates. how does lind work with multiple thetas?

thetas <- 1

my_Lambdat <- parsed_data$reTrms$Lambdat

my_Lambdat@x <- thetas[parsed_data$reTrms$Lind]

penalty <- t(b) %*% solve(my_Lambdat) %*% b

# penalised likelihood
sum(stat * (risk_score - log(at_risk))) - penalty


lp(parms = c(beta, b),
   X = c("X1", "X2", "X3"), stat_time = stat_time, dij = stat, theta = thetas, cluster = "M1", data = ds)

new_lp <- function(params, formula, data, theta) {

  # This sort var could be pulled out of the formula.
  ds_sorted <- sortAndIndex(data, sort_vars = stat_time)

  parsed_data <- lme4::lFormula(formula, data = ds_sorted)

  beta_index <- seq_len(ncol(parsed_data$X) - 1)

  beta <- params[beta_index]
  b <- params[-beta_index]

  # drop the intercept column from the X model.matrix
  risk_score <- parsed_data$X[, -1] %*% beta + crossprod(parsed_data$reTrms$Zt, b)

  exp_risk_score <- exp(risk_score)

  # calculate risk sets
  # I do it this way to preserve the class and other slots. rev(risk_score) coerces to a numeric vector.
  # probably no point, as cumsum convert it anyway, and I need to remake the matrix with Matrix.
  rev_exp_risk_score <- exp_risk_score
  rev_exp_risk_score@x <- rev(exp_risk_score@x)

  at_risk <- Matrix(rev(cumsum(rev_exp_risk_score)), ncol = 1)

  # would be good to peel stat out of the Surv() obeject in the parsed formula...
  stat <- Matrix(ds_sorted$stat, ncol = 1)

  # penalty

  parsed_data$reTrms$Lambdat@x <- theta[parsed_data$reTrms$Lind]

  penalty <- 0.5 * t(b) %*% solve(parsed_data$reTrms$Lambdat) %*% b

  penalised_likelihood <- sum(stat * (risk_score - log(at_risk))) - penalty

  penalised_likelihood@x

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

new_lp(params = c(beta, b), survival::Surv(stat_time, stat) ~ X1 + X2 + X3 + (1 | M1), data = ds, theta = theta)

# In a speed test, the new function is ~4 times faster.

microbenchmark::microbenchmark(
  lp(parms = c(beta, b),
     X = c("X1", "X2", "X3"), stat_time = stat_time, dij = stat, theta = theta, cluster = "M1", data = ds),
  new_lp(params = c(beta, b), survival::Surv(stat_time, stat) ~ X1 + X2 + X3 + (1 | M1), data = ds, theta = theta)
)


fit0 <- survival::coxph(survival::Surv(stat_time, stat) ~ X1 + X2 + X3, data = ds)

my_start_parameters <- c(coef(fit0), rep(0, 10))
names(my_start_parameters) <- c(c("X1", "X2", "X3"), paste0("Z", seq_len(10)))

microbenchmark::microbenchmark(
  fit_beta_b <- optim(par = my_start_parameters,
                      fn = lp,
                      gr = lp_grd,
                      X = c("X1", "X2", "X3"),
                      stat_time = stat_time ,
                      cluster = "M1",
                      dij = stat,
                      theta = theta,
                      data = ds,
                      method = "BFGS",
                      control = list(fnscale = -1)),
  new_fit_beta_b <- optim(par = my_start_parameters,
                          fn = new_lp,
                          gr = NULL,
                          formula = survival::Surv(stat_time, stat) ~ X1 + X2 + X3 + (1 | M1),
                          data = ds,
                          theta = theta,
                          method = "BFGS",
                          control = list(fnscale = -1))
)

cbind(
  fit_beta_b$par,
  new_fit_beta_b$par)


##### Look at other random effects structures. does it work? compare to coxme.


