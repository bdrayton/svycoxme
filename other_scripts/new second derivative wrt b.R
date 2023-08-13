

# Matrix version of second derviative wrt b. Used to estimate theta.

ds <- one_dataset(~X1 + X2 + X3 + (1 | M1),
                  dists = list(X1 = ~rnorm(n),
                               X2 = ~rnorm(k * nk),
                               X3 = ~rbinom(n, 1, 0.5),
                               M1 = ~rep(1:k, each = nk),
                               error = ~rexp(n, 10),
                               stat = ~sample(rep(c(0, 1), round(n * c(0.2, 0.8))), n)),
                  dist_args = list(k = 10, nk = 4, n = 10*4),
                  coefficients = c(1, 1, 1),
                  random_effect_variance = list(M1 = 0.5)
)

beta = c(1, 1, 1)
b <- c(rnorm(10, mean = 0, sd = 1))
theta <- 1

lp_hess_bb <- function(params, formula, data, theta) {

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

  at_risk_Z <- fast_risk_sets(exp_risk_score_Z)

  # This is slow, but it does what I want. Could be a target for speed increases.
  at_risk_ZtZ <- lapply(seq_along(ZtZ_exp_risk_score), function(start){
    Reduce("+", ZtZ_exp_risk_score[start:Zt_ncol])
  })

  term1 <- lapply(seq_len(nrow(at_risk_Z)), function(i){
    crossprod(at_risk_Z[i, , drop = FALSE] ) / at_risk[i]^2
  })

  term2 <- lapply(seq_along(at_risk_ZtZ), function(i){
    at_risk_ZtZ[[i]] / at_risk[i]
  })

  ll_part <- function(m1, m2, stat){
    stat * (m1 - m2)
  }

  ll_parts <- mapply(ll_part, stat = stat@x, m1 = term1, m2 = term2, SIMPLIFY = FALSE)

  bb_unpenalised <- Reduce("+", ll_parts)

  parsed_data$reTrms$Lambdat@x <- theta[parsed_data$reTrms$Lind]

  penalty <- solve(parsed_data$reTrms$Lambdat)

  bb = bb_unpenalised - penalty

  bb

}

lp_hess_bb(params = c(beta, b), survival::Surv(stat_time, stat) ~ X1 + X2 + X3 + (1 | M1), data = ds, theta = theta)

my_beta = c(1, 1, 1)
my_theta = 1
my_k = 10
my_nk = 4

my_X = c("X1", "X2", "X3")

my_params <- c(my_beta[seq_along(my_X)], b)

bb(parms = my_params, X = my_X, stat_time =  stat_time, cluster = "M1", dij = stat,
   data = ds, theta = my_theta, return_matrix = TRUE)

microbenchmark::microbenchmark(
  lp_hess_bb(params = c(beta, b), survival::Surv(stat_time, stat) ~ X1 + X2 + X3 + (1 | M1), data = ds, theta = theta)
  ,
  bb(parms = my_params, X = my_X, stat_time =  stat_time, cluster = "M1", dij = stat,
     data = ds, theta = my_theta, return_matrix = TRUE)
)


# The ZtZ matricies are hard.
# I need one per subject, which is lots. They should all be sparse. I think this is best served by using a list of sparse matricies.
# Hoping I can construct that from parsed_data.

params = c(beta, b)
formula = survival::Surv(stat_time, stat) ~ X1 + X2 + X3 + (1 | M1)
data = ds
theta = theta






  stat <- Matrix(ds_sorted$stat, ncol = 1)

  likelihood_gradients_unpenalised <- colSums(stat * (t(parsed_data$reTrms$Zt) - at_risk_Z/at_risk))

  penalty <- t(b) %*% solve(parsed_data$reTrms$Lambdat)

  # this accessing @x is probably poor practice?
  likelihood_gradients_unpenalised - penalty@x

parsed_data


M = matrix(c(1, 0, 1), nr = 1)

t(M) %*% M

Zt <- parsed_data$reTrms$Zt

Zt_ncol <- ncol(Zt)

ZtZ_exp_risk_score <- vector(mode = "list", length = Zt_ncol)

for (i in seq_len(Zt_ncol)) {

  ZtZ_exp_risk_score[[i]] <- tcrossprod(Zt[,i, drop = FALSE]) * exp_risk_score[i]

}

rm(i)

head(ZtZ_exp_risk_score)

tcrossprod(Zt[, 1, drop = FALSE])

# multiply each list matrix by the corresponding exponentiated risk score.

# This is slow, but it does what I want. Could be a target for speed increases.
at_risk_ZtZ <- lapply(seq_along(ZtZ_exp_risk_score), function(start){
  Reduce("+", ZtZ_exp_risk_score[start:Zt_ncol])
})

term1 <- lapply(seq_len(nrow(at_risk_Z)), function(i){
  crossprod(at_risk_Z[i, , drop = FALSE] ) / at_risk[i]^2
})

term2 <- lapply(seq_along(at_risk_ZtZ), function(i){
  at_risk_ZtZ[[i]] / at_risk[i]
})

ll_part <- function(m1, m2, stat){
  stat * (m1 - m2)
}

ll_parts <- mapply(ll_part, stat = stat@x, m1 = term1, m2 = term2, SIMPLIFY = FALSE)

bb_unpenalised <- Reduce("+", ll_parts)

# pretty sure theta needs to be updated in this. Currently doesn't matter because theta is fixed at 1
penalty = solve(parsed_data$reTrms$Lambdat)

bb_matrix = bb_unpenalised - penalty





bb_matrix -bb(parms = my_params, X = my_X, stat_time =  stat_time, cluster = "M1", dij = stat,
   data = ds, theta = my_theta, return_matrix = TRUE)




# ok which parts are wrong?

# save ds somewhere so I can load it and run bb in a separate session.
readr::write_rds(ds, file = "C:/Users/Bradley/OneDrive - The University of Auckland/PhD/data/ds.rds")


# risk sets are all the same
# penalty is the same.
# term2 is the same.










