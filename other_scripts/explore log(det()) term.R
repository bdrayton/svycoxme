

# There is a term in the full laplace approximation of the log likelihood.
# log(det(sum(cumulative_hazard_i * risk_score_i * Z_it(Z_i) - solve(D) )
# if the determinant is negative, the log is undefined.
# here I test if it is possible to get negative determinants by testing
# combinations of b, beta, and theta


theta = 0.5
beta = c(1, -1, 1)

nk = 10
k = 50
n = nk * k

formula <- survival::Surv(stat_time, stat) ~ X1 + X2 + X3 + (1 | M)


one_rep <- function(b_zero = FALSE) {

ds <- one_dataset(formula,
            dists = list(X1 = ~rnorm(n),
                         X2 = ~rnorm(k * nk),
                         X3 = ~rbinom(n, 1, 0.5),
                         M = ~rep(1:k, each = nk),
                         error = ~rexp(n, 10),
                         stat = ~sample(rep(c(0, 1), round(n * c(0.2, 0.8))), n)),
            dist_args = list(k = k, nk = nk, n = n),
            coefficients = beta,
            random_effect_variance = list(M = theta)
)

# make the Z %*% t(Z) matricies

ds_sorted <- sortAndIndex(ds, sort_vars = stat_time)

parsed_data <- lme4::lFormula(formula, data = ds_sorted)

nc <- ncol(parsed_data$reTrms$Zt)

ZZt <- lapply(seq_len(nc),
              function(i){
                Matrix::tcrossprod(parsed_data$reTrms$Zt[, i, drop = FALSE])
              })

b <- dplyr::bind_cols(cluster = ds$M, attr(ds, "random_effects")) |>
  dplyr::distinct() |> dplyr::pull(M)

# what if all b = 0, i.e. no random effects, or starting guess for random effects
# b <- rep(0.01, length(b)) * sign(runif(length(b), min = -1, max = 1))

if(b_zero) {

  b <- rep(0.0, length(b))

}


# drop the intercept column from the X model.matrix
risk_score <- parsed_data$X[, -1] %*% beta + Matrix::crossprod(parsed_data$reTrms$Zt, b)

exp_risk_score <- exp(risk_score)

rev_exp_risk_score <- exp_risk_score
rev_exp_risk_score@x <- rev(exp_risk_score@x)

at_risk <- Matrix::Matrix(rev(cumsum(rev_exp_risk_score)), ncol = 1)

parsed_data$reTrms$Lambdat@x <- theta[parsed_data$reTrms$Lind]

cumhaz <- cumsum(ds_sorted$stat / at_risk)

# plot(ds_sorted$stat_time, exp(-cumhaz))


D_inv <- Matrix::solve(parsed_data$reTrms$Lambdat)

# calculate NA cumhaz

# NA_cumhaz <- ds_sorted|>
#   dplyr::mutate(at_risk = dplyr::n() + 1 - dplyr::row_number(),
#                 NA_i = stat/at_risk,
#                 NA_cumhaz = cumsum(NA_i)) |>
#   dplyr::pull(NA_cumhaz)
#
# plot(cumhaz[-(400:500)], NA_cumhaz[-(400:500)])


term2_parts <- lapply(seq_along(ZZt), function(i){

  cumhaz[i] * exp_risk_score[i] * ZZt[[i]]

})

log(det(Reduce("+", term2_parts) - D_inv))

}

# seems to fail for "bad" guesses of random effects. won't fail at true param.
# i think it fails when all(b == 0), which is when theta = 0. this likelihood is
# undefined at that point.
correct_b <- replicate(1000, one_rep())

table(is.nan(correct_b))

# sometimes fails (15/1000) at the true values of beta and b.
zero_b <- replicate(1000, one_rep(b_zero = TRUE))

# about half are NAN.
table(is.nan(zero_b))















































