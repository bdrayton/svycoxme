


my_formula <- survival::Surv(stat_time, stat)~X1 + X2 + X3 + (1 | M1) + (1 | M2)
my_k = 50
my_nk = 10
my_theta = c(M1 = 2, M2 = 1)
my_beta = c(1, -0.7, 0.5)

ds <- one_dataset(my_formula,
                  dists = list(X1 = ~rnorm(n),
                               X2 = ~rnorm(k * nk),
                               X3 = ~rbinom(n, 1, 0.5),
                               M1 = ~rep(1:k, each = nk),
                               M2 = ~rep(c("l", "r"), ceiling(n/2))[seq_len(n)],
                               error = ~rexp(n, 10),
                               stat = ~sample(rep(c(0, 1), round(n * c(0.2, 0.8))), n)),
                  dist_args = list(k = my_k, nk = my_nk, n = my_k * my_nk),
                  coefficients = my_beta,
                  random_effect_variance = my_theta
)

coxme_fit <- coxme::coxme(my_formula, data = ds)

ests <- est_parameters(my_formula, ds, control = control.list(grad = TRUE,
                                                              ndeps = c(1e-2, 1e-2)))

ds_sorted <- sortAndIndex(ds, sort_vars = stat_time)

parsed_data <- lme4::lFormula(my_formula, data = ds_sorted)

# set up D
D <- parsed_data$reTrms$Lambdat

D@x <- ests$theta[parsed_data$reTrms$Lind]

D_inv <- solve(D)

# the derivatives of D
d_D_d_theta <- vector(mode = "list", length = length(parsed_data$reTrms$theta))

for (i in seq_along(parsed_data$reTrms$theta)) {

  d_D_d_theta[[i]] <- parsed_data$reTrms$Lambdat

  d_D_d_theta[[i]]@x <- c(1, 0)[(parsed_data$reTrms$Lind != i) + 1]

}

# In Ripatti Palmgren they only use this part of the Hessian matrix.
Kbb <- attr(ests, "beta_b_est")$hessian[-(1:n_fixed), -(1:n_fixed)]

b <- Matrix(ests$b, ncol = 1)

K_inv <- solve(Kbb)

term3_matrix <- term2_matrix <- term1_matrix <- matrix(NA, nrow = length(ests$theta), ncol = length(ests$theta))

for (i in seq_along(ests$theta)) {

  for (j in seq_along(ests$theta)) {

    term1_matrix[i, j] <- -sum(diag(D_inv %*% d_D_d_theta[[j]] %*% D_inv %*% d_D_d_theta[[i]]))

    term2_matrix[i, j] <- sum(diag(K_inv %*% D_inv %*% d_D_d_theta[[j]] %*% D_inv %*% K_inv %*% D_inv %*% d_D_d_theta[[i]] %*% D_inv
                                   -K_inv %*% (D_inv %*% d_D_d_theta[[j]] %*% D_inv %*% d_D_d_theta[[i]] %*% D_inv +
                                                 D_inv %*% d_D_d_theta[[i]] %*% D_inv %*% d_D_d_theta[[j]] %*% D_inv)))

    term3_temp <- t(b) %*% (D_inv %*% d_D_d_theta[[j]] %*% D_inv %*% d_D_d_theta[[i]] %*% D_inv +
                              D_inv %*% d_D_d_theta[[i]] %*% D_inv %*% d_D_d_theta[[j]] %*% D_inv) %*% b

    term3_matrix[i, j] <- term3_temp@x

    rm(term3_temp)

  }

}

term1_matrix
term2_matrix
term3_matrix

diag(-solve(-0.5 * (-term1_matrix + term2_matrix + term3_matrix)))




