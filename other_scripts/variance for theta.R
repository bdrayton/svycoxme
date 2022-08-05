


# analytical hessian for theta

my_formula <- survival::Surv(stat_time, stat)~X1 + X2 + X3 + (1 | M1) + (1 | M2)
my_k = 10
my_nk = 4
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

ds_sorted <- sortAndIndex(ds, sort_vars = stat_time)

parsed_data <- lme4::lFormula(my_formula, data = ds_sorted)

stat <- Matrix(unclass(parsed_data$fr[,1])[, "status"], ncol = 1)

theta_start <- get_start_theta(length(parsed_data$reTrms$flist))

n_fixed <- length(attr(terms(coxme:::formula1(my_formula)$fixed), "order"))

fixed_formula <- lme4:::getFixedFormula(my_formula)

fit0 <- survival::coxph(fixed_formula, data = ds)


start_params <- c(coef(fit0),
                  rnorm(n = length(parsed_data$reTrms$Lind),
                        mean = 0,
                        sd = sqrt(theta_start[parsed_data$reTrms$Lind])))

control.list <- function(re_only = TRUE, convergence_threshold = 0.00001, max_iter = 10){
  list(re_only = re_only, convergence_threshold = convergence_threshold, max_iter = max_iter)
}

parameters <- est_parameters(my_formula, ds, method = "ppl",
                             control = control.list(max_iter = 100))


parameters2 <- est_parameters(my_formula, ds, method = "ppl",
                              start_params = c(parameters$beta, parameters$b),
                              theta_start = parameters$theta,
                              control = control.list())

attr(parameters2, "iterations")
attr(parameters, "iterations")


# compare estimates and hessians. In particular, hessians 1 and 4, which should
# be the same in parameters2.

max(abs(parameters$beta - parameters2$beta))

max(abs(parameters$b - parameters2$b))

max(abs(parameters$theta - parameters2$theta))


max(abs(attr(parameters,  "theta_ests")$hessian - attr(parameters2, "theta_ests")$hessian))

max(abs(attr(parameters,  "beta_b_est")$hessian - attr(parameters2, "beta_b_est")$hessian))


(hessian_estimates <- theta_ests[grepl("hessian", names(theta_ests))])
(hessian_estimates2 <- theta_ests2[grepl("hessian", names(theta_ests2))])





theta <- parameters$theta

ests <- optim(par = attr(parameters, "beta_b_est")$par,
              fn = lp,
              gr = lp_gr,
              formula = formula,
              parsed_data = make_ppl(parsed_data),
              other_args = list(n_fixed = n_fixed,
                                stat = stat,
                                theta = theta),
              method = "BFGS",
              control = list(fnscale = -1),
              hessian = TRUE)

# set up D
parsed_data$reTrms$Lambdat@x <- theta[parsed_data$reTrms$Lind]

D <- parsed_data$reTrms$Lambdat

D_inv <- solve(D)

# the derivatives of D
d_D_d_theta <- vector(mode = "list", length = length(parsed_data$reTrms$theta))

for (i in seq_along(parsed_data$reTrms$theta)) {

  d_D_d_theta[[i]] <- parsed_data$reTrms$Lambdat

  d_D_d_theta[[i]]@x <- c(theta[i], 0)[(parsed_data$reTrms$Lind != i) + 1]

}
#




# In Ripatti Palmgren they only use this part of the Hessian matrix.
Kbb <- ests$hessian[-(1:n_fixed), -(1:n_fixed)]

b <- Matrix(ests$par[-(1:n_fixed)], ncol = 1)

K_inv <- solve(Kbb)

# now, for each theta, calculate the terms.

# some common terms
#
# D_inv_dD <- lapply(d_D_d_theta, function(dD) {
#   D_inv %*% dD
# })
#
# D_inv_dD_D_inv <- lapply(d_D_d_theta, function(dD) {
#   D_inv %*% dD %*% D_inv
# })
#
# D_inv_dD_D_inv_dD_D_inv <- lapply(d_D_d_theta, function(dD) {
#   D_inv %*% dD %*% D_inv %*% dD %*% D_inv
# })

term3_matrix <- term2_matrix <- term1_matrix <- matrix(NA, nrow = length(theta), ncol = length(theta))


for (i in seq_along(theta)) {

  for (j in seq_along(theta)) {

    term1_matrix[i, j] <- -sum(diag(D_inv %*% d_D_d_theta[[j]] %*% D_inv %*% d_D_d_theta[[i]]))

    term2_matrix[i, j] <- sum(diag(K_inv %*% D_inv %*% d_D_d_theta[[j]] %*% D_inv %*% K_inv %*% D_inv %*% d_D_d_theta[[i]] %*% D_inv
                                      +K_inv %*% (D_inv %*% d_D_d_theta[[j]] %*% D_inv %*% d_D_d_theta[[i]] %*% D_inv +
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



#
# # first term
# term1 <- lapply(seq_along(theta), function(i) {
#
#   -0.5 * sum(diag(D_inv_dD[[i]] %*% D_inv_dD[[i]]))
# })

#
# # second term
# term2 <- lapply(seq_along(theta), function(i) {
#
#   -0.5 * sum(diag(K_inv %*% D_inv_dD_D_inv[[i]] %*% K_inv %*% D_inv_dD_D_inv[[i]] +
#                     2 * K_inv %*% D_inv_dD_D_inv_dD_D_inv[[i]]))
#
# })
#
# term2
#
# # third term
#
# term3 <- lapply(seq_along(theta), function(i) {
#
#     -t(b) %*% (2 * D_inv_dD_D_inv_dD_D_inv[[i]]) %*% b
#
# })
#
# term3

# lapply(seq_along(theta), function(i){
#   term1[[i]] + term2[[i]] + term3[[i]]
# })

# hessian estimates
theta_ests <- attr(parameters, "theta_ests")
hessian_estimates <- theta_ests[grepl("hessian", names(theta_ests))]

analytical_hessian <- -0.5 * (-term1_matrix + term2_matrix + term3_matrix)

I <- solve(-analytical_hessian)

sqrt(diag(I))

lapply(hessian_estimates, function(x){
  sqrt(diag(-solve(x)))
})





