


# analytical hessian for theta

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

ests <- est_parameters(my_formula, ds, control = control.list(grad = FALSE))

ests_w_gr  <- est_parameters(my_formula, ds, control = control.list(grad = TRUE))


ests2 <- est_parameters(my_formula, ds, start_params = c(ests$beta, ests$b), theta_start = ests$theta)
ests3 <- est_parameters(my_formula, ds, start_params = c(ests2$beta, ests2$b), theta_start = ests2$theta)

# thetas
cbind(
  unlist(coxme::VarCorr(coxme_fit)),
  ests$theta
  # ests2$theta,
  # ests3$theta
)

# betas
cbind(
  coefficients(coxme_fit),
  ests$beta,
  ests2$beta,
  ests3$beta)

# theta hessian
# more stable when sample size is a bit bigger. (50 clusters of 10)
# massively divergent when sample size is smaller (10 clusters of 4).
cbind(
  sqrt(diag(-solve(attr(ests, "theta_est")$hessian))),
  sqrt(diag(-solve(attr(ests2, "theta_est")$hessian))),
  sqrt(diag(-solve(attr(ests3, "theta_est")$hessian))))


# beta hessian
cbind(
  sqrt(diag(-solve(attr(ests, "beta_b_est")$hessian)))[1:3],
  sqrt(diag(-solve(attr(ests2, "beta_b_est")$hessian)))[1:3],
  sqrt(diag(-solve(attr(ests3, "beta_b_est")$hessian)))[1:3],
  sqrt(diag(vcov(coxme_fit))))

ds_sorted <- sortAndIndex(ds, sort_vars = stat_time)

parsed_data <- lme4::lFormula(my_formula, data = ds_sorted)

stat <- Matrix(unclass(parsed_data$fr[,1])[, "status"], ncol = 1)


theta_ipl(ests$theta, formula = my_formula, parsed_data = make_ppl(parsed_data),
          other_args = list(start_params = attr(ests, "beta_b_est")$par,
                            re_only = TRUE,
                            n_fixed = 3,
                            stat = stat))
# debugonce(theta_ipl_gr)
numDeriv::grad(theta_ipl, x = ests$theta, parsed_data = make_ppl(parsed_data),
                  other_args = list(start_params = attr(ests, "beta_b_est")$par,
                                    re_only = TRUE,
                                    n_fixed = 3,
                                    stat = stat))

theta_ipl_gr(ests$theta, formula = my_formula, parsed_data = make_ppl(parsed_data),
             other_args = list(start_params = attr(ests, "beta_b_est")$par,
                               re_only = TRUE,
                               n_fixed = 3,
                               stat = stat))



test_thetas <- expand.grid(c(ests$theta[1], seq(from = max(ests$theta[1]-1, 0.001), to = ests$theta[1]+0.5, length.out = 21)),
                           c(ests$theta[2], seq(from = max(ests$theta[2]-0.5, 0.001), to = ests$theta[2]+0.5, length.out = 21)))

theta_likelihoods <- apply(test_thetas, 1, function(theta){
  theta_ipl(theta, formula = my_formula, parsed_data = make_ppl(parsed_data),
            other_args = list(start_params = attr(ests, "beta_b_est")$par,
                              re_only = TRUE,
                              n_fixed = 3,
                              stat = stat))
})

theta_gradients <- apply(test_thetas, 1, function(theta){
  theta_ipl_gr(theta, formula = my_formula, parsed_data = make_ppl(parsed_data),
            other_args = list(start_params = attr(ests, "beta_b_est")$par,
                              re_only = TRUE,
                              n_fixed = 3,
                              stat = stat))
})




plot_data <- dplyr::filter(cbind(test_thetas, theta_likelihoods), Var1 > 0, Var2 > 0) |>
  dplyr::arrange(theta_likelihoods) |>
  as.matrix()

threejs::scatterplot3js(plot_data[,1], plot_data[,2], plot_data[,3],
                        color = c("blue", "red")[(plot_data[,3] == max(plot_data[,3])) + 1],
                        size = 0.1)


fig <- plotly::plot_ly(x = test_thetas[,1], y = test_thetas[,2], z = theta_likelihoods, type = 'mesh3d')

fig

grad_data <- cbind(test_thetas, t(theta_gradients), theta_likelihoods)

names(grad_data)

threejs::scatterplot3js(grad_data[,3], grad_data[,4], grad_data[,5],
                        color = c("blue", "red")[(grad_data[,5] == max(grad_data[,5])) + 1],
                        size = 0.1)


# evaluating the gradient function.

# set one theta, and vary the other by a small amount.

# calculate likelihoods and change over that interval. calculate gradient in the middle of the interval.



test_thetas2 <- data.frame(changing = ests$theta[1] + seq(-0.5, 0.5, by = 0.1),
                           constant = ests$theta[2])

likelihoods <- apply(test_thetas2, 1, function(theta){
                     theta_ipl(theta, formula = my_formula, parsed_data = make_ppl(parsed_data),
                               other_args = list(start_params = attr(ests, "beta_b_est")$par,
                                                 re_only = TRUE,
                                                 n_fixed = 3,
                                                 stat = stat))
  })


# gradients at each point
gradients <- apply(test_thetas2, 1, function(theta){
  theta_ipl_gr(theta, formula = my_formula, parsed_data = make_ppl(parsed_data),
            other_args = list(start_params = attr(ests, "beta_b_est")$par,
                              re_only = TRUE,
                              n_fixed = 3,
                              stat = stat))
})

gradients2 <- apply(test_thetas2, 1, function(theta){
  numDeriv::grad(theta_ipl, x = theta, parsed_data = make_ppl(parsed_data),
                 other_args = list(start_params = attr(ests, "beta_b_est")$par,
                                   re_only = TRUE,
                                   n_fixed = 3,
                                   stat = stat))
})

plot(test_thetas2$changing, likelihoods)
text(test_thetas2$changing, likelihoods, round(gradients[1,], 2), adj = c(0, 0))
text(test_thetas2$changing, likelihoods, round(gradients2[1,], 2), adj = c(0, 1))


ests$theta
coxme_fit$vcoef


theta_ipl_gr_num(ests$theta,
                 formula = my_formula,
                 parsed_data = make_ppl(parsed_data),
          other_args = list(start_params = attr(ests, "beta_b_est")$par,
                            re_only = TRUE,
                            n_fixed = 3,
                            stat = stat))




test_thetas[max(theta_likelihoods) == theta_likelihoods,]

coxme_fit$vcoef

coxme_fit$loglik

theta_ipl(unlist(coxme_fit$vcoef), formula = my_formula, parsed_data = make_ppl(parsed_data),
          other_args = list(start_params = attr(ests, "beta_b_est")$par,
                            re_only = TRUE,
                            n_fixed = 3,
                            stat = stat))


attr(ests, "theta_est")$hessian
numDeriv::hessian(theta_ipl, x = ests$theta, parsed_data = make_ppl(parsed_data),
               other_args = list(start_params = attr(ests, "beta_b_est")$par,
                                 re_only = TRUE,
                                 n_fixed = 3,
                                 stat = stat))

pracma::hessian(theta_ipl, x0 = ests2$theta, parsed_data = make_ppl(parsed_data),
                other_args = list(start_params = attr(ests2, "beta_b_est")$par,
                                  re_only = TRUE,
                                  n_fixed = 3,
                                  stat = stat))



microbenchmark::microbenchmark(
optim(par = ests$theta,
                   fn = theta_ipl,
                   gr = NULL,
                   formula = my_formula,
                   parsed_data = make_ppl(parsed_data),
                   other_args = list(n_fixed = 3,
                                     start_params = attr(ests, "beta_b_est")$par,
                                     stat = stat,
                                     re_only = TRUE),
                   method = "L-BFGS-B",
                   control = list(fnscale = -1),
                   lower = 0.00001, upper = Inf,
                   hessian = TRUE),
optim(par = ests$theta,
      fn = theta_ipl,
      gr = NULL,
      formula = my_formula,
      parsed_data = make_ppl(parsed_data),
      other_args = list(n_fixed = 3,
                        start_params = attr(ests, "beta_b_est")$par,
                        stat = stat,
                        re_only = TRUE),
      method = "L-BFGS-B",
      control = list(fnscale = -1),
      lower = 0.00001, upper = Inf))





theta_est_test$hessian

optimHess(ests$theta, fn = theta_ipl, gr = NULL, parsed_data = make_ppl(parsed_data),
          other_args = list(n_fixed = 3,
                            start_params = attr(ests, "beta_b_est")$par,
                            stat = stat,
                            re_only = TRUE),
          control = list(fnscale = -1))


theta_start <- get_start_theta(length(parsed_data$reTrms$flist))

n_fixed <- length(attr(terms(coxme:::formula1(my_formula)$fixed), "order"))

fixed_formula <- lme4:::getFixedFormula(my_formula)

fit0 <- survival::coxph(fixed_formula, data = ds)


start_params <- c(coef(fit0),
                  rnorm(n = length(parsed_data$reTrms$Lind),
                        mean = 0,
                        sd = sqrt(theta_start[parsed_data$reTrms$Lind])))

# debugonce(est_parameters)
parameters <- est_parameters(my_formula, ds, method = "ppl",
                             control = control.list(max_iter = 100))

# pull out estimates
thetas <- Reduce("rbind", lapply(parameters$theta_est, "[[", "par"))

theta_est <- data.frame(thetas)
names(theta_est) <- paste("theta", seq_len(ncol(theta_est)), sep = "_")

library(tidyverse)
tidyr::pivot_longer(theta_est, cols = everything()) %>%
  group_by(name) %>%
  mutate(iteration = row_number()) %>%
  ggplot(aes(iteration, value)) + geom_point() + facet_grid(vars(name))

# betas
betas <- Reduce("rbind", lapply(parameters$beta_b_est, "[[", "par"))

beta_est <- data.frame(betas) %>% select(starts_with("X"))
# names(theta_est) <- paste("theta", seq_len(ncol(theta_est)), sep = "_")

tidyr::pivot_longer(beta_est, cols = everything()) %>%
  group_by(name) %>%
  mutate(iteration = row_number()) %>%
  ggplot(aes(iteration, value)) + geom_point() + facet_grid(vars(name))




parameters2 <- est_parameters(my_formula, ds, method = "ppl",
                              start_params = tail(parameters$beta_b_est,1)[[1]]$par,
                              theta_start = tail(parameters$theta_est,1)[[1]]$par,
                              control = control.list(max_iter = 0))

parameters$theta_est[[1]]$par
parameters$theta_est[[101]]$par
parameters$theta_est[[2]]$par



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





