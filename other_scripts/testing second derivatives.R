

my_beta = c(1, -0.7, 0.5)
my_theta = 0.2
my_k = 10
my_nk = 50

sample_data <- one_dataset(control = list(k = my_k, nk = my_nk, beta = my_beta, theta = my_theta))

b <- attr(sample_data, "random_effects")

my_params <- c(my_beta[c(-2, -3)], b)


BB(parms = my_params, X = c("X1"), t = t, cluster = "M", dij = stat,
   data = sample_data)

bb(parms = my_params, X = c("X1"), t = t, cluster = "M", dij = stat,
   data = sample_data, theta = my_theta)

Bb(parms = my_params, X = c("X1"), t = t, cluster = "M", dij = stat,
   data = sample_data)

my_hessian <- ppl_hessian(parms = my_params, X = c("X1"), t = t,
                          cluster = "M", dij = stat, data = sample_data,
                          theta = my_theta)


fit <- survival::coxph(survival::Surv(t, stat) ~ X1, data = sample_data)

start_parameters = c(coef(fit), rep(0, length(b)))

names(start_parameters) <- c("X1", paste0("Z", seq_len(length(b))))

fit_optim <- optim(par = start_parameters,
                   fn = lp,
                   gr = lp_grd,
                   X = c("X1"),
                   t = t,
                   cluster = "M",
                   dij = stat,
                   D = my_theta * diag(length(b)),
                   data = sample_data,
                   method = "BFGS",
                   control = list(fnscale = -1),
                   hessian = TRUE)


fit_optim$hessian
my_hessian


diag(fit_optim$hessian)
diag(my_hessian)

sqrt(diag(solve(-my_hessian)))
sqrt(diag(solve(-fit_optim$hessian)))

coxme::coxme(survival::Surv(t, stat) ~ X1 + (1|M), data = sample_data)












