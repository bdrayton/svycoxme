

my_beta = c(1, -0.7, 0.5)
my_theta = 0.2
my_k = 10
my_nk = 10

my_X = c("X1", "X2", "X3")

sample_data <- one_dataset(control = list(k = my_k, nk = my_nk, beta = my_beta, theta = my_theta))

b <- attr(sample_data, "random_effects")

my_params <- c(my_beta[seq_along(my_X)], b)

BB(parms = my_params, X = my_X, t = t, cluster = "M", dij = stat,
   data = sample_data)

bb(parms = my_params, X = my_X, t = t, cluster = "M", dij = stat,
   data = sample_data, theta = my_theta, return_matrix = TRUE)

Bb(parms = my_params, X = my_X, t = t, cluster = "M", dij = stat,
   data = sample_data)

my_hessian <- ppl_hessian(parms = my_params, X = my_X, t = t,
                          cluster = "M", dij = stat, data = sample_data,
                          theta = my_theta)

fit <- survival::coxph(survival::Surv(t, stat) ~ X1 + X2 + X3, data = sample_data)

start_parameters = c(coef(fit), rep(0, length(b)))

names(start_parameters) <- c(my_X, paste0("Z", seq_len(length(b))))

fit_optim <- optim(par = start_parameters,
                   fn = lp,
                   gr = lp_grd,
                   X = my_X,
                   t = "t",
                   cluster = "M",
                   dij = stat,
                   theta = my_theta,
                   data = sample_data,
                   method = "BFGS",
                   control = list(fnscale = -1),
                   hessian = TRUE)

my_hessian <- ppl_hessian(parms = fit_optim$par, X = my_X, t = t,
                          cluster = "M", dij = stat, data = sample_data,
                          theta = my_theta)

# any meaningful differences between numeric and analytical hessian?

mean(abs(my_hessian))

# What is coxme doing?

fit <- coxme::coxme(survival::Surv(t, stat) ~ X1 + X2 + X3 + (1|M), data = sample_data)


hmat <- as.matrix(fit$hmat)

class(as.matrix(fit$hmat)) * t(as.matrix(fit$hmat))

str(fit$hmat)

L <- as.matrix(fit$hmat)
D <- diag(fit$hmat)

methods("diag")
search()

library(bdsmatrix)
L %*% diag(D) %*% t(L)

diag(my_hessian)


# test the k_ppl function
debugonce(bb2)
bb2(parms = my_params, X = my_X, t = t, cluster = "M", dij = stat,
    data = sample_data, theta = my_theta, return_matrix = TRUE)





