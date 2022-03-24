

my_beta = c(1, -0.7, 0.5)
my_theta = 1

sample_data <- one_dataset(control = list(k = 50, nk = 4, beta = my_beta, theta = my_theta))

my_parms <- c(my_beta, b <- attr(sample_data, "random_effects"))

D = my_theta * diag(length(b))

Z_matrix <- model.matrix( ~ as.factor(M) - 1, data = sample_data)

colnames(Z_matrix) <- paste0("Z", seq(ncol(Z_matrix)))

data_with_Z <- dplyr::bind_cols(sample_data, data.frame(Z_matrix))


d2 <- sortAndIndex(data_with_Z, t)

d3 <- calcLinearPredictor(data = d2, X = c("X1", "X2", "X3"), Z = colnames(Z_matrix), parms = my_parms)


# checking the linear predictor
# X <- d2[, c("X1", "X2", "X3")] |> as.matrix()
#
# lp <- X %*% myParms[1:3] + myParms[-(1:3)][d2$M]
#
# cbind(d3$lp, lp, d3$lp - lp)
#
# exp_lp <- exp(X %*% myParms[1:3] + myParms[-(1:3)][d2$M])
#
# cbind(d3$A, exp_lp, d3$A - exp_lp)

d4 <- calcRiskSets(d3)

fit <- survival::coxph(survival::Surv(t, stat) ~ X1 + X2 + X3,data = sample_data)

start_parameters = c(coef(fit), rep(0, length(b)))

fit_optim <- optim(par = start_parameters,
                   fn = minuslp,
                   X = c("X1", "X2", "X3"),
                   t = t,
                   cluster = M,
                   dij = stat,
                   D = D,
                   data = sample_data)



fit_coxph <-
  tryCatch(
    survival::coxph(survival::Surv(t, stat) ~ X1 + X2 + X3 +
                      survival::frailty(M,
                                        distribution = "gaussian",
                                        method = "reml"),
                    control = survival::coxph.control(iter.max = 5000,
                                                      outer.max = 20,
                                                      timefix = FALSE),
                    data = sample_data),
    warning = function(w) w,
    error = function(e) e)

fit_coxph$penalty

fit_coxme <- coxme::coxme(survival::Surv(t, stat) ~ X1 + X2 + X3 + (1 | M), data = sample_data)

b_fit <- coxme::ranef(fit_coxme)$M

theta_fit <- coxme::VarCorr(fit_coxme)$M

D_theta_fit <- theta_fit * diag(length(b_fit))

lp(c(coxme::fixef(fit_coxme), coxme::ranef(fit_coxme)$M),
   X = c("X1", "X2", "X3"), t = t, cluster = M, dij = stat, D = D_theta_fit, data = sample_data)

fit_coxme$penalty
fit_coxme$loglik

fit_coxph$penalty

fit_coxph$loglik


D_theta_fit <- 1.216438 * diag(length(b_fit))


lp(c(fit_coxph$coefficients, fit_coxph$frail),
   X = c("X1", "X2", "X3"), t = t, cluster = M, dij = stat, D = D_theta_fit, data = sample_data)

my_beta = c(0.5, -0.5, 1)
my_theta = 0.2

sample_data <- one_dataset(control = list(k = 10, nk = 10, beta = my_beta, theta = my_theta))

d2 <- sortAndIndex(sample_data, t)

debugonce(calcCrossProducts)

calcCrossProducts(data = d2, f1 = c("X1", "X2", "X3"), n1 = "Xr", n2 = "Xs")




Z_matrix <- model.matrix( ~ as.factor(M) - 1, data = sample_data)

colnames(Z_matrix) <- paste0("Z", seq(ncol(Z_matrix)))

data_with_Z <- dplyr::bind_cols(sample_data, data.frame(Z_matrix))

b <- parms[-seq_len(length.out = length(X))]

D <- my_theta * diag(dplyr::n_distinct(sample_data$M))

b <- attr(sample_data, "random_effects")

b %*% solve(D) %*% b



my_beta = c(X1 = 1, X2 = -0.7, X3 = 0.5)
my_theta = 0.2

sample_data <- one_dataset(control = list(k = 10, nk = 4, beta = my_beta, theta = my_theta))

my_parms <- c(my_beta, b <- attr(sample_data, "random_effects"))

names(my_parms) <- c(paste0("X", 1:3), paste0("Z", seq_len(length(b))))
names(b) <- paste0("Z", seq_len(length(b)))

D = my_theta * diag(length(b))

rownames(solve(D) %*% b)


debugonce(dlp_beta)

dlp_beta(my_parms, X = c("X1", "X2", "X3"), t = t, cluster = "M", dij = stat, data = sample_data)

debugonce(dlp_b)

dlp_b(my_parms, X = c("X1", "X2", "X3"), t = t, cluster = "M", dij = stat, data = sample_data, D = D)



lp_grd(my_parms, X = c("X1", "X2", "X3"), t = t, cluster = "M", dij = stat, data = sample_data, D = D)


fit <- survival::coxph(survival::Surv(t, stat) ~ X1 + X2 + X3, data = sample_data)

start_parameters = c(coef(fit), rep(0, length(b)))

names(start_parameters) <- c(paste0("X", 1:3), paste0("Z", seq_len(length(b))))

fit_optim <- optim(par = start_parameters,
                   fn = lp,
                   gr = lp_grd,
                   X = c("X1", "X2", "X3"),
                   t = t,
                   cluster = "M",
                   dij = stat,
                   D = D,
                   data = sample_data,
                   method = "BFGS",
                   control = list(fnscale = -1))



lp_grd(my_parms, X = c("X1", "X2", "X3"), t = t, cluster = "M", dij = stat, data = sample_data, D = D)
debugonce(dlp_b)
dlp_b(start_parameters, X = c("X1", "X2", "X3"), t = t, cluster = "M", dij = stat, data = sample_data, D = D)




fit_coxph <-
  tryCatch(
    survival::coxph(survival::Surv(t, stat) ~ X1 + X2 + X3 +
                      survival::frailty(M,
                                        distribution = "gaussian",
                                        method = "reml"),
                    control = survival::coxph.control(iter.max = 5000,
                                                      outer.max = 20,
                                                      timefix = FALSE),
                    data = sample_data),
    warning = function(w) w,
    error = function(e) e)

fit_coxph$penalty

fit_coxme <- coxme::coxme(survival::Surv(t, stat) ~ X1 + X2 + X3 + (1 | M), data = sample_data)

b_fit <- coxme::ranef(fit_coxme)$M

theta_fit <- coxme::VarCorr(fit_coxme)$M

D_theta_fit <- theta_fit * diag(length(b_fit))

lp(c(coxme::fixef(fit_coxme), coxme::ranef(fit_coxme)$M),
   X = c("X1", "X2", "X3"), t = t, cluster = M, dij = stat, D = D_theta_fit, data = sample_data)

fit_coxme$penalty
fit_coxme$loglik

fit_coxph$penalty

fit_coxph$loglik

cbind(fit_optim$par[1:3],
      coxme::fixef(fit_coxme),
      coef(fit_coxph))



#######################



fit_optim <- optim(par = start_parameters,
                   fn = lp,
                   gr = lp_grd,
                   X = c("X1", "X2", "X3"),
                   t = t,
                   cluster = "M",
                   dij = stat,
                   D = D,
                   data = sample_data,
                   method = "BFGS",
                   control = list(fnscale = -1))


fit_coxme <- coxme::coxme(survival::Surv(t, stat) ~ X1 + X2 + X3 + (1 | M), data = sample_data,
                          vfixed = list(M = my_theta))


cbind(fit_optim$par[1:3],
      coxme::fixef(fit_coxme))



######################

















my_beta <- c(1, 1, 1)
my_theta <- 1

sample_data <- one_dataset(control = list(k = 30, nk = 5, beta = my_beta, theta = my_theta))

myParms <- c(my_beta, b <- attr(sample_data, "random_effects"))

D = theta * diag(b)

debugonce(lp)

lp(my_parms, X = c("X1", "X2", "X3"), t = t, cluster = M, dij = stat, D = D, data = sample_data)
minuslp(myParms, X = c("X1", "X2", "X3"), t = t, cluster = M, dij = stat, D = D, data = sample_data)



my_beta <- c(1, 1, 1)
my_theta <- 1

library(coxme)
fit_coxme <- coxme::coxme(survival::Surv(t, stat) ~ X1 + X2 + X3 + (1 | M), data = sample_data)

lp(c(fixef(fit_coxme), ranef(fit_coxme)$M),
   X = c("X1", "X2", "X3"), t = t, cluster = M, dij = stat, D = D_theta_fit, data = sample_data)

# These should be the same as returned by coxme
fit_coxme$loglik
fit_coxme$penalty

b_fit <- ranef(fit_coxme)$M

theta_fit <- VarCorr(fit_coxme)$M

D_theta_fit <- theta_fit * diag(length(b_fit))



0.5 * b_fit %*% D_theta_fit %*% b_fit

fit_coxme$penalty


one_rep <- function(){
sample_data <- one_dataset(control = list(k = 30, nk = 5, beta = my_beta, theta = my_theta))

myParms <- c(my_beta, b <- attr(sample_data, "random_effects"))

D = my_theta * diag(b)


fit <- survival::coxph(survival::Surv(t, stat) ~ X1 + X2 + X3,data = sample_data)

start_parameters = c(coef(fit), rep(0, length(b)))

fit_optim <- optim(par = start_parameters,
      fn = minuslp,
      X = c("X1", "X2", "X3"),
      t = t,
      cluster = M,
      dij = stat,
      D = D,
      data = sample_data)

fit_coxph <-
  tryCatch(
    survival::coxph(survival::Surv(t, stat) ~ X1 + X2 + X3 +
                      survival::frailty(M,
                                        distribution = "gaussian",
                                        method = "reml"),
                    control = survival::coxph.control(iter.max = 5000,
                                                      outer.max = 20,
                                                      timefix = FALSE),
                    data = sample_data),
    warning = function(w) w,
    error = function(e) e)

list(
  optim_pars = fit_optim$par[1:3],
  coxph_pars = coef(fit_coxph))

}

one_rep()



