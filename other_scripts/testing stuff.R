


sample_data <- one_dataset(control = list(k = 50, nk = 4, beta = c(1, -0.7, 0.5), theta = 1))

d1 <- dplyr::bind_cols(sample_data, data.frame(Z))

theta = 1

D = theta * diag(50)





d2 <- sortAndIndex(d1, t)

d3 <- calcLinearPredictor(data = d2, X = c("X1", "X2", "X3"), Z = colnames(Z), parms = myParms)


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







debugonce(minuslp)

my_beta <- c(1, 1, 1)
my_theta <- 1

sample_data <- one_dataset(control = list(k = 30, nk = 5, beta = my_beta, theta = my_theta))

myParms <- c(my_beta, b <- attr(sample_data, "random_effects"))

D = theta * diag(b)

debugonce(lp)

     lp(myParms, X = c("X1", "X2", "X3"), t = t, cluster = M, dij = stat, D = D, data = sample_data)
minuslp(myParms, X = c("X1", "X2", "X3"), t = t, cluster = M, dij = stat, D = D, data = sample_data)

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




my_beta <- c(1, 1, 1)
my_theta <- 1


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



