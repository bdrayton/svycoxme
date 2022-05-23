




my_beta = c(1, -0.7, 0.5)
my_theta = 1
my_k = 10
my_nk = 10
my_X = c("X1", "X2", "X3")

ds <- one_dataset(list(k = my_k, nk = my_nk, beta = my_beta, theta = my_theta))

my_b <- attr(ds, "random_effects")

kpp <- K_prime_prime(parms = c(my_beta, my_b),
              X = my_X,
              t = t,
              dij = stat,
              theta = my_theta,
              cluster = "M",
              data = ds)

kbb <- bb(parms = c(my_beta, my_b),
   X = my_X,
   t = t,
   dij = stat,
   theta = my_theta,
   cluster = "M",
   data = ds,
   return_matrix = TRUE)

-0.5 * log(det(kpp))
-0.5 * log(det(kbb))


# for a given theta, maximise beta and b, and then calculate the full likelihood.

test_theta <- 0.8

my_D <-  test_theta * diag(length(my_b))

fitph <- survival::coxph(survival::Surv(t, stat) ~ X1 + X2 + X3, data = ds)

start_parms <- c(coef(fitph), rep(0, length(my_b)))

names(start_parms) <- c(my_X, paste0("Z", seq_len(my_nk)))

fit_optim <- optim(par = start_parms,
                   fn = lp,
                   gr = lp_grd,
                   X = my_X,
                   t = "t" ,
                   cluster = "M",
                   dij = stat,
                   D = my_D,
                   theta = test_theta,
                   data = ds,
                   method = "BFGS",
                   control = list(fnscale = -1))


# return likelihood for given theta, and estimated beta and b

kpp <- K_prime_prime(parms = fit_optim$par,
                     X = my_X,
                     t = t,
                     dij = stat,
                     theta = test_theta,
                     cluster = "M",
                     data = ds)

kbb <- bb(parms = fit_optim$par,
          X = my_X,
          t = t,
          dij = stat,
          theta = test_theta,
          cluster = "M",
          data = ds,
          return_matrix = TRUE)


b <- fit_optim$par[-seq_len(length(my_X))]


-0.5 * ( log(det(my_D)) + log(det(kpp)) +  inner(b) * test_theta)


-0.5 * ( log(det(my_D)) + log(det(kbb)) +  inner(b) * test_theta)


debugonce(optim_ipl)
optim_ipl(theta = test_theta, start_parms = start_parms,
          X = my_X, t = "t", cluster = "M", dij = stat, data = ds)



optim(c(0.8), fn = optim_ipl, gr = NULL, method = "L-BFGS-B", control = list(fnscale = 1),
      start_parms = start_parms, X = my_X, t = "t", cluster = "M", dij = stat, data = ds,
      lower = 0.01, upper = 100)

coxfit <- coxme::coxme(survival::Surv(t, stat) ~ X1 + X2 + X3 + (1|M), data = ds)

coxme::VarCorr(coxfit)$M



