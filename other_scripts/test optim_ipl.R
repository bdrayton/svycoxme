



my_beta = c(1, -0.7, 0.5)
my_theta = 3
my_k = 50
my_nk = 10

max_iter = 100
convergence_threshold = 0.00001 # same as coxme

my_X = c("X1", "X2", "X3")


sample_data <- one_dataset(control = list(k = my_k, nk = my_nk, beta = my_beta, theta = my_theta))

fit <- survival::coxph(survival::Surv(t, stat) ~ X1 + X2 + X3, data = sample_data)

nb <- dplyr::n_distinct(sample_data$M)

start_parameters = c(coef(fit), rep(0, nb))

names(start_parameters) <- c(my_X, paste0("Z", seq_len(nb)))

debugonce(optim_ipl)
optim_ipl(theta = 1,
          start_parms = start_parameters,
          X = my_X,
          t = "t",
          cluster = "M",
          dij = stat,
          data = sample_data)

optim(par = c(1),
      fn = optim_ipl,
      gr = NULL,
      method = "Brent",
      lower = 0.01, upper = 100,
      start_parms = start_parameters,
      X = my_X,
      t = "t",
      cluster = "M",
      dij = stat,
      data = sample_data,
      control = list(fnscale = -1))


lapply(X = list(0.01, 0.5, 1, 10, 50), function(val){
  optim_ipl(theta = val,
            start_parms = start_parameters,
            X = my_X,
            t = "t",
            cluster = "M",
            dij = stat,
            data = sample_data)
})


coxme_fit <- coxme::coxme(survival::Surv(t, stat) ~ X1 + X2 + X3 + (1|M), data = sample_data)

coxme::VarCorr(coxme_fit)$M











