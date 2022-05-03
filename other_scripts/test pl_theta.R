




my_beta = c(1, -0.7, 0.5)
my_theta = 0.2
my_k = 10
my_nk = 10

my_X = c("X1", "X2", "X3")

sample_data <- one_dataset(control = list(k = my_k, nk = my_nk, beta = my_beta, theta = my_theta))

b <- attr(sample_data, "random_effects")

my_params <- c(my_beta[seq_along(my_X)], b)


K_ppl <- bb(parms = my_params, X = my_X, t = t, cluster = "M", dij = stat,
   data = sample_data, theta = my_theta, return_matrix = TRUE)

debugonce(pl_theta)

pl_theta(my_theta, b, K_ppl)

optim(par = 0, fn = pl_theta, b = b, K_ppl = K_ppl, method = "Brent", lower = 0, upper = 100,
      control = list(fnscale = -1))

