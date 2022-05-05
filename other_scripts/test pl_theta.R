
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


# compare with theta from coxme

coxme_fit <- coxme::coxme(survival::Surv(t, stat) ~ X1 + X2 + X3 + (1|M), data = sample_data)

(coxme_theta <- coxme::VarCorr(coxme_fit)$M)

b <- coxme::ranef(coxme_fit)$M

K_ppl <- bb(parms = c(coxme::fixef(coxme_fit), b), X = my_X, t = t, cluster = "M", dij = stat,
            data = sample_data, theta = coxme_theta, return_matrix = TRUE)

hmat <- coxme_fit$hmat

L <- as.matrix(hmat)
D <- diag(hmat)

# check that they're the same(ish)
(L %*% diag(D) %*% t(L))[seq_len(my_nk), seq_len(my_nk)] - (-1 * K_ppl)

coxme_K_ppl <- L %*% diag(D) %*% t(L)

optim(par = 0, fn = pl_theta, b = b, K_ppl = coxme_K_ppl, method = "Brent", lower = 0, upper = 100,
      control = list(fnscale = -1))$par

theta <- seq(0.001, 0.5, by = 0.00001)

likes <- sapply(theta, pl_theta, b = b, K_ppl = coxme_K_ppl)

plot(theta, likes, type = "l")

theta[which(likes == max(likes))]

coxme_theta




D <- my_theta * diag(length(b))

c((-1/2) * ( log(det(D)) + log(det(K_ppl)) + t(b) %*% solve(D) %*% b ))















