

my_beta = c(1, -0.7, 0.5)
my_theta = 1
my_k = 50
my_nk = 10

ds <- one_dataset(list(k = my_k, nk = my_nk, beta = my_beta, theta = my_theta))

fit <- coxme::coxme(survival::Surv(t, stat) ~ X1 + X2 + X3 + (1|M), data = ds)

coxme_est_theta <- coxme::VarCorr(fit)$M

coxme_est_parms <- c(coxme::fixef(fit), coxme::ranef(fit)$M)

names(coxme_est_parms) <- paste0(rep(c("X", "Z"), c(3, my_k)), c(1:3, seq_len(my_k)))

my_loglik <- lp(parms = coxme_est_parms,
   X = c("X1", "X2", "X3"),
   cluster = "M",
   t = t, dij = stat,
   theta = coxme_est_theta,
   data = ds)

fit$loglik["Penalized"] - my_loglik

fit$penalty - attr(my_loglik, "penalty")


my_u <- lp_grd(parms = coxme_est_parms,
               X = c("X1", "X2", "X3"),
               cluster = "M",
               t = t,
               dij = stat,
               theta = coxme_est_theta,
               data = ds)

# the gradient in coxme have a mixed up order

cbind(my_u[order(my_u)], fit$u[order(fit$u)])









