


my_k = 50
my_theta = 1

ds <- one_dataset(list(k = my_k, nk = 10, beta = c(1, -0.7, 0.5), theta = my_theta))

fit <- coxme::coxme(survival::Surv(t, stat) ~ X1 + X2 + X3 + (1|M), data = ds)

est_theta <- coxme::VarCorr(fit)$M

D = est_theta * diag(my_k)

est_parms <- c(coxme::fixef(fit), coxme::ranef(fit)$M)

names(est_parms) <- paste0(rep(c("X", "Z"), c(3, my_k)), c(1:3, seq_len(my_k)))

my_loglik <- lp(parms = est_parms,
   X = c("X1", "X2", "X3"),
   t = t, dij = stat, D =  D, data = ds)

(fit$loglik["Penalized"] - my_loglik) < 1e-10

fit$penalty == attr(my_loglik, "penalty")


my_u <- lp_grd(parms = est_parms,
               X = c("X1", "X2", "X3"),
               cluster = "M",
               t = t,
               dij = stat,
               D =  D,
               data = ds)

cbind(my_u, fit$u)







