my_beta = c(1, -0.7, 0.5)
my_theta = 0.2
my_k = 10
my_nk = 10
my_X = c("X1", "X2", "X3")

ds <- one_dataset(list(k = my_k, nk = my_nk, beta = my_beta, theta = my_theta))


fit0 <- survival::coxph(survival::Surv(stat_time, stat) ~ X1 + X2 + X3, data = ds)

my_start_parameters <- c(coef(fit0), rep(0, my_k))
names(my_start_parameters) <- c(my_X, paste0("Z", seq_len(my_k)))

test_loop <- estimate_parameters_loop(start_theta = 0.5,
                                      start_parms = my_start_parameters,
                                      X = my_X,
                                      stat_time = stat_time,
                                      cluster = "M",
                                      dij = stat,
                                      data = ds, max_iter = 200)

ests <- tail(test_loop$estimate_history, 1)[[1]]

theta_ipl_gr(one_theta = ests$new_theta,
             parms = ests$new_beta_b,
             X = my_X,
             stat_time = stat_time,
             dij = stat,
             cluster = "M",
             data = ds)

coxme_fit <- coxme::coxme(survival::Surv(stat_time, stat) ~ X1 + X2 + X3 + (1|M), data = ds)

library(bdsmatrix)
D <- diag(coxme_fit$hmat)
L <- as.matrix(coxme_fit$hmat)

hess <- L %*% diag(D) %*% t(L)

hess_22 <- hess[seq_len(my_k), seq_len(my_k)]

inv_hess_22 <- solve(hess_22)



theta_inv <- 1/coxme::VarCorr(coxme_fit)$M

0.5 * (inner(b) * (theta_inv)^2 - sum(my_k * theta_inv) - sum(diag(inv_hess_22)))

debugonce(theta_ipl_gr)

theta_ipl_gr(one_theta = coxme::VarCorr(coxme_fit)$M,
             parms = c(coxme::fixef(coxme_fit), coxme::ranef(coxme_fit)$M),
             X = my_X,
             stat_time = stat_time,
             dij = stat,
             cluster = "M",
             data = ds)


theta_ipl_gr(one_theta = my_theta,
             parms = c(coxme::fixef(coxme_fit), coxme::ranef(coxme_fit)$M),
             X = my_X,
             stat_time = stat_time,
             dij = stat,
             cluster = "M",
             data = ds)



# library(coxme)

# coxfit <- coxme::coxme(survival::Surv(stat_time, stat) ~ X1 + X2 + X3 + (1|M), data = ds)

fit0 <- survival::coxph(survival::Surv(stat_time, stat) ~ X1 + X2 + X3, data = ds)

my_start_parameters <- c(coef(fit0), rep(0, my_k))
names(my_start_parameters) <- c(my_X, paste0("Z", seq_len(my_k)))

fit2 <- estimate_parameters2(start_theta = 0.5,
                    start_parms = my_start_parameters,
                    X = my_X,
                    stat_time = stat_time,
                    cluster = "M",
                    dij = stat,
                    data = ds)

coxme_fit <- coxme::coxme(survival::Surv(stat_time, stat) ~ X1 + X2 + X3 + (1|M), data = ds)

theta_est <- coxme::VarCorr(coxme_fit)$M
b_est <- coxme::ranef(coxme_fit)$M

library(bdsmatrix)

L <- as.matrix(coxme_fit$hmat)
D <- diag(coxme_fit$hmat)

hess <- L %*% diag(D) %*% t(L)

hess_22 <- hess[1:50, 1:50]

inv_hess_22 <- solve(hess_22)

var_theta <- 2 * theta_est^2 * (length(b_est) + sum(diag(inv_hess_22 %*% inv_hess_22)) -
                     (2/theta_est) * sum(diag(inv_hess_22)))^(-1)

theta_est + c(-1.96, 1.96) * var_theta








test_loop <- estimate_parameters_loop(start_theta = 0.5,
                    start_parms = my_start_parameters,
                    X = my_X,
                    stat_time = stat_time,
                    cluster = "M",
                    dij = stat,
                    data = ds, max_iter = 200)

test_loop$converged

ests <- tail(test_loop$estimate_history, 1)[[1]]

coxme_fit <- coxme::coxme(survival::Surv(stat_time, stat) ~ X1 + X2 + X3 + (1|M), data = ds)

coxme::VarCorr(coxme_fit)
ests$new_theta

cbind(c(ests$new_beta_b, theta = ests$new_theta),
      c(coxme::fixef(coxme_fit),
        coxme::ranef(coxme_fit)$M, theta = coxme::VarCorr(coxme_fit)$M))




# fit2 <- optim_ipl(theta = VarCorr(coxfit)$M, start_parms = my_start_parameters, X = my_X,
                  # stat_time = stat_time, cluster = "M", dij = stat, data = ds)

# optim_ipl(theta = VarCorr(coxfit)$M, start_parms = my_start_parameters, X = my_X,
#           stat_time = stat_time, cluster = "M", dij = stat, data = ds, likelihood_only = TRUE)

# cbind(
#   my_start_parameters,
#   fit2[[1]]$par)

optimfit <- optim(c(0.5), optim_ipl, start_parms = my_start_parameters, X = my_X,
      stat_time = stat_time, cluster = "M", dij = stat, data = ds,
      likelihood_only = TRUE, method = "BFGS", control = list(fnscale = -1))

VarCorr(coxfit)

coxfit$loglik

optimfit

# Alright that doesn't work. Lets go back to the back and forth method, but with these newer likelihoods.

# estimate beta and b for fixed theta


# start shit
optim_loop <- function(start_theta,
                       start_parms,
                       X,
                       stat_time,
                       cluster,
                       dij,
                       data = ds) {

  fit_beta_b <- optim(par = start_parms,
                     fn = lp,
                     gr = lp_grd,
                     X = X,
                     stat_time = {{ stat_time }},
                     cluster = cluster,
                     dij = {{ dij }},
                     theta = start_theta,
                     data = data,
                     method = "BFGS",
                     control = list(fnscale = -1))

  fit_theta <- optim(par = c(start_theta),
                     fn = theta_ipl,
                     gr = NULL,
                     parms = fit_beta_b$par,
                     X = X,
                     stat_time = {{ stat_time }},
                     dij = {{ dij }},
                     cluster = cluster,
                     data = data,
                     method = "BFGS",
                     control = list(fnscale = -1))

  # check convergence
  if(fit_beta_b$convergence != 0 | fit_theta$convergence != 0 ) stop("failed to converge")

  list(new_theta = fit_theta$par,
       new_beta_b = fit_beta_b$par)

}


theta_ipl <-   function(one_theta, parms, X, stat_time, dij, cluster, data ){

  D <- one_theta * diag(length(parms) - length(X))

  kbb <- bb(parms = parms,
            X = X,
            stat_time = {{ stat_time }},
            dij = {{ dij }},
            theta = one_theta,
            cluster = cluster,
            data = data,
            return_matrix = TRUE)

  penalised <- lp(parms = parms,
                  X = X,
                  stat_time = {{ stat_time }},
                  dij = {{ dij }},
                  theta = one_theta,
                  cluster = cluster,
                  data = data)

  ipl <- -0.5 * log(det(D)) - 0.5 * log(det(kbb)) + penalised

  attr(ipl, "penalty") <- NULL

  ipl

}

theta_ipl(one_theta = 0.5,
          parms = my_start_parameters,
          X = my_X,
          stat_time = stat_time,
          dij = stat,
          cluster = "M",
          data = ds)

debugonce(optim_loop)

optim_loop(start_theta = 0.5,
           start_parms = my_start_parameters,
           X = my_X,
           stat_time = stat_time,
           cluster = "M",
           dij = stat,
           data = ds)




(inner(kbb) + sum(diag(solve(kbb))))/length(kbb)


coxme_b <- ranef(coxfit)$M

coxme_beta <- fixef(coxfit)

coxme_theta <- VarCorr(coxfit)$M

coxme_D <- coxme_theta * diag(length(coxme_b))

coxme_loglik <- coxfit$loglik

pl_theta(theta = coxme_theta, b = coxme_b,
         K_ppl = K_prime_prime(parms = c(coxme_beta, coxme_b),
                               X = my_X,
                               t = stat_time,
                               dij = stat,
                               theta = coxme_theta,
                               cluster = "M",
                               data = ds))





kbb <- bb(parms = c(coxme_beta, coxme_b),
          X = my_X,
          t = stat_time,
          dij = stat,
          theta = coxme_theta,
          cluster = "M",
          data = ds,
          return_matrix = TRUE)

debugonce(lp)

lp(parms = c(coxme_beta, coxme_b),
   X = my_X,
   stat_time = stat_time,
   dij = stat,
   theta = coxme_theta,
   cluster = "M",
   data = ds)

my_start_params <- c(coxme_beta, coxme_b)
names(my_start_params) <- c(my_X, paste0("Z", 1:10))

lp_grd(parms = my_start_params,
   X = my_X,
   stat_time = stat_time,
   dij = stat,
   theta = coxme_theta,
   cluster = "M",
   data = ds)




my_int_lik <- -0.5 * log(det(coxme_D)) - 0.5 * log(det(kbb)) + penalised

my_int_lik - coxme_loglik[2]


# try and get a parabola for theta likelihoods.

coxph_fit <- coxph(Surv(stat_time, stat) ~ X1 + X2 + X3, data = ds)

my_start_params <- c(coef(coxph_fit), rep(0, my_nk))

names(my_start_params) <- c(my_X, paste0("Z", seq_len(my_nk)))

optim_ipl(theta = coxme_theta, start_parms = my_start_params, X = my_X, stat_time = stat_time, cluster = "M", dij = stat, data = ds)

test_thetas <- seq(0.01, 5, by = 1)

test_ipl <- sapply(test_thetas, function(one_theta){
  optim_ipl(theta = one_theta,
            start_parms = my_start_params,
            X = my_X,
            stat_time = stat_time,
            cluster = "M",
            dij = stat,
            data = ds)
})

# clearly no maximum.
plot(test_thetas, test_ipl)


# what about plotting lp(theta) for fixed b and beta ?

library(coxme)

coxfit <- coxme::coxme(survival::Surv(stat_time, stat) ~ X1 + X2 + X3 + (1|M), data = ds)

coxme_b <- ranef(coxfit)$M

coxme_beta <- fixef(coxfit)


which(max(test_ipl_fixed) == test_ipl_fixed)

test_thetas <- seq(0.5, .71, length.out = 1000)

theta_ipl <-   function(one_theta, parms){
  D <- one_theta * diag(length(coxme_b))

  kbb <- bb(parms = parms,
            X = my_X,
            t = stat_time,
            dij = stat,
            theta = one_theta,
            cluster = "M",
            data = ds,
            return_matrix = TRUE)

  penalised <- lp(parms = parms,
                  X = my_X,
                  stat_time = stat_time,
                  dij = stat,
                  theta = one_theta,
                  cluster = "M",
                  data = ds)


  -0.5 * log(det(D)) - 0.5 * log(det(kbb)) + penalised

}


test_ipl_fixed <- sapply(test_thetas, theta_ipl, parms = c(coxme_beta, coxme_b))

plot(test_thetas, test_ipl_fixed, type = "l")

test_thetas[max(test_ipl_fixed) == test_ipl_fixed]

theta_ipl(test_thetas[max(test_ipl_fixed) == test_ipl_fixed], parms = c(coxme_beta, coxme_b)) - theta_ipl(VarCorr(coxfit)$M, parms = c(coxme_beta, coxme_b))

coxfit$loglik

VarCorr(coxfit)$M

debug(optim_ipl)
sapply(c(0.5, 0.6, 0.7), function(one_theta){
  optim_ipl(theta = one_theta,
            start_parms = my_start_params,
            X = my_X,
            stat_time = stat_time,
            cluster = "M",
            dij = stat,
            data = ds)
})


# penalty = c(0.5 * t(b)%*%solve(D)%*%b)
penalty = inner(a = b) * 1/theta

theta <- 0.1

D <- theta * diag(length(coxme_b))

c(0.5 * t(coxme_b)%*%solve(D)%*%coxme_b) - 0.5 * inner(a = coxme_b) * 1/theta

















