
my_formula <- survival::Surv(stat_time, stat)~X1 + X2 + X3 + (1 | M1) + (1| M1:M2)

ds <- one_dataset(my_formula,
                  dists = list(X1 = ~rnorm(n),
                               X2 = ~rnorm(k * nk),
                               X3 = ~rbinom(n, 1, 0.5),
                               M1 = ~rep(1:k, each = nk),
                               M2 = ~rep(c("l", "r"), ceiling(n/2))[seq_len(n)],
                               error = ~rexp(n, 10),
                               stat = ~sample(rep(c(0, 1), round(n * c(0.2, 0.8))), n)),
                  dist_args = list(k = 50, nk = 4, n = 50*4),
                  coefficients = c(1, -0.7, 0.5),
                  random_effect_variance = c(M1 = 1, `M1:M2` = 1.5)
)

ds$`M1:M2` <- with(ds, interaction(M1, M2))

coxfit <- coxme::coxme(survival::Surv(stat_time, stat)~X1 + X2 + X3 + (1 | M1) + (1| `M1:M2`), data = ds)

coxme::VarCorr(coxfit)
ds_sorted <- sortAndIndex(ds, sort_vars = stat_time)

parsed_data <- lme4::lFormula(my_formula, data = ds_sorted)

stat <- Matrix(unclass(parsed_data$fr[,1])[, "status"], ncol = 1)

theta_start <- get_start_theta(length(parsed_data$reTrms$flist))

n_fixed <- length(attr(terms(coxme:::formula1(my_formula)$fixed), "order"))

fixed_formula <- lme4:::getFixedFormula(my_formula)

fit0 <- survival::coxph(fixed_formula, data = ds)

start_params <- c(coef(fit0),
                  rnorm(n = length(parsed_data$reTrms$Lind),
                        mean = 0,
                        sd = sqrt(theta_start[parsed_data$reTrms$Lind])))

theta_ipl_gr(theta = unlist(coxme::VarCorr(coxfit)),
             formula = my_formula,
             parsed_data = make_ppl(parsed_data),
             other_args = list(n_fixed = n_fixed,
                               start_params = start_params,
                               stat = stat,
                               re_only = TRUE,
                               reltol = 1e-13))

# calculate the gradient for many theta combinations.

thetas <- expand.grid(theta_1 = seq(1, 4, by = 0.5),
                      theta_2 = seq(1, 4, by = 0.5))

res <- apply(thetas, 1,
      theta_ipl_gr,
      formula = my_formula,
      parsed_data = make_ppl(parsed_data),
      other_args = list(n_fixed = n_fixed,
                        start_params = start_params,
                        stat = stat,
                        re_only = TRUE,
                        reltol = 1e-13))

library(tidyverse)
data.frame(thetas, t(res)) %>%
ggplot(aes(theta_1, X1, colour = factor(theta_2))) + geom_line()

data.frame(thetas, t(res)) %>%
  ggplot(aes(theta_2, X2, colour = factor(theta_1))) + geom_line()

#######

my_formula <- survival::Surv(stat_time, stat)~X1 + X2 + X3 + (1 | M1)

ds <- one_dataset(my_formula,
                  dists = list(X1 = ~rnorm(n),
                               X2 = ~rnorm(k * nk),
                               X3 = ~rbinom(n, 1, 0.5),
                               M1 = ~rep(1:k, each = nk),
                               M2 = ~rep(c("l", "r"), ceiling(n/2))[seq_len(n)],
                               error = ~rexp(n, 10),
                               stat = ~sample(rep(c(0, 1), round(n * c(0.2, 0.8))), n)),
                  dist_args = list(k = 50, nk = 4, n = 50*4),
                  coefficients = c(1, -0.7, 0.5),
                  random_effect_variance = c(M1 = 1)
)

coxfit <- coxme::coxme(survival::Surv(stat_time, stat)~X1 + X2 + X3 + (1 | M1), data = ds)

coxme::VarCorr(coxfit)
ds_sorted <- sortAndIndex(ds, sort_vars = stat_time)

parsed_data <- lme4::lFormula(my_formula, data = ds_sorted)

stat <- Matrix(unclass(parsed_data$fr[,1])[, "status"], ncol = 1)

theta_start <- get_start_theta(length(parsed_data$reTrms$flist))

n_fixed <- length(attr(terms(coxme:::formula1(my_formula)$fixed), "order"))

fixed_formula <- lme4:::getFixedFormula(my_formula)

fit0 <- survival::coxph(fixed_formula, data = ds)

start_params <- c(coef(fit0),
                  rnorm(n = length(parsed_data$reTrms$Lind),
                        mean = 0,
                        sd = sqrt(theta_start[parsed_data$reTrms$Lind])))

theta_ipl_gr(theta = unlist(coxme::VarCorr(coxfit)),
             formula = my_formula,
             parsed_data = make_ppl(parsed_data),
             other_args = list(n_fixed = n_fixed,
                               start_params = start_params,
                               stat = stat,
                               re_only = TRUE,
                               reltol = 1e-13))

# calculate the gradient for many theta combinations.

thetas <- seq(0.5, 3, by = 0.1)

res <- sapply(thetas,
             theta_ipl_gr,
             formula = my_formula,
             parsed_data = make_ppl(parsed_data),
             other_args = list(n_fixed = n_fixed,
                               start_params = start_params,
                               stat = stat,
                               re_only = TRUE,
                               reltol = 1e-13))

data.frame(thetas, res) %>%
  ggplot(aes(thetas, res)) + geom_line()

######################

# try without local optimisation.

ests <- optim(par = start_params,
fn = lp,
gr = lp_gr,
formula = my_formula,
parsed_data = make_ppl(parsed_data),
other_args = c(list(theta = unlist(coxme::VarCorr(coxfit))), list(n_fixed = n_fixed,
                                                                  start_params = start_params,
                                                                  stat = stat,
                                                                  re_only = TRUE)),
method = "BFGS",
control = list(fnscale = -1, reltol = 1e-13),
hessian = TRUE)


theta_ipl_gr2(theta = unlist(coxme::VarCorr(coxfit)),
             formula = my_formula,
             parsed_data = make_ppl(parsed_data),
             other_args = list(n_fixed = n_fixed,
                               start_params = start_params,
                               stat = stat,
                               re_only = TRUE,
                               reltol = 1e-13),
             ests = ests)


res <- sapply(thetas,
              theta_ipl_gr2,
              formula = my_formula,
              parsed_data = make_ppl(parsed_data),
              other_args = list(n_fixed = n_fixed,
                                start_params = start_params,
                                stat = stat,
                                re_only = TRUE,
                                reltol = 1e-13),
              ests = ests)

thetas[res == min(res)]


##### calculate the integrated partial likelihood using theta_ipl,
##### return the beta and b, hessian etc at that point and use to calculate the gradient.

my_formula <- survival::Surv(stat_time, stat)~X1 + X2 + X3 + (1 | M1)

ds <- one_dataset(my_formula,
                  dists = list(X1 = ~rnorm(n),
                               X2 = ~rnorm(k * nk),
                               X3 = ~rbinom(n, 1, 0.5),
                               M1 = ~rep(1:k, each = nk),
                               M2 = ~rep(c("l", "r"), ceiling(n/2))[seq_len(n)],
                               error = ~rexp(n, 10),
                               stat = ~sample(rep(c(0, 1), round(n * c(0.2, 0.8))), n)),
                  dist_args = list(k = 50, nk = 4, n = 50*4),
                  coefficients = c(1, -0.7, 0.5),
                  random_effect_variance = c(M1 = 2)
)

coxfit <- coxme::coxme(survival::Surv(stat_time, stat)~X1 + X2 + X3 + (1 | M1), data = ds)

coxme::VarCorr(coxfit)
ds_sorted <- sortAndIndex(ds, sort_vars = stat_time)

parsed_data <- lme4::lFormula(my_formula, data = ds_sorted)

stat <- Matrix(unclass(parsed_data$fr[,1])[, "status"], ncol = 1)

theta_start <- get_start_theta(length(parsed_data$reTrms$flist))

n_fixed <- length(attr(terms(coxme:::formula1(my_formula)$fixed), "order"))

fixed_formula <- lme4:::getFixedFormula(my_formula)

fit0 <- survival::coxph(fixed_formula, data = ds)

start_params <- c(coef(fit0),
                  rnorm(n = length(parsed_data$reTrms$Lind),
                        mean = 0,
                        sd = sqrt(theta_start[parsed_data$reTrms$Lind])))

theta_max_coxme <- unlist(coxme::VarCorr(coxfit))

my_ests <- est_parameters(my_formula, data = ds,
                          control = control.list(factr = 1e3,
                                                 reltol = 1e-13,
                                                 ndeps = 1e-2))
# use my maximum.
theta_max <- my_ests$theta

theta_ipl_est <- theta_ipl(theta_max, formula = my_formula, parsed_data = make_ppl(parsed_data),
          other_args = list(n_fixed = n_fixed,
                            start_params = start_params,
                            stat = stat,
                            re_only = TRUE,
                            reltol = 1e-13))

# first derivatives of D
d_D_d_theta <- vector(mode = "list", length = length(parsed_data$reTrms$theta))

for (i in seq_along(parsed_data$reTrms$theta)) {

  d_D_d_theta[[i]] <- parsed_data$reTrms$Lambdat

  d_D_d_theta[[i]]@x <- c(1, 0)[(parsed_data$reTrms$Lind != i) + 1]

}



b <- Matrix::Matrix(attr(theta_ipl_est, "beta_b")[-seq(n_fixed)], ncol = 1)

Kbb <- attr(theta_ipl_est, "Kbb")

Kbb_inv <- solve(Kbb)

D <- attr(theta_ipl_est, "D")

D_inv <- Matrix::solve(D)

calc_grad <- function(dD){

  grad <- -0.5 * ( sum(diag( D_inv %*% dD ))
                   + sum(diag( Kbb_inv %*% D_inv %*% dD %*% D_inv ))
                   - t(b) %*% D_inv %*% dD %*% D_inv %*% b )

  grad@x

}

# debugonce(calc_grad)

calc_grad(diag(length(b)))

D[1] == theta_max

# first term
50/theta_max

sum(diag(D_inv %*% diag(length(b))))

# second term

m_temp <- D_inv %*% diag(length(b)) %*% D_inv

m_temp@x[1] - 1/theta_max^2

sum(diag(Kbb_inv %*% m_temp))

sum(diag( (1/theta_max^2) * Kbb_inv ))

sum(diag(Kbb_inv %*% D_inv %*% diag(length(b)) %*% D_inv))

# third term
t(b) %*% m_temp %*% b
sum(b^2 * (- 1/theta_max^2))

# what does the gradient look like around the theta ipl maximum. It should
# be positive and decreasing, then 0 at the maximum, then negative.

# select seq of thetas around the maximum

test_thetas <- seq(theta_max-0.1, theta_max+0.1, by = 0.01)



theta_lls <- vector(mode = "numeric", length = length(test_thetas))
grads <- vector(mode = "numeric", length = length(test_thetas))

for (i in seq_along(test_thetas)) {

  theta_ipl_est <- theta_ipl(test_thetas[i], formula = my_formula, parsed_data = make_ppl(parsed_data),
                             other_args = list(n_fixed = n_fixed,
                                               start_params = start_params,
                                               stat = stat,
                                               re_only = TRUE,
                                               reltol = 1e-13))

  theta_lls[i] <- c(theta_ipl_est)

  # first derivatives of D
  d_D_d_theta <- vector(mode = "list", length = length(parsed_data$reTrms$theta))

  for (j in seq_along(parsed_data$reTrms$theta)) {

    d_D_d_theta[[j]] <- parsed_data$reTrms$Lambdat

    d_D_d_theta[[j]]@x <- c(1, 0)[(parsed_data$reTrms$Lind != j) + 1]

  }

  b <- attr(theta_ipl_est, "beta_b")[-seq(n_fixed)]
  Kbb <- attr(theta_ipl_est, "Kbb")
  Kbb_inv <- solve(Kbb)
  D <- attr(theta_ipl_est, "D")
  D_inv <- Matrix::solve(D)

  calc_grad <- function(dD){

    grad <- -0.5 * ( sum(diag( D_inv %*% dD ))
                     + sum(diag( Kbb_inv %*% D_inv %*% dD %*% D_inv ))
                     - t(b) %*% D_inv %*% dD %*% D_inv %*% b )

    grad@x

  }

  grads[i] <- calc_grad(d_D_d_theta[[1]])

}

plot(test_thetas, theta_lls)
abline(v = theta_max)

plot(test_thetas, grads)
abline(v = theta_max, h = 0)






theta_ipl_gr_test <- function(theta){

theta_ipl_est <- theta_ipl(theta, formula = my_formula, parsed_data = make_ppl(parsed_data),
                           other_args = list(n_fixed = n_fixed,
                                             start_params = start_params,
                                             stat = stat,
                                             re_only = TRUE,
                                             reltol = 1e-13))

# first derivatives of D
d_D_d_theta <- vector(mode = "list", length = length(parsed_data$reTrms$theta))

for (j in seq_along(parsed_data$reTrms$theta)) {

  d_D_d_theta[[j]] <- parsed_data$reTrms$Lambdat

  d_D_d_theta[[j]]@x <- c(1, 0)[(parsed_data$reTrms$Lind != j) + 1]

}

b <- attr(theta_ipl_est, "beta_b")[-seq(n_fixed)]
Kbb <- attr(theta_ipl_est, "Kbb")
Kbb_inv <- solve(Kbb)
D <- attr(theta_ipl_est, "D")
D_inv <- Matrix::solve(D)

calc_grad <- function(dD){

  grad <- -0.5 * ( sum(diag( D_inv %*% dD ))
                   + sum(diag( Kbb_inv %*% D_inv %*% dD %*% D_inv ))
                   - t(b) %*% D_inv %*% dD %*% D_inv %*% b )

  grad@x

}

calc_grad(d_D_d_theta[[1]])

}


uniroot(f = theta_ipl_gr_test, interval = theta_max + c(-0.15, 0.15))
theta_max_coxme







grad <- -0.5 * ( sum(diag( D_inv %*% dD )) +
                  sum(diag( Kbb_inv %*% D_inv %*% dD %*% D_inv )) -
                  t(b) %*% D_inv %*% dD %*% D_inv %*% b )


numDeriv::grad(func = theta_ipl, x = theta_max,
               formula = my_formula, parsed_data = make_ppl(parsed_data),
               other_args = list(n_fixed = n_fixed,
                                 start_params = start_params,
                                 stat = stat,
                                 re_only = TRUE,
                                 reltol = 1e-13))




















