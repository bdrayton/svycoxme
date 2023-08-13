# The algorithm for fittin beta b and theta in coxme is to have a likelihood function for theta that
# 1) optimises beta and b for a given theta,
# 2) caluculates the ipl for estimated beta and b, and given theta.
# Then use this function to optimise theta, beta, and b simulataneously.
# Previous attempts by me to do this didn't work, so I reverted to the back and forth method.
# Here I manage to get it working with svycoxme, and compare the ipl curve from coxme to points
# along the curve calculated using modified svycoxme functions.
#

# for theta fixed ie given, estimate b and  beta, and return the integrated log likelihood.

ds <- one_dataset(~X1 + X2 + X3 + (1 | M),
            dists = list(X1 = ~rnorm(n),
                         X2 = ~rnorm(k * nk),
                         X3 = ~rbinom(n, 1, 0.5),
                         M = ~rep(1:k, each = nk),
                         error = ~rexp(n, 10),
                         stat = ~sample(rep(c(0, 1), round(n * c(0.2, 0.8))), n)),
            dist_args = list(k = 50, nk = 4, n = 200),
            coefficients = c(1, 1, 1),
            random_effect_variance = list(M = 1.5)
)


thetas <- seq(0.001, 3, by = 0.01)

ll <- sapply(thetas, function(theta) {

  fit <- coxme::coxme(survival::Surv(stat_time, stat) ~ X1 + X2 + X3 + (1 | M), data = ds, vfixed = list(M = theta))

  fit$loglik[2]

})



# does this work with svycoxme?

# For a given theta, estimate beta and b, calculate the ill.

fit0 <- survival::coxph(survival::Surv(stat_time, stat) ~ X1 + X2 + X3, data = ds)

my_start_parameters <- c(coef(fit0), rep(0, 50))
names(my_start_parameters) <- c(c("X1", "X2", "X3"), paste0("Z", seq_len(50)))


estimate_parameters_temp <- function(start_theta,
                                    start_parms,
                                    X,
                                    stat_time,
                                    cluster,
                                    dij,
                                    data,
                                    theta_hessian = FALSE) {

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

  fit_beta_b$par

}

# run this for some thetas. fewer cos it's slow.

some_thetas <- seq(0.001, 3, by = 0.5)

ll2 <- sapply(some_thetas, function(one_theta){

  r1 <- estimate_parameters_temp(start_theta = one_theta,
                           start_parms = my_start_parameters,
                           X = c("X1", "X2", "X3"),
                           stat_time = stat_time,
                           cluster = "M",
                           dij = stat,
                           data = ds)

  theta_ipl(one_theta, r1, c("X1", "X2", "X3"), stat_time, stat, "M", ds)

})



# try and do the nested optim thingo

theta_ipl_temp <- function(one_theta, start_parms, data){

  r1 <- estimate_parameters_temp(start_theta = one_theta,
                                 start_parms = start_parms,
                                 X = c("X1", "X2", "X3"),
                                 stat_time = stat_time,
                                 cluster = "M",
                                 dij = stat,
                                 data = data)

  theta_ipl(one_theta, r1, c("X1", "X2", "X3"), stat_time, stat, "M", data)

}

theta_ipl_temp(0.5, start_parms = my_start_parameters, data = ds)

theta_est <- optim(par = 0.5,
      fn = theta_ipl_temp,
      gr = NULL,
      start_parms = my_start_parameters,
      data = ds,
      method = "L-BFGS-B",
      control = list(fnscale = -1),
      lower = 0.00001,
      upper = 3)

# coxme integrated partial log likelihood, by theta.
plot(thetas, ll, type = "l")
max_theta <- thetas[ll == max(ll)]
# coxme theta hat
abline(v = max_theta)

# some points along the ipl using my functions (slow so not the whole line)
points(some_thetas, ll2, col = "red")

# theta hat, using the nested optim function above
abline(v = theta_est$par, col = "red")

# this method is know to underestimate theta on average. Additionally, my
# method is (always?) estimates a slightly bigger theta, so is less biased on average.










