

# 1000 random data draws with unchanging random effects

ds <- one_dataset(~X1 + X2 + X3 + (1 | M1) + (1| M2),
                  dists = list(X1 = ~rnorm(n),
                               X2 = ~rnorm(k * nk),
                               X3 = ~rbinom(n, 1, 0.5),
                               M1 = ~rep(1:k, each = nk),
                               M2 = ~rep(c("l", "r"), ceiling(n/2))[seq_len(n)],
                               error = ~rexp(n, 10),
                               stat = ~sample(rep(c(0, 1), round(n * c(0.2, 0.8))), n)),
                  dist_args = list(k = 50, nk = 10, n = 500),
                  coefficients = c(1, 1, 1),
                  random_effect_variance = c(M1 = 1, M2 = 2),
                  random_effect_seed = c(M1 = 125, M2 = 457)
)



res_list <- lapply(1:1000, function(i){

  ds <- one_dataset(~X1 + X2 + X3 + (1 | M1) + (1| M2),
                  dists = list(X1 = ~rnorm(n),
                               X2 = ~rnorm(k * nk),
                               X3 = ~rbinom(n, 1, 0.5),
                               M1 = ~rep(1:k, each = nk),
                               M2 = ~rep(c("l", "r"), ceiling(n/2))[seq_len(n)],
                               error = ~rexp(n, 10),
                               stat = ~sample(rep(c(0, 1), round(n * c(0.2, 0.8))), n)),
                  dist_args = list(k = 50, nk = 10, n = 500),
                  coefficients = c(1, 1, 1),
                  random_effect_variance = c(M1 = 1, M2 = 2),
                  random_effect_seed = c(M1 = 125, M2 = 457)
  )

  fit <- coxme::coxme(Surv(stat_time, stat) ~ (1 | M1) + (1 | M2), data = ds)

  # attr(ds, "random_effects")
  coxme::ranef(fit)

})

re1 <- sapply(res_list, "[[", "M1", simplify = "matrix")

re2 <- sapply(res_list, "[[", "M2", simplify = "matrix")



# should be about 0
sum(attr(ds, "random_effects")$M1)
sum(rowMeans(x = re1))


cbind(
  attr(ds, "random_effects")$M1,
  rowMeans(x = re1))

cbind(
  attr(ds, "random_effects")$M2,
  rowMeans(x = re2))


# bootstrap a theta CI
ds <- one_dataset(~X1 + X2 + X3 + (1 | M1) + (1| M2),
                  dists = list(X1 = ~rnorm(n),
                               X2 = ~rnorm(k * nk),
                               X3 = ~rbinom(n, 1, 0.5),
                               M1 = ~rep(1:k, each = nk),
                               M2 = ~rep(c("l", "r"), ceiling(n/2))[seq_len(n)],
                               error = ~rexp(n, 10),
                               stat = ~sample(rep(c(0, 1), round(n * c(0.2, 0.8))), n)),
                  dist_args = list(k = 50, nk = 10, n = 500),
                  coefficients = c(1, 1, 1),
                  random_effect_variance = c(M1 = 1, M2 = 2),
                  random_effect_seed = c(M1 = 125, M2 = 457)
)


theta.fun <- function(dat, original_fit){

  n <- nrow(dat)

  resampled_data <- dat[sample(n, size = n, replace = TRUE), ]

  new_fit <- update(original_fit, data = resampled_data)

  unlist(coxme::VarCorr(new_fit))

}

fit <- coxme::coxme(Surv(stat_time, stat) ~ (1 | M1) + (1 | M2), data = dat)

theta_boot_reps <- replicate(1000, theta.fun(dat = ds, original_fit = fit), simplify = "matrix")

data.frame(mean_est=rowMeans(theta_boot_reps),
           t(apply(theta_boot_reps, 1, quantile, c(0.025, 0.975))))






