
# look at how coxme interprets random effects formulas

debugonce(one_dataset)
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
                  random_effect_seed = c(M1 = 125, M2 = 457),
                  seed = 10
)

ds2 <- one_dataset(~X1 + X2 + X3 + (1 | M1) + (1| M2),
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
                  random_effect_seed = c(M1 = 123, M2 = 456),
                  seed = 10
)

all(attr(ds2, "random_effects")$M1 == attr(ds, "random_effects")$M1)

all(attr(ds2, "random_effects")$M2 == attr(ds, "random_effects")$M2)

table(ds$X1 == ds2$X1)

table(ds$stat == ds2$stat)

table(ds$stat_time == ds2$stat_time)



fit <- coxme(Surv(stat_time, stat) ~ (1 | M1) + (1 | M2), data = ds)

VarCorr(fit)


cbind(
  attr(ds, "random_effects")$M1,
  ranef(fit)$M1
)





head(ds)

ds$left <- (ds$M2 == "l") + 0
ds$right <- (ds$M2 == "r") + 0

library(coxme)

coxme(Surv(stat_time, stat) ~ (1 | M1), data = ds)

coxme(Surv(stat_time, stat) ~ (1 | M1) + (1 | left), data = ds)

coxme(Surv(stat_time, stat) ~ (1 | M1) + (1 | left) + (1 | right), data = ds)

fit <- coxme(Surv(stat_time, stat) ~ (1 | M1) + (1 | left/M1) + (1 | right/M1), data = ds)

ranef(fit)

my_k = 50
my_nk = 4
my_n = my_nk * my_k


survey:::svycoxph.survey.design









