
# look at how coxme interprets random effects formulas

ds <- one_dataset(~X1 + X2 + X3 + (1 | M1) + (1| M2),
                  dists = list(X1 = ~rnorm(n),
           /
             X2 = ~rnorm(k * nk),
                               X3 = ~rbinom(n, 1, 0.5),
                               M1 = ~rep(1:k, each = nk),
                               M2 = ~rep(c("l", "r"), ceiling(n/2))[seq_len(n)],
                               error = ~rexp(n, 10),
                               stat = ~sample(rep(c(0, 1), round(n * c(0.2, 0.8))), n)),
                  dist_args = list(k = 50, nk = 4, n = 200),
                  coefficients = c(1, 1, 1),
                  random_effect_variance = list(M1 = 1, M2 = 2)
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









