

# test new estimation setup.

library(Matrix)

my_formula <- survival::Surv(stat_time, stat)~X1 + X2 + X3 + (1 | M)
my_k = 20
my_nk = 4
my_theta = 1
my_beta = c(1, -0.7, 0.5)

ds <- one_dataset(my_formula,
                  dists = list(X1 = ~rnorm(n),
                               X2 = ~rnorm(k * nk),
                               X3 = ~rbinom(n, 1, 0.5),
                               M = ~rep(1:k, each = nk),
                               error = ~rexp(n, 10),
                               stat = ~sample(rep(c(0, 1), round(n * c(0.2, 0.8))), n)),
                  dist_args = list(k = my_k, nk = my_nk, n = my_k*my_nk),
                  coefficients = my_beta,
                  random_effect_variance = list(M = my_theta)
)


est_parameters(my_formula, ds)




coxme_fit <- coxme::coxme(my_formula, data = ds)

coxme::VarCorr(coxme_fit)

coxme::random.effects(coxme_fit)
coxme::fixed.effects(coxme_fit)



