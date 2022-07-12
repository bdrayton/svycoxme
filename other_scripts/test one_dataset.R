# test one_dataset

# test with three fixed and one random effects
debugonce(one_dataset)
one_dataset(~X1 + X2 + X3 + (1 | M),
            dists = list(X1 = ~rnorm(n),
                         X2 = ~rnorm(k * nk),
                         X3 = ~rbinom(n, 1, 0.5),
                         M = ~rep(1:k, each = nk),
                         error = ~rexp(n, 10),
                         stat = ~sample(rep(c(0, 1), round(n * c(0.2, 0.8))), n)),
            dist_args = list(k = 50, nk = 4, n = 200),
            coefficients = c(1, 1, 1),
            random_effect_variance = list(M = 1)
)

# test with three fixed and two random effects
one_dataset(~X1 + X2 + X3 + (1 | M1) + (1| M2),
            dists = list(X1 = ~rnorm(n),
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

# test with no random effects
# debugonce(one_dataset)
one_dataset(~X1 + X2 + X3,
            dists = list(X1 = ~rnorm(n),
                         X2 = ~rnorm(k * nk),
                         X3 = ~rbinom(n, 1, 0.5),
                         error = ~rexp(n, 10),
                         stat = ~sample(rep(c(0, 1), round(n * c(0.2, 0.8))), n)),
            dist_args = list(k = 50, nk = 4, n = 200),
            coefficients = c(1, 1, 1)
)


# test with no fixed effects
# debugonce(one_dataset)
one_dataset(~(1 | M),
            dists = list(M = ~rep(1:k, each = nk),
                         error = ~rexp(n, 10),
                         stat = ~sample(rep(c(0, 1), round(n * c(0.2, 0.8))), n)),
            dist_args = list(k = 50, nk = 4, n = 200),
            random_effect_variance = list(M = 1)
)

# test with no fixed or random effects

# debugonce(one_dataset)
one_dataset(~1,
            dists = list(error = ~rexp(n, 10),
                         stat = ~sample(rep(c(0, 1), round(n * c(0.2, 0.8))), n)),
            dist_args = list(k = 50, nk = 4, n = 200)
)


