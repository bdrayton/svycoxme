

k <- 50
nk <- 10

theta = 1

true_coefs = c(X1 = 1, X2 = 0.5, X3 = 1.5, Z1 = 2)

debugonce(one_dataset)

the_data <- one_dataset(~X1 + X2 + X3 + Z1 + (1 | M),
                        dists = list(Z1 = ~rnorm(n),
                                     X1 = ~rnorm(n),
                                     X1Z = ~rnorm(n, mean = Z1),
                                     X2 = ~rep(rnorm(k), each = nk),
                                     X3 = ~rep(rbinom(k, 1, 0.5), each = nk),
                                     M = ~rep(1:k, each = nk)),
                        error = ~rexp(n, 10),
                        stat = ~sample(rep(c(0, 1), round(n * c(0.2, 0.8))), n),
                        dist_args = list(k = k, nk = nk,
                                         n = k * nk),
                        coefficients = true_coefs,
                        random_effect_variance = c(M=theta)
)

head(data)
