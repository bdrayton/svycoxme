
# generate population
k <- 2000
nk <- 10
specs = list(theta = 1)
true_coefs = c(X1 = 1, X2 = 2, X3 = -1.5, Z1 = 0.5)
the_data <- one_dataset(~X1 + X2 + X3 + Z1 + (1 | M),
                        dists = list(X1 = ~rnorm(n),
                                     X2 = ~rep(rnorm(k), each = nk),
                                     X3 = ~rep(rbinom(k, 1, 0.5), each = nk),
                                     Z1 = ~rep(rbinom(k, 1, 0.5), each = nk),
                                     M = ~rep(1:k, each = nk)),
                                     error = ~rexp(n, 10),
                                     stat = ~sample(rep(c(0, 1), round(n * c(0.2, 0.8))), n),
                        dist_args = list(k = k, nk = nk,
                                         n = k * nk),
                        coefficients = true_coefs,
                        random_effect_variance = c(M=specs$theta)
)

pop <- dplyr::mutate(the_data, id = M)

# sample from it
# one cluster sample
pop_clusters <- dplyr::distinct(pop, id, Z1) %>%
  dplyr::mutate(pr_sel = 2/75 * (1 + Z1))

sample_of_cluster = dplyr::slice_sample(pop_clusters, n = 200, weight_by = pr_sel)

my_samp <- pop[pop$id %in% sample_of_cluster$id, ]

my_samp <- dplyr::left_join(my_samp, pop_clusters, by = c("Z1", "id"))

my_samp$weights <- my_samp$pr_sel^-1

my_samp$rweights <- (1/my_samp$pr_sel)/mean(1/my_samp$pr_sel)

my_samp <- my_samp[order(my_samp$stat_time), ]

coxme_fit <- coxme::coxme(survival::Surv(stat_time, stat)~ X1 + X2 + X3 + Z1 + (1|id),
                          data = my_samp, weights = rweights)

parts <- make_parts(coxme_fit, my_samp)

ui <- calc_ui(parts)

matrix_parts <- lapply(parts, as.matrix)

ui2 <- C_calc_ui(time_start     = matrix_parts$time_start,
                 time_stop      = matrix_parts$time_stop,
                 stat           = matrix_parts$stat,
                 weights        = matrix_parts$weights,
                 exp_risk_score = matrix_parts$exp_risk_score,
                 S0             = matrix_parts$S0,
                 X              = matrix_parts$X,
                 S1_X           = matrix_parts$S1_X,
                 weighted = TRUE)

all.equal(as.matrix(ui[,1:4]), ui2, check.attributes = FALSE)

