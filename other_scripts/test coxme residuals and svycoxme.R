
# test the residuals.coxme and svycoxme.survey.design2


# generate population
nk <- 10
specs <- list()
specs$theta <- 1
specs$n_clusters = 5000
k = specs$n_clusters
specs$n_clusters_in_sample = 100
true_coefs <- c(X1 = 1, X2 = -0.7, X3 = 0.5)

the_data <- one_dataset(~X1 + X2 + X3 + (1 | M),
                        dists = list(X1 = ~rnorm(n),
                                     X2 = ~rep(rnorm(k), each = nk),
                                     X3 = ~rep(rbinom(k, 1, 0.5), each = nk),
                                     M = ~rep(1:k, each = nk)),
                        error = ~rexp(n, 10),
                        censoring_time = ~runif(n, 0, 0.8),
                        dist_args = list(k = k, nk = nk,
                                         n = k * nk),
                        coefficients = true_coefs,
                        random_effect_variance = c(M=specs$theta)
)

pop <- dplyr::mutate(the_data, id = M)

# sample from it
# one cluster sample
samp_cluster_ids <- unique(pop$id)[sample.int(specs$n_clusters, specs$n_clusters_in_samp)]

my_samp <- pop[pop$id %in% samp_cluster_ids, ]

# for memory reasons
# rm(list = c('pop', 'the_data'))

my_samp$prob <- (specs$n_clusters_in_samp/specs$n_clusters)
my_samp$weights <- my_samp$prob^-1

# rescale_weights
my_samp$rweights <- (1/my_samp$prob)/mean(1/my_samp$prob)

my_samp <- my_samp[order(my_samp$stat_time), ]

# the regular fit.
coxme_fit <- coxme::coxme(survival::Surv(stat_time, stat)~ X1 + X2 + X3 + (1|id),
                          data = my_samp, weights = rweights)

# calculate ui and get information
# parts <- make_parts(coxme_fit, my_samp)

# debugonce(residuals.coxme)
debugonce(make_parts.coxme)
my_resids = resid(coxme_fit, data = my_samp, type = 'dfbeta')

coxme_fit_pop  <- coxme::coxme(survival::Surv(stat_time, stat)~ X1 + X2 + X3 + (1|id),
                               data = pop)

debugonce(residuals.coxme)
debugonce(make_parts.coxme)
my_resids_pop = resid(coxme_fit_pop, data = pop, type = 'dfbeta')

ui <- calc_ui(parts)
vv <- get_information(coxme_fit)

uivv <- ui %*% vv


# define design, add ui
my_des <- svydesign(~id, weights = ~weights, data = my_samp)
my_des_jackknife <- as.svrepdesign(my_des, type = "JK1")
my_des_bootstrap <- as.svrepdesign(my_des, type = "bootstrap")
# my_des_subbootstrap <- as.svrepdesign(my_des, type = "subbootstrap")

svycoxme_fit_jackknife <- svycoxme(survival::Surv(stat_time, stat) ~ X1 + X2 + X3 + (1 | id), des = my_des_jackknife)
debugonce(svycoxme)
svycoxme_fit <- svycoxme(survival::Surv(stat_time, stat) ~ X1 + X2 + X3 + (1 | id), des = my_des)

svycoxme_fit_jackknife$var - svycoxme_fit$var

svycoxme_fit_bootstrap <- svycoxme(survival::Surv(stat_time, stat) ~ X1 + X2 + X3 + (1 | id), des = my_des_bootstrap)

svycoxme_fit_bootstrap$var










