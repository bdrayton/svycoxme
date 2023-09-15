# test the syvcoxme methods


# these are the parameter combinations
true_coefs = c(X1 = 1, X2 = -0.7, X3 = 0.5, Z1 = log(2))

the_data <- one_dataset(~X1 + X2 + X3 + Z1 + (1 | M),
                        dists = list(X1 = ~rnorm(n),
                                     X2 = ~rep(rnorm(k), each = nk),
                                     X3 = ~rep(rbinom(k, 1, 0.5), each = nk),
                                     Z1 = ~rep(rbinom(k, 1, 0.5), each = nk),
                                     M = ~rep(1:k, each = nk),
                                     error = ~rexp(n, 10),
                                     stat = ~sample(rep(c(0, 1), round(n * c(0.2, 0.8))), n)),
                        dist_args = list(k = 5000, nk = 10,
                                         n = 5000 * 10),
                        coefficients = true_coefs,
                        random_effect_variance = c(M=1)
)

pop <- dplyr::mutate(the_data, id = M)

# sample from it
# one cluster sample
pop_clusters <- dplyr::select(pop, id, Z1) %>%
  dplyr::mutate(pr_sel = 2/75 * (1 + Z1))

sample_of_cluster = dplyr::slice_sample(pop_clusters, n = 100, weight_by = pr_sel)

my_samp <- pop[pop$id %in% sample_of_cluster$id, ]

my_samp <- dplyr::left_join(my_samp, pop_clusters, by = c("Z1", "id"))

my_samp$weights <- my_samp$pr_sel^-1

# rescale_weights
my_samp$rweights <- (1/my_samp$pr_sel)/mean(1/my_samp$pr_sel)

my_samp <- my_samp[order(my_samp$stat_time), ]

form <- survival::Surv(stat_time, stat)~ X1 + X2 + X3 + Z1 + (1|id)

coxme_fit <- coxme::coxme(form, data = my_samp, weights = rweights)

summary(coxme_fit)

my_des <- svydesign(~id, weights = ~weights, data = my_samp)
my_des_jackknife <- as.svrepdesign(my_des, type = "JK1")
my_des_bootstrap <- as.svrepdesign(my_des, type = "bootstrap")

svycoxme_fit <- eval(bquote( svycoxme(.(form), des = my_des) ))
svycoxme_fit_jackknife <- eval(bquote( svycoxme(form, des = my_des_jackknife) ))
svycoxme_fit_bootstrap <- eval(bquote( svycoxme(form, des = my_des_bootstrap) ))

debugonce(svycoxme::summary.svycoxme)

test_fit <- svycoxph(survival::Surv(stat_time, stat)~ X1 + X2 + X3 + Z1, design = my_des)


survey:::print.svycoxph

class(test_fit)

debugonce(coxme:::summary.coxme)
debugonce(coxme:::print.coxme)
summary(svycoxme_fit)


summary(test_fit)

svycoxme_fit$survey.design

survey:::summary.svycoxph

summary.svycoxme

methods(summary)



summary(svycoxme_fit)

summary(svycoxme_fit_jackknife)
summary(svycoxme_fit_jackknife)










