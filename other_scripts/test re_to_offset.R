


#
specs <- list(my_formula = survival::Surv(stat_time, stat)~X1 + X2 + X3 + (1 | M.1) + (1 | M.2) + (1 | M.2:M.1),
              my_k = 50,
              my_nk = 10,
              my_theta = c(M.1 = 2, M.2 = 1, `M.2:M.1` = 0.5),
              my_beta = c(1, -0.7, 0.5),
              my_ndeps = rep(0.001, 2))



the_data <- one_dataset(specs$my_formula,
                        dists = list(X1 = ~rnorm(n),
                                     X2 = ~rep(rnorm(k), each = nk),
                                     X3 = ~rep(rbinom(k, 1, 0.5), each = nk),
                                     M.1 = ~rep(1:k, each = nk),
                                     M.2 = ~rep(c("l", "r"), ceiling(n/2))[seq_len(n)],
                                     error = ~rexp(n, 10),
                                     stat = ~sample(rep(c(0, 1), round(n * c(0.2, 0.8))), n)),
                        dist_args = list(k = specs$my_k, nk = specs$my_nk,
                                         n = specs$my_k * specs$my_nk),
                        coefficients = specs$my_beta,
                        random_effect_variance = specs$my_theta)

the_data$`M.1:M.2` <- with(the_data, interaction(M.1, M.2))



coxfit <- coxme::coxme(survival::Surv(stat_time, stat)~X1 + X2 + X3 + (1 | M.1) + (1 | M.2) + (1 | `M.1:M.2`),
                       data = the_data)

the_data$offset <- re_to_offset(the_data, coxfit)

coxfit2 <- survival::coxph(survival::Surv(stat_time, stat)~X1 + X2 + X3 + offset(offset), data = the_data, robust = TRUE, cluster = `M.1:M.2`)

# coefficients are about the same
all.equal(coef(coxfit),coef(coxfit2))

# standard errors are different
all.equal(confint(coxfit2),
          confint(coxfit))

sqrt(diag(vcov(coxfit)))
coef(coxfit) + 1.96 * sqrt(diag(vcov(coxfit2)))

#
specs <- list(my_formula = survival::Surv(stat_time, stat)~X1 + X2 + X3 + (1 | M2) + (1 | M2/M1),
              my_k = 50,
              my_nk = 10,
              my_theta = c(M1 = 2, M2 = 1, `M1:M2` = 0.5),
              my_beta = c(1, -0.7, 0.5),
              my_ndeps = rep(0.001, 2))


the_data <- one_dataset(specs$my_formula,
                        dists = list(X1 = ~rnorm(n),
                                     X2 = ~rep(rnorm(k), each = nk),
                                     X3 = ~rep(rbinom(k, 1, 0.5), each = nk),
                                     M1 = ~rep(1:k, each = nk),
                                     M2 = ~rep(c("l", "r"), ceiling(n/2))[seq_len(n)],
                                     error = ~rexp(n, 10),
                                     stat = ~sample(rep(c(0, 1), round(n * c(0.2, 0.8))), n)),
                        dist_args = list(k = specs$my_k, nk = specs$my_nk,
                                         n = specs$my_k * specs$my_nk),
                        coefficients = specs$my_beta,
                        random_effect_variance = specs$my_theta)

# the_data$`M.1:M.2` <- with(the_data, interaction(M.1, M.2))
#

coxfit <- coxme::coxme(survival::Surv(stat_time, stat)~X1 + X2 + X3 + (1 | M1) + (1 | M2/M1),
                       data = the_data)

the_data$offset <- re_to_offset(the_data, coxfit)

coxfit2 <- survival::coxph(survival::Surv(stat_time, stat)~X1 + X2 + X3 + offset(offset), data = the_data)

all.equal(coxfit$coefficients,coxfit2$coefficients)


#
specs <- list(my_formula = survival::Surv(stat_time, stat)~X1 + X2 + X3 + (1 | M2) + (1 | M2/M1) + (1 | M3),
              my_k = 50,
              my_nk = 10,
              my_theta = c(M1 = 2, M2 = 1, `M1:M2` = 0.5, M3 = 0.5),
              my_beta = c(1, -0.7, 0.5),
              my_ndeps = rep(0.001, 2))


the_data <- one_dataset(specs$my_formula,
                        dists = list(X1 = ~rnorm(n),
                                     X2 = ~rep(rnorm(k), each = nk),
                                     X3 = ~rep(rbinom(k, 1, 0.5), each = nk),
                                     M1 = ~rep(1:k, each = nk),
                                     M2 = ~rep(c("l", "r"), ceiling(n/2))[seq_len(n)],
                                     M3 = ~rep(1:nk, each = k),
                                     error = ~rexp(n, 10),
                                     stat = ~sample(rep(c(0, 1), round(n * c(0.2, 0.8))), n)),
                        dist_args = list(k = specs$my_k, nk = specs$my_nk,
                                         n = specs$my_k * specs$my_nk),
                        coefficients = specs$my_beta,
                        random_effect_variance = specs$my_theta)

# the_data$`M.1:M.2` <- with(the_data, interaction(M.1, M.2))
#

coxfit <- coxme::coxme(survival::Surv(stat_time, stat)~X1 + X2 + X3 + (1 | M1) + (1 | M2/M1) + (1 | M3),
                       data = the_data)

the_data$offset <- re_to_offset(the_data, coxfit)

coxfit2 <- survival::coxph(survival::Surv(stat_time, stat)~X1 + X2 + X3 + offset(offset), data = the_data)

all.equal(coxfit$coefficients,coxfit2$coefficients)





