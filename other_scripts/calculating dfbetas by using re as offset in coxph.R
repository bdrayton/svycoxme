

coxme_fit_d2 <- coxme(Surv(stat_time, stat) ~ X1 + X2 + X3 + (1 | M), data = sample_data)

dfbetas(coxme_fit_d2)

# get df betas by fitting a coxph model with clustering and offsets.

library(coxme)

fit2 <- coxme(Surv(time, status) ~ ph.ecog + age + (1|inst), lung)

offsets <- data.frame(value = ranef(fit2)$inst) |>
  tibble::rownames_to_column(var = "inst")

lung2 <- dplyr::mutate(lung, inst = as.character(inst))

lung3 <- dplyr::left_join(lung2, offsets, by = "inst")

fit1 <- coxph(Surv(time, status) ~ ph.ecog + age + offset(value), lung3)

coef(fit1) - coef(fit2)

dfbetas <- resid(fit1, type = "dfbetas")

# try with clustering

fit3 <- coxph(Surv(time, status) ~ ph.ecog + age + offset(value), lung3, cluster = inst)

dfbetas2 <- resid(fit3, type = "dfbetas")

all.equal(dfbetas, dfbetas2)

# with clustering and weighting

weights <- tibble::tibble(inst = as.character(unique(lung$inst))) |>
  dplyr::mutate(weight = round(runif(dplyr::n(), 5, 10)))

lung4 <- dplyr::left_join(lung3, weights, by = "inst")

fit4 <- coxph(Surv(time, status) ~ ph.ecog + age + offset(value), data = lung4, cluster = inst, weights = weight)

dfbetas3 <- resid(fit4, type = "dfbetas", weighted = TRUE)

all.equal(dfbetas, dfbetas3)


# need to turn random effects into an offset term

# need the variables used for random effects from the dataset.
# flatten out each effect.
# sum together appropriately.


ranefs <- ranef(fit2)
ranef_names <- names(ranefs)

# interactions are the difficult bit

re_df <- lapply(ranefs, function(x){

  data.frame(level = names(x), value = x)

})



# need to do something about names like v1:v2




data.frame(ranef(fit2))


test <- one_dataset(~X1 + X2 + X3 + (1 | M1) + (1| M1/M2),
                    dists = list(X1 = ~rnorm(n),
                                 X2 = ~rnorm(k * nk),
                                 X3 = ~rbinom(n, 1, 0.5),
                                 M1 = ~rep(1:k, each = nk),
                                 M2 = ~rep(c("l", "r"), ceiling(n/2))[seq_len(n)],
                                 error = ~rexp(n, 10),
                                 stat = ~sample(rep(c(0, 1), round(n * c(0.2, 0.8))), n)),
                    dist_args = list(k = 50, nk = 4, n = 200),
                    coefficients = c(1, 1, 1),
                    random_effect_variance = list(M1 = 1, M2 = 0.5, `M2:M1` = 2)
)

test$`M2.M1` = with(test, interaction(M1, M2))

fit_test <- coxme(Surv(stat_time, stat) ~ X1 + X2 + X3 + (1 | M1) + (1 | M2) + (1 | M2.M1), data = test)

res <- ranef(fit_test)
re_names <- names(res)

test4 <- test
test4$offset = 0

for(i in seq_along(res)){

  re <- data.frame(new_offset = res[[i]])

  re_class <- class(test4[,re_names[i]])

  if (re_class == "factor"){
    re[,re_names[i]] <- factor(rownames(re))
  } else {
    re[,re_names[i]] <- as(rownames(re), Class = re_class)
  }

  test4 <- dplyr::left_join(test4, re, by = re_names[i])

  test4$offset = test4$offset + test4$new_offset

  test4$new_offset <- NULL

}

# the coxph with offset should estimate the same coefs as coxme with the random effects.
fit_test4 <- coxph(Surv(stat_time, stat) ~ X1 + X2 + X3 + offset(offset), data = test4)


cbind(
  coef(fit_test),
  coef(fit_test4))

# try by hand
re1 <- data.frame(offset_1 = res[[1]])
re1$M1 <- as.integer(rownames(re1))

re2 <- data.frame(offset_2 = res[[2]])
re2$M2 <- rownames(re2)

test3 <- dplyr::left_join(
           dplyr::left_join(test, re1, by = "M1"),
           re2, by = "M2")

fit_test3 <- coxph(Surv(stat_time, stat) ~ X1 + X2 + X3 + offset(offset_1 + offset_2), data = test3)

# They are the same.
cbind(
  coef(fit_test),
  coef(fit_test2),
  coef(fit_test3))
















