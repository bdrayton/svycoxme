

my_beta = c(1, -0.7, 0.5)
my_theta = 0.2
my_k = 5
my_nk = 10

my_X = c("X1", "X2", "X3")

sample_data <- one_dataset(control = list(k = my_k, nk = my_nk, beta = my_beta, theta = my_theta))

b <- attr(sample_data, "random_effects")

my_params <- c(my_beta[seq_along(my_X)], b)

BB(parms = my_params, X = my_X, t = t, cluster = "M", dij = stat,
   data = sample_data)

bb(parms = my_params, X = my_X, t = t, cluster = "M", dij = stat,
   data = sample_data, theta = my_theta)

Bb(parms = my_params, X = my_X, t = t, cluster = "M", dij = stat,
   data = sample_data)

my_hessian <- ppl_hessian(parms = my_params, X = my_X, t = t,
                          cluster = "M", dij = stat, data = sample_data,
                          theta = my_theta)


fit <- survival::coxph(survival::Surv(t, stat) ~ X1 + X2 + X3, data = sample_data)

start_parameters = c(coef(fit), rep(0, length(b)))

names(start_parameters) <- c(my_X, paste0("Z", seq_len(length(b))))

fit_optim <- optim(par = start_parameters,
                   fn = lp,
                   gr = lp_grd,
                   X = my_X,
                   t = t,
                   cluster = "M",
                   dij = stat,
                   D = my_theta * diag(length(b)),
                   data = sample_data,
                   method = "BFGS",
                   control = list(fnscale = -1),
                   hessian = TRUE)

fit_optim$hessian

my_hessian <- ppl_hessian(parms = fit_optim$par, X = my_X, t = t,
                          cluster = "M", dij = stat, data = sample_data,
                          theta = my_theta)

my_hessian - fit_optim$hessian

BB(parms = fit_optim$par, X = my_X, t = t, cluster = "M", dij = stat,
   data = sample_data)


bb(parms = fit_optim$par, X = my_X, t = t, cluster = "M", dij = stat,
   data = sample_data, theta = my_theta)


Bb(parms = fit_optim$par, X = my_X, t = t, cluster = "M", dij = stat,
   data = sample_data)





diag(fit_optim$hessian)
diag(my_hessian)

sqrt(diag(solve(-my_hessian)))
sqrt(diag(solve(-fit_optim$hessian)))

coxme::coxme(survival::Surv(t, stat) ~ X1 + X2 + X3 + (1|M), data = sample_data)

# ok so the question is how to get the correct off diagonal terms in the second derivatives.
# I might be doing the wrong matrix operation, the dot product, so I will try others to see if
# something else works.


sorted_data <- sortAndIndex(sample_data, t)

data_with_Z <- add_Z(sorted_data, "M")

data_with_lp <- calcLinearPredictor(data = data_with_Z, X = my_X, Z = attr(data_with_Z, "Z_names"), parms = fit_optim$par)

data_with_risksets <- calcRiskSets(data_with_lp)

X_risksets <- calcRiskSets(data_with_lp, vars = my_X, varCol = "Xr")

head(X_risksets)

X_matrix <- data_with_lp[,my_X, drop = TRUE] %>% as.matrix()

XtX <- lapply(seq_len(nrow(X_matrix)), function(i){

  X_matrix[i, ] %*% t(X_matrix[i, ])

})

A <- data_with_lp$A

XtX_A <- mapply("*", XtX, A, SIMPLIFY = FALSE)

cumsum_XtX_A <- lapply(seq_along(XtX_A), function(i){
  Reduce('+', XtX_A[i:length(XtX_A)])
})

cumsum_A <- data_with_risksets$cumsum_A

term2 <- mapply('/', cumsum_XtX_A, cumsum_A, SIMPLIFY = FALSE)

cumsum_A_squared = cumsum_A^2

split_risksets <- split(X_risksets$cumsum_Xr_A, X_risksets$index)

term1 <- lapply(seq_along(split_risksets), function(i){

  (split_risksets[[i]] %*% t(split_risksets[[i]])) / cumsum_A_squared[i]

})

parts <- mapply('-', term1, term2, SIMPLIFY = FALSE)

Reduce('+', parts[data_with_risksets$stat==1])

fit_optim$hessian[my_X, my_X]



term1_list[[1]] - term2[[1]]

parts[[1]]

rbind(
  term1_list_t[[1]],
  term1_list_t[[1]],
  term1_list_t[[1]]) - term2[[1]]




