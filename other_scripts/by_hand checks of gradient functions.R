

#' check the parts of the gradient likelihood functions. Trying to determine if there is a mistake,
#' as my gradients are different to coxme.

my_k = 50
my_theta = 1

ds <- one_dataset(list(k = my_k, nk = 10, beta = c(1, -0.7, 0.5), theta = my_theta))

fit <- coxme::coxme(survival::Surv(t, stat) ~ X1 + X2 + X3 + (1|M), data = ds)

est_theta <- coxme::VarCorr(fit)$M

D = est_theta * diag(my_k)

est_parms <- c(coxme::fixef(fit), coxme::ranef(fit)$M)

names(est_parms) <- paste0(rep(c("X", "Z"), c(3, my_k)), c(1:3, seq_len(my_k)))

my_loglik <- lp(parms = est_parms,
                X = c("X1", "X2", "X3"),
                t = t, dij = stat, D =  D, data = ds)

(fit$loglik["Penalized"] - my_loglik) < 1e-10

fit$penalty == attr(my_loglik, "penalty")


my_u <- lp_grd(parms = est_parms,
               X = c("X1", "X2", "X3"),
               cluster = "M",
               t = t,
               dij = stat,
               D =  D,
               data = ds)

cbind(my_u, fit$u)

sum(my_u) - sum(fit$u)

X <- c("X1", "X2", "X3")

Z_formula <- formula(glue::glue(" ~ as.factor( M ) - 1"))

Z_matrix <- model.matrix(Z_formula, data = ds)

Z_names <- paste0("Z", seq(ncol(Z_matrix)))

colnames(Z_matrix) <- Z_names

ds_with_Z <- dplyr::bind_cols(ds, data.frame(Z_matrix))

b <- est_parms[-seq_len(length.out = length(X))]

names(b) <- Z_names

B <- est_parms[seq_len(length.out = length(X))]

sortedIndexedData <- sortAndIndex(ds_with_Z, sort_vars = t )

addedLP <- calcLinearPredictor(sortedIndexedData, X = X, Z = colnames(Z_matrix), parms = est_parms)

dlp_beta(est_parms,
         X = c("X1", "X2", "X3"),
         cluster = "M",
         t = t,
         dij = stat,
         D =  D,
         data = ds)


by_hand_lp <- as.matrix(sortedIndexedData[, c("X1", "X2", "X3")]) %*% B + b[sortedIndexedData$M]

cumsum_exp_lp = rev(cumsum(exp(rev(by_hand_lp))))


cumsum_X1_exp_lp = rev(cumsum(rev(sortedIndexedData$X1 * exp(by_hand_lp))))
cumsum_X2_exp_lp = rev(cumsum(rev(sortedIndexedData$X2 * exp(by_hand_lp))))
cumsum_X3_exp_lp = rev(cumsum(rev(sortedIndexedData$X3 * exp(by_hand_lp))))

sum(sortedIndexedData$stat * (sortedIndexedData$X1 - cumsum_X1_exp_lp/cumsum_exp_lp))
sum(sortedIndexedData$stat * (sortedIndexedData$X2 - cumsum_X2_exp_lp/cumsum_exp_lp))
sum(sortedIndexedData$stat * (sortedIndexedData$X3 - cumsum_X3_exp_lp/cumsum_exp_lp))


cumsum_Z1_exp_lp = rev(cumsum(rev(sortedIndexedData$Z1 * exp(by_hand_lp))))

Z_penalties <- b %*% solve(D)

sum(sortedIndexedData$stat * (sortedIndexedData$Z1 - cumsum_Z1_exp_lp/cumsum_exp_lp)) - Z_penalties[1]


sort(fit$u[4:53])

Z_penalties <- b %*% solve(D)

names(Z_penalties) <- Z_names

byhand_grads <- sapply(Z_names, function(one_Z){

  cumsum_Z_exp_lp = rev(cumsum(rev(sortedIndexedData[,one_Z] * exp(by_hand_lp))))

  sum(sortedIndexedData$stat * (sortedIndexedData[one_Z] - cumsum_Z_exp_lp/cumsum_exp_lp)) - Z_penalties[one_Z]

})


cbind(byhand_grads, dlp_b(est_parms,
                          X = c("X1", "X2", "X3"),
                          cluster = "M",
                          t = t,
                          dij = stat,
                          D =  D,
                          data = ds),
      fit$u[4:53])

plot(sort(byhand_grads), sort(fit$u[4:53]))


plot(sort(byhand_grads) - sort(fit$u[4:53]))

plot(sort(fit$u[4:53]))
points(sort(byhand_grads), col = "red")


order(fit$u[4:53])
order(byhand_grads)


grad_order = ds %>%
  dplyr::distinct(M) %>%
  dplyr::pull(M)

plot(fit$u[4:53], byhand_grads[grad_order])

str(fit)


dlp_b(est_parms,
         X = c("X1", "X2", "X3"),
         cluster = "M",
         t = t,
         dij = stat,
         D =  D,
         data = ds)




max(addedLP$lp - by_hand_lp )

max(exp(by_hand_lp) - exp(addedLP$lp))

fit$loglik


addedCumsums <- calcRiskSets(addedLP, X, "Xr")



ll <- addedCumsums %>%
  dplyr::mutate(li =  {{ dij }} * (value - cumsum_Xr_A / cumsum_A)) %>%
  dplyr::summarise(ll = sum(li), .groups = "drop")

ll
















