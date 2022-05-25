

my_beta = c(1, -0.7, 0.5)
my_theta = 1
my_k = 10
my_nk = 10

my_X = c("X1", "X2", "X3")

ds <- one_dataset(list(k = my_k, nk = my_nk, beta = my_beta, theta = my_theta))

fit <- coxme::coxme(survival::Surv(stat_time, stat) ~ X1 + X2 + X3 + (1|M), data = ds)

coxme_est_theta <- coxme::VarCorr(fit)$M

coxme_est_parms <- c(coxme::fixef(fit), coxme::ranef(fit)$M)

parm_names <- c(my_X, paste0("Z", seq_len(my_k)))

names(coxme_est_parms) <- parm_names

my_loglik <- lp(parms = coxme_est_parms,
   X = my_X,
   cluster = "M",
   t = stat_time,
   dij = stat,
   theta = coxme_est_theta,
   data = ds)

# what is the difference between our log likelihoods?
fit$loglik["Penalized"] - my_loglik

# and penalties?
fit$penalty - attr(my_loglik, "penalty")


my_u <- lp_grd(parms = coxme_est_parms,
               X = c("X1", "X2", "X3"),
               cluster = "M",
               t = stat_time,
               dij = stat,
               theta = coxme_est_theta,
               data = ds)

# the gradient in coxme have a mixed up order


# find the parameters that gradients refer to. lp_grd gives pretty much the
# same gradients as coxme, but i know the parametres. this code gives the parameter
# order for the coxme gradients, which I'm hoping is inhereted by the coxme hessian.

ordered_coxme_u <-
   tibble::tibble(coxme_u = fit$u,
                  coxme_u_order = seq_along(coxme_u)) %>%
   dplyr::arrange(coxme_u)


ordered_my_u <-
   tibble::tibble(parm = parm_names,
                  my_u = my_u) %>%
   dplyr::arrange(my_u)

# it looks like the only reordering is to put the fixed effects after the random effects.
dplyr::bind_cols(
   ordered_coxme_u,
   ordered_my_u) %>%
   dplyr::arrange( coxme_u_order) %>%
   dplyr::pull(parm)

my_hessian <- svycoxme::ppl_hessian(parms = coxme_est_parms,
                                    X = c("X1", "X2", "X3"),
                                    cluster = "M",
                                    t = t,
                                    dij = stat,
                                    theta = coxme_est_theta,
                                    data = ds)

# reorder my hessian to be the same as coxme (at least what I think coxme is)

reordered_names <- colnames(my_hessian)[c(seq_len(my_k) + 3, 1:3)]

reord_hessian <- my_hessian[reordered_names, reordered_names]

<<<<<<< Updated upstream
# decompose my hessian
=======
library(bdsmatrix)

>>>>>>> Stashed changes
gchol_reord_hessian <- gchol(reord_hessian)

gchol_my_hessian <- gchol(my_hessian)

# I can get the original matrix back, but not perfectly. introduces errors of magnitude less than e-14

L <- as.matrix(gchol_reord_hessian)

D <- diag(gchol_reord_hessian)

# maximum error
max(L %*% diag(D) %*% t(L) - reord_hessian)

# Is fit_coxme$hmat the same as my_hessian after the appropriate transformations and reordering?

L <- as.matrix(fit$hmat)

D <- diag(fit$hmat)

back_trans_hmat <- L %*% diag(D) %*% t(L)

back_trans_hmat - (-1 * reord_hessian)

<<<<<<< Updated upstream
# for the fixed effects, estimates are of the same (to almost machine error).
=======
# for the fixed effects, estimates are the same.
>>>>>>> Stashed changes
vcov(fit)
solve(-my_hessian)[1:3, 1:3]

# looking at the variances and ignoring non-diagonal terms, they are the same.
cbind(
   myvar = diag(solve(-reord_hessian)),
   coxmevar = diag(fit$variance))

# solve(hmat) == coxme$variance
# More precisely, the inverse of coxme$hmat differs from coxme$variance by less
# than machine error.
hmat_inv <- solve(fit$hmat, full = TRUE)
all.equal(fit$variance, hmat_inv)




