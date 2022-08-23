


# can optim work when the graident function is not 0 at the maximum?
# call a costly function (theta_ipl) once per iteration.
# https://stackoverflow.com/questions/50897490/r-optimization-pass-value-from-function-to-gradient-with-each-iteration


theta_ipl_est <- theta_ipl(test_thetas[i], formula = my_formula, parsed_data = make_ppl(parsed_data),
                           other_args = list(n_fixed = n_fixed,
                                             start_params = start_params,
                                             stat = stat,
                                             re_only = TRUE,
                                             reltol = 1e-13))
#
# theta_ipl_fn_gr <- function(theta, formula, parsed_data, other_args) {
#
#     quant1 <- NULL
#   prev_param <- NULL
#
#   update_theta_ipl <- function(param_vec) {
#     if (!identical(param_vec, prev_param)) {
#       quant1 <<- theta_ipl(theta = theta, formula = formula,
#                            parsed_data = parsed_data, other_args = other_args)
#       prev_param <<- param_vec
#     }
#   }
#
#   optfunc<-function(param_vec){
#     update_theta_ipl(param_vec)
#     loglikelihood<-sum(quant1)**2
#     return(loglikelihood)
#   }
#
#   optgr<-function(param_vec){
#     update_theta_ipl(param_vec)
#     mygrad<-sum(quant1)
#     return(mygrad)
#   }
#   return(list(fn = optfunc, gr = optgr))
#
# }


theta_ipl_fn_gr <- function(formula, parsed_data, other_args) {

  quant1 <- NULL
  prev_param <- NULL

  update_theta_ipl <- function(param_vec) {
    if (!identical(param_vec, prev_param)) {
      theta_ipl_est <<- theta_ipl(theta = param_vec, formula = formula,
                           parsed_data = parsed_data, other_args = other_args)
      prev_param <<- param_vec
    }
  }

  optfunc<-function(param_vec){
    update_theta_ipl(param_vec)
    ll<-c(theta_ipl_est)
    return(ll)
  }

  optgr<-function(param_vec){
    update_theta_ipl(param_vec)

    d_D_d_theta <- vector(mode = "list", length = length(parsed_data$reTrms$theta))

    for (j in seq_along(parsed_data$reTrms$theta)) {

      d_D_d_theta[[j]] <- parsed_data$reTrms$Lambdat

      d_D_d_theta[[j]]@x <- c(1, 0)[(parsed_data$reTrms$Lind != j) + 1]

    }

    b <- attr(theta_ipl_est, "b")
    Kbb <- attr(theta_ipl_est, "Kbb")
    Kbb_inv <- solve(Kbb)
    D <- attr(theta_ipl_est, "D")
    D_inv <- Matrix::solve(D)

    calc_grad <- function(dD){

      grad <- -0.5 * ( sum(diag( D_inv %*% dD ))
                       + sum(diag( Kbb_inv %*% D_inv %*% dD %*% D_inv ))
                       - t(b) %*% D_inv %*% dD %*% D_inv %*% b )

      grad@x

    }

    grads <- unlist(lapply(d_D_d_theta, calc_grad))

    return(grads)
  }

  return(list(fn = optfunc, gr = optgr))

}




# test this method

my_formula <- survival::Surv(stat_time, stat)~X1 + X2 + X3 + (1 | M1)

ds <- one_dataset(my_formula,
                  dists = list(X1 = ~rnorm(n),
                               X2 = ~rnorm(k * nk),
                               X3 = ~rbinom(n, 1, 0.5),
                               M1 = ~rep(1:k, each = nk),
                               M2 = ~rep(c("l", "r"), ceiling(n/2))[seq_len(n)],
                               error = ~rexp(n, 10),
                               stat = ~sample(rep(c(0, 1), round(n * c(0.2, 0.8))), n)),
                  dist_args = list(k = 50, nk = 4, n = 50*4),
                  coefficients = c(1, -0.7, 0.5),
                  random_effect_variance = c(M1 = 2)
)

coxfit <- coxme::coxme(survival::Surv(stat_time, stat)~X1 + X2 + X3 + (1 | M1), data = ds)

coxme::VarCorr(coxfit)
ds_sorted <- sortAndIndex(ds, sort_vars = stat_time)

parsed_data <- lme4::lFormula(my_formula, data = ds_sorted)

stat <- Matrix(unclass(parsed_data$fr[,1])[, "status"], ncol = 1)

theta_start <- get_start_theta(length(parsed_data$reTrms$flist))

n_fixed <- length(attr(terms(coxme:::formula1(my_formula)$fixed), "order"))

fixed_formula <- lme4:::getFixedFormula(my_formula)

fit0 <- survival::coxph(fixed_formula, data = ds)

start_params <- c(coef(fit0),
                  rnorm(n = length(parsed_data$reTrms$Lind),
                        mean = 0,
                        sd = sqrt(theta_start[parsed_data$reTrms$Lind])))

theta_max_coxme <- unlist(coxme::VarCorr(coxfit))

my_ests <- est_parameters(my_formula, data = ds,
                          control = control.list(factr = 1e3,
                                                 reltol = 1e-13,
                                                 ndeps = 1e-2))
# use my maximum.
theta_max <- my_ests$theta

theta_ipl_est <- theta_ipl(theta_max, formula = my_formula, parsed_data = make_ppl(parsed_data),
                           other_args = list(n_fixed = n_fixed,
                                             start_params = start_params,
                                             stat = stat,
                                             re_only = TRUE,
                                             reltol = 1e-13))



do.call(optim, c(list(par=paramvec,method="BFGS"),optfngr() ))

theta_ipl_fn_gr(formula = my_formula,
               parsed_data = make_ppl(parsed_data),
               other_args = list(n_fixed = n_fixed,
                                 start_params = start_params,
                                 stat = stat,
                                 re_only = TRUE,
                                 reltol = 1e-13))


do.call(optim, c(list(par = theta_start,
                    method = "L-BFGS-B",
                    control = list(fnscale = -1, factr = 1e3),
                    lower = 0.00001, upper = Inf,
                    hessian = TRUE),
                 theta_ipl_fn_gr(
                    formula = my_formula,
                    parsed_data = make_ppl(parsed_data),
                    other_args = list(n_fixed = n_fixed,
                                      start_params = start_params,
                                      stat = stat,
                                      re_only = TRUE,
                                      reltol = 1e-13))))


# try without the memoisation.

do.call(optim, c(list(par = theta_start,
                      method = "L-BFGS-B",
                      control = list(fnscale = -1, factr = 1e3),
                      lower = 0.00001, upper = Inf,
                      hessian = TRUE,
                      fn = theta_ipl,
                      gr = NULL,
                   formula = my_formula,
                   parsed_data = make_ppl(parsed_data),
                   other_args = list(n_fixed = n_fixed,
                                     start_params = start_params,
                                     stat = stat,
                                     re_only = TRUE,
                                     reltol = 1e-13))))


# timing:

microbenchmark::microbenchmark(

do.call(optim, c(list(par = theta_start,
                      method = "L-BFGS-B",
                      control = list(fnscale = -1, factr = 1e3),
                      lower = 0.00001, upper = Inf,
                      hessian = TRUE),
                 theta_ipl_fn_gr(
                   formula = my_formula,
                   parsed_data = make_ppl(parsed_data),
                   other_args = list(n_fixed = n_fixed,
                                     start_params = start_params,
                                     stat = stat,
                                     re_only = TRUE,
                                     reltol = 1e-13))))
,
do.call(optim, c(list(par = theta_start,
                      method = "L-BFGS-B",
                      control = list(fnscale = -1, factr = 1e3),
                      lower = 0.00001, upper = Inf,
                      hessian = TRUE,
                      fn = theta_ipl,
                      gr = NULL,
                      formula = my_formula,
                      parsed_data = make_ppl(parsed_data),
                      other_args = list(n_fixed = n_fixed,
                                        start_params = start_params,
                                        stat = stat,
                                        re_only = TRUE,
                                        reltol = 1e-13))))
, times = 10)


