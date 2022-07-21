
# trials with generic function and class-specific methods.

lp <- function(params, formula, data, theta){
  UseMethod("lp", data)
}


lp.ppl <- function(params, formula, data, theta) {

  # This sort var could be pulled out of the formula.
  ds_sorted <- sortAndIndex(data, sort_vars = stat_time)

  parsed_data <- lme4::lFormula(formula, data = ds_sorted)

  beta_index <- seq_len(ncol(parsed_data$X) - 1)

  beta <- params[beta_index]
  b <- params[-beta_index]

  # drop the intercept column from the X model.matrix
  risk_score <- parsed_data$X[, -1] %*% beta + Matrix::crossprod(parsed_data$reTrms$Zt, b)

  exp_risk_score <- exp(risk_score)

  # calculate risk sets
  # I do it this way to preserve the class and other slots. rev(risk_score) coerces to a numeric vector.
  # probably no point, as cumsum convert it anyway, and I need to remake the matrix with Matrix.
  rev_exp_risk_score <- exp_risk_score
  rev_exp_risk_score@x <- rev(exp_risk_score@x)

  at_risk <- Matrix::Matrix(rev(cumsum(rev_exp_risk_score)), ncol = 1)

  # would be good to peel stat out of the Surv() object in the parsed formula.
  # This is awkward because Surv() doesn't preserve the variable names, so I may
  # need to parse the formula myself, and including consideration of the "type"
  # attribute from the evaluated call.

  stat <- Matrix::Matrix(ds_sorted$stat, ncol = 1)

  # penalty

  parsed_data$reTrms$Lambdat@x <- theta[parsed_data$reTrms$Lind]

  penalty <- 0.5 * t(b) %*% Matrix::solve(parsed_data$reTrms$Lambdat) %*% b

  penalised_likelihood <- sum(stat * (risk_score - log(at_risk))) - penalty

  penalised_likelihood@x

}

lp.ppl_extra <- function(params, formula, data, theta) {

  # This sort var could be pulled out of the formula.
  ds_sorted <- sortAndIndex(data, sort_vars = stat_time)

  parsed_data <- lme4::lFormula(formula, data = ds_sorted)

  beta_index <- seq_len(ncol(parsed_data$X) - 1)

  beta <- params[beta_index]
  b <- params[-beta_index]

  # drop the intercept column from the X model.matrix
  risk_score <- parsed_data$X[, -1] %*% beta + Matrix::crossprod(parsed_data$reTrms$Zt, b)

  exp_risk_score <- exp(risk_score)

  # calculate risk sets
  # I do it this way to preserve the class and other slots. rev(risk_score) coerces to a numeric vector.
  # probably no point, as cumsum convert it anyway, and I need to remake the matrix with Matrix.
  rev_exp_risk_score <- exp_risk_score
  rev_exp_risk_score@x <- rev(exp_risk_score@x)

  at_risk <- Matrix::Matrix(rev(cumsum(rev_exp_risk_score)), ncol = 1)

  # would be good to peel stat out of the Surv() object in the parsed formula.
  # This is awkward because Surv() doesn't preserve the variable names, so I may
  # need to parse the formula myself, and including consideration of the "type"
  # attribute from the evaluated call.

  stat <- Matrix::Matrix(ds_sorted$stat, ncol = 1)

  # penalty

  parsed_data$reTrms$Lambdat@x <- theta[parsed_data$reTrms$Lind]

  D <- parsed_data$reTrms$Lambdat

  D_inverse <- Matrix::solve(D)

  penalty <- 0.5 * t(b) %*% D_inverse %*% b

  penalised_likelihood <- sum(stat * (risk_score - log(at_risk))) - penalty

  # terms 1 and 2 from ripatti and palmgren, which they drop.

  # determinant returns the logarithm by default. Matrix just calls the base
  # method by default. IDK if there speed advantage from sparseness.

  term1 <- -0.5 * log(det(D))

  # second term is more complicated. Need cumulative hazard (breslow), and Zt(Z)

  cumulative_hazard <- Matrix::Matrix(cumsum(stat/at_risk), ncol = 1)


  Zt <- parsed_data$reTrms$Zt

  Zt_ncol <- ncol(Zt)

  ZtZ_exp_risk_score <- vector(mode = "list", length = Zt_ncol)

  for (i in seq_len(Zt_ncol)) {

    ZtZ_exp_risk_score[[i]] <-  cumulative_hazard[i] * exp_risk_score[i] * tcrossprod(Zt[,i, drop = FALSE])

  }

  term2 <- -0.5 * determinant(Reduce("+", ZtZ_exp_risk_score) - D_inverse)$modulus

  term1 + term2 + penalised_likelihood@x

}

debugonce(lp.ppl_extra)

lp.ppl_extra(rnorm(n = 5), formula = survival::Surv(stat_time, stat)~X1 + X2 + X3 + (1 | M), data = ds, theta = 1)



make_ppl <- function(data){

  stopifnot(is.data.frame(data))

  class(data) <- c("ppl", oldClass(data))

  data

}

make_pplextra <- function(data){

  stopifnot(is.data.frame(data))

  class(data) <- c("ppl_extra", oldClass(data))

  data

}



# make objects with classes ppl, pplextra

my_formula <- survival::Surv(stat_time, stat)~X1 + X2 + X3 + (1 | M)

ds <- one_dataset(my_formula,
            dists = list(X1 = ~rnorm(n),
                         X2 = ~rnorm(k * nk),
                         X3 = ~rbinom(n, 1, 0.5),
                         M = ~rep(1:k, each = nk),
                         error = ~rexp(n, 10),
                         stat = ~sample(rep(c(0, 1), round(n * c(0.2, 0.8))), n)),
            dist_args = list(k = 2, nk = 2, n = 4),
            coefficients = c(1, 1, 1),
            random_effect_variance = list(M = 1)
)

ds.ppl <- make_ppl(ds)
ds.pplextra <- make_pplextra(ds)

fixed_formula <- lme4:::getFixedFormula(my_formula)

fit0 <- survival::coxph(fixed_formula, data = ds)

start_params <- c(coef(fit0), rep(0, 50))

lp(start_params, formula = my_formula, data = ds.ppl, theta = 1)
lp(start_params, formula = my_formula, data = ds.pplextra, theta = 1)

# can I use optim with a generic function?

ests <- optim(par = start_params,
              fn = lp,
              gr = NULL,
              formula = my_formula,
              data = ds.ppl,
              theta = 1,
              method = "Nelder-Mead",
              control = list(fnscale = -1))

ests2 <- optim(par = start_params,
              fn = lp,
              gr = NULL,
              formula = my_formula,
              data = ds.pplextra,
              theta = 1,
              method = "Nelder-Mead",
              control = list(fnscale = -1))

cbind(ests$par[1:3], ests2$par[1:3])




























