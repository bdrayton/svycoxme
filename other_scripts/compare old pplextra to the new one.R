
library(Matrix)

my_formula <- survival::Surv(stat_time, stat)~X1 + X2 + X3 + (1 | M)
my_k = 20
my_nk = 4
my_theta = 1

ds <- one_dataset(my_formula,
                  dists = list(X1 = ~rnorm(n),
                               X2 = ~rnorm(k * nk),
                               X3 = ~rbinom(n, 1, 0.5),
                               M = ~rep(1:k, each = nk),
                               error = ~rexp(n, 10),
                               stat = ~sample(rep(c(0, 1), round(n * c(0.2, 0.8))), n)),
                  dist_args = list(k = my_k, nk = my_nk, n = my_k*my_nk),
                  coefficients = c(1, 1, 1),
                  random_effect_variance = list(M = 1)
)

ds_sorted <- sortAndIndex(ds, sort_vars = stat_time)

parsed_data <- lme4::lFormula(my_formula, data = ds_sorted)

Zt <- parsed_data$reTrms$Zt

Matrix::tcrossprod(Zt)


pplextra <- function(d, X, Z, B, b, D) {

  D_inverse = solve(D)

  lp = X%*%B + Z%*%b

  elp = exp(X%*%B + Z%*%b)

  risksets = calcRiskSetsMatrix(elp)

  hazard = cumsum(d / risksets)

  ZZ_list = NULL

  for (i in 1:nrow(Z)) {

    zt = Z[i,]

    ZZ_list[[i]] = zt %*% t(zt)

  }

  bparts <- lapply(seq_along(ZZ_list), function(i){
    hazard[i] * elp[i] * ZZ_list[[i]]
  })

  as.numeric(
    -0.5 * determinant(D)$modulus +
      -0.5 * determinant(Reduce(`+`, bparts) - D_inverse )$modulus +
      sum(d * (lp - log(risksets))) - t(b) %*% D_inverse %*% b)

}






# compare old ppl and ppl extra with new versions

ds.ppl <- make_ppl(ds)
ds.pplextra <- make_pplextra(ds)

fixed_formula <- lme4:::getFixedFormula(my_formula)

fit0 <- survival::coxph(fixed_formula, data = ds)

start_params <- c(coef(fit0), rep(0, my_k))

sortedData <- add_Z(ds_sorted, cluster = "M")

xColumns <- paste0("X", 1:3)
zColumns <- paste0("Z", 1:my_k)

B = coef(fit0)
b = rep(0, my_k)

D_theta = my_theta * diag(my_k)

##

ppl(d = sortedData$stat, X = as.matrix(sortedData[, xColumns]),
    Z = as.matrix(sortedData[, zColumns]), B = B, b = b, D = D_theta)

lp(start_params, formula = my_formula, data = ds.ppl, theta = my_theta)


##
pplextra(d = sortedData$stat, X = as.matrix(sortedData[, xColumns]),
         Z = as.matrix(sortedData[, zColumns]), B = B, b = b, D = D_theta)

lp(start_params, formula = my_formula, data = ds.pplextra, theta = my_theta)

##


microbenchmark::microbenchmark(
  ppl(d = sortedData$stat, X = as.matrix(sortedData[, xColumns]),
      Z = as.matrix(sortedData[, zColumns]), B = B, b = b, D = D_theta),
  pplextra(d = sortedData$stat, X = as.matrix(sortedData[, xColumns]),
           Z = as.matrix(sortedData[, zColumns]), B = B, b = b, D = D_theta),
  lp(start_params, formula = my_formula, data = ds.ppl, theta = my_theta),
  lp(start_params, formula = my_formula, data = ds.pplextra, theta = my_theta),
  lp.ppl(start_params, formula = my_formula, data = ds.ppl, theta = my_theta),
  lp.ppl_extra(start_params, formula = my_formula, data = ds.pplextra, theta = my_theta)
)














