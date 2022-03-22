


one_dataset <- function(control) {

  n = control$k*control$nk

  M = rep(1:control$k, each = control$nk)

  X1 = rnorm(n, 0, 1)
  X2 = rep(rnorm(control$k, 0, 1), each = control$nk)

  # cluster level binary treatment allocation
  X3 = rep(rep(c(1, 0), ceiling(control$k/2))[1:control$k], each = control$nk)

  X = cbind(X1, X2, X3)

  b = rep(rnorm(control$k, 0, sqrt(control$theta)), each = control$nk)

  error = rexp(n, 10)

  t = exp(-X%*%control$beta - b) * error

  stat =  sample(rep(c(0, 1), round(n*c(0.2, 0.8))), n)

  data.frame(X, t, stat, M)

}

sample_data <- one_dataset(control = list(k = 50, nk = 4, beta = c(1, -0.7, 0.5), theta = 1))

d1 <- sortAndIndex(sample_data, t)

calcLinearPredictor














































