#' @export

make_parts <- function(x, data, ...){

    UseMethod("make_parts", x)

}


#' @export

make_parts.coxme <- function(coxme.object, data){

  # extract the formula from the object call.
  # eval in case a symbol was passed in e.g.
  ## my_form = Surv(time, status) ~ X + (1 + M)
  ## coxme(my_form, data = my_data)

  form <- eval(coxme.object$call[[2]])

  # response is actually a (masked) matrix with time and status.
  # but, see trello card.

  # need to use getFixedFormula, otherwise model.frame may complain about random effects terms.
  response <- model.response(model.frame(lme4:::getFixedFormula(form), data))

  time <- response[,"time"]

  stat <- response[,"status"]

  weights <- weights(coxme.object)

  if(is.null(weights)){
    # weights <- rep(1, length(stat)) # this is wrong. if the weights are 1, the ui2s are on the wrong scale.
    # better to throw an error than give the wrong answer.
    stop("I need weights. Add weights = your_weights to your model call.")
  }

  # reorder things
  time_order <- order(time)

  time <- time[time_order]
  stat <- stat[time_order]

  weights <- weights[time_order]

  ds_sorted <- data[time_order, ]

  parsed_data <- lme4::lFormula(form, data = ds_sorted)

  # drop the intercept column from the X model.matrix
  # I use X a few times, so extracting. but could I use this each time?
  # does it make a difference anyone would care about in terms of memory or speed?
  X <- parsed_data$X[ , -1, drop = FALSE]
  Zt <- parsed_data$reTrms$Zt
  Z <- t(Zt)


  # I think these need to be p*1 and k*1 matrices
  beta <- Matrix(coxme::fixef(coxme.object), ncol = 1)
  # need to unlist these and line them up with Zt
  # they should be lined up by default, but this will need some testing
  b <- Matrix(unlist(coxme::ranef(coxme.object)), ncol = 1)

  risk_score <- X %*% beta + Matrix::crossprod(Zt, b)

  # weighted
  exp_risk_score <- weights * exp(risk_score)

  # this is S0_hat in binder, a n * 1 matrix.
  at_risk <- fast_risk_sets(exp_risk_score)

  # this is S1_hat, a n * p matrix
  exp_risk_score_X <- exp_risk_score * X
  n <- nrow(exp_risk_score)
  exp_risk_score_Z <- Matrix(as.numeric(exp_risk_score) * as.numeric(Z), nrow = n)

  # also an n * p matrix
  at_risk_X <- fast_risk_sets(exp_risk_score_X)

  # the equivalent for Z is an n * q matrix
  # if I could not transpose Zt that would be nice (faster).
  at_risk_Z <- fast_risk_sets(exp_risk_score_Z)

  # here we gain another dimension. n * p * p
  # we want to sum over i, so splitting over i into a list or array seems sensible.
  # however, I'll then need to recombine using Reduce, which is slow. All these matricies are dense,
  # so an array might work.
  # solution: change the dimensions of XtX from p * p to 1 * (p^2).
  # the matrix will then be n * p^2, and any operations using it will need to know this.

  # this is a p * n matrix. the ith column is t(X[i, ]) %*% X[i, ]
  XtX_i <- apply(X, 1, tcrossprod)

  # I need XtZ and ZtZ matricies too.
  # I think the column indexing here... apply will pass tcrossprod a vector, not a n * 1 matrix, so it's the same results as
  # apply(t(parsed_data$reTrms$Zt), 1, tcrossprod), but I don't need to transpose.

  ZtZ_i <- apply(Zt, 2, tcrossprod)

  exp_risk_score_XtX <- Matrix(as.numeric(t(XtX_i)) * as.numeric(exp_risk_score), nrow = n)

  exp_risk_score_ZtZ <- Matrix(as.numeric(t(ZtZ_i)) * as.numeric(exp_risk_score), nrow = n)

  XtZ_i <- mapply(tcrossprod, split(X, row(X)), split(Zt, col(Zt)))

  exp_risk_score_XtZ <- Matrix(as.numeric(t(XtZ_i)) * as.numeric(exp_risk_score), nrow = n)

  # I think I can give these to fast_risk_sets.
  at_risk_XtX <- fast_risk_sets(exp_risk_score_XtX)

  at_risk_ZtZ <- fast_risk_sets(exp_risk_score_ZtZ)

  at_risk_XtZ <- fast_risk_sets(exp_risk_score_XtZ)

  # I also nee the penalty, divided into n parts. I'm going to ignore it for now.
  # will need D(theta), divided up appropritely for the calculation of D_i.
  # will need b * D(theta), divided up appropriately for the calculation of U_i
  # Will need to weight it appropriately.

  # this is the penalty matrix.
  # assume that each term has one theta. Must be shared frailty of some sort.

  # theta <- unlist(coxme::VarCorr(coxme.object))
  # parsed_data$reTrms$Lambdat@x <- theta[parsed_data$reTrms$Lind]
  # D <- parsed_data$reTrms$Lambdat

  # should I add in the Z equivalents, or can I treat them as X components?
  # add them in. there are cross product terms too.

  list(stat = stat,
       time = time,
       weights = weights,
       S0 = at_risk,
       S1_X = at_risk_X,
       S1_Z = at_risk_Z,
       exp_risk_score = exp(risk_score),
       weighted_exp_risk_score = exp_risk_score,
       S2_XtX = at_risk_XtX,
       S2_ZtZ = at_risk_ZtZ,
       S2_XtZ = at_risk_XtZ,
       X = X,
       Z = Z)

}


#' @export

make_parts.coxph <- function(coxph.object, data){

  form <- eval(coxph.object$call[[2]])

  model_frame <- model.frame(form, data)

  response <- model.response(model_frame)

  time <- response[,"time"]

  stat <- response[,"status"]

  weights <- weights(coxph.object)

  if(is.null(weights)){
    weights <- rep(1, length(stat))
  }

  # reorder things
  time_order <- order(time)

  time <- time[time_order]
  stat <- stat[time_order]

  weights <- weights[time_order]

  X <- model.matrix(coxph.object)[time_order, , drop = FALSE]

  # I think these need to be p*1 and k*1 matrices
  beta <- Matrix(coef(coxph.object), ncol = 1)

  risk_score <- X %*% beta

  # weighted
  exp_risk_score <- weights * exp(risk_score)

  # this is S0_hat in binder, a n * 1 matrix.
  at_risk <- fast_risk_sets(exp_risk_score)

  # this is S1_hat, a n * p matrix
  exp_risk_score_X <- exp_risk_score * X

  # also an n * p matrix
  at_risk_X <- fast_risk_sets(exp_risk_score_X)


  # here we gain another dimension. n * p * p
  # we want to sum over i, so splitting over i into a list or array seems sensible.
  # however, I'll then need to recombine using Reduce, which is slow. All these matricies are dense,
  # so an array might work.
  # solution: change the dimensions of XtX from p * p to 1 * (p^2).
  # the matrix will then be n * p^2, and any operations using it will need to know this.


  # this is a p * n matrix. the ith column is t(X[i, ]) %*% X[i, ]
  XtX_i <- matrix(apply(X, 1, tcrossprod), ncol = nrow(X))




  exp_risk_score_XtX <- t(XtX_i) * exp_risk_score

  # I think I can give this to fast_risk_sets.
  at_risk_XtX <- fast_risk_sets(exp_risk_score_XtX)

  # should I add in the Z equivalents, or can I treat them as X components?
  # add them in.

  list(stat = stat,
       time = time,
       weights = weights,
       S0 = at_risk,
       S1 = at_risk_X,
       exp_risk_score = exp(risk_score),
       weighted_exp_risk_score = exp_risk_score,
       S2 = at_risk_XtX,
       X = X)

}

#' @export

make_Di <- function(parts){

  with(parts, {

    weights * stat * (t(matrix(apply(S1/S0, 1, tcrossprod), nrow = ncol(S2))) - S2/S0)
  })

}

#'
#' @export

calc_ui <- function(parts){

  with(parts, {

    weights * stat * (X - S1/S0)

  })

}


#'
#' @export

calc_ui2 <- function(parts){

  n <- length(parts$stat)
  N_hat <- sum(parts$weights)
  p <- ncol(parts$X)

  # first term is the same

  parts <- within(parts, {
    S1_S0 = S1/S0
  })

  term1 <- with(parts, {

    weights * stat * (X - S1_S0)

  })

  # divide each exp(beta*X) by the series of risk sets.
  # this needs to be an array with the third dimension = p (number of covariates.)
  p1 <- with(parts, {

    array(tcrossprod(exp_risk_score, 1/S0), dim = c(n, n, p))

  })

  # replicate X n times to make a n*n matrix from n*1 X.

  # to update this for multiple X, need to replicate each X

  # minus S0/S1, which is also n*1 from each row
  # Matrix(as.numeric(X), ncol = n, nrow = n) - Matrix(as.numeric(S1/S0), ncol = n, nrow = n, byrow = TRUE)

  X_mat <- apply(parts$X, 2, function(x){

    matrix(x, nrow = n, ncol = n)

  })

  dim(X_mat) <- c(n, n, p)

  S1_S0_mat  <- apply(parts$S1_S0, 2, function(x){

    matrix(x, nrow = n, ncol = n, byrow = TRUE)

  })

  dim(S1_S0_mat) <- c(n, n, p)

  p2 <- X_mat - S1_S0_mat

  # yi is the indicator I(ti >= tj). Each ti needs to be compared to all tj
  # because ti is sorted, this should be an upper triangle matrix (if no ties)

  # to update to array, this matrix needs to be stacked p times.
  # array will need this to be a matrix, not Matrix


  yi_temp <- with(parts, {

    matrix(rep(time, each = n) <= rep(time, n), ncol = n, nrow = n, byrow = TRUE)

  })

  yi <- array(yi_temp, dim = c(n, n, p))


  # matrix with n reps of stat, as rows.
  # to update to array, this matrix needs to be stacked p times.
  # array will need this to be a matrix, not Matrix

  dj_temp <- with(parts, {

    matrix(stat, nrow = n, ncol = n, byrow = TRUE)

  })

  dj <- array(dj_temp, dim = c(n, n, p))


  # matrix with n reps of weights, as rows.
  # to update to array, this matrix needs to be stacked p times.
  # array will need this to be a matrix, not Matrix

  wj_temp <- with(parts, {

    matrix(weights, nrow = n, ncol = n, byrow = TRUE)

  })

  wj <- array(wj_temp, dim = c(n, n, p))

  # the second term is the row sums after taking the product of these matricies

  # sum across j
  # division by N_hat is a lin thing. Not in the binder formula.
  # what does rowSums return with an array?
  # I think I need to split dimension 3 (layer?) and do row sums over those matricies
  # what I want returned is an n * p matrix

  term2 <- apply(-dj * wj * yi * p1 * p2, 3, rowSums)/N_hat

  r <- term1 + term2

  ## if you need to examine the parts
  attr(r, "terms") <- list(term1 = term1, term2 = term2, dj = dj, wj = wj, yi = yi, p1 = p1, p2 = p2, N_hat = N_hat)

  r

}

#'
#' @export
#'
#' should only work for p = 1, and be the same as the new version.
#' update: it is the same.
#' update: it was, but I change a denominator from n to N_hat

old_calc_ui2 <- function(parts){

  n <- length(parts$stat)

  # first term is the same

  term1 <- with(parts, {

    weights * stat * (X - S1/S0)

  })

  # divide each exp(beta*X) by the series of risk sets.
  p1 <- with(parts, {

    tcrossprod(exp_risk_score, 1/S0)

  })

  # replicate X n times to make a n*n matrix from n*1 X.

  # to update this for multiple X, need to replicate each X

  # minus S0/S1, which is also n*1 from each row
  p2 <- with(parts, {

    Matrix(as.numeric(X), ncol = n, nrow = n) - Matrix(as.numeric(S1/S0), ncol = n, nrow = n, byrow = TRUE)

  })

  # yi is the indicator I(ti >= tj). Each ti needs to be compared to all tj
  # because ti is sorted, this should be an upper triangle matrix (if no ties)

  yi <- with(parts, {

    Matrix(rep(time, each = n) <= rep(time, n), ncol = n, nrow = n, byrow = TRUE)

  })

  # matrix with n reps of stat, as rows.
  dj <- with(parts, {

    Matrix(stat, nrow = n, ncol = n, byrow = TRUE)

  })

  # matrix with n reps of weights, as rows.
  wj <- with(parts, {

    Matrix(weights, nrow = n, ncol = n, byrow = TRUE)

  })

  # the second term is the row sums after taking the product of these matricies

  # sum across j
  # division by n is a lin thing. Not in the binder formula.
  term2 <- rowSums(-dj * wj * yi * p1 * p2)/n

  r <- term1 + term2

  attr(r, "terms") <- list(term1, dj, wj, yi, p1, p2)

  r

}







