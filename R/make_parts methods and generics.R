#' @export

make_parts <- function(x, data, ...){

    UseMethod("make_parts", x)

}


#' @export

make_parts.coxme <- function(coxme.object, data, weights){

  # extract the formula from the object call.
  # eval in case a symbol was passed in e.g.
  ## my_form = Surv(time, status) ~ X + (1 + M)
  ## coxme(my_form, data = my_data)

  form <- eval(coxme.object$call[[2]])

  # response is actually a (masked) matrix with time and status.
  # but, see trello card.

  # need to use getFixedFormula, otherwise model.frame may complain about random effects terms.
  response <- model.response(model.frame(lme4:::getFixedFormula(form), data))

  # this is fine for some Surv objects, but depends on attr(Surv_object, "type").
  # If type is "counting", then there are three columns, start, stop, and status.
  time <- response[,"time"]

  stat <- response[,"status"]

  # I need the actual weights, but to get the correct point estimates for coxme, I need to use
  # rescaled weights, so this will be wrong. Weights are now passed in by the user.

  #
  # weights <- weights(coxme.object)
  #
  # if(is.null(weights)){
  #   # weights <- rep(1, length(stat)) # this is wrong. if the weights are 1, the ui2s are on the wrong scale.
  #   # better to throw an error than give the wrong answer.
  #   # stop("I need weights. Add weights = your_weights to your model call.")
  #
  #   # i rescale the weights now, and if they are 1 then coxme drops them, but i need the
  #   # rest of this to run when this happens. plus I don't use ui2 anymore.
  #   weights <- rep(1, length(stat))
  #
  # }

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

  # check here that the results has the same X Z ordering as this line from calculating Di.
  # all_rows <- t(mapply(Matrix::tcrossprod, split(parts$S1_X, row(parts$S1_X)), split(parts$S1_Z, row(parts$S1_Z)))) |> Matrix()
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

  theta <- unlist(coxme::VarCorr(coxme.object))
  parsed_data$reTrms$Lambdat@x <- theta[parsed_data$reTrms$Lind]
  D <- parsed_data$reTrms$Lambdat

  penalty <- t(b) %*% D

  # divided by n and repeated
  ui_penalty <- (penalty/n)[rep(1, n),]

  # should I add in the Z equivalents, or can I treat them as X components?
  # add them in. there are cross product terms too.

  r <- list( stat = Matrix(stat, ncol = 1),
             time = Matrix(time, ncol = 1),
             weights = Matrix(weights, ncol = 1),
             S0 = at_risk,
             S1_X = at_risk_X,
             S1_Z = at_risk_Z,
             exp_risk_score = exp(risk_score),
             weighted_exp_risk_score = exp_risk_score,
             S2_XtX = at_risk_XtX,
             S2_ZtZ = at_risk_ZtZ,
             S2_XtZ = at_risk_XtZ,
             X = X,
             Z = Z,
             ui_penalty = ui_penalty,
             di_penalty = D)

  class(r) <- c("coxme_parts", class(r))

  return(r)

}


#' @export

make_parts.coxph <- function(coxph.object, data, weights){

  form <- eval(coxph.object$call[[2]])

  model_frame <- model.frame(form, data)

  response <- model.response(model_frame)


  response_type = attr(response, 'type')

  if(response_type == "right"){

    time_stop = response[,"time"]
    time_start = rep(0, length(time_stop))
    stat = response[,'status']

    # time <- response[,"time"]
    # stat <- response[,"status"]

  } else if(response_type == "counting") {

    time_start = response[,"start"]
    time_stop = response[,"stop"]
    stat = response[,"status"]

  } else {

    stop("response type is not supported")

  }

  ## weights now passed in by user.
  # weights <- weights(coxph.object)
  #
  # if(is.null(weights)){
  #   weights <- rep(1, length(stat))
  # }

  # reorder things

  # if(response_type == "right"){
  #   time_order <- order(time)
  #
  #   time <- time[time_order]
  #   stat <- stat[time_order]
  # }
  # else if(response_type == "counting"){

    time_order = order(time_stop, time_start)

    time_start = time_start[time_order]
    time_stop = time_stop[time_order]

  # }

  weights <- weights[time_order]

  X <- model.matrix(coxph.object)[time_order, , drop = FALSE]

  # I think these need to be p*1 and k*1 matrices
  beta <- Matrix(coef(coxph.object), ncol = 1)

  risk_score <- X %*% beta

  # weighted
  exp_risk_score <- weights * exp(risk_score)

  # this is S0_hat in binder, a n * 1 matrix.
  # this is the faster way (i think, not tested) to do it, if you don't have counting time.
  # at_risk <- fast_risk_sets(exp_risk_score)

  # for counting time, need an n * n matrix of in_risk_set indicators
  # for each start time, check if an observation's start time is greater or equal
  # for each end time, check if an observation's end time is less or equal
  # in_risk_set if both are true
  # this way also works for 'right' data (if start times are set to 0).

  n = nrow(response)

  start_test <- time_start <= rep(time_start, each = n)
  stop_test <- time_stop >= rep(time_stop, each = n)

  in_risk_set_matrix <- matrix(start_test & stop_test, nrow = n)

  at_risk <- colSums(in_risk_set_matrix * exp_risk_score[,rep(1,n)])

  # this is S1_hat, a n * p matrix
  exp_risk_score_X <- exp_risk_score * X

  # also an n * p matrix
  # at_risk_X <- fast_risk_sets(exp_risk_score_X)
  at_risk_X <- apply(exp_risk_score_X, 2, function(X_j){

    colSums(in_risk_set_matrix * X_j)

  })

  # here we gain another dimension. n * p * p
  # we want to sum over i, so splitting over i into a list or array seems sensible.
  # however, I'll then need to recombine using Reduce, which is slow. All these matricies are dense,
  # so an array might work.
  # solution: change the dimensions of XtX from p * p to 1 * (p^2).
  # the matrix will then be n * p^2, and any operations using it will need to know this.


  # this is a p * n matrix. the ith column is t(X[i, ]) %*% X[i, ]
  # I think this stays the same for counting time.
  XtX_i <- matrix(apply(X, 1, tcrossprod), ncol = nrow(X))

  exp_risk_score_XtX <- t(XtX_i) * exp_risk_score

  # I think I can give this to fast_risk_sets.
  # at_risk_XtX <- fast_risk_sets(exp_risk_score_XtX)
  # this changes with counting time to:

  at_risk_XtX <- apply(exp_risk_score_XtX, 2, function(X_j){

    colSums(in_risk_set_matrix * X_j)

  })

  r <- list(stat = Matrix(stat, ncol = 1, sparse = FALSE),
       time_start = Matrix(time_start, ncol = 1),
       time_stop = Matrix(time_stop, ncol = 1),
       weights = Matrix(weights, ncol = 1),
       S0 = at_risk,
       S1 = at_risk_X,
       exp_risk_score = exp(risk_score),
       weighted_exp_risk_score = exp_risk_score,
       S2 = at_risk_XtX,
       X = X)

  # maybe coxph_parts would be a better class...
  class(r) <- c("coxph_parts", class(r))

  r

}


#' @export

calc_Di <- function(x, ...){

  UseMethod("calc_Di", x)

}


#' @export

calc_Di.coxph_parts <- function(parts){

  D_beta_beta <- with(parts, {

    S1S0_squared <- apply(S1/S0, 1, tcrossprod)

    if(is.matrix(S1S0_squared)) {
      t_S1S0_squared <- t(S1S0_squared)
    } else {
      t_S1S0_squared <- Matrix(S1S0_squared, ncol = 1)
    }

    weights * stat * (t_S1S0_squared - S2/S0)

  })

  # assemble the Hessian

  n <- nrow(D_beta_beta)
  N_hat <- sum(parts$weights)

  n_fixef <- sqrt(ncol(D_beta_beta))

  # with division by n.
  # top_left <- matrix(colSums(D_beta_beta)/n, nrow = n_fixef)

  top_left <- matrix(colSums(D_beta_beta)/N_hat, nrow = n_fixef)

  top_left

}

#' @export

calc_Di.coxme_parts <- function(parts){

  D_beta_beta <- with(parts, {

    weights * stat * (t(apply(S1_X/S0, 1, tcrossprod)) - S2_XtX/S0)

  })

  D_b_b <- with(parts, {

    weights * stat * (t(apply(S1_Z/S0, 1, tcrossprod)) - as(S2_ZtZ, "unpackedMatrix")/S0)

  })

  D_beta_b <- with(parts, {

    # numerator for first term
    S1_XtS1_Z <- Matrix(t(mapply(Matrix::tcrossprod, split(S1_X, row(S1_X)), split(S1_Z, row(S1_Z)))))

    # denominator for first term
    ncols <- ncol(X) * ncol(Z)
    S0tS0 <- Matrix(rep(apply(S0, 1, tcrossprod), ncols), nrow = nrow(S2_XtZ))

    stat * weights * (S1_XtS1_Z/S0tS0 - as(S2_XtZ, "unpackedMatrix") / S0)


  })

  # assemble the Hessian

  # used to divide by n, now use Nhat. see eqn 2.5 in Lin 2000
  n <- nrow(D_beta_beta)
  Nhat <- sum(parts$weights)

  n_fixef <- sqrt(ncol(D_beta_beta))

  n_ranef <- sqrt(ncol(D_b_b))

  top_left <- matrix(colSums(D_beta_beta), nrow = n_fixef)

  # minus penalty
  bottom_right <- matrix(colSums(D_b_b), nrow = n_ranef) - parts$di_penalty

  top_right <- matrix(colSums(D_beta_b), nrow = n_fixef)

  bottom_left <- t(top_right)

  ret <- rbind(
    cbind(top_left, top_right),
    cbind(bottom_left, bottom_right)
  )

  # division by Nhat as in Lin.
  ret/Nhat

}


#'
#' @export

calc_ui <- function(x, ...){

  UseMethod("calc_ui", x)

}


#'
#' @export

calc_ui.coxph_parts <- function(parts){

  n = sum(parts$weights)

  lin_term2 <- parts$X

  lin_term2[] <- NA

  for(i in seq_len(nrow(lin_term2))) {

    lin_term2[i, ] <- with(parts,{
      irep <- rep(i, length(parts$stat))
      # no division by n
      # colSums((stat * weights * (time[irep]>=time) * exp_risk_score[irep] * (1/S0)) * (X[irep, ] - S1/S0))
      Yi_at_tj = (time_stop[irep]>=time_stop & time_start <= time_start)
      colSums((stat * weights * Yi_at_tj * exp_risk_score[irep] * (1/S0)) * (X[irep, ] - S1/S0))
    })
  }

  lin_term1 <- with(parts,{
    stat * (X - S1/S0)
  })

  lin_score <- lin_term1 - lin_term2

  lin_score

}


#'
#' @export

calc_ui.coxme_parts <- function(parts){

  n = sum(parts$weights)

  lin_X_term2 <- parts$X
  lin_Z_term2 <- parts$Z

  lin_X_term2[] <- NA
  lin_Z_term2[] <- NA

  for(i in seq_len(nrow(lin_X_term2))) {

    irep <- rep(i, length(parts$stat))

    lin_X_term2[i, ] <- with(parts,{
      # no division by n
      colSums((stat * weights * (time[irep]>=time) * exp_risk_score[irep] * (1/S0)) * (X[irep, ] - S1_X/S0))
    })

    lin_Z_term2[i, ] <- with(parts,{
      colSums((stat * weights * (time[irep]>=time) * exp_risk_score[irep] * (1/S0)) * (Z[irep, ] - S1_Z/S0))
    })

  }

  lin_term2 <- cbind(lin_X_term2, lin_Z_term2)

  lin_term1_beta <- with(parts, {

    # weights * stat * (X - S1_X/S0)
    stat * (X - S1_X/S0)

  })

  lin_term1_b <- with(parts, {

    stat * ((Z - S1_Z/S0) - ui_penalty)

  })

  lin_term1 <- cbind(lin_term1_beta, lin_term1_b)

  lin_score <- lin_term1 - lin_term2

  lin_score

}

#'
#' # this is the binder ui and should agree with Claudia's ui
#' #'
#' #' @export
#'
#' calc_ui_binder.coxph(parts){
#'
#'   # this is the estimated ui, arrived at by subbing in estimates to formula 3.7 in binder.
#'
#'   # there are three terms
#'
#'   term_1 <- with(parts,{
#'
#'     stat * (X - S1/S0)
#'
#'   })
#'
#'
#'
#'
#'
#'   term_2 <- with(parts, {
#'
#'    stat * weight * (X * exp_risk_score / S0)
#'
#'
#'   })
#'
#'
#'
#'
#'
#' }
#'
#' # need one x observation.
#' # need all the observations from that x onwards.
#' # need weights, stat, exp risk score
#'
#' make_term2 <- function() {
#'
#'
#'
#' }
#'



#
# #'
# # @export
#'
#' calc_ui2 <- function(parts){
#'
#'   n <- length(parts$stat)
#'   N_hat <- sum(parts$weights)
#'   p <- ncol(parts$X)
#'
#'   # first term is the same
#'
#'   parts <- within(parts, {
#'     S1_S0 = S1/S0
#'   })
#'
#'   term1 <- with(parts, {
#'
#'     weights * stat * (X - S1_S0)
#'
#'   })
#'
#'   # divide each exp(beta*X) by the series of risk sets.
#'   # this needs to be an array with the third dimension = p (number of covariates.)
#'   p1 <- with(parts, {
#'
#'     array(tcrossprod(exp_risk_score, 1/S0), dim = c(n, n, p))
#'
#'   })
#'
#'   # replicate X n times to make a n*n matrix from n*1 X.
#'
#'   # to update this for multiple X, need to replicate each X
#'
#'   # minus S0/S1, which is also n*1 from each row
#'   # Matrix(as.numeric(X), ncol = n, nrow = n) - Matrix(as.numeric(S1/S0), ncol = n, nrow = n, byrow = TRUE)
#'
#'   X_mat <- apply(parts$X, 2, function(x){
#'
#'     matrix(x, nrow = n, ncol = n)
#'
#'   })
#'
#'   dim(X_mat) <- c(n, n, p)
#'
#'   S1_S0_mat  <- apply(parts$S1_S0, 2, function(x){
#'
#'     matrix(x, nrow = n, ncol = n, byrow = TRUE)
#'
#'   })
#'
#'   dim(S1_S0_mat) <- c(n, n, p)
#'
#'   p2 <- X_mat - S1_S0_mat
#'
#'   # yi is the indicator I(ti >= tj). Each ti needs to be compared to all tj
#'   # because ti is sorted, this should be an upper triangle matrix (if no ties)
#'
#'   # to update to array, this matrix needs to be stacked p times.
#'   # array will need this to be a matrix, not Matrix
#'
#'
#'   yi_temp <- with(parts, {
#'
#'     matrix(rep(time, each = n) <= rep(time, n), ncol = n, nrow = n, byrow = TRUE)
#'
#'   })
#'
#'   yi <- array(yi_temp, dim = c(n, n, p))
#'
#'
#'   # matrix with n reps of stat, as rows.
#'   # to update to array, this matrix needs to be stacked p times.
#'   # array will need this to be a matrix, not Matrix
#'
#'   dj_temp <- with(parts, {
#'
#'     matrix(stat, nrow = n, ncol = n, byrow = TRUE)
#'
#'   })
#'
#'   dj <- array(dj_temp, dim = c(n, n, p))
#'
#'
#'   # matrix with n reps of weights, as rows.
#'   # to update to array, this matrix needs to be stacked p times.
#'   # array will need this to be a matrix, not Matrix
#'
#'   wj_temp <- with(parts, {
#'
#'     matrix(weights, nrow = n, ncol = n, byrow = TRUE)
#'
#'   })
#'
#'   wj <- array(wj_temp, dim = c(n, n, p))
#'
#'   # the second term is the row sums after taking the product of these matricies
#'
#'   # sum across j
#'   # division by N_hat is a lin thing. Not in the binder formula.
#'   # what does rowSums return with an array?
#'   # I think I need to split dimension 3 (layer?) and do row sums over those matricies
#'   # what I want returned is an n * p matrix
#'
#'   term2 <- apply(-dj * wj * yi * p1 * p2, 3, rowSums)/N_hat
#'
#'   r <- term1 + term2
#'
#'   ## if you need to examine the parts
#'   attr(r, "terms") <- list(term1 = term1, term2 = term2, dj = dj, wj = wj, yi = yi, p1 = p1, p2 = p2, N_hat = N_hat)
#'
#'   r
#'
#' }
#'
#' #'
# # @export
# #'
# #' @details
# #' @details
#' #' should only work for p = 1, and be the same as the new version.
#' #' update: it is the same.
#' #' update: it was, but I change a denominator from n to N_hat
#'
#' old_calc_ui2 <- function(parts){
#'
#'   n <- length(parts$stat)
#'
#'   # first term is the same
#'
#'   term1 <- with(parts, {
#'
#'     weights * stat * (X - S1/S0)
#'
#'   })
#'
#'   # divide each exp(beta*X) by the series of risk sets.
#'   p1 <- with(parts, {
#'
#'     tcrossprod(exp_risk_score, 1/S0)
#'
#'   })
#'
#'   # replicate X n times to make a n*n matrix from n*1 X.
#'
#'   # to update this for multiple X, need to replicate each X
#'
#'   # minus S0/S1, which is also n*1 from each row
#'   p2 <- with(parts, {
#'
#'     Matrix(as.numeric(X), ncol = n, nrow = n) - Matrix(as.numeric(S1/S0), ncol = n, nrow = n, byrow = TRUE)
#'
#'   })
#'
#'   # yi is the indicator I(ti >= tj). Each ti needs to be compared to all tj
#'   # because ti is sorted, this should be an upper triangle matrix (if no ties)
#'
#'   yi <- with(parts, {
#'
#'     Matrix(rep(time, each = n) <= rep(time, n), ncol = n, nrow = n, byrow = TRUE)
#'
#'   })
#'
#'   # matrix with n reps of stat, as rows.
#'   dj <- with(parts, {
#'
#'     Matrix(stat, nrow = n, ncol = n, byrow = TRUE)
#'
#'   })
#'
#'   # matrix with n reps of weights, as rows.
#'   wj <- with(parts, {
#'
#'     Matrix(weights, nrow = n, ncol = n, byrow = TRUE)
#'
#'   })
#'
#'   # the second term is the row sums after taking the product of these matricies
#'
#'   # sum across j
#'   # division by n is a lin thing. Not in the binder formula.
#'   term2 <- rowSums(-dj * wj * yi * p1 * p2)/n
#'
#'   r <- term1 + term2
#'
#'   attr(r, "terms") <- list(term1, dj, wj, yi, p1, p2)
#'
#'   r
#'
#' }
#'
#'
#'
#'



