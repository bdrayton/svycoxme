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

  # need to use getFixedFormula, otherwise model.frame may complain about random effects terms.
  response <- model.response(model.frame(lme4:::getFixedFormula(form), data))

  n = nrow(response)

  response_type = attr(response, 'type')

  if(response_type == "right"){

    time_stop = response[,"time"]
    time_start = rep(0, length(time_stop))
    stat = response[,'status']

  } else if(response_type == "counting") {

    time_start = response[,"start"]
    time_stop = response[,"stop"]
    stat = response[,"status"]

  } else {

    stop("response type is not supported")

  }

  # reorder things
  time_order = order(time_stop, time_start)

  time_start = time_start[time_order]
  time_stop = time_stop[time_order]

  weights <- weights(coxme.object)
  if(is.null(weights)){
    weights <- rep(1, length(time_order))
  }
  weights <- weights[time_order]

  ds_sorted <- data[time_order, ]

  parsed_data <- lme4::lFormula(form, data = ds_sorted)

  # drop the intercept column from the X model.matrix
  # I use X a few times, so extracting. but could I use this each time?
  # does it make a difference anyone would care about in terms of memory or speed?
  X <- parsed_data$X[ , -1, drop = FALSE]
  Zt <- parsed_data$reTrms$Zt
  Z <- t(Zt)

  nX <- ncol(X)

  # I think these need to be p*1 and k*1 matrices
  beta <- Matrix::Matrix(coxme::fixef(coxme.object), ncol = 1)
  # need to unlist these and line them up with Zt
  # they should be lined up by default, but this will need some testing
  b <- Matrix::Matrix(unlist(coxme::ranef(coxme.object)), ncol = 1)

  risk_score <- X %*% beta + Matrix::crossprod(Zt, b)

  # weighted
  exp_risk_score <- weights * exp(risk_score)

  start_test <- time_stop > rep(time_start, each = n)
  stop_test  <- time_stop <= rep(time_stop, each = n)

  in_risk_set_matrix <- matrix(start_test & stop_test, nrow = n, byrow = TRUE)

  # this is S0_hat in binder, a n * 1 matrix.
  S0 = Matrix::crossprod(in_risk_set_matrix, exp_risk_score)

  # this is S1_hat, n * p matrix and n * q matrix
  exp_risk_score_X <- exp_risk_score[, rep(1, nX)] * X
  exp_risk_score_Z <- Matrix::Matrix(as.numeric(exp_risk_score) * as.numeric(Z), nrow = n)

  S1_X <- Matrix::crossprod(in_risk_set_matrix, exp_risk_score_X)

  S1_Z <- Matrix::crossprod(in_risk_set_matrix, exp_risk_score_Z)

  # I also need the penalty, divided into n parts.
  # will need D(theta), divided up appropriately for the calculation of D_i.
  # will need b * D(theta), divided up appropriately for the calculation of U_i
  # Will need to weight it appropriately.

  # this is the penalty matrix.
  # assume that each term has one theta. Must be shared frailty of some sort.
  # need to weight it.

  # get cluster level weight as... what?
  # mean of all weights in the cluster? only works if all obs in the cluster are sampled. it which case, weight for any obs
  # is the cluster level weight. this is our situation, so I just take one weight per cluster.

  # cluster_weights = weights %*% Z / colSums(Z)

  theta <- unlist(coxme::VarCorr(coxme.object))
  parsed_data$reTrms$Lambdat@x <- theta[parsed_data$reTrms$Lind]
  D <- solve(parsed_data$reTrms$Lambdat)

  # the slowest way of doing it...
  # wD = diag(cluster_weights@x) %*% D

  penalty <- Matrix::crossprod(b, D)

  # # divided by n and repeated
  # ui_penalty <- (penalty/n)[rep(1, n),]
  # # split penalty among the relevant cluster.

  ui_penalty = Z * (penalty/colSums(Z))[rep(1, n),]

  r <- list( stat = Matrix::Matrix(stat, ncol = 1),
             time_start = Matrix::Matrix(time_start, ncol = 1),
             time_stop = Matrix::Matrix(time_stop, ncol = 1),
             weights = Matrix::Matrix(weights, ncol = 1),
             S0 = S0,
             S1_X = S1_X,
             S1_Z = S1_Z,
             exp_risk_score = exp(risk_score),
             weighted_exp_risk_score = exp_risk_score,
             X = X,
             Z = Z,
             ui_penalty = ui_penalty,
             di_penalty = D)

  class(r) <- c("coxme_parts", class(r))

  return(r)

}


#' @export

make_parts.coxph <- function(coxph.object, data){

  form <- eval(coxph.object$call[[2]])

  model_frame <- model.frame(form, data)

  response <- model.response(model_frame)

  n = nrow(response)

  response_type = attr(response, 'type')

  if(response_type == "right"){

    time_stop = response[,"time"]
    time_start = rep(0, length(time_stop))
    stat = response[,'status']

  } else if(response_type == "counting") {

    time_start = response[,"start"]
    time_stop = response[,"stop"]
    stat = response[,"status"]

  } else {

    stop("response type is not supported")

  }

  time_order = order(time_stop, time_start)

  time_start = time_start[time_order]
  time_stop = time_stop[time_order]

  # get weights from the model object
  weights <- weights(coxph.object)
  if(is.null(weights)){
    weights <- rep(1, length(time_order))
  }
  weights <- weights[time_order]

  X <- model.matrix(coxph.object)[time_order, , drop = FALSE]
  nX <- ncol(X)

  # I think these need to be p*1 and k*1 matrices
  beta <- Matrix::Matrix(coef(coxph.object), ncol = 1)

  risk_score <- X %*% beta

  # weighted
  exp_risk_score <- weights * exp(risk_score)

  start_test <- time_stop > rep(time_start, each = n)
  stop_test  <- time_stop <= rep(time_stop, each = n)

  in_risk_set_matrix <- matrix(start_test & stop_test, nrow = n, byrow = TRUE)

  S0 = Matrix::crossprod(in_risk_set_matrix, exp_risk_score)

  exp_risk_score_X <- exp_risk_score[,rep(1, nX)] * X

  S1 <- Matrix::crossprod(in_risk_set_matrix, exp_risk_score_X)

  r <- list(stat = Matrix::Matrix(stat, ncol = 1, sparse = FALSE),
       time_start = Matrix::Matrix(time_start, ncol = 1),
       time_stop = Matrix::Matrix(time_stop, ncol = 1),
       weights = Matrix::Matrix(weights, ncol = 1),
       S0 = S0,
       S1 = S1,
       exp_risk_score = exp(risk_score),
       weighted_exp_risk_score = exp_risk_score,
       # S2 = at_risk_XtX,
       X = X,
       original_order = time_order)

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
      t_S1S0_squared <- Matrix::Matrix(S1S0_squared, ncol = 1)
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
    S1_XtS1_Z <- Matrix::Matrix(t(mapply(Matrix::tcrossprod, split(S1_X, row(S1_X)), split(S1_Z, row(S1_Z)))))

    # denominator for first term
    ncols <- ncol(X) * ncol(Z)
    S0tS0 <- Matrix::Matrix(rep(apply(S0, 1, tcrossprod), ncols), nrow = nrow(S2_XtZ))

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

get_information <- function(x, ...){

    UseMethod("get_information", x)

}

#'
#' @export

get_information.coxph <- function(coxph.object){

    vv <- coxph.object$naive.var

  if(is.null(vv)){
    vv <- vcov(coxph.object)
  }

  vv

}

#'
#' @export

get_information.coxme <- function(coxme.object){

  # this is the information matrix, with the dense fixed effects part in the bottom right corner.
  # converting to a regular matrix. Probably better to retain the sparse structure. I don't know how to do that
  # and still do the re-ordering.

  inf1 <- as.matrix(coxme.object$variance)

  # rearrange to put it it the top left corner.
  nfixed <- length(coxme::fixed.effects(coxme.object))
  nrandom <- length(unlist(coxme::random.effects(coxme.object)))

  ntotal <- nfixed + nrandom

  inf2 <- inf1[c((nrandom+1):(ntotal), 1:nrandom) , c((nrandom+1):(ntotal), 1:nrandom)]

  inf2

}





#'
#' @export

calc_ui <- function(x, ...){

  UseMethod("calc_ui", x)

}


#'
#' @export

calc_ui.coxph_parts <- function(parts, weighted = TRUE){

  # it's easier to debug if you can access the parts without using parts$ or with()
  # env <- environment()
  #
  # lapply(names(parts), function(part){
  #   assign(part, parts[[part]], pos = env)
  # })

  lin_term2 <- parts$X

  lin_term2[] <- NA

  for(i in seq_len(nrow(lin_term2))) {

      irep <- rep(i, length(parts$stat))
      # Yi_at_tj = (time_start[irep] < time_stop & time_stop[irep] >= time_stop)
      # lin_term2[i, ] <- Matrix::colSums((stat * weights * Yi_at_tj * exp_risk_score[irep] * (1/S0)) * (X[irep, ] - S1/S0))

    lin_term2[i, ] <-
      with(parts,{
      Yi_at_tj = (time_start[irep] < time_stop & time_stop[irep] >= time_stop)
      Matrix::colSums((stat * weights * Yi_at_tj * exp_risk_score[irep] * (1/S0)) * (X[irep, ] - S1/S0))
    })

  }

  lin_term1 <- with(parts,{
    stat * (X - S1/S0)
  })

  lin_score <- lin_term1 - lin_term2

  if (weighted) {
    lin_score = parts$weights * lin_score
  }

  lin_score

}


#'
#' @export

calc_ui.coxme_parts <- function(parts, weighted = TRUE){

  lin_X_term2 <- parts$X
  lin_Z_term2 <- parts$Z

  lin_X_term2[] <- NA
  lin_Z_term2[] <- NA

  nX <- ncol(lin_X_term2)
  nZ <- ncol(lin_Z_term2)

  # for column reps
  nXreps = rep(1, nX)
  nZreps = rep(1, nZ)

  for(i in seq_len(nrow(lin_X_term2))) {

    irep <- rep(i, length(parts$stat))
    Yi_at_tj <- with(parts, time_start[irep] < time_stop & time_stop[irep] >= time_stop)

    temp = with(parts, stat * weights * Yi_at_tj * exp_risk_score[irep] * (1/S0))

    lin_X_term2[i, ] <- with(parts,{
      Matrix::colSums(temp[, nXreps] * (X[irep, ] - S1_X/S0))
    })

    lin_Z_term2[i, ] <- with(parts,{
      Matrix::colSums(temp[, nZreps] * (Z[irep, ] - S1_Z/S0))
    })

  }

  lin_term2 <- cbind(lin_X_term2, lin_Z_term2)

  lin_term1_beta <- with(parts, {

    # weights * stat * (X - S1_X/S0)
    stat[, nXreps] * (X - S1_X/S0)

  })

  lin_term1_b <- with(parts, {
    # don't ignore penalty
    if (weighted){
      stat[, nZreps] * (Z - S1_Z/S0) - ui_penalty/weights[,nZreps]
    } else {
      stat[, nZreps] * (Z - S1_Z/S0) - ui_penalty
    }
    # stat[, nZreps] * (Z - S1_Z/S0)
  })

  lin_term1 <- cbind(lin_term1_beta, lin_term1_b)

  lin_score <- lin_term1 - lin_term2

  if (weighted) {
    lin_score = parts$weights[,c(nXreps, nZreps)] * lin_score
  }

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



