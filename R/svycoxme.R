




#' Survey-weighted mixed-effects Cox models
#'
#' Fit a mixed-effect proportional hazards model to data from a complex design.
#'
#' @param formula Model formula.
#' @param design `survey.design` object. It must contain all variables in the formula.
#' @param subset Expression to select a subpopulation.
#' @param rescale Rescale weights to improve numerical stability.
#' @param ... Other arguments passed to `coxme`.
#'
#' @return An object of class `svycoxme`.
#'
#' @export
#'
#' @useDynLib svycoxme, .registration=TRUE
#' @importFrom Rcpp evalCpp

svycoxme <- function(formula, design, subset = NULL, ...) {
  survey:::.svycheck(design)
  UseMethod("svycoxme", design)
}


#' @method svycoxme DBIsvydesign
#' @export

svycoxme.DBIsvydesign <- function(formula, design, subset, ...) {
  warning("coxme has not been implemented for class = \"DBIsvydesign\"")

  call = match.call()

  invisible(call)

}

#' @method svycoxme survey.design
#' @export


svycoxme.survey.design <-
  function(formula,
           design,
           subset = NULL,
           rescale = TRUE,
           ...) {
    subset <- substitute(subset)
    subset <- eval(subset, model.frame(design), parent.frame())
    if (!is.null(subset))
      design <- design[subset, ]

    if (any(weights(design) < 0))
      stop("weights must be non-negative")

    data <- model.frame(design)

    g <- match.call()
    g$formula <- eval.parent(g$formula)
    g$design <- NULL
    g$var <- NULL
    g$rescale <- NULL
    if (is.null(g$weights))
      g$weights <- quote(.survey.prob.weights)
    else
      g$weights <- bquote(.survey.prob.weights * .(g$weights))
    ## g[[1]]<-quote(coxph)
    g[[1]] <- quote(coxme::coxme)

    # g$data<-quote(data)
    g$data <- as.name("data")
    g$subset <- quote(.survey.prob.weights > 0)
    # g$model <- TRUE
    # the equivalent to this would be to set x = TRUE and y = TRUE
    # that may be useful if I define how residuals.coxme works.
    # currently data is passed to residuals.coxme, so I don't need x and y.
    # This is inefficient if x + y is a small subset of the data.

    ##need to rescale weights for stability
    ## unless the user doesn't want to
    if (rescale)
      data$.survey.prob.weights <- (1 / design$prob) / mean(1 / design$prob)
    else
      data$.survey.prob.weights <- 1 / design$prob

    if (!all(all.vars(formula) %in% names(data)))
      stop("all variables must be in design= argument")

    # g<-with(list(data=data), eval(g))
    g <- eval(g)

    ## not needed.
    # if (inherits(g, "coxph.penal"))
    #   warning("svycoxph does not support penalised terms")

    g$call <- match.call()
    g$call[[1]] <- as.name(.Generic)
    g$printcall <- sys.call(-1)
    g$printcall[[1]] <- as.name(.Generic)
    # class(g)<-c("svycoxph", class(g))
    class(g) <- c("svycoxme", class(g))
    g$survey.design <- design

    nas <- g$na.action
    if (length(nas))
      design <- design[-nas, ]
    # subset data here?
    dbeta.subset <-
      resid(g,
            data = data,
            weighted = TRUE,
            type = "dfbeta")
    if (nrow(design) == NROW(dbeta.subset)) {
      dbeta <- as.matrix(dbeta.subset)
    } else {
      dbeta <- matrix(0,
                      ncol = NCOL(dbeta.subset),
                      nrow = nrow(design))
      dbeta[is.finite(design$prob), ] <- dbeta.subset
    }
    g$inv.info <- g$var


    if (inherits(design, "survey.design2")) {
      g$variance <- svyrecvar(dbeta,
                              design$cluster,
                              design$strata,
                              design$fpc,
                              postStrata = design$postStrata)
    }
    # I'm not sure if these will work correctly. Needs testing.
    else if (inherits(design, "twophase")) {
      warning('twophase design has not been tested')
      g$variance <- twophasevar(dbeta, design)
    }
    else if (inherits(design, "twophase2")) {
      warning('twophase2 design has not been tested')
      g$variance <- twophase2var(dbeta, design)
    }
    else if (inherits(design, "pps")) {
      warning('pps design has not been tested')
      g$variance <- ppsvar(dbeta, design)
    }
    else {
      g$variance <- svyCprod(
        dbeta,
        design$strata,
        design$cluster[[1]],
        design$fpc,
        design$nPSU,
        design$certainty,
        design$postStrata
      )
    }

    n_coef <- length(coef(g))

    # fixed effects only
    g$var <- g$variance[seq(n_coef), seq(n_coef)]

    g$wald.test <- coef(g) %*% solve(g$var, coef(g))
    g$ll <- g$loglik
    g$loglik <- rep(NA_real_, 3)
    g$rscore <- NULL
    g$score <- NA
    g$degf.resid <- degf(design) - length(coef(g)[!is.na(coef(g))]) + 1

    g
  }



#' @method svycoxme svyrep.design
#' @export

svycoxme.svyrep.design <-
  function (formula,
            design,
            subset = NULL,
            rescale = NULL,
            ...,
            control = coxme::coxme.control(),
            starts = "mean",
            return.replicates = FALSE,
            vfixed = NULL,
            na.action,
            multicore = getOption("survey.multicore"),
            cores = 2) {
    subset <- substitute(subset)
    subset <- eval(subset, design$variables, parent.frame())
    if (!is.null(subset))
      design <- design[subset,]
    if (multicore && !requireNamespace("future.apply", quietly = TRUE))
      multicore <- FALSE
    if (multicore) {
      message("future.apply is used for parallel processing")
      if (cores > parallelly::availableCores()) {
        cores = parallelly::availableCores()
        warning("cores is greater than parallelly::availableCores()")
      }
      message("replicate fits will be processed on ", cores, " cores")
    }
    data <- design$variables
    g <- match.call()
    g$design <- NULL
    g$return.replicates <- NULL
    g$weights <- quote(.survey.prob.weights)
    ## change to coxme
    # g[[1]] <- quote(coxph)
    g[[1]] <- quote(coxme::coxme)
    g$x <- TRUE
    scale <- design$scale
    rscales <- design$rscales
    if (is.null(rescale))
      pwts <- design$pweights / mean(design$pweights)
    else if (rescale)
      pwts <- design$pweights / sum(design$pweights)
    if (is.data.frame(pwts))
      pwts <- pwts[[1]]
    if (!all(all.vars(formula) %in% names(data)))
      stop("all variables must be in design= argument")
    .survey.prob.weights <- pwts
    g$control = control
    g$vfixed = vfixed
    full <- with(data, eval(g))
    # full <- eval(g)

    # Not needed with coxme
    # if (inherits(full, "coxph.penal"))
    #   warning("svycoxph does not support penalised terms")
    nas <- attr(full$model, "na.action")
    nreps <- ncol(design$repweights)
    betas <- matrix(ncol = length(coef(full)), nrow = nreps)
    thetas <-
      matrix(ncol = length(coxme::VarCorr(full)), nrow = nreps)
    full_frails <- unlist(coxme::random.effects(full))
    full_frails_names <- names(full_frails)
    frails <- matrix(ncol = length(full_frails), nrow = nreps)

    wts <- design$repweights
    if (!design$combined.weights) {
      pw1 <- pwts
      rwt <- pw1 / mean(pw1)
    }
    else {
      rwt <- 1 / mean(as.vector(wts[, 1]))
      pw1 <- rwt
    }
    if (length(nas))
      wts <- wts[-nas,]
    beta0 <- coef(full)
    # vinit needs an unnamed list of start values, matched by position.
    theta0 <- unname(lapply(coxme::VarCorr(full), unname))
    EPSILON <- 1e-10

    # if (full$method %in% c("efron", "breslow")) {
    #   if (attr(full$y, "type") == "right")
    #     fitter <- coxph.fit
    #   else if (attr(full$y, "type") == "counting")
    #     fitter <- survival::agreg.fit
    #   else stop("invalid survival type")
    # }
    # else fitter <- survival::agexact.fit
    ## Make fitter coxme, always. Would be interesting to test with different Surv() types to see what happens. I think this is done now. only two surv types work.
    ## Or would it be better to use coxph.fit, with an offset term? no

    g$init <- beta0
    g$vinit <- theta0

    # will fix theta, which should decrease the perturbation is bootstrap fixed effects...
    # edit: g$vfixed needs to be set to the true value of theta.
    # if(fix.v){
    #   g$vfixed <- theta0
    # }



    ## multicore
    if (multicore) {

      old_plan = future::plan()

      future::plan(future::multisession, workers = cores)

      replicate_fit_function <- function(i){

        weights_temp = as.vector(wts[, i]) * pw1

        .survey.prob.weights <- weights_temp[which(weights_temp != 0)]

        data_temp = data[which(weights_temp != 0),]

        # handle errors here, but for future, consider coxme wrapper with error
        # handling that gets called instead of coxme.

        fit <- try(with(data_temp, eval(g)))

        if (inherits(fit, "try-error")) {
          list(
             beta   = rep(NA, ncol(betas))
            ,theta  = rep(NA, ncol(thetas))
            ,frails = rep(NA, ncol(frails))
          )
        } else {
          list(
             beta   = coef(fit)
            ,theta  = unlist(coxme::VarCorr(fit))
            ,frails = unlist(coxme::random.effects(fit))
          )
        }
      }

      replicate_fits <- future.apply::future_lapply(1:ncol(wts), replicate_fit_function)

      # unpack the list of replicate fit results into the matrices.
      # leaving this as a separate step, as indexing into the same matrix from multiple workers
      # seem like it could lead to problems (I'm not sure).

      for(i in 1:ncol(wts)){
        betas[i,]  <- replicate_fits[[i]][["beta"]]
        thetas[i,] <- replicate_fits[[i]][["theta"]]
        new_frails <- replicate_fits[[i]][["frails"]]
        frails[i, which(names(full_frails) %in% names(new_frails))] <- new_frails
      }

      future::plan(old_plan)

    }
    else {
      for (i in 1:ncol(wts)) {
        # .survey.prob.weights <- as.vector(wts[,i]) * pw1 + EPSILON

        weights_temp = as.vector(wts[, i]) * pw1

        .survey.prob.weights <- weights_temp[which(weights_temp != 0)]

        data_temp = data[which(weights_temp != 0),]

        # handle errors here, but for future, consider coxme wrapper with error
        # handling that gets called instead of coxme.

        fit <- try(with(data_temp, eval(g)))

        if (inherits(fit, "try-error")) {
          ## the cells are already NA by default
          # betas[i, ] <- rep(NA, ncol(betas))
          # thetas[i, ] <- rep(NA, ncol(thetas))
          # frails[i, ] <- rep(NA, ncol(frails))
        } else {
          betas[i,] <- coef(fit)
          thetas[i,] <- unlist(coxme::VarCorr(fit))
          new_frails <- unlist(coxme::random.effects(fit))
          frails[i, which(names(full_frails) %in% names(new_frails))] <- new_frails
        }

        # updating initial betas and thetas may improve computation time,
        # particularly in later iterations. nah, it doesn't. It's not slower either.
        if (starts == "mean") {
          g$init <- colMeans(betas, na.rm = TRUE)
          g$vinit <- colMeans(thetas, na.rm = TRUE)
        }


      }
    }
    if (length(nas))
      design <- design[-nas,]
    v <- svrVar(betas, scale, rscales, mse = design$mse, coef = beta0)
    full$var <- v
    # add in stuff for bootstrapping theta
    v <- svrVar(thetas, scale, rscales, mse = FALSE, coef = theta0)
    full$vvar <- v

    if (return.replicates) {
      attr(betas, "scale") <- design$scale
      attr(betas, "rscales") <- design$rscales
      attr(betas, "mse") <- design$mse
      full$replicates <- betas
      full$replicates_theta <- thetas
      full$replicates_frail <- frails
    }

    full$naive.var <- NULL
    full$wald.test <- coef(full) %*% solve(full$var, coef(full))
    full$loglik <- rep(NA_real_, 3)
    full$rscore <- NULL
    full$score <- NA
    full$degf.residual <-
      degf(design) + 1 - length(coef(full)[!is.na(coef(full))])
    class(full) <- c("svrepcoxme", "svycoxme", class(full))
    full$call <- match.call()
    full$printcall <- sys.call(-1)
    full$survey.design <- design
    full
  }

#' @exportS3Method survey::svycontrast

svycontrast.svycoxme <- function(stat, contrasts, add = FALSE, ...) {
  stop("svycontrast has not been implemented for \"class = svycoxme\" ")

  # return something?

}

#' @exportS3Method survey::svycontrast

svycontrast.svrepcoxme <- function(stat, contrasts, add = FALSE, ...) {
  stop("svycontrast has not been implemented for \"class = svrepcoxme\" ")

  # return something?

}


#' @exportS3Method stats::AIC
AIC.svycoxme <- function(object, ...) {
  stop("No AIC for survey models")

  # return something?

}


#' Calculate residuals for a 'coxme' fit
#'
#' Calculates score, dfbeta, or dfbetas residuals for a mixed-effects proportional hazards model. Only fixed-effect components are calculated; see Details.
#'
#' An observation's contribution to the score vector includes values for every fixed and random effect in the fitted model. In many cases, the number of random effects will be large, and most residuals will be zero. Until efficient sparse computation is implemented, it is too expensive computationally and on memory to calculate the random effect residual terms, so they are excluded. This is likely to change, and the parameter \code{include_re} is include for future expansion.
#'
#' @param object an object inheriting from class \code{coxme}. This includes the output from \code{coxme} and \code{svycoxme} functions.
#' @param data the data used to generate \code{object}.
#' @param type character string indicating the type of residual desired. Possible values are "score", "dfbeta"', "dfbetas".
#' @param weighted	if TRUE and the model was fit with case weights, then the weighted residuals are returned.
#' @param include_re logical flag indicating if residuals for random effects should be returned. This flag is currently ignored; see Details.
#' @param ...	other unused arguments.
#'
#' @return The score residuals are each observation's contribution to the score vector. Two transformations of this are often more useful: dfbeta is the approximate change in the coefficient vector if that observation were dropped, and dfbetas is the approximate change in the coefficients, scaled by the standard error for the coefficients.
#'
#'
#' @method residuals coxme
#' @export

residuals.coxme <- function (object,
                             data,
                             type = c("score", "dfbeta", "dfbetas"),
                             weighted = (type %in% c("dfbeta", "dfbetas")),
                             include_re = FALSE,
                             ...) {
  type <- match.arg(type)
  otype <- type

  if (!any(type == c("score", "dfbeta", "dfbetas"))) {
    stop(paste(type, " residuals have not been implemented."))
  }

  if (type == "dfbeta" || type == "dfbetas") {
    otype <- type
    type <- "score"
    if (missing(weighted))
      weighted <- TRUE
  }

  form <- eval(object$call[[2]])

  # need to use getFixedFormula, otherwise model.frame may complain about random effects terms.
  response <-
    model.response(model.frame(lme4:::getFixedFormula(form), data))

  n = as.integer(nrow(response))

  response_type = attr(response, 'type')

  if (response_type == "right") {
    time_stop = response[, "time"] |> unname()
    time_start = rep(0, length(time_stop)) |> unname()
    stat = response[, 'status'] |> unname()
  } else if (response_type == "counting") {
    time_start = response[, "start"] |> unname()
    time_stop = response[, "stop"] |> unname()
    stat = response[, "status"] |> unname()
  } else {
    stop("response type is not supported")
  }

  # reorder things later. depends on strata too.
  # time_order = order(time_stop, time_start)

  # time_start = time_start[time_order]
  # time_stop = time_stop[time_order]

  weights <- weights(object)

  if (is.null(weights)) {
    weights <- rep(1, n)
  }

  # weights <- weights[time_order]

  # ds_sorted <- data[time_order, ]

  # parsed_data <- lme4::lFormula(form, data = ds_sorted)
  parsed_data <- lme4::lFormula(form, data = data)

  # drop the intercept column from the X model.matrix
  X <- parsed_data$X[,-1, drop = FALSE]
  Zt <- parsed_data$reTrms$Zt
  Z <- Matrix::t(Zt)

  nX <- as.integer(ncol(X))

  # use object$linear.predictor.
  #
  #   beta <- Matrix::Matrix(coxme::fixef(object), ncol = 1)
  #   b <- Matrix::Matrix(unlist(coxme::ranef(object)), ncol = 1)
  #
  #   risk_score <- X %*% beta + Matrix::crossprod(Zt, b)
  #


  vv <- vcov(object)

  strat <- object$strata
  method <- object$ties

  Terms <- object$terms
  strats <- attr(Terms, "specials")$strata

  # set up strata, order.
  if (is.null(strat)) {
    ord <- order(time_stop,-stat)
    newstrat <- integer(n)
    istrat <- integer(n)
  } else {
    istrat <- as.integer(strat)
    ord <- order(istrat, time_stop,-stat)
    newstrat <- c(diff(as.numeric(istrat[ord])) != 0, 1)
  }
  newstrat[n] <- 1
  X <- X[ord,]
  time_start <- time_start[ord]
  time_stop <- time_stop[ord]
  stat <- stat[ord]

  exp_risk_score = exp(object$linear.predictor)[ord]
  # exp_risk_score <- exp(risk_score)[ord]

  istrat <- istrat[ord]

  if (is.null(strat)) {
    sort1 <- order(time_start)
  } else {
    sort1 <- order(istrat, time_start)
  }

  storage.mode(time_start) <- "double"
  storage.mode(time_stop) <- "double"
  storage.mode(stat) <- "double"
  storage.mode(X) <- "double"

  storage.mode(newstrat) <- "integer"
  storage.mode(exp_risk_score) <- storage.mode(weights) <- "double"
  if (type == "score") {
    resid = agscore3(
      time_start,
      time_stop,
      stat,
      covar = X,
      strata = istrat,
      score = exp_risk_score,
      weights = weights[ord],
      sort1 = sort1 - 1L,
      method = as.integer(method == "efron")
    )

  }

  rr <- matrix(0, n, nX)

  if (nX > 1) {
    rr[ord, ] <- resid
    dimnames(rr) <- list(as.character(seq_len(n)),
                         names(object$coefficients))
  }
  else {
    rr[ord] <- resid
  }

  if (otype == "dfbeta") {
    rr <- rr %*% vv
  }

  else if (otype == "dfbetas") {
    rr <- (rr %*% vv) %*% diag(sqrt(1 / diag(vv)))
  }

  if (weighted) {
    rr <- rr * weights
  }

  if (!is.null(object$na.action)) {
    rr <- naresid(object$na.action, rr)
  }

  rr

}




# this is the old method. replaced with much faster method above.
#' @export
residuals2.coxme <-
  function (object,
            data,
            weighted = TRUE,
            include_re = FALSE,
            type = c("score", "dfbeta", "dfbetas"),
            ...) {
    type <- match.arg(type)
    otype <- type

    if (!any(type == c("score", "dfbeta", "dfbetas"))) {
      stop(paste(type, " residuals have not been implemented."))

    }


    if (type == "dfbeta" || type == "dfbetas") {
      otype <- type
      type <- "score"
    }

    # vv <- object$naive.var
    # if (is.null(vv)){

    # get_information get the covariance matrix for fixed and random effects
    # i only need fixed effects because i'm ignoring the random effects (too slow to compute the residuals)
    # actually, I may want, depends on include_re.

    if (include_re) {
      vv <- get_information.coxme(object)
    } else {
      vv <- vcov(object)
    }

    strat <- object$strata
    if (!is.null(strat))
      stop("Handling models with strata has not been implemented")

    if (type == "score") {
      parts <- make_parts.coxme(object, data)

      # I've modified make_parts.coxme to return class(matrix), so don't need this conversion.
      # the C++ method needs regular matrices
      # parts <- lapply(parts, as.matrix)

      # rr <- resid <- calc_ui(parts, weighted = weighted)
      if (!include_re) {
        rr <- with(parts, {
          C_calc_ui(
            time_start = time_start,
            time_stop = time_stop,
            stat = stat,
            weights = weights,
            exp_risk_score = exp_risk_score,
            S0 = S0,
            S1_X = S1_X,
            X = X,
            weighted = TRUE
          )
        })
      } else {
        stop("residuals for random effects are not implemented")
        rr <- with(parts, {
          C_calc_ui(
            time_start = time_start,
            time_stop = time_stop,
            stat = stat,
            weights = weights,
            exp_risk_score = exp_risk_score,
            S0 = S0,
            S1_X = cbind(S1_X, S1_Z),
            X = cbind(X, Z),
            weighted = TRUE
          )
        })

        rr = rr - cbind(matrix(0, nrow = nrow(rr), ncol = ncol(parts$X)), parts$ui_penalty)

      }

      if (otype == "dfbeta") {
        rr <- rr %*% vv
      }
      else if (otype == "dfbetas") {
        rr <- (rr %*% vv) %*% diag(sqrt(1 / diag(vv)))
      }
    }

    if (!is.null(object$na.action)) {
      rr <- naresid(object$na.action, rr)
    }

    rr

  }

#' @method summary svycoxme
#' @export
#'

summary.svycoxme <- function(object, ...) {
  print(object$survey.design,
        varnames = FALSE,
        design.summaries = FALSE,
        ...)
  NextMethod()

}


#' @method print svycoxme
#' @export
#'

print.svycoxme <- function (x, ...) {
  print(x$survey.design,
        varnames = FALSE,
        design.summaries = FALSE,
        ...)
  NextMethod()

}

#' @method logLik svycoxme
#' @export

logLik.svycoxme <- function(x, ...) {
  NextMethod()

}

#' @method anova svycoxme
#' @export

anova.svycoxme <- function(x, ...) {
  warning("anova has not been implemented for \"class = svycoxme\" ")

}

#' @method formula svycoxme
#' @export

formula.svycoxme <- function(x, ...) {
  NextMethod()

}

#' @method predict svycoxme
#' @export

predict.svycoxme <- function(x, ...) {
  warning("predict has not been implemented for \"class = svycoxme\"")


}

#' @method vcov svycoxme
#' @export

vcov.svycoxme <- function(object, ...) {
  NextMethod()

}
