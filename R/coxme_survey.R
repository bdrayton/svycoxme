
#' functions needed to make svycoxme go.
#'
#' svycoxme
#'
#'
#'
#'


#' Need a function that takes the random effect information from coxme and converts it into an offset.
#' This function will need to handle interaction terms from the random effects.


#' Generic method
#'
#' rescale should always be true in my opinion. Coxme often fails with unscaled weights. When it runs, the estimates it produces are really wrong.
#'

svycoxme <- function(formula, design, subset=NULL, ...){
  survey:::.svycheck(design)
  UseMethod("svycoxme", design)
}

svycoxme.DBIsvydesign <- function(formula, design, subset, ...){

  warning("coxme has not been implemented for class = \"DBIsvydesign\"")

  # return something?

}

## implemented
# svycoxme.survey.design <- function(formula, design, subset, ...){
#
#   warning("coxme has not been implemented for class = \"survey.design\"")
#
#   # return something?
#
# }

## implemented
# svycoxme.svyrep.design <- function(formula, design, subset, ...){
#
#   warning("coxme has not been implemented for \"class = svyrep.design\" ")
#
#   # return something?
#
# }

svycontrast <- function(stat, contrasts, add = FALSE, ...){

  UseMethod("svycontrast")

}

svycontrast.svycoxme <- function(){

  warning("svycontrast has not been implemented for \"class = svycoxme\" ")

  # return something?

}

svycontrast.svrepcoxme <- function(){

  warning("svycontrast has not been implemented for \"class = svrepcoxme\" ")

  # return something?

}

#' svycoxme.survey.design
#'
#' check next release for the fix to the rescale issues.
#'
#'
#' @export

svycoxme.survey.design <- function(formula, design, subset=NULL, rescale=TRUE, control = coxme::coxme.control(), ...){

  subset<-substitute(subset)
  subset<-eval(subset, model.frame(design),parent.frame())
  if (!is.null(subset))
    design<-design[subset,]

  if(any(weights(design)<0)) stop("weights must be non-negative")

  data<-model.frame(design)

  g<-match.call()
  g$formula<-eval.parent(g$formula)
  g$design<-NULL
  g$var<-NULL
  g$rescale <- NULL
  if (is.null(g$weights))
    g$weights<-quote(.survey.prob.weights)
  else
    g$weights<-bquote(.survey.prob.weights*.(g$weights))
  ## g[[1]]<-quote(coxph)
  g[[1]] <- quote(coxme::coxme)
  # g$data<-quote(data)
  g$data <- as.name("data")
  g$subset<-quote(.survey.prob.weights>0)
  g$control = control
  # g$model <- TRUE
  # the equivalent to this would be to set x = TRUE and y = TRUE
  # that may be useful if I define how residuals.coxme works.
  # currently data is passed to residuals.coxme, so I don't need x and y.
  # This is inefficient if x + y is a small subset of the data.

  ##need to rescale weights for stability
  ## unless the user doesn't want to
  if (rescale)
    data$.survey.prob.weights<-(1/design$prob)/mean(1/design$prob)
  else
    data$.survey.prob.weights <- 1/design$prob

  if (!all(all.vars(formula) %in% names(data)))
    stop("all variables must be in design= argument")

  # g<-with(list(data=data), eval(g))
  g<-eval(g)

  ## not needed.
  # if (inherits(g, "coxph.penal"))
  #   warning("svycoxph does not support penalised terms")

  g$call<-match.call()
  g$call[[1]]<-as.name(.Generic)
  g$printcall<-sys.call(-1)
  g$printcall[[1]]<-as.name(.Generic)
  # class(g)<-c("svycoxph", class(g))
  class(g)<-c("svycoxme", class(g))
  g$survey.design<-design

  nas<-g$na.action
  if (length(nas))
    design<-design[-nas,]

  # subset data here?
  dbeta.subset <- resid(object = g, data = data, type = "dfbeta")
  if (nrow(design)==NROW(dbeta.subset)){
    dbeta<-as.matrix(dbeta.subset)
  } else {
    dbeta<-matrix(0,ncol=NCOL(dbeta.subset),nrow=nrow(design))
    dbeta[is.finite(design$prob),]<-dbeta.subset
  }
  g$inv.info<-g$var

  if (inherits(design,"survey.design2"))
    g$variance<-svyrecvar(dbeta, design$cluster,
                     design$strata, design$fpc,
                     postStrata=design$postStrata)
  # I'm not sure if these will work correctly. Needs testing.
  else if (inherits(design, "twophase"))
    g$variance<-twophasevar(dbeta, design)
  else if(inherits(design, "twophase2"))
    g$variance<-twophase2var(dbeta, design)
  else if(inherits(design, "pps"))
    g$variance<-ppsvar(dbeta,design)
  else
    g$variance<-svyCprod(dbeta, design$strata,
                    design$cluster[[1]], design$fpc,design$nPSU,
                    design$certainty,design$postStrata)

  n_coef <- length(coef(g))

  g$var <- g$variance[seq(n_coef), seq(n_coef)]

  g$wald.test<-coef(g)%*%solve(g$var, coef(g))
  g$ll<-g$loglik
  # g$loglik<-NULL
  g$rscore<-NULL
  g$score<-NA
  g$degf.resid<-degf(design)-length(coef(g)[!is.na(coef(g))])+1

  g
}


#' svycoxme
#'
#' function for rep weight designs
#'
#' @import survey
#'
#' @export

svycoxme.svyrep.design <- function (formula, design, subset = NULL, rescale = NULL, ...,
                                    control = coxme::coxme.control(),
          return.replicates = FALSE, na.action, multicore = getOption("survey.multicore")){
  subset <- substitute(subset)
  subset <- eval(subset, design$variables, parent.frame())
  if (!is.null(subset))
    design <- design[subset, ]
  if (multicore && !requireNamespace(parallel, quietly = TRUE))
    multicore <- FALSE
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
    pwts <- design$pweights/mean(design$pweights)
  else if (rescale)
    pwts <- design$pweights/sum(design$pweights)
  if (is.data.frame(pwts))
    pwts <- pwts[[1]]
  if (!all(all.vars(formula) %in% names(data)))
    stop("all variables must be in design= argument")
  .survey.prob.weights <- pwts
  g$control = control
  full <- with(data, eval(g))
  # full <- eval(g)

  # Not needed with coxme
  # if (inherits(full, "coxph.penal"))
  #   warning("svycoxph does not support penalised terms")
  nas <- attr(full$model, "na.action")
  betas <- matrix(ncol = length(coef(full)), nrow = ncol(design$repweights))
  thetas <- matrix(ncol = length(coxme::VarCorr(full)), nrow = ncol(design$repweights))
  wts <- design$repweights
  if (!design$combined.weights) {
    pw1 <- pwts
    rwt <- pw1/mean(pw1)
  }
  else {
    rwt <- 1/mean(as.vector(wts[, 1]))
    pw1 <- rwt
  }
  if (length(nas))
    wts <- wts[-nas, ]
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
  ## Make fitter coxme, always. Would be interesting to test with different Surv() types to see what happens.
  ## Or would it be better to use coxph.fit, with an offset term?

  g$init <- beta0
  g$vinit <- theta0

## Ignore multicore for now
  # if (multicore) {
  #   betas <- do.call(rbind, parallel::mclapply(1:ncol(wts),
  #                                              function(i) {
  #                                                fitter(full$x, full$y, full$strata, full$offset,
  #                                                       coef(full), coxph.control(), as.vector(wts[,
  #                                                                                                  i]) * pw1 + EPSILON, full$method, names(full$resid))$coef
  #                                              }))
  # }
  # else {
    for (i in 1:ncol(wts)) {

      .survey.prob.weights <- as.vector(wts[,i]) * pw1 + EPSILON

      # handle errors here, but for future, consider coxme wrapper with error
      # handling that gets called instead of coxme.

      fit <- try(with(data, eval(g)))
      # fit <- with(data, eval(g))

      # assuming the method is the problem
      # if(inherits(fit, "try-error")) {
      #
      #   g$optpar <- list(method = "Brent",
      #                    control=list(reltol = 1e-5))
      #
      #   fit <- try(with(data, eval(g)))
      #
      # }

      # assuming the start values are the problem
      # if(inherits(fit, "try-error")) {
      #
      #   g$init <- NULL
      #   g$vinit <- NULL
      #
      #   fit <- try(with(data, eval(g)))
      #
      # }
      if(inherits(fit, "try-error")) {
        betas[i, ] <- rep(NA, ncol(betas))
        thetas[i, ] <- rep(NA, ncol(thetas))
      } else {
        betas[i, ] <- coef(fit)
        thetas[i, ] <- unlist(coxme::VarCorr(fit))
      }

    }
  # }
  if (length(nas))
    design <- design[-nas, ]
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
  }
  full$naive.var <- NULL
  full$wald.test <- coef(full) %*% solve(full$var, coef(full))
  full$loglik <- c(NA, NA)
  full$rscore <- NULL
  full$score <- NA
  full$degf.residual <- degf(design) + 1 - length(coef(full)[!is.na(coef(full))])
  class(full) <- c("svrepcoxme", "svycoxme", class(full))
  full$call <- match.call()
  full$printcall <- sys.call(-1)
  full$survey.design <- design
  full
}

#' modifies a formula.
#'
#' removes random effect terms
#' adds  + offset(offset)
#' intented to modify the formula slot in coxme call when fitting
#' coxph model with random effects as offsets.
#'
#' @export

fix_formula <- function(formula){

  # need to retain the formula as a call, but need a formula for formula1 to work.
  if(inherits(formula, "call")) formula <- eval(formula)

  # keep the bit with fixed effects.
  fixed_part <- coxme:::formula1(formula)$fixed

  # add offset(offset)
  new_form <- update.formula(fixed_part, ~ . + offset(offset))

  # convert it back to a call
  attributes(new_form) <- NULL

  return(new_form)

}


#' depends on survival coxph to do all the work.
#' converts the coxme random effects into an offset,
#' fits a coxph model with that offset, and then uses the
#' coxph machinery. Warning that residuals may be affected in weird ways that
#' we haven't thought through yet.
#' need the data set used for coxme, which is unlike most residual functions (that calculate residuals as part of the fitting function)
#'
#' @export


residuals.coxme <- function (object, data, weighted = TRUE,
                             type = c("score", "dfbeta", "dfbetas"), ...){

  type <- match.arg(type)
  otype <- type

  if(!any(type == c("score", "dfbeta", "dfbetas"))) {

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
    # vv <- get_information.coxme(object)

    vv <- vcov(object)

  # }

  strat <- object$strata
  if(!is.null(strat)) stop("Handling models with strata has not been implemented")

  if (type == "score") {

    parts <- make_parts.coxme(object, data)

    # the C++ method needs regular matrices
    parts <- lapply(parts, as.matrix)

    # rr <- resid <- calc_ui(parts, weighted = weighted)
    rr <- with(parts,{
          C_calc_ui(time_start = time_start,
                    time_stop = time_stop,
                    stat = stat,
                    weights = weights,
                    exp_risk_score = exp_risk_score,
                    S0 = S0,
                    S1_X = S1_X,
                    X = X,
                    weighted = TRUE)
          })

    if (otype == "dfbeta") {
        rr <- rr %*% vv
    }
    else if (otype == "dfbetas") {
        rr <- (rr %*% vv) %*% diag(sqrt(1/diag(vv)))
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

  print(object$survey.design, varnames = FALSE, design.summaries = FALSE, ...)
  NextMethod()

}


#' @method print svycoxme
#' @export
#'

print.svycoxme <- function (x, ...) {

  print(x$survey.design, varnames = FALSE, design.summaries = FALSE, ...)
  NextMethod()

}

#' @method logLik svycoxme
#' @export

logLik.svycoxme <- function(x, ...){

  NextMethod()

}

#' @method anova svycoxme
#' @export

anova.svycoxme <- function(x, ...){

  warning("anova has not been implemented for \"class = svycoxme\" ")

}

#' @method formula svycoxme
#' @export

formula.svycoxme <- function(x, ...){

  NextMethod()

}

#' @method predict svycoxme
#' @export

predict.svycoxme <- function(x, ...){

  warning("predict has not been implemented for \"class = svycoxme\"" )

  # this could be the method to use for data generation.
  # or something similar, that takes:
  # a data set
  # a model formula
  # fixed and random effects
  # further specification of the event scenario (recurrent, clustered, etc)
  # returns:
  # a matrix of responses, approriately identified, structure defined by the
  # Surv() call, plus any other specs required.



}

#' @method vcov svycoxme
#' @export

vcov.svycoxme <- function(object, ...){

NextMethod()

}


