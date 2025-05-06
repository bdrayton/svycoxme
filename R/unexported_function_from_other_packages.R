# define methods than need to be accessed with :::
# This means I don't depend on methods that aren't intended to be exported.


# survey:::.svycheck

.svycheck <- function (object)
{
  if (inherits(object, "survey.design") && !is.null(object$nPSU))
    warning("This is an old-style design object. Please use as.svydesign2 to update it.")
}


#lme4:::RHSForm

RHSForm <- function (form, as.form = FALSE)
{
  rhsf <- form[[length(form)]]
  if (as.form)
    stats::reformulate(deparse(rhsf))
  else rhsf
}


# lme4:::`RHSForm<-`

`RHSForm<-` <- function (formula, value)
{
  formula[[length(formula)]] <- value
  formula
}


# lme4:::getFixedFormula

getFixedFormula <- function (form)
{
  RHSForm(form) <- lme4::nobars(RHSForm(form))
  form
}

# survey:::ppsvar

ppsvar <- function (x, design)
{
  postStrata <- design$postStrata
  est <- design$variance
  if (!is.null(postStrata)) {
    for (psvar in postStrata) {
      if (inherits(psvar, "greg_calibration")) {
        if (psvar$stage == 0) {
          x <- qr.resid(psvar$qr, x/psvar$w) * psvar$w
        }
        else {
          stop("calibration within clusters not yet available for PPS designs")
        }
      }
      else {
        psw <- attr(psvar, "weights")
        postStrata <- as.factor(psvar)
        psmeans <- rowsum(x/psw, psvar, reorder = TRUE)/as.vector(table(factor(psvar)))
        x <- x - psmeans[match(psvar, sort(unique(psvar))),
        ] * psw
      }
    }
  }
  dcheck <- design$dcheck
  if (length(dcheck) != 1)
    stop("Multistage not implemented yet")
  rval <- switch(est, HT = htvar.matrix(rowsum(x, dcheck[[1]]$id,
                                               reorder = FALSE), dcheck[[1]]$dcheck), YG = ygvar.matrix(rowsum(x,
                                                                                                               dcheck[[1]]$id, reorder = FALSE), dcheck[[1]]$dcheck),
                 stop("can't happen"))
  rval
}

# survey:::htvar.matrix

htvar.matrix <- function (xcheck, Dcheck)
{
  if (is.null(dim(xcheck)))
    xcheck <- as.matrix(xcheck)
  rval <- apply(xcheck, 2, function(xicheck) apply(xcheck,
                                                   2, function(xjcheck) as.matrix(Matrix::crossprod(xicheck,
                                                                                                    Dcheck %*% xjcheck))))
  if (is.null(dim(rval)))
    dim(rval) <- c(1, 1)
  rval
}

# survey:::ygvar.matrix

ygvar.matrix <- function (xcheck, Dcheck)
{
  ht <- htvar.matrix(xcheck, Dcheck)
  if (is.null(dim(xcheck))) {
    corr <- sum(Dcheck %*% (xcheck * xcheck))
  }
  else {
    corr <- apply(xcheck, 2, function(xicheck) apply(xcheck,
                                                     2, function(xjcheck) sum(Dcheck %*% (xicheck * xjcheck))))
  }
  rval <- ht - corr
}

