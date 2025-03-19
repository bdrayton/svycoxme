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





