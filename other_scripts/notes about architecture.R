#' I want to test the difference between
#'   - Ripatti Likelihoods
#'   - REML likelihoods
#'   - Ripatti but without dropping the math likelihoods.
#'
#' For each of these likelihoods, I want to look at various models
#'   - random effects only
#'   - fixed and random effects
#'   - nested random effects.
#'   - crossed random effects.
#'
#' For each model, I want to look at a range of
#'   - sample sizes,
#'   - cluster sizes,
#'   - number of clusters
#'   - varying thetas
#'
#' Any for each combination of these things, I want several repetitions.
#'
#'
#' Set up:
#' the onerep() function could take a control list that specifies everything above,
#' plus all other (constant, but specified) parameters.
#'
#' alternatively, I could have a bespoke definition of onerep() for each parameter combination.
#' this would require a function that parses the control list and returns... an environment? with
#' everything in it. Like lme4.
#'
#' I think I would like to dispatch this stuff using a generic method call lp() or something,
#' and writing lp.ripatti lp.reml lp.full etc. This is probably overkill but I'd like practice
#' at doing that. This approach is really amenable to extention.
#'
#'




