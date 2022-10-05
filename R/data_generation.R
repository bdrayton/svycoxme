
#' extract term.labels attribute
#'
#' \link{\code{terms}} returns the formula it is given with a decomposition of attributes attached.
#' This function is a wrapper for `attr(., "term.labels")`
#'
#'
#'


get_group_label <- function(random_term){

  term_list <- coxme:::formula2(random_term)

  as.character(term_list$group)

}


#' Check formula against dists names
#'
#' Checks that the variables in the right hand side of the formula all have matching
#' expressions is the dists list.
#'
#' @param var_labels Names of the variables in the right hand side of the formula
#' @param dist_labels Names of the expressions in dists
#'
#' @export


check_var_labels <- function(var_labels, dist_labels){

  test <- var_labels %in% dist_labels

  not_in_dist_labels <- paste0(var_labels[!test], collapse = ", ")

  if(!all(test))
    stop(paste0("All variables in formula must be in dists.\n  These variables are in your formula but not in dists: ", not_in_dist_labels))

}



#'  create a dataset for testing cox model stuff.
#'
#' @param formula Model equation for the model you plan to fit to the data. Only the right hand side is used at the moment.
#' @param dists expressions for generating the variables in the formula
#' @param dist_args any values that are needed for the expressions in dists
#' @param coefficients vector of fixed effect coefficients. If there are more coefficients than fixed effects, only those up to the number of fixed effects will be used and the rest ignored.
#' @param random_effect_variance named list of thetas. Names to match formula terms.
#' @param seed passed to set.seed. Use if you need the random processes to be reproducible.
#'
#'
#'
#' @export

one_dataset <- function(formula, dists, dist_args, coefficients = c(), random_effect_variance = list(), seed = NULL,
                        random_effect_seed = NULL){

  arbitraty_seed <- runif(1, -99999, 99999)

  if(!is.null(random_effect_seed) & (length(random_effect_variance) != length(random_effect_seed)))
    stop("If random_effect_seed is not NULL, it must be the same length as random_effect_variance")

  original_formula <- formula

  if(!is.null(seed)) set.seed(seed)

  fixed_terms <- stats::terms(lme4::nobars(formula))

  # These are the fixed term variables in the formula, even if there are interactions.
  fixed_term_variables <- colnames(attr(fixed_terms, "factors"))

  # checking coefficients needs to account for interactions, which depends on levels of the variables.
  # im going to remove this for now, and reinstate a more sophisticated version if that becomes necessary.
  # if(length(fixed_term_variables) > length(coefficients))
  #   stop("Each fixed effect variable in the formula must have a corresponding coefficient")

  # all_dists_are_formula <- all(sapply(dists, class) == "formula")

  # if a user forgets the ~ on one of there expressions, it will fail with an error before now.
  # if(!all_dists_are_formula)
  #   stop("All terms in dists must be left hand formulas e.g. ~rnorm(n)")


  bars = lme4::findbars(formula)

  barnames = lme4:::barnames(bars)

  random_terms = terms(as.formula(paste("~", paste(c(1, barnames), collapse = "+"))))

  random_term_variables <- rownames(attr(random_terms, "factors"))

  variables <- c(fixed_term_variables, random_term_variables)

  check_var_labels(variables, names(dists))

  vars <- lapply(dists[variables], lazyeval::f_eval, data = dist_args)

  vars_df <- list2DF(vars)

  error <- lazyeval::f_eval(dists$error, data = dist_args)

  vars_df$stat <- lazyeval::f_eval(dists$stat, data = dist_args)

  lp_random_effects <- list()

  if(length(barnames)) {

    vars_df$ytemp <- 1

    formula <- update(formula, ytemp ~ .)

    parsed_data <- lme4::lFormula(formula, data = vars_df)

    vars_df <- dplyr::select(vars_df, -ytemp)

    # for each random effect term, need to generate the random effect terms, and join them to the data.
    # if terms don't exist in the data, they will be in flist.

    re_terms <- names(parsed_data$reTrms$cnms)
    names(re_terms) <- re_terms

    if(length(re_terms) != length(random_effect_variance))
      stop("I need one variance for each random effect term. Note that (1 | M1/M2) breaks down into M1, M2 and M1:M2")

# save the seed. It will be restored after this.

    lp_random_effects_list <- lapply(re_terms, function(re){

      re_factor <- parsed_data$reTrms$flist[[re]]

      group_counts <- table(re_factor)

      # set a seed
      if(!is.null(random_effect_seed)) set.seed(random_effect_seed[[re]])

      b <- rnorm(length(group_counts), mean = 0, sd = sqrt(random_effect_variance[[re]]))

      b_rep <- dplyr::left_join(data.frame(re_name = re_factor),
                                data.frame(re_name = names(group_counts), re = b), by = "re_name")

      # This is the random effect repeated for each observation
      re <- b_rep$re

      # this is the unique random effects
      names(b) <- names(group_counts)
      attr(re, "random_effects") <- b

      re

    })

    lp_random_effects <- Reduce("cbind", lp_random_effects_list)

    random_effects <- lapply(lp_random_effects_list, attr, "random_effects")

    # reset randomness
    set.seed(arbitraty_seed)

  } else {
    random_effects <- NULL
  }



  # transpose non-zero length coefficients
  if(length(coefficients))
    coefficients <- t(coefficients)

  risk_score_parts_list <- list(fixed = tcrossprod(as.matrix(vars_df[, fixed_term_variables]), coefficients),
                                random = lp_random_effects)

  risk_score_parts_list <- risk_score_parts_list[lapply(risk_score_parts_list, length)>0]

  risk_score_parts_matrix <- Reduce(cbind, risk_score_parts_list)

  if(!is.null(risk_score_parts_matrix)) {
    risk_score <- .rowSums(risk_score_parts_matrix,
                           m = nrow(risk_score_parts_matrix), n = ncol(risk_score_parts_matrix))
  } else {
    risk_score <- rep(0, length(error))
    vars_df <- data.frame(stat_time = numeric(length(error)))
  }

  vars_df$stat_time <- exp(-risk_score) * error


  # add random effects as an attribute

  attr(vars_df, "random_effects") <- random_effects

  vars_df

}
