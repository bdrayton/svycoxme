

# #' Check formula against dists names
# #'
# #' Checks that the variables in the right hand side of the formula all have matching
# #' expressions is the dists list.
# #'
# #' @param var_labels Names of the variables in the right hand side of the formula
# #' @param dist_labels Names of the expressions in dists
# #'
# #' @export
#
#
# check_var_labels <- function(var_labels, dist_labels){
#
#   test <- var_labels %in% dist_labels
#
#   not_in_dist_labels <- paste0(var_labels[!test], collapse = ", ")
#
#   if(!all(test))
#     stop(paste0("All variables in formula must be in dists.\n  These variables are in your formula but not in dists: ", not_in_dist_labels))
#
# }



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

one_dataset <- function(formula, dists, dist_args, error, censoring_time, coefficients = c(), random_effect_variance = list(), seed = NULL,
                        random_effect_seed = NULL){

  # arbitraty_seed <- runif(1, -99999, 99999)

  # if(!is.null(random_effect_seed) & (length(random_effect_variance) != length(random_effect_seed)))
  #   stop("If random_effect_seed is not NULL, it must be the same length as random_effect_variance")

  original_formula <- formula

  # if(!is.null(seed)) set.seed(seed)

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

  # add in sequential evaluation of the dist arguments, and an evironment to put the
  # resulting variables. allows subsequent dists to refererence previously defined variables.
  # vars <- lapply(dists[variables], lazyeval::f_eval, data = dist_args)
  #
  # vars_df <- list2DF(vars)

  var_env <- new.env()

  invisible(lapply(names(dist_args), function(arg){
    assign(arg, dist_args[[arg]], pos = var_env)
  }))

  # need variables in the order or dist names.
  # easier to just pull dist names.
  dist_names <- names(dists)

  for (i in seq(length(dist_names))) {

    assign(dist_names[i],
           do.call("eval", list(expr = dists[[i]][[2]], envir = var_env), envir = var_env),
           envir = var_env)

  }

  vars_df <- eval(parse(text = paste("data.frame(", paste(dist_names, collapse = ","), ")")), envir = var_env)

  error <- lazyeval::f_eval(error, data = dist_args)

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
      # if(!is.null(random_effect_seed)) set.seed(random_effect_seed[[re]])

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
    # set.seed(arbitraty_seed)

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
    vars_df <- data.frame(event_time = numeric(length(error)))
  }

  vars_df$event_time <- exp(-risk_score) * error

  vars_df$censoring_time <- lazyeval::f_eval(censoring_time, data = dist_args)

  # add stat_time and stat, based on event and censoring times.
  vars_df$stat = vars_df$event_time <= vars_df$censoring_time

  vars_df$stat_time = with(vars_df, ifelse(stat, event_time, censoring_time))

  # add random effects as variables
  # for each item in lp_random_effects_list, join it to vars_df.

  for (i in seq_along(lp_random_effects_list)){

    # much simpler
    var_name <- re_terms[i]

    vars_df[, paste0(var_name, "_re")] <- lp_random_effects_list[[i]]
    #
    #
    # join_by <- var_name
    # var_re <- lp_random_effects_list[[i]]
    # var_levels <- names(var_re)
    #
    # re_temp <- setNames(data.frame(var_levels, var_re), c(var_name, paste(var_name, "ranef", sep = "_")))
    #
    # # any names like M1:M2, M1.M2, M1/M2 need to be split
    # name_get_split <- grepl("\\.|:|/", var_name)
    #
    # if (name_get_split) {
    #
    #   sep <- gsub(".*([\\.|:|/]).*", "\\1")
    #
    #   new_cols <- unlist(strsplit(var_name, sep, fixed = TRUE))
    #   join_by <- var_name
    #
    #   # don't need sep it seems. Good because '.' would need to be masked
    #   re_temp <- tidyr::separate(re_temp,
    #                              col = dplyr::all_of(var_name),
    #                              into = new_cols)

    # }
    #
    # vars_df <- dplyr::left_join(vars_df, by = join_by)

  }

  # add random effects as an attribute

  attr(vars_df, "random_effects") <- random_effects

  vars_df

}



#' Create a dataset of single or recurrent events for subjects with possibly time-varying covariate data.
#'
#' @param formula Model equation for the data generation.
#' @param data Data frame containing the covariates in the \code{formula}, an \code{id} column identifying subjects, and start and end times for time-varying covariates.
#' @param coefficients Named vector of fixed effect coefficients.
#' @param random_effect_variance Named vector of random effects. Names must match random effect formula terms.
#' @param id The data column that identifies observations belonging to one subject.
#' @param baseline_hazard scalar or vector with stepwise hazards over time
#' @param t change points for \code{baseline_hazard}
#' @param event option for a single event or recurrent events per subject
#' @param origin when the time scale starts
#' @param max_events for recurrent events, maximum events per subject. See details.
#'
#' @details Covariate data is passed to \code{draw_event_times} via the \code{data} parameter, which
#' accepts a data.frame. If data is time varying, covariates will need start and end times. In either single
#' or recurrent event generation, providing an end time to \code{data} will provide an upper limit on event times, that is,
#' administrative censoring at the end time provided. The default, when an end time is not provided, is to set it to \code{Inf}. For recurrent events, this results in \code{max_events} being generated. Known bugs: max_events is being applied per row, not per subject.
#'
#'
#' @export


draw_event_times <- function(formula, data, coefficients = c(), random_effect_variance = list(), id,
                             baseline_hazard = 1, t = 0, event = c("single", "recurrent"),
                             origin = 0, max_events = 1000) {

  Call <- match.call()
  event <- match.arg(event)

  Call_id <- Call[["id"]]
  if (is.null(Call_id)){
    stop("id required")
  } else if (is.name(Call_id)) {
    idx <- as.character(Call_id)
  }

  fixed_terms <- stats::terms(lme4::nobars(formula))
  # some sort of checks for formula and coefficients. match names? length?
  # trying for informative warning, because if they don't align, a cryptic
  # error is thrown later: Error in base::tcrossprod(x, y) : non-conformable arguments
  # if(length(fixed_terms) != length(coefficients())) {
  #   stop("Coefficients and fixed terms don't match ")
  # }

  # These are the fixed term variables in the formula, even if there are interactions.
  fixed_term_variables <- colnames(attr(fixed_terms, "factors"))

  bars = lme4::findbars(formula)

  barnames = lme4:::barnames(bars)

  random_terms = terms(as.formula(paste("~", paste(c(1, barnames), collapse = "+"))))

  random_term_variables <- rownames(attr(random_terms, "factors"))

  variables <- c(fixed_term_variables, random_term_variables)

  if(attr(fixed_terms, "response") == 0) stop("formula must have a Surv response")

  # strata
  # find the strata spec in the form.
  # make sure there is a baseline hazard for each stratum

  # clustering. dont allow a cluster() call in the formula. instead ask for (1|cluster_id).
  # Or allow it and convert it.
  # will need an associated random effect.
  # cluster in the sandwich means that the variances are cluster-specific. need to allow input of variance structure?

  # data will have start and stop times, as specified, but only if there is time-varying covariates.
  # otherwise there wont be times.
  # and there wont be a stat variable.

  Surv_Call <- Call$formula[[2]]
  # get time names out of the Call.
  # cant use the attribute, because the data need the variables in the call before calling model.frame.
  # but I'll use the length of the Surv_call instead. This will probably break if anything but one or two times
  # are passed to Surv. Maybe I need to evaluate Surv? Nope. need the vars first.
  # if a stop time has not been specified, set it to Inf.

  # Surv_evaluated <- eval(Surv_Call)

  # if(attr(model_response, "type") == "counting"){
  if(length(Surv_Call) == 4){
    start_name <- Surv_Call[[2]]
    stop_name <- Surv_Call[[3]]
    status_name <- Surv_Call[[4]]

    Surv_vars <- c(start_name, stop_name, status_name)

    # } else if(attr(model_response, "type") == "right"){
  } else if(length(Surv_Call) == 3) {
    start_name <- "origin"
    stop_name <- Surv_Call[[2]]    # actually called time, not stop
    status_name <- Surv_Call[[3]]

    if(start_name == stop_name) {stop("You need to call your event times something other than, 'origin'")}

    Surv_vars <- c(start_name, stop_name, status_name)

  }

  # convert Surv_vars from a list of names to a character vector
  Surv_vars <- sapply(Surv_vars, as.character)

  # check if variables in the Surv() Call are already in the data set.
  Surv_vars_in_data <- Surv_vars %in% colnames(data)
  # add the ones that aren't
  # need a message or warning? Probably not. Just explain in documentation
  # message(paste(Surv_vars[Surv_vars_in_data], collapse = ", "), paste(" will be modified"))

  if(any(!Surv_vars_in_data)){

    data[Surv_vars[!Surv_vars_in_data]] <- double(length = nrow(data))

    # If stop time was added to the data, set to Inf.
    if(!Surv_vars_in_data[2]){
      data[[stop_name]] <- Inf
    }

  }

  model_frame <- model.frame(lme4:::getFixedFormula(formula), data)
  model_response <- model.response(model_frame)
  X <- model_frame[, -1, drop = FALSE]

  lp_random_effects <- list()

  if(length(barnames)) {

    data$ytemp <- 1

    formula <- update(formula, ytemp ~ .)

    # control options stop lme4 complaining when there is only one obs per group.
    # this is probably only a problem when making tiny test data sets.
    # additionally, I might suppress the warning,
    # 'fixed-effect model matrix is rank deficient so dropping 1 column / coefficient'
    # it's not relevant here, because rank deficiency is a problem for estimation,
    # not data generation. And I already have that matrix.
    parsed_data <- lme4::lFormula(formula, data = data,
                                  control = lme4::lmerControl(check.nobs.vs.nlev = "ignore",
                                                              check.nobs.vs.rankZ = "ignore",
                                                              check.nobs.vs.nRE = "ignore"))

    data <- dplyr::select(data, -ytemp)

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
      # if(!is.null(random_effect_seed)) set.seed(random_effect_seed[[re]])

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
    # set.seed(arbitraty_seed)

  } else {
    random_effects <- NULL
  }


  # transpose non-zero length coefficients
  if(length(coefficients))
    coefficients <- t(coefficients)

  risk_score_parts_list <- list(fixed = tcrossprod(as.matrix(data[, fixed_term_variables]), coefficients),
                                random = lp_random_effects)

  risk_score_parts_list <- risk_score_parts_list[lapply(risk_score_parts_list, length)>0]

  risk_score_parts_matrix <- Reduce(cbind, risk_score_parts_list)

  if(!is.null(risk_score_parts_matrix)) {
    risk_score <- .rowSums(risk_score_parts_matrix,
                           m = nrow(risk_score_parts_matrix), n = ncol(risk_score_parts_matrix))
  } else {
    risk_score <- rep(0, length(error))
    data <- data.frame(stat_time = numeric(length(error)))
  }

  d_list <- C_draw_event_times(id = as.integer(data[[idx]]),
                               start_time = as.double(data[[Surv_vars[1]]]),
                               end_time = as.double(data[[Surv_vars[2]]]),
                               status = as.integer(data[[Surv_vars[3]]]),
                               X = as.matrix(X),
                               risk_score = risk_score,
                               baseline_hazard = baseline_hazard,
                               baseline_hazard_start = t,
                               origin = origin,
                               single = as.integer(event == "single"))

  data_with_events <- data.frame(d_list[[1]],
                                 d_list[[5]])

  names(data_with_events) <- c(idx, colnames(X))

  data_with_events[Surv_vars[1]] <- d_list[[2]]
  data_with_events[Surv_vars[2]] <- d_list[[3]]
  data_with_events[Surv_vars[3]] <- d_list[[4]]

  # lets see where we're at at this point.
  return(data_with_events)

}


