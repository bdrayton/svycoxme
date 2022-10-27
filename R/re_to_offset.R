


# function takes a data set, and coxme object that used that data and returns a
# numeric vector of length nrow(data) containing the sum of random effects for
# an observation.

#' turn a list of coxme random effects into a single offset column.
#' @export

re_to_offset <- function(data, model) {

  data_og_names <- names(data)

  # Identify the random effects
  f <- formula(model)
  re_terms <- lme4::findbars(lme4:::RHSForm(f))
  re_terms_nobar <- lapply(re_terms, get_random_bit)

  # handle interactions
  re_vars <- unique(unlist(lapply(re_terms_nobar, split_name)))

  # get dataset with just unique id and random effect variables

  unique_id_name <- paste0("id_", paste0(sample(letters, 10), collapse = ""))

  data[,unique_id_name] <- seq(nrow(the_data))

  re_data <- data[, c(unique_id_name, re_vars)]

  random_effects_list <- coxme::ranef(model)

  re_dfs <- mapply(process_re, re = random_effects_list, re_name = names(random_effects_list),
         MoreArgs = list(re_vars = re_vars,
                         re_data = re_data,
                         id_col = unique_id_name),
         SIMPLIFY = FALSE)

  merged.data.frame = suppressWarnings(Reduce(function(...) merge(..., all=TRUE, by = unique_id_name, sort = FALSE), re_dfs))

  # drop the id column, which is always the first column

  offset = rowSums(merged.data.frame[,-1])
  offset

}


process_re <- function(re, re_name, re_vars, re_data, id_col) {

  #process each random effect using re_vars, the re name, and the re_data

  is_interaction <- !(re_name %in% re_vars)

  if(is_interaction){

    # need to separate names
    sep <- gsub(paste0(re_vars, collapse = "|") , "", re_name)

    # need to know order of variables in the name.
    name_pos <- unlist(lapply(re_vars, gregexpr, re_name))

    # filter re_vars involved in this re term
    re_vars_this_term <- re_vars[name_pos != -1]

    name_pos_this_term <- name_pos[name_pos != -1]

    name_order <- order(name_pos_this_term)

    # break the value names into columns of re levels
    levels <- Reduce(rbind, strsplit(names(re), split = sep, fixed = TRUE)) |> data.frame()

    # label with the re variable name
    names(levels) <- re_vars_this_term[name_order]

    r_data <- merge(re_data, data.frame(levels, value = re),
                    by = re_vars_this_term, all.x = TRUE, sort = FALSE)

  } else {

    r_data <- data.frame(value = re)
    r_data[re_name] <- names(re)

    r_data <- merge(re_data, r_data, by = re_name, all.x = TRUE, sort = FALSE)

  }

  r_data[, c(id_col, "value")]

}



# extract random variables from the formula

get_random_bit <- function(fbit){
  # check it is in shape we expect
  if(is.call(fbit)){
    if(!all(fbit[[1]] == "|", fbit[[2]] == 1))
      stop("Formula error: Expected a random intercept term of the form (1 | x)")

    fbit[[3]]

  } else {
    stop("Formula error: Not a formula")
  }
}

split_name <- function(one_term){

  if(is.call(one_term)) return(as.character(one_term)[-1])

  if(is.name(one_term)) one_term <- as.character(one_term)

  if(grepl("/", one_term)) {

    sep = gsub(".*/).*", "\\1", one_term)
    vars = strsplit(one_term, split = sep)
    return(unlist(vars))

  }

  return(one_term)

}










