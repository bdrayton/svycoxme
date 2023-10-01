
# random number management


#' @export
assign_streams <- function(cl, random_streams){

  parallel::clusterApply(cl, x = random_streams, function(x){

    assign("thread_id", x$thread_id, envir = .GlobalEnv)

    assign(".Random.seed", x$seed, envir = .GlobalEnv)

  })

}

#' @export
update_streams <- function(cl, file){

  current_streams <- parallel::clusterApply(cl, x = seq(length(cl)), function(x){
    list(thread_id = thread_id,
         seed = .Random.seed)
  })

  saveRDS(current_streams, file = file)

}

