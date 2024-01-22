
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



# works like parallel::clusterSetRNGStream, but uses a file with seeds
# from independent streams

#' @export

clusterSetRNGStreamFromFile <- function(cl, file) {
  cl <- defaultCluster(cl)

  nc <- length(cl)

  seeds <- readRDS(file)

  if(nc > length(seeds)) stop("Not enough seeds in file")

  nc_seq <- seq_along(cl)

  seeds <- seeds[nc_seq]

  for (i in nc_seq) {
    expr <- substitute(assign(".Random.seed", seed, envir = .GlobalEnv),
                       list(seed = seeds[[i]]))
    sendCall(cl[[i]], eval, list(expr))
  }
  checkForRemoteErrors(lapply(cl, recvResult))
  invisible()
}


# takes the current seeds from a cluster and saves them to a file, only
# overwriting the first length(cl) seeds out of all the seeds on file.

#' @export

updateRNGStreamFileFromCluster <- function(cl, file) {

  cl <- defaultCluster(cl)

  nc_seq <- seq_along(cl)

  old_seeds <- readRDS(file)

  current_seeds <- parallel::clusterApply(cl, x = nc_seq, function(x){
    return(.Random.seed)
  })

  # checkForRemoteErrors(lapply(cl, recvResult))

  old_seeds[nc_seq] <- current_seeds

  saveRDS(old_seeds, file = file)

}


