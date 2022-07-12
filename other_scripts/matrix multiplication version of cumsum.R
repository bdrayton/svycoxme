


# faster than cumsum? matrix multiplication version of cumsums with sparse matricies.

library(Matrix)


make_test_matrix <- function() {

  cells <- sample.int(n = 2, 25, replace = TRUE, prob = c(0.9, 0.1)) - 1

  rnorm(1) * Matrix(cells, nrow = 5)

}


matricies <- replicate(4, make_test_matrix(), simplify = "list")

# convert matricies into one matrix with one per column
big_matrix <- Matrix(data = 0, nrow = 25, ncol = 4)
big_matrix[,1] <- matricies[[1]]
big_matrix[,2] <- matricies[[2]]
big_matrix[,3] <- matricies[[3]]
big_matrix[,4] <- matricies[[4]]



fun1 <- function(){

  big_matrix <- Matrix(data = 0, nrow = 25, ncol = 4)
  big_matrix[,1] <- matricies[[1]]
  big_matrix[,2] <- matricies[[2]]
  big_matrix[,3] <- matricies[[3]]
  big_matrix[,4] <- matricies[[4]]

list(
  Matrix(rowSums(big_matrix), nrow = 5),
  Matrix(rowSums(big_matrix[,2:4]), nrow = 5),
  Matrix(rowSums(big_matrix[,3:4]), nrow = 5),
  Matrix(big_matrix[,4], nrow = 5))
}



fun2 <- function(){
  lapply(seq_along(matricies), function(i) Reduce("+", matricies[i:4]))
}

microbenchmark::microbenchmark(
fun1(),
fun2()
)

fun3 <- function(big_matrix){

  ind_matrix <- Matrix(
                c(rep(c(1, 0), c(4, 0)),
                  rep(c(1, 0), c(3, 1)),
                  rep(c(1, 0), c(2, 2)),
                  rep(c(1, 0), c(1, 3))),
                nrow = 4, sparse = TRUE)

  # this is the cumsums

  big_matrix %*% ind_matrix

}

fun3(big_matrix)

microbenchmark::microbenchmark(fun2(), fun3(big_matrix))

# break up into the matrices

mtemp <- cumsums[, 1, drop = FALSE]

dim(mtemp) <- c(5, 5)













