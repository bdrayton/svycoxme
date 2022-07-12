

# speed test for cumsums

#

### If this is too slow, consider matrixStats::colCumsums(testMatrix)

# matrix with 10 columns, 200 rows

library(Matrix)



fun1 <- function(a_matrix){

  # flip it
  rev_index <- rev(seq_len(nrow(a_matrix)))
  m <- a_matrix[rev_index, ]

  # cumsums
  m <- matrixStats::colCumsums(as.matrix(m))

  # flip it back
  m <- m[rev_index,]

  # restore class
  Matrix(m)
}

fun2 <- function(a_matrix) {
  at_risk_list <- apply(a_matrix, 2, function(column){
    Matrix(rev(cumsum(rev(column))), ncol = 1)
  })

  Reduce(cbind, at_risk_list)
}


test_matrix <- Matrix(data = abs(rnorm(2 * 4)), nrow = 2, ncol = 4)

fun1(test_matrix)

fun2(test_matrix)


test_matrix <- Matrix(data = abs(rnorm(200 * 40)), nrow = 200, ncol = 40)

microbenchmark::microbenchmark(
  fun1(test_matrix),
  fun2(test_matrix)
)





