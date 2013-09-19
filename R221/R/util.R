
## Get fractions of a second of the current time
fracSec <- function() {
  as.numeric(format(Sys.time(), "%OS")) * 10^5
}

## Run all tests of the package
runTests <- function() {
  require(testthat)
  tdir <- system.file("tests", package="R221")
  test_dir(tdir)
}

## Inverse Wishart draw, from package MCMCpack
riwish <- function (v, S) {
  return(solve(rwish(v, solve(S))))
}

## Wishart draw, from package MCMCpack
rwish <- function (v, S) {
  if (!is.matrix(S))
    S <- matrix(S)
  if (nrow(S) != ncol(S)) {
    stop(message = "S not square in rwish().\n")
  }
  if (v < nrow(S)) {
    stop(message = "v is less than the dimension of S in rwish().\n")
  }
  p <- nrow(S)
  CC <- chol(S)
  Z <- matrix(0, p, p)
  diag(Z) <- sqrt(rchisq(p, v:(v - p + 1)))
  if (p > 1) {
    pseq <- 1:(p - 1)
    Z[rep(p * pseq, pseq) + unlist(lapply(pseq, seq))] <-
      rnorm(p * (p - 1)/2)
  }
  return(crossprod(Z %*% CC))
}
