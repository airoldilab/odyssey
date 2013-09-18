
context("nnm")

test_that("nnm in silico generator works", {
  library(R221)
  set.seed(42)
  psi <- cbind(c(1,.9), c(.9,1))
  theta <- .1 * diag(4)
  nu <- c(0, 0, 0, 0, -Inf, -Inf, -Inf, -Inf)
  data <- nnmSim(psi, theta, nu, nGenes=100)

  res <- structure(c(0.984998163452428, 0.83564528905883, 0.762651890157275, 
                     0.789790398822123, 0.83564528905883, 0.914469344405947,
                     0.742630110437584, 0.786413476587148, 0.762651890157275,
                     0.742630110437584, 0.955479790749899, 0.902174818849571,
                     0.789790398822123, 0.786413476587148, 0.902174818849571, 
                     1.04620332867223), .Dim = c(4L, 4L))

  expect_that(cov(data), is_equivalent_to(res))
})

