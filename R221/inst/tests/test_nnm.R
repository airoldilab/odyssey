
context("nnm")

toydata <- function(nGenes=100) {
  psi <- cbind(c(1,.9), c(.9,1))
  theta <- .1 * diag(4)
  nu <- c(0, 0, 0, 0, -Inf, -Inf, -Inf, -Inf)
  nnmSim(psi, theta, nu, nGenes=nGenes)
}

test_that("nnm in silico generator works", {
  library(R221)
  set.seed(42)
  data <- toydata()$data

  res <- structure(c(0.984998163452428, 0.83564528905883, 0.762651890157275, 
                     0.789790398822123, 0.83564528905883, 0.914469344405947,
                     0.742630110437584, 0.786413476587148, 0.762651890157275,
                     0.742630110437584, 0.955479790749899, 0.902174818849571,
                     0.789790398822123, 0.786413476587148, 0.902174818849571, 
                     1.04620332867223), .Dim = c(4L, 4L))

  expect_that(cov(data), is_equivalent_to(res))
})

test_that("nnm L draws work", {
  library(R221)
  set.seed(42)
  data <- toydata()

  newL <- nnmDrawL(data)$L
  res <- structure(c(0.871834434126793, 0.778325042076193,
                     0.778325042076193, 0.882075457999555),
                   .Dim = c(2L, 2L))
  expect_that(cov(newL), is_equivalent_to(res))
})

test_that("nnm psi draw works", {
  library(R221)
  set.seed(42)
  data <- toydata(nGenes=1000)

  newPsi <- nnmDrawPsi(data)$psi
  res <- structure(c(0.962955449928442, 0.859851093804473,
                     0.859851093804473, 0.959728678885832),
                   .Dim = c(2L, 2L))
  expect_that(newPsi, is_equivalent_to(res))
})
