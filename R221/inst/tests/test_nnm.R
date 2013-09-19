
context("nnm")

toydata <- function(nGenes=100) {
  psi <- cbind(c(1,.9), c(.9,1))
  theta <- .1 * diag(4)
  nu <- c(0, 0, 0, 0, -Inf, -Inf, -Inf, -Inf)
  nnmSim(psi, theta, nu, nGenes=nGenes)
}

toydataMiss <- function(nGenes=100) {
  psi <- cbind(c(1,.9), c(.9,1))
  theta <- .1 * diag(4)
  nu <- c(-2, -2, -1, 0, -2, -3, -4, -Inf)
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

test_that("nnm nu draws work", {
  library(R221)
  set.seed(42)
  data <- toydataMiss(nGenes=1000)

  newNu <- nnmDrawNu(data)$nu
  res <- c(-2.10819031907204, -2.76605673346096, -3.6214431722642, 0,
           -1.97858084563688, -2.07328990552117, -1.24537573857442, -Inf)

  expect_that(newNu, is_equivalent_to(res))
})

test_that("nnm theta draws work", {
  library(R221)
  set.seed(42)
  data <- toydata()

  newTheta <- nnmDrawTheta(data)$theta
  res <- structure(c(0.105223780533292, 0, 0, 0, 0, 0.0885303890153848,
                     0, 0, 0, 0, 0.0925352784573173, 0, 0, 0, 0,
                     0.0923744659251901), .Dim = c(4L, 4L))
  expect_that(newTheta, is_equivalent_to(res))
})
