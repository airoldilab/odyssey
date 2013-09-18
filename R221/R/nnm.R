
## state:
## - L
## - data
## - psi
## - theta
## - nu, first four elements are the slope, next 4 the intercept

nnmDrawL <- function(state) {
  ## Conditional normal draw
  TODO
}

nnmDrawMissing <- function(state) {
  ## Imputation, MH step, TODO
  TODO
}

nnmDrawNu <- function(state) {
  ## Logistic regression
  TODO
}

nnmDrawPsi <- function(state) {
  ## Inverse-Wishart
  TODO
}

nnmDrawTheta <- function(state) {
  ## Linear regression
  TODO
}

nnmFit <- function(data, start, steps=50) {
  ## MCMC
  samples <- vector("list", steps)
  samples[[1]] <- state <- start
  for (i in (1:steps)[-1]) {
    state <- nnmDrawL(state)
    state <- nnmDrawMissing(state)
    state <- nnmDrawNu(state)
    state <- nnmDrawPsi(state)
    state <- nnmDrawTheta(state)
    samples[[i]] <- state
  }
  samples
}

nnmStart <- function(data) {
  ## Start state for MCMC
  TODO
}

