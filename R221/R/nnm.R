
## state:
## - L
## - data
## - psi
## - theta
## - nu, first four elements are the slope, next 4 the intercept

nnmDrawL <- function(state) {
  ## Conditional normal draw
  B <- cbind(c(1,0), c(1,0), c(0,1), c(0,1))
  V.X <- t(B) %*% state$psi %*% B + state$theta
  C.L.X <- state$psi %*% B
  V.LX <- rbind(cbind(state$psi, C.L.X), cbind(t(C.L.X), V.X))

  for (i in 1:2) {
    S11 <- V.LX[i, i, drop=FALSE]
    S12 <- V.LX[i, 3:6, drop=FALSE]
    S21 <- t(S12)
    S22 <- V.LX[3:6, 3:6]
    mu <- S12 %*% solve(S22) %*% t(state$data)
    sigma <- S11 - S12 %*% solve(S22) %*% S21
    state$L[,i] <- rnorm(length(mu), mean=mu, sd=sqrt(sigma))
  }
  state
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
  state$psi <- riwish(nrow(state$L), t(state$L) %*% state$L)
  state
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

nnmSim <- function(psi, theta, nu, nGenes=5000) {
  ## Simulate artificial data
  L <- rmvnorm(nGenes, sigma=psi)
  R <- rmvnorm(nGenes, sigma=theta)
  X <- L[,c(1,1,2,2)] + R

  for (i in 1:4) {
    cens.prob <- 1/(1+exp(-(nu[i+4] + nu[i] * X[,i])))
    X[,i][runif(length(cens.prob)) < cens.prob] <- NA
  }

  list(data=X, L=L, R=R, psi=psi, theta=theta, nu=nu)
}
