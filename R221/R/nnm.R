
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
  lr <- function(x, y) {
    if (all(y==1)) { return( c(0, -Inf) ) }
    lfit <- glm(y ~ x, data=data.frame(x=x, y=y),
                family=binomial(link="logit"))
    lfitCoef <- as.matrix(coef(summary(lfit)))
    c(-rnorm(1, lfitCoef[1,1], lfitCoef[1,2]),
      -rnorm(1, lfitCoef[2,1], lfitCoef[2,2]))
  }

  for (i in 1:ncol(state$data)) {
    x <- state$data[,i]
    y <- state$obs[,i] + 0
    newnu <- lr(x, y)
    state$nu[i] <- newnu[1]
    state$nu[i+4] <- newnu[2]
  }
  state
}

nnmDrawPsi <- function(state) {
  ## Inverse-Wishart
  state$psi <- riwish(nrow(state$L), t(state$L) %*% state$L)
  state
}

nnmDrawTheta <- function(state) {
  ## Linear regression
  for (i in 1:ncol(state$data)) {
    Y <- state$data[, i, drop=FALSE]
    Y <- Y - state$L[, (i+1) %/% 2]
    df <- nrow(state$L)
    scale <- (t(Y) %*% Y) / nrow(state$L)
    state$theta[i,i] <- as.numeric(df * scale / rchisq(1, df))
  }
  state
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
  X <- X2 <- L[,c(1,1,2,2)] + R

  for (i in 1:4) {
    cens.prob <- 1/(1+exp(-(nu[i+4] + nu[i] * X2[,i])))
    X2[,i][runif(length(cens.prob)) < cens.prob] <- NA
  }

  list(data=X, cdata=X2,
       obs=!is.na(X2), L=L, R=R, psi=psi, theta=theta, nu=nu)
}
