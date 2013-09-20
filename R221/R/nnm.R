
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

invlogit <- function(x, nu0, nu1) {
  1/(1+exp(-(nu0+nu1*x)))
}

targetDens <- function(x, nu0, nu1, mu, sig2) {
  invlogit(x, nu0, nu1) * dnorm(x, mu, sqrt(sig2))
}

targetHessian <- function(x, nu0, nu1, sig2) {
  -nu1^2/(exp(nu0+nu1*x)+1) + nu1^2/(1+exp(nu0+nu1*x))^2 - 1/sig2
}

nnmDrawMissing <- function(state) {
  ## Imputation, MH step
  B <- cbind(c(1,0), c(1,0), c(0,1), c(0,1))
  V.X <- t(B) %*% state$psi %*% B + state$theta
  C.L.X <- state$psi %*% B
  V.LX <- rbind(cbind(state$psi, C.L.X), cbind(t(C.L.X), V.X))

  missIdx <- which(!state$obs, arr.ind=TRUE)
  for (i in seq_len(nrow(missIdx))) {
    mr <- missIdx[i, 1]
    mc <- missIdx[i, 2]
    nu0 <- state$nu[mc+4]
    nu1 <- state$nu[mc]

    S11 <- V.LX[mc+2, mc+2, drop=FALSE]
    S12 <- V.LX[mc+2, -(mc+2), drop=FALSE]
    S21 <- t(S12)
    S22 <- V.LX[-(mc+2), -(mc+2), drop=FALSE]

    muCond <- S12 %*% solve(S22) %*% c(state$L[mr,], state$data[mr,-mc])
    sig2Cond <- S11 - S12 %*% solve(S22) %*% S21

    sdMu <- sqrt(-1 / targetHessian(muCond, nu0, nu1, sig2Cond))

    xMax <- optimize(f=function(x) -targetDens(x, nu0, nu1, muCond,
                       sig2Cond),
                     interval=c(muCond-sdMu*3, muCond+sdMu*3))$minimum
    V <- -1 / targetHessian(xMax, nu0, nu1, sig2Cond)

    proposal <- rnorm(1, mean=xMax, sd=sqrt(V))

    r1 <- targetDens(proposal, nu0, nu1, muCond, sig2Cond) /
      targetDens(state$data[mr, mc], nu0, nu1, muCond, sig2Cond)
    r2 <- dnorm(state$data[mr, mc], mean=xMax, sd=sqrt(V)) /
      dnorm(proposal, mean=xMax, sd=sqrt(V))
    r <- r1 * r2
    if (r>=1 || runif(1) < r) {
      state$data[mr, mc] <- proposal
    }
  }
  state
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
  state <- list()

  state$obs  <- !is.na(data)
  state$data <- data
  state$data[is.na(state$data)] <- 0

  state$L <- cbind(rowMeans(state$data[,1:2]),
                   rowMeans(state$data[,3:4]))
  state$psi <- cov(state$L)
  state$theta <- 0.1 * diag(4)

  if (any(is.na(data))) {
    state$nu <- nnmDrawNu(state)$nu
  } else {
    state$nu <- c(0, 0, 0, 0, -Inf, -Inf, -Inf, -Inf)
  }

  state
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

nnmRun <- function(data, start, noSamples, thin=50, outputdir) {
  ## Run a simuation, save samples to files
  for (i in 1:noSamples) {
    samples <- nnmFit(data, start, steps=thin)
    start <- sample <- samples[[thin]]
    save(sample, file=sprintf("%s/sample-%i.Rdata", outputdir, i))
  }
  invisible(TRUE)
}

##############################################################
# JOBS
##############################################################

nnmRunInsilico <- function(ID) {
  ## Get my ID
  if (missing(ID)) {
    args <- commandArgs(TRUE)
    ID <- as.numeric(args[1])
  }
  outdir <- sprintf("out-%i", ID)
  dir.create(outdir)

  pars <- seq(.1, .9, by=.1)
  psi12 <- pars[ID]

  ## Random seed
  seed <- fracSec()
  set.seed(seed)

  ## Session information
  session <- sessionInfo()

  ## In-silico data
  psi <- cbind(c(1, psi12), c(psi12, 1))
  theta <- .1 * diag(4)
  nu <- c(0, 0, 0, 0, -Inf, -Inf, -Inf, -Inf)
  truth <- nnmSim(psi, theta, nu, nGenes=1000)
  start <- nnmStart(truth$cdata)

  save(truth, start, seed, session, file=sprintf("%s/simdata.Rdata", outdir))
  nnmRun(data$cdata, start, noSamples=100, thin=50, outputdir=outdir)
}
