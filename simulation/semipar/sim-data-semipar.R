library(MASS)
library(spam)

simulate.data <- function(seed, n,
                          tau.u, tau.z, phi.u, phi.z, rho,
                          beta = c(0, 1, 1), gamma = c(1),
                          misspecified = FALSE, nonnorm = FALSE,
                          VERBOSE = FALSE) {
  
  set.seed(seed)
  
  dexpit <- function(x) {
    1 / (1 + exp(-x)) * 2 - 1
  }
  
  
  CAR.W <- matrix(0, nrow = n, ncol = n)
  CAR.W[cbind(2 : n, 1 : (n - 1))] <- 1
  CAR.W[cbind(1 : (n - 1), 2 : n)] <- 1
  CAR.D <- diag(colSums(CAR.W))
  NB <- list()
  for (i in 1 : nrow(CAR.W)) {
    NB[[i]] <- which(CAR.W[i, ] == 1)
  }
  
  n.short <- 4
  CAR.W.shortcut <- matrix(0, nrow = n.short, ncol = n.short)
  CAR.W.shortcut[cbind(2 : n.short, 1 : (n.short - 1))] <- 1
  CAR.W.shortcut[cbind(1 : (n.short - 1), 2 : n.short)] <- 1
  CAR.W.shortcut[1, n.short] <- 1
  CAR.W.shortcut[n.short, 1] <- 1
  CAR.D.shortcut <- diag(colSums(CAR.W.shortcut))
  
  
  V.from.params <- function(vars, CAR.W, CAR.D) {
    tau.u <- exp(vars[1])
    tau.z <- exp(vars[3])
    tau.e <- exp(vars[5])
    rho.u <- dexpit(vars[2])
    rho.z <- dexpit(vars[4])
    beta.conf <- vars[6]
    
    n <- nrow(CAR.W)
    R.0 <- diag(n)
    
    G.0 <- (CAR.D - rho.u * CAR.W)
    H.0 <- (CAR.D - rho.z * CAR.W)
    
    Q.0 <- - diag(sqrt(diag(G.0) * diag(H.0)))
    
    G.0.inv <- solve(G.0)
    B.0 <- - G.0.inv %*% Q.0
    
    R <- tau.e * R.0
    G <- tau.u * G.0
    H <- tau.z * H.0
    Q <- beta.conf * tau.u * Q.0
    B <- beta.conf * B.0
    A.inv <- tau.z * H.0 - beta.conf ^ 2 * tau.u * t(Q.0) %*% G.0.inv %*% Q.0
    
    V.inv <-  solve(G.0.inv / tau.u + diag(n) / tau.e)
    
    return(list(V.inv = V.inv, A.inv = A.inv, B = B, R = R, G = G, H = H, Q = Q))
  }
  
  
  
  
  #####################
  ### simulate data ###
  #####################
  
  G <- (CAR.D - phi.u * CAR.W) * tau.u
  H <- (CAR.D - phi.z * CAR.W) * tau.z
  
  Q <- -diag(rho * sqrt(diag(G) * diag(H)))
  
  if (VERBOSE) {
    print("solving")
  }
  P <- as.spam(cbind(rbind(G, t(Q)), rbind(Q, H)))
  LP <- chol(P)
  if (VERBOSE) {
    print("cholesky")
  }
  if (VERBOSE) {
    print("done")
  }
  
  x <- runif(n, -1 / 2, 1 / 2)
  
  if (misspecified) {
    u <- backsolve(chol(G), rnorm(n))
    z <- gamma[1] * x + u + rnorm(n)
  } else {
    v <- backsolve(LP, rnorm(2 * n))
    u <- v[1 : n]
    z <- gamma[1] * x + v[n + (1 : n)]
  }
  
  if (nonnorm) {
    y <- rpois(n, exp((2 / (1 + exp(-6 * z)) - 1) + atan(u)))
  } else {
    y <- rpois(n, exp((2 / (1 + exp(-6 * z)) - 1) + u))
  }
  
  data <- data.frame(y = y, x = x, z = z)
  attr(data, "spatial") <- list(NB = NB, CAR.D = CAR.D, CAR.W = CAR.W,
                                CAR.D.shortcut = CAR.D.shortcut,
                                CAR.W.shortcut = CAR.W.shortcut)
  
  attr(data, "truth") <- list(
    u = u,
    tau.u = tau.u,
    tau.z = tau.z,
    phi.u = phi.u,
    phi.z = phi.z,
    rho = rho,
    beta = beta,
    gamma = gamma
  )
  
  return(data)
}

