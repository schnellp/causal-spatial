

library(arm)
library(spam)
library(Matrix)
library(lme4)

fit.affine <- function(y, z, Xz, Xy,
                       offset = rep(0, length(y)),
                       NB, CAR.D, CAR.W,
                       seed = 0,
                       init = list(),
                       print.iters = 100,
                       burn.in = 2000,
                       n.iters = 5000,
                       CONSTRAIN_SCALES = TRUE,
                       ESTIMATE_RHO = TRUE,
                       SPATIAL = TRUE,
                       USE_PRIOR = FALSE,
                       COND_PRIOR = FALSE,
                       CAR.D.shortcut = NULL, CAR.W.shortcut = NULL,
                       debug.control = list(),
                       VERBOSE = FALSE,
                       coeff.prior = 1 / 100,
                       tau.prior = 5,
                       rho.prior = 1,
                       prop.scale = 1 / 4) {
  
  start.time = Sys.time()
  
  missing.y <- which(is.na(y))
  
  if (is.null(Xz)) {
    Xz <- matrix(0, nrow = length(y), ncol = 1)
  }
  
  # VERBOSE <- FALSE
  # CONSTRAIN_SCALES <- TRUE
  # ESTIMATE_RHO <- TRUE
  # SPATIAL <- TRUE
  # USE_PRIOR <- TRUE
  
  if (length(coeff.prior) == 1) {
    coeff.prior <- rep(coeff.prior, ncol(Xy))
  }
  
  spline.penalty <- 1
  
  if (!is.null(debug.control$ESTIMATE_U)) {
    ESTIMATE_U <- debug.control$ESTIMATE_U
  } else {
    ESTIMATE_U <- TRUE
  }
  
  
  
  set.seed(seed)
  
  expit <- function(x) {
    1 / (1 + exp(-x))
  }
  
  logit <- function(p) {
    log(p) - log(1 - p)
  }
  
  dexpit <- function(x) {
    1 / (1 + exp(-x)) * 2 - 1
  }
  
  dlogit <- function(p) {
    log(1 + p) - log(1 - p)
  }
  
  ldmvt <- function(x, df, mu = rep(0, length(x)), P = diag(length(x))) {
    - (df + length(x)) / 2 * log1p(t(x - mu) %*% P %*% (x - mu) / df)
  }
  
  pirls.poisson.u <- function(y, offset, prec.0, mu.0,
                              n.iters = 10, init = 0) {
    
    
    u <- numeric(n.iters)
    u[1] <- log(y + 0.5) - offset
    
    for (i in 2 : n.iters) {
      eta <- u[i - 1] + offset
      mu <- exp(eta)
      W <- mu
      z <- eta - offset + (y - mu) / W
      u[i] <- 1 / (W + prec.0) * (prec.0 * mu.0 + W * z)
    }
    
    mu.post <- u[n.iters]
    P.post <- W + prec.0
    
    return(list(mu = mu.post, prec = P.post))
  }
  
  ldfu <- function(u, y, offset, prec, mu) {
    y * (offset + u) - exp(offset + u) - prec * (u - mu) ^ 2 / 2
  }
  
  ldfn <- function(u, prec, mu) {
    dnorm(u, mu, sqrt(2 / prec), log = TRUE)
  }
  
  pirls.poisson <- function(Xy, offset, y, P.0, mu.0,
                            n.iters = 10, init = rep(0, ncol(Xy)),
                            max.w = 1e8, max.eta.abs = 1e4) {
    g <- function(mu) {
      log(mu)
    }
    
    g.inv <- function(eta) {
      exp(eta)
    }
    
    a <- function(mu) {
      mu
    }
    
    beta <- matrix(NA, nrow = n.iters, ncol = ncol(Xy))
    
    beta[1, ] <- init
    
    
    for (i in 2 : n.iters) {
      eta <- Xy %*% beta[i - 1, ] + offset
      mu <- g.inv(eta)
      W <- diag(as.vector(pmin(a(mu), max.w)))
      z <- eta - offset + (y - mu) / diag(W)
      beta[i, ] <- as.vector(solve(t(Xy) %*% W %*% Xy + P.0) %*%
                               (P.0 %*% mu.0 + t(Xy) %*% W %*% z))
    }
    
    mu.post <- beta[n.iters, ]
    P.post <- t(Xy) %*% W %*% Xy + P.0
    
    return(list(mu = mu.post, P = P.post))
  }
  
  CAR.from.par <- function(tau, phi, CAR.D, CAR.W) {
    # (tau * (CAR.D - phi * CAR.W))
    tau * as.spam(CAR.D - phi * CAR.W)
  }
  
  Q.from.par <- function(rho, G, H) {
    # (diag(-rho * sqrt(diag(G) * diag(H))))
    -rho * diag.spam(sqrt(diag(G) * diag(H)))
  }
  
  chol.ev.approx <- function(A, kmax = 10) {
    J <- t(A) %*% A
    
    for (k in 1 : kmax) {
      R <- as.spam(chol(J))
      J <- R %*% t(R)
    }
    
    sqrt(diag(J))
  }
  
  log.prior.dep <- function(tau.u, tau.z, phi.u, phi.z, rho,
                            e.val) {
    if (CONSTRAIN_SCALES && phi.z > phi.u) {
      return(-Inf)
    }
    
    if (abs(rho) >= 0.5) {
      return(-Inf)
    }
    
    if (USE_PRIOR) {
      
      if (COND_PRIOR) {
        if (any(e.val <= 0)) {
          return(-Inf)
        } else {
          kappa <- max(e.val) / min(e.val)
          
          return(dexp(kappa, rate = 1 / 10, log = TRUE))
        }
      } else {
        return(sum(dgamma(c(tau.u, tau.z), tau.prior, tau.prior, log = TRUE)))
      }
      
      
      
    } else {
      return(0)
    }
  }
  
  log.prior.coeff <- function(coeff, coeff.prior) {
    sum(dnorm(coeff, 0, sqrt(1 / coeff.prior), log = TRUE))
  }
  
  log.lik.y <- function(beta, u, Xy, offset, y, sum = TRUE) {
    eta <- Xy %*% beta + u + offset
    ll <- dpois(y, exp(eta), log = TRUE)
    if (sum) {
      return(sum(ll))
    } else {
      return(ll)
    }
  }
  
  full.cond.dep <- function(tau.u, tau.z, phi.u, phi.z, rho,
                            G, H, Q,
                            u, z.res,
                            P, cP, e.val) {
    
    
    if (is.null(cP)) {
      return(-Inf)
    }
    
    (1 / 2) * (
      2 * determinant(cP)$modulus
    ) -
      (1 / 2) * drop(as.matrix( # all mult by tau.u?
        t(u) %*% G %*% u +
          2 * t(u) %*% Q %*% z.res +
          t(z.res) %*% H %*% z.res
      )) +
      log.prior.dep(tau.u = tau.u, tau.z = tau.z,
                    phi.u = phi.u, phi.z = phi.z,
                    rho = rho,
                    e.val)
  }
  
  mh.accept.log.prob <- function(log.dens.pro, log.dens.cur,
                                 log.prop.pro = 0, log.prop.cur = 0) {
    (log.dens.pro - log.prop.pro) - (log.dens.cur - log.prop.cur)
  }
  
  mh.accept <- function(log.dens.pro, log.dens.cur,
                        log.prop.pro = 0, log.prop.cur = 0) {
    as.logical(rbinom(1, 1,
                      min(1, exp(mh.accept.log.prob(
                        log.dens.pro = log.dens.pro,
                        log.dens.cur = log.dens.cur,
                        log.prop.pro = log.prop.pro,
                        log.prop.cur = log.prop.cur
                      )))))
  }
  
  post.par.gaussian <- function(prior.mu, prior.V, lik.mu, lik.V) {
    prior.V.inv <- solve(prior.V)
    lik.V.inv <- solve(lik.V) # ginv?
    post.V <- solve(prior.V.inv + lik.V.inv)
    post.mu <- as.vector(post.V %*% (prior.V.inv %*% prior.mu + lik.V.inv %*% lik.mu))
    
    return(list(mu = post.mu, V = post.V))
  }
  
  post.par.gaussian.prec <- function(prior.mu, prior.prec, lik.mu, lik.prec) {
    post.Prec <- prior.prec + lik.prec
    post.mu <- as.vector(solve(post.Prec) %*% (prior.prec %*% prior.mu + lik.prec %*% lik.mu))
    
    return(list(mu = post.mu, Prec = post.Prec))
  }
  
  
  ### config ###
  
  total.iters <- burn.in + n.iters
  
  ### init ###
  
  n <- nrow(Xy)
  py <- ncol(Xy)
  pz <- ncol(Xz)
  
  y[missing.y] <- 0
  
  if (is.null(init[["beta"]])) {
    coeff.prior[penalized] <- spline.penalty
    beta <- coefficients(bayesglm(y ~ 0 + offset(offset) + Xy,
                                  family = "poisson",
                                  prior.scale = c(1 / coeff.prior),
                                  prior.df = Inf,
                                  scaled = FALSE))
    for (i in missing.y) {
      log.rate.y <- offset[i] + Xy[i, ] %*% beta
      rate.y <- exp(log.rate.y)
      
      y.poss <- 0 : 9
      log.probs <- dpois(y.poss, rate.y, log = TRUE)
      y[i] <- sample(y.poss, 1, prob = exp(log.probs - max(log.probs)))
    }
    beta <- coefficients(bayesglm(y ~ 0 + offset(offset) + Xy,
                                  family = "poisson",
                                  prior.scale = c(1 / coeff.prior),
                                  prior.df = Inf,
                                  scaled = FALSE))
  } else {
    beta <- init[["beta"]]
  }
  
  if (is.null(init[["gamma"]])) {
    if (all(Xz == 0)) {
      gamma <- rep(0, ncol(Xz))
    } else {
      gamma <- coefficients(lm(z ~ 0 + Xz))
    }
    
  } else {
    gamma <- init[["gamma"]]
  }
  
  if (is.null(init[["tau.u"]])) {
    tau.u <- 1
  } else {
    tau.u <- init[["tau.u"]]
  }
  
  if (is.null(init[["tau.z"]])) {
    tau.z <- 1
  } else {
    tau.z <- init[["tau.z"]]
  }
  
  if (is.null(init[["phi.u"]])) {
    phi.u <- 0
  } else {
    phi.u <- init[["phi.u"]]
  }
  
  if (is.null(init[["phi.z"]])) {
    phi.z <- 0
  } else {
    phi.z <- init[["phi.z"]]
  }
  
  if (is.null(init[["rho"]])) {
    rho <- 0
  } else {
    rho <- init[["rho"]]
  }
  
  if (SPATIAL) {
    if (is.null(init[["u"]])) {
      # u <- log(y + 0.5)
      u <- rep(0, length(y))
    } else {
      u <- init[["u"]]
    }
  } else {
    u <- rep(0, length(y))
  }
  
  
  
  accept.u <- 0
  accept.z <- 0
  accept.prec <- 0
  accept.beta <- 0
  
  CAR.D <- as.spam(CAR.D)
  CAR.W <- as.spam(CAR.W)
  
  
  G <- CAR.from.par(tau = tau.u, phi = phi.u, CAR.D = CAR.D, CAR.W = CAR.W)
  G.inv <- solve(G)
  H <- CAR.from.par(tau = tau.z, phi = phi.z, CAR.D = CAR.D, CAR.W = CAR.W)
  H.inv <- solve(H)
  Q <- Q.from.par(rho = rho, G = G, H = H)
  
  P <- as.spam(rbind(cbind(G, Q),
                     cbind(Q, H)))
  
  cP <- tryCatch({
    suppressWarnings(chol(P))
  }, error = function(e) {
    NULL
  })
  
  # e.val <- range(chol.ev.approx(P, kmax = 10))
  if (COND_PRIOR) {
    G.pro.short <- CAR.from.par(tau = tau.u, phi = phi.u,
                                CAR.D = CAR.D.shortcut,
                                CAR.W = CAR.W.shortcut)
    H.pro.short <- CAR.from.par(tau = tau.z, phi = phi.z,
                                CAR.D = CAR.D.shortcut,
                                CAR.W = CAR.W.shortcut)
    Q.pro.short <- Q.from.par(rho = rho,
                              G = G.pro.short,
                              H = H.pro.short)
    P.pro.short <- as.spam(rbind(
      cbind(G.pro.short, Q.pro.short),
      cbind(Q.pro.short, H.pro.short)
    ))
    e.val <- range(eigen(P.pro.short, only.values = TRUE)$values)
  } else {
    e.val <- 1
  }
  
  
  trace <- matrix(NA,
                  nrow = total.iters,
                  ncol = length(beta) + length(gamma) +
                    length(tau.u) + length(tau.z) +
                    length(phi.u) + length(phi.z) +
                    length(rho) +
                    1 + 1)
  
  colnames(trace) <- c(
    paste0("beta.", colnames(Xy)),
    paste0("gamma.", colnames(Xz)),
    "tau.u", "tau.z", "phi.u", "phi.z", "rho",
    "spline.penalty",
    "acs.scalar"
  )
  
  z.coord <- which(str_detect(colnames(Xy), "^z."))
  log.acs.prog <- Xy[, -z.coord] %*% beta[-z.coord] + u
  log.acs.prog.offset <- max(log.acs.prog)
  acs.scalar <- exp(log.acs.prog.offset) * colMeans(exp(log.acs.prog - log.acs.prog.offset))
  
  trace[1, ] <- c(beta, gamma,
                  tau.u, tau.z, phi.u, phi.z, rho,
                  spline.penalty,
                  acs.scalar)
  
  u.mean <- rep(0, length(u))
  
  ### MCMC ###
  
  for (iter in 2 : total.iters) {
    if (iter %% print.iters == 0) {
      print(paste(iter,
                  round(accept.beta / iter, 2),
                  round(accept.prec / iter, 2),
                  round(accept.u / (iter * length(u)), 2),
                  round(accept.z / (iter * length(z)), 2),
                  ifelse(iter > burn.in + 1,
                         round(mean(trace[(burn.in + 1) : (iter - 1), "rho"]), 2),
                         NA),
                  ifelse(iter > burn.in + 1,
                         round(100 * mean(trace[(burn.in + 1) : (iter - 1), "beta.z.1"]), 2),
                         NA),
                  ifelse(iter > burn.in + 1,
                         round(100 * quantile(trace[(burn.in + 1) : (iter - 1), "beta.z.1"],
                                        prob = 0.025), 2),
                         NA),
                  ifelse(iter > burn.in + 1,
                         round(100 * quantile(trace[(burn.in + 1) : (iter - 1), "beta.z.1"],
                                        prob = 0.975), 2),
                         NA)
                  )
            )
    }
    
    ### imputation ###
    
    for (i in missing.y) {
      log.rate.y <- offset[i] + Xy[i, ] %*% beta + u[i]
      rate.y <- exp(log.rate.y)
      
      y.poss <- 0 : 9
      log.probs <- dpois(y.poss, rate.y, log = TRUE)
      y[i] <- sample(y.poss, 1, prob = exp(log.probs - max(log.probs)))
    }
    
    
    
    ### beta ###
    
    if (VERBOSE) {
      print("beta")
    }
    
    coeff.prior[penalized] <- spline.penalty
    
    fit.b <- bayesglm(y ~ 0 + Xy + offset(offset + u), family = "poisson",
                      prior.scale = c(1 / coeff.prior), prior.df = Inf,
                      scaled = FALSE)
    post.mu <- coefficients(fit.b)
    if (rcond(vcov(fit.b)) < 1e-15) {
      post.Prec <- diag(ncol(vcov(fit.b)))
    } else {
      post.Prec <- solve(vcov(fit.b))
    }
    

    # beta.pro <- as.vector(backsolve(chol(post.Prec), rnorm(py))) + post.mu
    # 
    # if (mh.accept(
    #   log.dens.cur = log.lik.y(beta = beta, u = u, Xy = Xy, offset = offset, y = y) +
    #   log.prior.coeff(beta, coeff.prior = coeff.prior),
    #   log.dens.pro = log.lik.y(beta = beta.pro, u = u, Xy = Xy, offset = offset, y = y) +
    #   log.prior.coeff(beta.pro, coeff.prior = coeff.prior),
    #   log.prop.cur = - t(beta - post.mu) %*% post.Prec %*% (beta - post.mu) / 2,
    #   log.prop.pro = - t(beta.pro - post.mu) %*% post.Prec %*% (beta.pro - post.mu) / 2
    #   )) {
    # 
    #   beta <- beta.pro
    #   accept.beta <- accept.beta + 1
    # }
    
    beta.pro <- beta + as.vector(backsolve(chol(post.Prec * 10), rnorm(py)))
    
    if (mh.accept(
      log.dens.cur = log.lik.y(beta = beta, u = u, Xy = Xy, offset = offset, y = y) +
      log.prior.coeff(beta, coeff.prior = coeff.prior),
      log.dens.pro = log.lik.y(beta = beta.pro, u = u, Xy = Xy, offset = offset, y = y) +
      log.prior.coeff(beta.pro, coeff.prior = coeff.prior),
      log.prop.cur = 0,
      log.prop.pro = 0
    )) {
      
      beta <- beta.pro
      accept.beta <- accept.beta + 1
    }
    
    ### spline penalty ###
    
    if (any(penalized)) {
      spline.penalty <- rgamma(1, 0.001 + sum(penalized) / 2,
                               0.001 + sum(beta[penalized] ^ 2) / 2 )
    }
    
    ### gamma ###
    
    if (VERBOSE) {
      print("gamma")
    }
    
    if (!all(Xz == 0)) {
      mu.gamma.prior <- rep(0, pz)
      Prec.gamma.prior <- diag(coeff.prior, pz)
      
      if (VERBOSE) {
        print("prec")
      }
      
      Prec.gamma.lik <- as.spam(t(Xz) %*% H %*% Xz)
      
      if (VERBOSE) {
        print("mu")
      }
      
      if (rho == 0) {
        mu.gamma.lik <- solve(Prec.gamma.lik) %*% t(Xz) %*% H %*% (z)
      } else {
        # mu.gamma.lik <- solve(Prec.gamma.lik) %*% t(Xz) %*% H %*% (z + H.inv %*% Q %*% u)
        mu.gamma.lik <- solve(Prec.gamma.lik) %*% t(Xz) %*% H %*% (z + solve(H, Q %*% u))
      }
      
      
      if (VERBOSE) {
        print("post")
      }
      
      post.par <- post.par.gaussian.prec(prior.mu = mu.gamma.prior,
                                         prior.prec = Prec.gamma.prior,
                                         lik.mu = mu.gamma.lik,
                                         lik.prec = Prec.gamma.lik)
      
      if (VERBOSE) {
        print("backsolve")
      }
      
      gamma <- as.vector(backsolve.spam(chol(as.spam(post.par$Prec)), rnorm(pz))) + post.par$mu
    }
    
    
    ### u ###
    
    if (VERBOSE) {
      print("u")
    }
    
    if (SPATIAL && ESTIMATE_U) {
      for (i in 1 : length(u)) {
        if (length(NB[[i]]) > 0) {
          mu.u.prior <- drop(phi.u / length(NB[[i]]) * sum(u[NB[[i]]]) +
                               rho * sqrt(tau.z / tau.u) * (z[i] - Xz[i, , drop = FALSE] %*% gamma))
        } else {
          mu.u.prior <- 0
        }
        
        
        prec.u.prior <- tau.u * max(length(NB[[i]]), 1)
        
        offset.u <- offset[i] + drop(Xy[i, , drop = FALSE] %*% beta)
        
        post.approx <- pirls.poisson.u(y = y[i], offset = offset.u,
                                       prec.0 = prec.u.prior, mu.0 = mu.u.prior,
                                       n.iters = 5)
        
        u.pro <- rnorm(1, post.approx$mu, sqrt(1 / post.approx$prec))
        
        if (mh.accept(
          log.dens.cur = ldfu(u = u[i], offset = offset.u, y = y[i],
                              mu = mu.u.prior, prec = prec.u.prior),
          log.dens.pro = ldfu(u = u.pro, offset = offset.u, y = y[i],
                              mu = mu.u.prior, prec = prec.u.prior),
          log.prop.cur = dnorm(u[i], post.approx$mu, sqrt(1 / post.approx$prec),
                               log = TRUE),
          log.prop.pro = dnorm(u.pro, post.approx$mu, sqrt(1 / post.approx$prec),
                               log = TRUE))) {
          
          u[i] <- u.pro
          accept.u <- accept.u + 1
        }
      }
    }
    
    if (iter > burn.in) {
      u.mean <- (u.mean * (iter - burn.in - 1) + u) / (iter - burn.in)
    }
    
    
    
    ### dependence ###
    
    if (SPATIAL) {
      
      if (VERBOSE) {
        print("dependence")
      }
      
      # prop.scale <- 1 / 4
      phi.u.pro <- dexpit(dlogit(phi.u) + rnorm(1, 0, (0.35 * prop.scale)))
      phi.z.pro <- dexpit(dlogit(phi.z) + rnorm(1, 0, (0.35 * prop.scale)))
      tau.u.pro <- exp(log(tau.u) + rnorm(1, 0, 0.20 * prop.scale))
      # tau.u.pro <- tau.u
      tau.z.pro <- exp(log(tau.z) + rnorm(1, 0, 0.20 * prop.scale))
      # tau.z.pro <- tau.z
      
      while (CONSTRAIN_SCALES && (phi.u.pro < phi.z.pro)) {
        phi.u.pro <- dexpit(dlogit(phi.u) + rnorm(1, 0, (0.35 * prop.scale)))
        phi.z.pro <- dexpit(dlogit(phi.z) + rnorm(1, 0, (0.35 * prop.scale)))
      }

      if (ESTIMATE_RHO) {
        rho.pro <- dexpit(dlogit(rho) + rnorm(1, 0, 0.05 * prop.scale))
      } else {
        rho.pro <- rho
      }

      # rho.pro <- rho

      if (VERBOSE) {
        print("construction")
      }

      G.pro <- as.spam(CAR.from.par(tau = tau.u.pro, phi = phi.u.pro,
                                 CAR.D = CAR.D, CAR.W = CAR.W))
      H.pro <- as.spam(CAR.from.par(tau = tau.z.pro, phi = phi.z.pro,
                                 CAR.D = CAR.D, CAR.W = CAR.W))
      Q.pro <- as.spam(Q.from.par(rho = rho.pro, G = G.pro, H = H.pro))

      if (VERBOSE) {
        print("log posterior densities...")
      }
      z.res <- as.vector(z - Xz %*% gamma)

      P.pro <- as.spam(rbind(cbind(G.pro, Q.pro),
                         cbind(Q.pro, H.pro)))

      if (VERBOSE) {
        print("cholesky")
      }
      cP.pro <- tryCatch({
        suppressWarnings(chol(P.pro))
      }, error = function(e) {
        NULL
      })

      if (USE_PRIOR && COND_PRIOR) {
        if (VERBOSE) {
          print("prior")
        }
        if (!is.null(CAR.D.shortcut)) {
          G.pro.short <- CAR.from.par(tau = tau.u.pro, phi = phi.u.pro,
                                      CAR.D = CAR.D.shortcut,
                                      CAR.W = CAR.W.shortcut)
          H.pro.short <- CAR.from.par(tau = tau.z.pro, phi = phi.z.pro,
                                      CAR.D = CAR.D.shortcut,
                                      CAR.W = CAR.W.shortcut)
          Q.pro.short <- Q.from.par(rho = rho.pro,
                                    G = G.pro.short,
                                    H = H.pro.short)
          P.pro.short <- as.spam(rbind(
            cbind(G.pro.short, Q.pro.short),
            cbind(Q.pro.short, H.pro.short)
          ))
          e.val.pro <- range(eigen(P.pro.short, only.values = TRUE)$values)
        } else {
          e.val.pro <- range(chol.ev.approx(P.pro, kmax = 10))
        }
      } else {
        e.val.pro <- e.val
      }
      
      # e.val.pro <- e.val

      if (VERBOSE) {
        print("log density eval")
      }
      log.post.cur <- full.cond.dep(tau.u = tau.u, tau.z = tau.z,
                                    phi.u = phi.u, phi.z = phi.z,
                                    rho = rho,
                                    G = G,
                                    H = H, Q = Q,
                                    u = u, z.res = z.res,
                                    P = P, cP = cP, e.val = e.val)
      log.post.pro <- full.cond.dep(tau.u = tau.u.pro, tau.z = tau.z.pro,
                                    phi.u = phi.u.pro, phi.z = phi.z.pro,
                                    rho = rho.pro,
                                    G = G.pro,
                                    H = H.pro, Q = Q.pro,
                                    u = u, z.res = z.res,
                                    P = P.pro, cP = cP.pro, e.val = e.val.pro)
      if (VERBOSE) {
        print("done")
      }

      if (mh.accept(log.dens.pro = log.post.pro,
                    log.dens.cur = log.post.cur)) {
        tau.u <- tau.u.pro
        tau.z <- tau.z.pro
        phi.u <- phi.u.pro
        phi.z <- phi.z.pro
        rho <- rho.pro

        G <- G.pro
        H <- H.pro
        Q <- Q.pro

        P <- P.pro
        cP <- cP.pro
        e.val <- e.val.pro

        accept.prec <- accept.prec + 1
      }
    } else { # non-spatial
      z.res <- z - Xz %*% gamma
      tau.z <- rgamma(1,
                      length(z) / 2,
                      sum(z.res ^ 2) / 2)
    }
    
    
    
    
    ### acs ###
    
    z.coord <- which(str_detect(colnames(Xy), "^z."))
    log.acs.prog <- Xy[, -z.coord] %*% beta[-z.coord] + u
    log.acs.prog.offset <- max(log.acs.prog)
    acs.scalar <- exp(log.acs.prog.offset) * colMeans(exp(log.acs.prog - log.acs.prog.offset))
    
    
    ### book-keeping ###
    
    trace[iter, ] <- c(beta, gamma,
                       tau.u, tau.z, phi.u, phi.z, rho,
                       spline.penalty,
                       acs.scalar)
  }
  
  return(list(trace = trace,
              final = list(
                beta = beta,
                gamma = gamma,
                tau.u = tau.u,
                tau.z = tau.z,
                phi.u = phi.u,
                phi.z = phi.z,
                rho = rho,
                u = u
              ),
              u.mean = u.mean,
              time = Sys.time() - start.time))
}
