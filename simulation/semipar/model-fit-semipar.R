
library(arm)
library(spam)
library(Matrix)
library(lme4)

fit.affine <- function(y, z, Xz, Xy,
                       NB, CAR.D, CAR.W,
                       seed = 0,
                       init = list(),
                       print.iters = 100,
                       burn.in = 1000,
                       n.iters = 10000,
                       penalized = rep(FALSE, ncol(Xy)),
                       CONSTRAIN_SCALES = TRUE,
                       ESTIMATE_RHO = TRUE,
                       SPATIAL = TRUE,
                       USE_PRIOR = FALSE,
                       USE_KAPPA_PRIOR = USE_PRIOR,
                       CAR.D.shortcut = NULL, CAR.W.shortcut = NULL,
                       debug.control = list(),
                       VERBOSE = FALSE,
                       coeff.prior = 1 / 100,
                       tau.prior = 5,
                       rho.prior = 1,
                       ESTIMATE_U = TRUE,
                       ESTIMATE_GAMMA = TRUE,
                       ESTIMATE_DEP = TRUE,
                       z.test.spline = NULL) {
  
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
      # W <- diag(as.vector(pmin(a(mu), max.w)))
      W <- diag(as.vector(a(mu)))
      z <- eta - offset + (y - mu) / diag(W)
      M <- t(Xy) %*% W %*% Xy + P.0
      if (rcond(M) < 1e-8) {
        print(range(Xy))
        print(range(M))
        print(range(diag(P.0)))
        print(range(diag(W)))
        print(beta)
        # print(diag(W))
        # print(t(Xy) %*% W %*% Xy + P.0)
      }
      # beta[i, ] <- as.vector(ginv(t(Xy) %*% W %*% Xy + P.0) %*%
      #                          (P.0 %*% mu.0 + t(Xy) %*% W %*% z))
      beta[i, ] <- solve(M, P.0 %*% mu.0 + t(Xy) %*% W %*% z)
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
    -rho * diag.spam(sqrt(diag(G) * diag(H)))
  }
  
  log.prior.dep <- function(tau.u, tau.z, phi.u, phi.z, rho,
                            e.val) {
    if (CONSTRAIN_SCALES && phi.z > phi.u) {
      return(-Inf)
    }
    
    # if (abs(rho) >= 0.5) {
    #   return(-Inf)
    # }
    
    if (USE_PRIOR) {
      if (USE_KAPPA_PRIOR) {
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
  
  log.lik.y <- function(beta, u, Xy, y, sum = TRUE) {
    eta <- Xy %*% beta + u
    ll <- dpois(y, exp(eta), log = TRUE)
    if (sum) {
      return(sum(ll))
    } else {
      return(ll)
    }
  }
  
  full.cond.u <- function(u,
                          G, G.inv, Q, sigma.sq.e, beta, gamma,
                          Xz, Xy, z, y) {
    u.res <- u + G.inv %*% Q %*% (z - Xz %*% gamma)
    y.res <- y - Xy %*% beta - u
    
    return(
      -(1 / 2) * (
        t(u.res) %*% G %*% u.res +
          t(y.res) %*% y.res / sigma.sq.e
      )
    )
  }
  
  full.cond.beta <- function(beta,
                             sigma.sq.e, Xy, y) {
    y.res <- y - Xy %*% beta - u
    
    return(
      -(1 / 2) * t(y.res) %*% y.res / sigma.sq.e +
        log.prior.beta(beta)
    )
  }
  
  full.cond.gamma <- function(gamma,
                              H, H.inv, Q, Xz, u, z) {
    z.res <- z - zX %*% gamma + H.inv %*% Q %*% u
    
    return(
      -(1 / 2) * t(z.res) %*% H %*% z.res +
        log.prior.gamma(gamma)
    )
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
      (1 / 2) * drop(as.matrix(
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
    lik.V.inv <- ginv(lik.V)
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
  
  if (is.null(init[["beta"]])) {
    coeff.prior[penalized] <- spline.penalty
    P.0 <- diag(coeff.prior)
    fit <- pirls.poisson(Xy = Xy, offset = rep(0, n), y = y,
                         P.0 = P.0, mu.0 = rep(0, py),
                         init = coefficients(glm(y ~ 0 + Xy, family = "poisson")))
    beta <- fit$mu
  } else {
    beta <- init[["beta"]]
  }
  
  if (is.null(init[["gamma"]])) {
    gamma <- coefficients(lm(z ~ 0 + Xz))
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
  
  e.val <- eigen(P, only.values = TRUE)$values
  
  trace <- matrix(NA,
                  nrow = total.iters,
                  ncol = length(beta) + length(gamma) +
                    length(tau.u) + length(tau.z) +
                    length(phi.u) + length(phi.z) +
                    length(rho) +
                    1 + 1 +
                    length(z.test))
  
  colnames(trace) <- c(
    paste0("beta.", colnames(Xy)),
    paste0("gamma.", colnames(Xz)),
    "tau.u", "tau.z", "phi.u", "phi.z", "rho",
    "spline.penalty",
    "u.1",
    paste0("acs_", round(z.test, 2))
  )
  
  z.coord <- which(str_detect(colnames(Xy), "^z."))
  if (is.null(z.test.spline)) {
    z.test.spline <- Xy[, z.coord]
  }
  log.acs.prog <- Xy[, -z.coord] %*% beta[-z.coord] + u
  log.acs.test <- outer(drop(log.acs.prog), z.test.spline %*% beta[z.coord], `+`)
  log.acs.test.offset <- max(log.acs.test)
  acs <- exp(log.acs.test.offset) * colMeans(exp(log.acs.test - log.acs.test.offset))
  
  trace[1, ] <- c(beta, gamma,
                  tau.u, tau.z, phi.u, phi.z, rho,
                  u[1], spline.penalty,
                  acs)
  
  ### MCMC ###
  
  for (iter in 2 : total.iters) {
    if (iter %% print.iters == 0) {
      print(paste(iter,
                  round(accept.beta / iter, 2),
                  round(accept.prec / iter, 2),
                  round(accept.u / (iter * length(u)), 2),
                  ifelse(iter > burn.in + 1,
                         round(mean(trace[(burn.in + 1) : (iter - 1), "spline.penalty"]), 2),
                         NA)
                  )
            )
    }
    
    ### beta ###
    
    if (VERBOSE) {
      print("beta")
    }
    
    
    coeff.prior[penalized] <- spline.penalty
    
    fit.b <- bayesglm(y ~ 0 + Xy + offset(u), family = "poisson",
                      prior.scale = c(1 / coeff.prior), prior.df = Inf,
                      scaled = FALSE)
    post.mu <- coefficients(fit.b)
    post.Prec <- solve(vcov(fit.b))
    
    
    tryCatch(
      beta.pro <- as.vector(backsolve(chol(post.Prec), rnorm(py))) + post.mu,
      error = function(e) {print(post.Prec)})
    

    if (mh.accept(
      log.dens.cur = log.lik.y(beta = beta, u = u, Xy = Xy, y = y) +
      log.prior.coeff(beta, coeff.prior = coeff.prior),
      log.dens.pro = log.lik.y(beta = beta.pro, u = u, Xy = Xy, y = y) +
      log.prior.coeff(beta.pro, coeff.prior = coeff.prior),
      log.prop.cur = - t(beta - post.mu) %*% post.Prec %*% (beta - post.mu) / 2,
      log.prop.pro = - t(beta.pro - post.mu) %*% post.Prec %*% (beta.pro - post.mu) / 2)) {

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
    
    if (ESTIMATE_GAMMA) {
      mu.gamma.prior <- rep(0, pz)
      Prec.gamma.prior <- diag(coeff.prior, pz)
      
      Prec.gamma.lik <- as.spam(t(Xz) %*% H %*% Xz)
      # mu.gamma.lik <- solve(Prec.gamma.lik) %*% t(Xz) %*% H %*% (z + H.inv %*% Q %*% u)
      mu.gamma.lik <- solve(Prec.gamma.lik) %*% t(Xz) %*% H %*% (z + solve(H, Q %*% u))
      
      post.par <- post.par.gaussian.prec(prior.mu = mu.gamma.prior,
                                         prior.prec = Prec.gamma.prior,
                                         lik.mu = mu.gamma.lik,
                                         lik.prec = Prec.gamma.lik)
      gamma <- as.vector(backsolve.spam(chol(as.spam(post.par$Prec)), rnorm(pz))) + post.par$mu
    }
    
    
    
    
    
    ### u ###
    
    if (VERBOSE) {
      print("u")
    }
    
    if (SPATIAL && ESTIMATE_U) {
      for (i in 1 : length(u)) {
        mu.u.prior <- drop(phi.u / length(NB[[i]]) * sum(u[NB[[i]]]) +
                             rho * sqrt(tau.z / tau.u) * (z[i] - Xz[i, , drop = FALSE] %*% gamma))

        prec.u.prior <- tau.u * length(NB[[i]])

        offset <- drop(Xy[i, , drop = FALSE] %*% beta)

        post.approx <- pirls.poisson.u(y = y[i], offset = offset,
                                       prec.0 = prec.u.prior, mu.0 = mu.u.prior,
                                       n.iters = 5)
        
        u.pro <- rnorm(1, post.approx$mu, sqrt(1 / post.approx$prec))
        
        if (mh.accept(
          log.dens.cur = ldfu(u = u[i], offset = offset, y = y[i],
                              mu = mu.u.prior, prec = prec.u.prior),
          log.dens.pro = ldfu(u = u.pro, offset = offset, y = y[i],
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
    
    
    ### dependence ###
    
    if (SPATIAL && ESTIMATE_DEP) {
      
      if (VERBOSE) {
        print("dependence")
      }
      
      prop.scale <- 1 / 4
      phi.u.pro <- dexpit(dlogit(phi.u) + rnorm(1, 0, (0.35 * prop.scale)))
      phi.z.pro <- dexpit(dlogit(phi.z) + rnorm(1, 0, (0.35 * prop.scale)))
      tau.u.pro <- exp(log(tau.u) + rnorm(1, 0, 0.20 * prop.scale))
      # tau.u.pro <- tau.u
      tau.z.pro <- exp(log(tau.z) + rnorm(1, 0, 0.20 * prop.scale))
      # tau.z.pro <- tau.z
      
      if (ESTIMATE_RHO) {
        rho.pro <- dexpit(dlogit(rho) + rnorm(1, 0, 0.5 * prop.scale))
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
      
      if (USE_KAPPA_PRIOR) {
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
        e.val.pro <- e.val
      }
      
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
        
        accept.prec <- accept.prec + 1
        
        e.val <- e.val.pro
      }
    }
    
    ### acs ###
    
    z.coord <- which(str_detect(colnames(Xy), "^z."))
    log.acs.prog <- Xy[, -z.coord] %*% beta[-z.coord] + u
    log.acs.test <- outer(drop(log.acs.prog), z.test.spline %*% beta[z.coord], `+`)
    log.acs.test.offset <- max(log.acs.test)
    acs <- exp(log.acs.test.offset) * colMeans(exp(log.acs.test - log.acs.test.offset))
    
    
    ### book-keeping ###
    
    trace[iter, ] <- c(beta, gamma,
                       tau.u, tau.z, phi.u, phi.z, rho,
                       spline.penalty,
                       u[1],
                       acs)
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
              )))
  
  
}





if (FALSE) {
  
  keep.plot <- seq(10, 11000, by = 10)
  plot((trace[, "tau.u"]))
  plot((trace[, "tau.z"]))
  plot((trace[, "phi.u"]))
  plot((trace[, "phi.z"]))
  # plot((trace[, "sigma.sq.e"]))
  plot(trace[, "rho"])
  plot(trace[, "beta.z"])
  plot(trace[, "u.1"])
  
  trace.trans <- matrix(NA, ncol = 5, nrow = nrow(trace))
  colnames(trace.trans) <- c("f", "g", "dl.phi.u", "dl.phi.z", "l.tau.u")
  trace.trans[, "l.tau.u"] <- log(trace[, "tau.u"])
  trace.trans[, "f"] <- log(trace[, "tau.z"]) - log(trace[, "tau.u"]) + log(1 - trace[, "rho"])
  trace.trans[, "g"] <- log(trace[, "tau.z"]) - log(trace[, "tau.u"]) - log(1 - trace[, "rho"])
  trace.trans[, "dl.phi.u"] <- dlogit(trace[, "phi.u"])
  trace.trans[, "dl.phi.z"] <- dlogit(trace[, "phi.z"])
  # trace.trans[, "rho"] <- dlogit(trace[, "rho"])
  
  
  cor(trace.trans[, c("f", "g", "dl.phi.u", "dl.phi.z")])
  round(var(trace.trans[, c("f", "g", "dl.phi.u", "dl.phi.z")]), 2)
  round(sqrt(diag(var(trace.trans[, c("f", "g", "dl.phi.u", "dl.phi.z")]))), 2)
  
  plot(log(trace[, "tau.z"]) - log(trace[, "tau.u"]), dlogit(trace[, "rho"]))
  
  plot(dlogit(trace[, c("phi.u", "phi.z")]))
  
  plot(log(trace[, "tau.z"] / trace[, "tau.u"]) + log(1 - trace[, "rho"]))
  plot(log(trace[, "tau.z"] / trace[, "tau.u"]) - log(1 - trace[, "rho"]))
  plot(log(trace[, "tau.z"] / trace[, "tau.u"]) + log(1 - trace[, "rho"]),
       log(trace[, "tau.z"] / trace[, "tau.u"]) - log(1 - trace[, "rho"]))
}
