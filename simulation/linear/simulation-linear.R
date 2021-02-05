library(coda)

source("model-fit-linear.R") # contains fit.affine() definition
source("sim-data-linear.R") # contains simulate.data() definition

simulate <- function(seed,
                     phi.u = 0.5, phi.z = 0.2, rho = 0.5,
                     misspecified = FALSE, nonnorm = FALSE, xucorr = FALSE,
                     model) {
  set.seed(seed)
  print(paste(seed, phi.u, phi.z, rho, misspecified, nonnorm))
  
  burn.in <- 1000
  n.iters <- 10000
  keep <- burn.in + (1 : n.iters)
  
  result <- array(NA, dim = c(10, 6),
                  dimnames = list(param = c("beta.", "beta.x", "beta.z",
                                            "gamma.x",
                                            "tau.u", "tau.z",
                                            "phi.u", "phi.z", "rho", "u.1"),
                                  summary = c("mean", "sd", "lower", "upper",
                                              "ess.burn", "ess.keep")))
  
  data <- simulate.data(seed = seed, n = 300,
                        tau.u = 1, tau.z = 1,
                        phi.u = phi.u, phi.z = phi.z,
                        rho = rho,
                        misspecified = misspecified,
                        nonnorm = nonnorm,
                        x.u.corr = xucorr)
  
  sp <- attributes(data)$spatial
  CAR.D <- sp$CAR.D
  CAR.W <- sp$CAR.W
  NB <- sp$NB
  
  
  y <- data$y
  z <- data$z - mean(data$z)
  Xy <- as.matrix(cbind(1, data$x - mean(data$x), z))
  colnames(Xy) <- c("", "x", "z")
  Xz <- as.matrix(cbind(data$x - mean(data$x)))
  colnames(Xz) <- c("x")
  
  print(model)
  fit <- fit.affine(y = y, z = z, Xz = Xz, Xy = Xy,
                    seed = seed,
                    CAR.D = CAR.D, CAR.W = CAR.W, NB = NB,
                    print.iters = 1000, burn.in = burn.in, n.iters = n.iters,
                    CAR.D.shortcut = sp$CAR.D.shortcut,
                    CAR.W.shortcut = sp$CAR.W.shortcut,
                    CONSTRAIN_SCALES = model["CONSTRAIN_SCALES"],
                    ESTIMATE_RHO = model["ESTIMATE_RHO"],
                    SPATIAL = model["SPATIAL"],
                    USE_PRIOR = model["USE_PRIOR"])
  trace <- fit$trace


  result[, ] <- cbind(colMeans(trace[keep, ]),
                      apply(trace[keep, ], 2, sd),
                      t(apply(trace[keep, ], 2, quantile,
                              prob = c(0.025, 0.975))),
                      effectiveSize(trace[-keep, ]),
                      effectiveSize(trace[keep, ]))
  
  return(result)
}

n.sims <- 500
dgm.names <- c("independent",
               "large-scale confounder",
               "large-scale exposure",
               "same scales",
               "misspecified",
               "non-normal",
               "xucorr")
phi.u <- c(0.5, 0.5, 0.2, 0.35, 0.5, 0.5, 0.5)
phi.z <- c(0.2, 0.2, 0.5, 0.35, 0.2, 0.2, 0.2)
rho <- c(0.0, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3)
misspecified <- c(FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE)
nonnorm <- c(FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE)
xucorr <- c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE)

result <- array(NA, dim = c(7, n.sims, 5, 10, 6))
dimnames(result) <- list(dgm = dgm.names,
                         seed = 1 : n.sims,
                         model = c("non-spatial",
                                   "spatial",
                                   "spatial-constrained",
                                   "affine",
                                   "affine-constrained"),
                         param = c("beta.", "beta.x", "beta.z",
                                   "gamma.x",
                                   "tau.u", "tau.z",
                                   "phi.u", "phi.z", "rho",
                                   "u.1"),
                         summary = c("mean", "sd",  "lower", "upper",
                                     "ess.burn", "ess.keep"))

models <- c("non-spatial",
            "spatial", "spatial-constrained",
            "affine", "affine-constrained")

model.args <- list(
  `non-spatial` = c(SPATIAL = FALSE, CONSTRAIN_SCALES = FALSE, ESTIMATE_RHO = FALSE, USE_PRIOR = FALSE),
  `spatial` = c(SPATIAL = TRUE, CONSTRAIN_SCALES = FALSE, ESTIMATE_RHO = FALSE, USE_PRIOR = FALSE),
  `spatial-constrained` = c(SPATIAL = TRUE, CONSTRAIN_SCALES = TRUE, ESTIMATE_RHO = FALSE, USE_PRIOR = TRUE),
  `affine` = c(SPATIAL = TRUE, CONSTRAIN_SCALES = FALSE, ESTIMATE_RHO = TRUE, USE_PRIOR = TRUE),
  `affine-constrained` = c(SPATIAL = TRUE, CONSTRAIN_SCALES = TRUE, ESTIMATE_RHO = TRUE, USE_PRIOR = TRUE)
)

for (dgm in 1 : 6) {
  for (seed in 1 : n.sims) {
    for (model in models) {
      if (is.na(result[dgm, seed, model, 1, 1])) {
        print(Sys.time())
        result[dgm.names[dgm], seed, model, , ] <- simulate(
          seed,
          phi.u = phi.u[dgm],
          phi.z = phi.z[dgm],
          rho = rho[dgm],
          misspecified = misspecified[dgm],
          nonnorm = nonnorm[dgm],
          xucorr = xucorr[dgm],
          model = model.args[[model]]
        )

        save(result, file = "result-linear-n300-b1k-r10k.RData")
      }
    }
  }
}

for (dgm in 7) {
  for (seed in 1 : 100) {
    for (model in models) {
      if (is.na(result[dgm, seed, model, 1, 1])) {
        print(Sys.time())
        result[dgm.names[dgm], seed, model, , ] <- simulate(
          seed,
          phi.u = phi.u[dgm],
          phi.z = phi.z[dgm],
          rho = rho[dgm],
          misspecified = misspecified[dgm],
          nonnorm = nonnorm[dgm],
          xucorr = xucorr[dgm],
          model = model.args[[model]]
        )
        
        save(result, file = "result-linear-n300-b1k-r10k.RData")
      }
    }
  }
}




round(colMeans(result["independent", , , , "mean"], na.rm = TRUE), 2)
round(colMeans((result["independent", , , "beta.z", "mean"] - 1) ^ 2, na.rm = TRUE), 2)
colMeans(result["independent", , , "beta.z", "lower"] < 1 &
           result["independent", , , "beta.z", "upper"] > 1, na.rm = TRUE)
boxplot(result["independent", , , "beta.z", "mean"])
abline(h = 1, lty = 3)

round(colMeans(result["large-scale confounder", , , , "mean"], na.rm = TRUE), 2)
round(colMeans((result["large-scale confounder", , , "beta.z", "mean"] - 1) ^ 2, na.rm = TRUE), 2)
colMeans(result["large-scale confounder", , , "beta.z", "lower"] < 1 &
           result["large-scale confounder", , , "beta.z", "upper"] > 1, na.rm = TRUE)
boxplot(result["large-scale confounder", , , "beta.z", "mean"])
abline(h = 1, lty = 3)

round(colMeans(result["large-scale exposure", , , , "mean"], na.rm = TRUE), 2)
round(colMeans((result["large-scale exposure", , , "beta.z", "mean"] - 1) ^ 2, na.rm = TRUE), 2)
colMeans(result["large-scale exposure", , , "beta.z", "lower"] < 1 &
           result["large-scale exposure", , , "beta.z", "upper"] > 1, na.rm = TRUE)
boxplot(result["large-scale exposure", , , "beta.z", "mean"])
abline(h = 1, lty = 3)

round(colMeans(result["same scales", , , , "mean"], na.rm = TRUE), 2)
round(colMeans((result["same scales", , , "beta.z", "mean"] - 1) ^ 2, na.rm = TRUE), 2)
colMeans(result["same scales", , , "beta.z", "lower"] < 1 &
           result["same scales", , , "beta.z", "upper"] > 1, na.rm = TRUE)
boxplot(result["same scales", , , "beta.z", "mean"])
abline(h = 1, lty = 3)


round(colMeans(result["misspecified", , , , "mean"], na.rm = TRUE), 2)
round(colMeans((result["misspecified", , , "beta.z", "mean"] - 1) ^ 2, na.rm = TRUE), 2)
colMeans(result["misspecified", , , "beta.z", "lower"] < 1 &
           result["misspecified", , , "beta.z", "upper"] > 1, na.rm = TRUE)
boxplot(result["misspecified", , , "beta.z", "mean"])
abline(h = 1, lty = 3)

round(colMeans(result["non-normal", , , , "mean"], na.rm = TRUE), 2)
round(colMeans((result["non-normal", , , "beta.z", "mean"] - 1) ^ 2, na.rm = TRUE), 2)
colMeans(result["non-normal", , , "beta.z", "lower"] < 1 &
           result["non-normal", , , "beta.z", "upper"] > 1, na.rm = TRUE)
boxplot(result["non-normal", , , "beta.z", "mean"])
abline(h = 1, lty = 3)

round(colMeans(result["xucorr", , , , "mean"], na.rm = TRUE), 2)
round(colMeans((result["xucorr", , , "beta.z", "mean"] - 1) ^ 2, na.rm = TRUE), 2)
colMeans(result["xucorr", , , "beta.z", "lower"] < 1 &
           result["xucorr", , , "beta.z", "upper"] > 1, na.rm = TRUE)
boxplot(result["xucorr", , , "beta.z", "mean"])
abline(h = 1, lty = 3)

########################
### output for paper ###
########################

print.result.table <- function(dgm) {
  result.table <- matrix(NA, ncol = 4, nrow = 5)
  colnames(result.table) <- c("Bias", "Std. Err.", "RMSE", "Coverage")
  rownames(result.table) <- c("OLS", "GLS", "GLS-RS", "Affine", "Affine-RS")
  
  result.table[, "Bias"] <- colMeans(result[dgm, , ,"beta.z" , "mean"] - 1, na.rm = TRUE)
  result.table[, "Std. Err."] <- apply(result[dgm, , ,"beta.z" , "mean"], 2, sd, na.rm = TRUE)
  result.table[, "RMSE"] <- sqrt(colMeans((result[dgm, , ,"beta.z" , "mean"] - 1) ^ 2, na.rm = TRUE))
  result.table[, "Coverage"] <- colMeans(result[dgm, , , "beta.z", "lower"] < 1 &
                                           result[dgm, , , "beta.z", "upper"] > 1, na.rm = TRUE)
  
  charmat <- matrix(paste0("& ", format(round(result.table, 2), nsmall = 2), ""),
                    nrow = nrow(result.table),
                    ncol = ncol(result.table))
  charmat <- cbind("& &", rownames(result.table), charmat, "\\")
  prmatrix(
    charmat,
    rowlab = rep("", nrow(charmat)),
    collab = rep("", ncol(charmat)),
    quote = FALSE)
}

print.result.table("independent")
print.result.table("large-scale confounder")
print.result.table("large-scale exposure")
print.result.table("same scales")
print.result.table("misspecified")
print.result.table("non-normal")
print.result.table("xucorr")
