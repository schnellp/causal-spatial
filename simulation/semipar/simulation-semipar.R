library(coda)
library(ggplot2)
library(stringr)
library(splines)

source("model-fit-semipar.R")
source("sim-data-semipar.R")

n.knots <- 5

construct.knots <- function(x, max.knots=20) {
  u <- unique(x)
  n.knots <- min(ceiling(length(u) / 4), max.knots)
  quantile(u, ((1:n.knots) + 1) / (n.knots + 2))
}

expand.spline <- function(x, degree, knots) {
  
  X <- sapply(1:degree, function(p, x) {x ^ p }, x=x)
  if (is.null(dim(X))) { # when x is a single number, X is a vector
    X <- t(X)
  }
  colnames(X)[1:degree] <- paste0(".", 1:degree)
  
  Z <- matrix(NA, nrow=length(x), ncol=length(knots))
  colnames(Z) <- paste0(".KNOT.", 1:length(knots))
  for (k in 1:length(knots)) {
    Z[, k] <- abs(x - knots[k]) ^ degree
  }
  
  spline.grid <- cbind(X, Z)
  attr(spline.grid, "penalized") <- c(rep(0, degree), rep(1, length(knots)))
  
  spline.grid
}

simulate <- function(seed,
                     phi.u = 0.5, phi.z = 0.2, rho = 0.3,
                     misspecified = FALSE, nonnorm = FALSE,
                     model,
                     degree = 3, z.test = round(seq(-2, 2, by = 0.1), 1)) {
  set.seed(seed)
  print(paste(seed, phi.u, phi.z, rho, misspecified, nonnorm))
  
  burn.in <- 5000
  n.iters <- 10000
  keep <- burn.in + (1 : n.iters)
  
  
  data <- simulate.data(seed = seed, n = 300,
                        tau.u = 1, tau.z = 1,
                        phi.u = phi.u, phi.z = phi.z,
                        rho = rho,
                        misspecified = misspecified,
                        nonnorm = nonnorm)
  
  sp <- attributes(data)$spatial
  CAR.D <- sp$CAR.D
  CAR.W <- sp$CAR.W
  NB <- sp$NB
  
  
  y <- data$y
  z <- data$z - mean(data$z)
  Xy <- as.matrix(cbind(1, data$x - mean(data$x), z))
  knots <- construct.knots(z, max.knots = n.knots)
  # knots <- seq(-3, 3, length.out = n.knots)
  Xy.spline <- expand.spline(z, degree, knots)
  z.test.spline <- expand.spline(z.test, degree, knots)
  # Xy.spline <- bs(z, degree = degree, df = degree + n.knots)
  # colnames(Xy.spline) <- paste0("z.KNOT.", colnames(Xy.spline))
  colnames(Xy.spline) <- paste0("z", colnames(Xy.spline))
  Xy <- cbind(1, x = data$x - mean(data$x), Xy.spline)
  penalized <- str_detect(colnames(Xy), "\\.KNOT\\.")
  Xz <- as.matrix(cbind(data$x - mean(data$x)))
  colnames(Xz) <- c("x")
  
  
  print(model)
  fit <- fit.affine(y = y, z = z, Xz = Xz, Xy = Xy,
                    seed = seed,
                    CAR.D = CAR.D, CAR.W = CAR.W, NB = NB,
                    print.iters = 1000, burn.in = burn.in, n.iters = n.iters,
                    CAR.D.shortcut = sp$CAR.D.shortcut,
                    CAR.W.shortcut = sp$CAR.W.shortcut,
                    penalized = penalized,
                    CONSTRAIN_SCALES = model["CONSTRAIN_SCALES"],
                    ESTIMATE_RHO = model["ESTIMATE_RHO"],
                    SPATIAL = model["SPATIAL"],
                    USE_PRIOR = model["USE_PRIOR"],
                    z.test.spline = z.test.spline)
  trace <- fit$trace
  
  ###
  
  D.test <- expand.spline(z.test, degree, knots)
  spline.coeff.trace <- trace[, str_detect(colnames(trace), "beta.z.")]
  spline.trace <- spline.coeff.trace %*% t(D.test)
  
  
  
  ###
  
  record <- trace
  
  result <- array(NA, dim = c(ncol(Xy) + 8 + length(z.test), 6),
                  dimnames = list(param = c(paste0("beta.", colnames(Xy)),
                                            "gamma.x",
                                            "tau.u", "tau.z",
                                            "phi.u", "phi.z", "rho", "spline.penalty", "u.1",
                                            paste0("acs_", z.test)),
                                  summary = c("mean", "sd", "lower", "upper",
                                              "ess.burn", "ess.keep")))


  result <- cbind(colMeans(record[keep, ]),
                      apply(record[keep, ], 2, sd),
                      t(apply(record[keep, ], 2, quantile,
                              prob = c(0.025, 0.975))),
                      effectiveSize(record[-keep, ]),
                      effectiveSize(record[keep, ]))
  
  return(result)
}

n.sims <- 500
dgm.names <- c("independent",
               "large-scale confounder",
               "large-scale exposure",
               "same scales",
               "misspecified",
               "non-normal")
phi.u <- c(0.5, 0.5, 0.2, 0.35, 0.5, 0.5)
phi.z <- c(0.2, 0.2, 0.5, 0.35, 0.2, 0.2)
rho <- c(0.0, 0.3, 0.3, 0.3, 0.3, 0.3)
misspecified <- c(FALSE, FALSE, FALSE, FALSE, TRUE, FALSE)
nonnorm <- c(FALSE, FALSE, FALSE, FALSE, FALSE, TRUE)

z.test <- round(seq(-2, 2, by = 0.1), 1)

result <- array(NA, dim = c(6, n.sims, 5, 13 + n.knots + length(z.test), 6))
dimnames(result) <- list(dgm = dgm.names,
                         seed = 1 : n.sims,
                         model = c("non-spatial",
                                   "spatial",
                                   "spatial-constrained",
                                   "affine",
                                   "affine-constrained"),
                         param = c("beta.", "beta.x",
                                   "beta.z.1", "beta.z.2", "beta.z.3",
                                   paste0("beta.z.KNOT.", 1 : n.knots),
                                   "gamma.x",
                                   "tau.u", "tau.z",
                                   "phi.u", "phi.z", "rho",
                                   "spline.penalty",
                                   "u.1",
                                   paste0("acs_", z.test)),
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

for (dgm in c(2)) {
  for (seed in 1 : 500) {
    for (model in c("spatial", "affine-constrained")) {
      if (is.na(result[dgm, seed, model, 1, 1])) {
        print(Sys.time())
        result[dgm.names[dgm], seed, model, , ] <- simulate(
          seed,
          phi.u = phi.u[dgm],
          phi.z = phi.z[dgm],
          rho = rho[dgm],
          misspecified = misspecified[dgm],
          nonnorm = nonnorm[dgm],
          model = model.args[[model]],
          z.test = z.test
        )
        
        save(result, file = "result-semipar-n300-b5k-r10k.RData")
      }
    }
  }
}


true.acs <- function(z) {
  (2 / (1 + exp(-6 * z))) - 1
}

set.seed(0)
data <- data.frame()
for (i in 1 : 100) {
  data.partial <- simulate.data(seed = i, n = 300,
                                tau.u = 1, tau.z = 1,
                                phi.u = 0.5, phi.z = 0.2,
                                rho = 0.3,
                                misspecified = FALSE,
                                nonnorm = FALSE)
  data <- rbind(data, cbind(data.partial$x, attributes(data.partial)$truth$u))
}
colnames(data) <- c("x", "u")

prod.truth <- mean(exp(data$u))

acs.idx <- which(stringr::str_detect(dimnames(result)$param, "acs_"))

log.acs.affine.constrained.mean <- colMeans(log(result["large-scale confounder",
                                                       ,
                                                       "affine-constrained",
                                                       acs.idx,
                                                       "mean"]),
                                            na.rm = TRUE)



log.acs.spatial.unconstrained.mean <- colMeans(log(result["large-scale confounder",
                                                          ,
                                                          "spatial",
                                                          acs.idx,
                                                          "mean"]),
                                               na.rm = TRUE)

log.acs.spatial.unconstrained <- log(result["large-scale confounder",
                                            ,
                                            "spatial",
                                            acs.idx,
                                            "mean"])

log.acs.affine.constrained <- log(result["large-scale confounder",
                                         ,
                                         "affine-constrained",
                                         acs.idx,
                                         "mean"])

plot.data <- rbind(
  data.frame(Exposure = z.test,
             Model = "Spatial unconstriained",
             Summary = "Posterior Mean",
             Quantity = "Mean",
             Value = log.acs.spatial.unconstrained.mean),
  data.frame(Exposure = z.test,
             Model = "Spatial unconstriained",
             Summary = "95% CI",
             Quantity = "Lower",
             Value = apply(log.acs.spatial.unconstrained, 2, quantile, prob = 0.025, na.rm = TRUE)),
  data.frame(Exposure = z.test,
             Model = "Spatial unconstriained",
             Summary = "95% CI",
             Quantity = "Upper",
             Value = apply(log.acs.spatial.unconstrained, 2, quantile, prob = 0.975, na.rm = TRUE)),
  data.frame(Exposure = z.test,
             Model = "Affine constrained",
             Summary = "Posterior Mean",
             Quantity = "Mean",
             Value = log.acs.affine.constrained.mean),
  data.frame(Exposure = z.test,
             Model = "Affine constrained",
             Summary = "95% CI",
             Quantity = "Lower",
             Value = apply(log.acs.affine.constrained, 2, quantile, prob = 0.025, na.rm = TRUE)),
  data.frame(Exposure = z.test,
             Model = "Affine constrained",
             Summary = "95% CI",
             Quantity = "Upper",
             Value = apply(log.acs.affine.constrained, 2, quantile, prob = 0.975, na.rm = TRUE))
)

plot.data <- as.data.table(plot.data)
plot.data[, key := paste0(Exposure, Model)]
setorder(plot.data, key)

plot_dta_mean <- subset(plot.data, Quantity == 'Mean')
plot_dta_LB <- subset(plot.data, Quantity == 'Lower')
plot_dta_UB <- subset(plot.data, Quantity == 'Upper')

sum(plot_dta_mean$key == plot_dta_LB$key)
sum(plot_dta_mean$key == plot_dta_UB$key)

plot_dta <- copy(plot_dta_mean)
plot_dta[, LB := plot_dta_LB$Value]
plot_dta[, UB := plot_dta_UB$Value]

true_dta <- data.frame(Exposure = unique(plot.data$Exposure))
true_dta$Value <- true.acs(true_dta$Exposure) + log(prod.truth)

theme_set(theme_light())

pdf("sim-semipar.pdf", width = 6, height = 4)
ggplot(as.data.frame(plot_dta), aes(x = Exposure, y = Value)) +
  geom_line(aes(linetype = Model), size = 1) +
  geom_ribbon(aes(fill = Model, color = Model, ymin = LB,
                  ymax = UB), size = 0.2, alpha = 0.1) +
  scale_color_manual(values = c('black', 'grey75')) +
  scale_fill_manual(values = c('black', 'grey75')) +
  scale_linetype_manual(values = c("dashed", "dotted")) +
  geom_line(data = true_dta, aes(x = Exposure, y = Value), linetype = 1,
            size = 1, col = 'black') +
  ylab("Log PAERC") +
  theme(legend.position = c(0.25, 0.8),
        legend.key.width = grid::unit(3, "lines"),
        legend.box.background = element_rect())
dev.off()


