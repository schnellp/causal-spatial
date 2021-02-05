source("data-compilation.R")
source("model-fit-semipar.R")



######################
### prior shortcut ###
######################

m.s <- 4
n.s <- 4
N.s <- m.s * n.s

coord.to.index <- function(row, col) {
  (row - 1) * m.s + col
}

CAR.W.s <- matrix(0, nrow = N.s, ncol = N.s)
for (row in 1 : m.s) {
  for (col in 1 : n.s) {
    i <- coord.to.index(row, col)
    if (row > 1) {
      CAR.W.s[i, coord.to.index(row - 1, col)] <- 1
    }
    if (row < m.s) {
      CAR.W.s[i, coord.to.index(row + 1, col)] <- 1
    }
    if (col > 1) {
      CAR.W.s[i, coord.to.index(row, col - 1)] <- 1
    }
    if (col < n.s) {
      CAR.W.s[i, coord.to.index(row, col + 1)] <- 1
    }
  }
}
CAR.D.s <- diag(colSums(CAR.W.s))

################
### analysis ###
################

n.knots <- 10
degree <- 3
z.test <- seq(0, 25, by = 0.1)

construct.knots <- function(x, max.knots=20) {
  u <- unique(x)
  n.knots <- min(ceiling(length(u) / 4), max.knots)
  quantile(u, ((1:n.knots) + 1) / (n.knots + 2))
}

expand.spline <- function(x, degree, knots, bound = Inf) {
  
  X <- sapply(1:degree, function(p, x) {x ^ p }, x=x)
  if (is.null(dim(X))) { # when x is a single number, X is a vector
    X <- t(X)
  }
  colnames(X)[1:degree] <- paste0(".", 1:degree)
  
  Z <- matrix(NA, nrow=length(x), ncol=length(knots))
  colnames(Z) <- paste0(".KNOT.", 1:length(knots))
  for (k in 1:length(knots)) {
    Z[, k] <- ifelse(x > knots[k],
                     ((x - knots[k]) / 100) ^ degree,
                     0)
    # q <- abs(x - knots[k]) ^ degree
    # Z[, k] <- sign(q) * pmin(abs(q), bound)
  }
  
  spline.grid <- cbind(X, Z)
  attr(spline.grid, "penalized") <- c(rep(0, degree), rep(1, length(knots)))
  
  spline.grid
}

y <- data$Deaths
offset <- log(data$Expected)
min.nzz <- min(data$PCT_HHNV1MI[data$PCT_HHNV1MI > 0])
z <- log(pmax(data$PCT_HHNV1MI, min.nzz)) - mean(log(pmax(data$PCT_HHNV1MI, min.nzz)))

saved.na.action <- options('na.action')
options(na.action = 'na.pass')
Xy <- model.matrix( ~ PCT_HHNV1MI +
                      SmokeRate +
                      log(TotPop) + log(PopPerSQM) +
                      log(MedianHValue) + log(MedianHHInc) +
                      PctUrban + PctHighSchool + PctPoor + PctFemale +
                      PctMovedIn5 + PctOccupied + PctWhite + PctBlack + PctHisp,
                    data = data)
colnames(Xy)[2] <- "z"
Xy <- scale(Xy)
Xy[, 1] <- 1

Xz <- model.matrix( ~ 1 +
                      SmokeRate +
                      log(TotPop) + log(PopPerSQM) +
                      log(MedianHValue) + log(MedianHHInc) +
                      PctUrban + PctHighSchool + PctPoor + PctFemale +
                      PctMovedIn5 + PctOccupied + PctWhite + PctBlack + PctHisp,
                    data = data)
Xz <- scale(Xz)
Xz[, 1] <- 1
options(na.action = saved.na.action$na.action)

knots <- construct.knots(data$PCT_HHNV1MI, max.knots = n.knots)
Xy.spline <- expand.spline(data$PCT_HHNV1MI, degree, knots, bound = 10)
z.test.spline <- expand.spline(z.test, degree, knots)
colnames(Xy.spline) <- paste0("z", colnames(Xy.spline))
Xy <- cbind(Xy[, !c(colnames(Xy) == "z")], Xy.spline)
penalized <- str_detect(colnames(Xy), "\\.KNOT\\.")

burn.in <- 5000
n.iters <- 10000
keep <- burn.in + (1 : n.iters)

fit.nonspatial <- fit.affine(y = y, z = z, Xz = Xz, Xy = Xy,
                             offset = offset,
                             seed = 0,
                             CAR.D = CAR.D, CAR.W = CAR.W, NB = NB,
                             print.iters = 100, burn.in = burn.in, n.iters = n.iters,
                             CONSTRAIN_SCALES = FALSE,
                             ESTIMATE_RHO = FALSE,
                             SPATIAL = FALSE)

fit.spatial.unconstrained <- fit.affine(y = y, z = z, Xz = Xz, Xy = Xy,
                             offset = offset,
                             seed = 0,
                             init = fit.nonspatial$final,
                             CAR.D = CAR.D, CAR.W = CAR.W, NB = NB,
                             print.iters = 100, burn.in = burn.in, n.iters = n.iters,
                             CONSTRAIN_SCALES = FALSE,
                             ESTIMATE_RHO = FALSE,
                             SPATIAL = TRUE,
                             USE_PRIOR = TRUE,
                             prop.scale = 1 / 5)

fit.spatial.constrained <- fit.affine(y = y, z = z, Xz = Xz, Xy = Xy,
                                      offset = offset,
                                      seed = 0,
                                      init = fit.nonspatial$final,
                                      CAR.D = CAR.D, CAR.W = CAR.W, NB = NB,
                                      # print.iters = 100, burn.in = 100, n.iters = 1000,
                                      print.iters = 100, burn.in = burn.in, n.iters = n.iters,
                                      CONSTRAIN_SCALES = TRUE,
                                      ESTIMATE_RHO = FALSE,
                                      SPATIAL = TRUE,
                                      USE_PRIOR = TRUE,
                                      COND_PRIOR = TRUE,
                                      CAR.D.shortcut = CAR.D.s,
                                      CAR.W.shortcut = CAR.W.s,
                                      prop.scale = 1 / 5)

fit.affine.unconstrained <- fit.affine(y = y, z = z, Xz = Xz, Xy = Xy,
                                     offset = offset,
                                     seed = 0,
                                     init = fit.spatial.unconstrained$final,
                                     CAR.D = CAR.D, CAR.W = CAR.W, NB = NB,
                                     print.iters = 100, burn.in = burn.in, n.iters = n.iters,
                                     CONSTRAIN_SCALES = FALSE,
                                     ESTIMATE_RHO = TRUE,
                                     SPATIAL = TRUE,
                                     USE_PRIOR = TRUE,
                                     COND_PRIOR = TRUE,
                                     CAR.D.shortcut = CAR.D.s,
                                     CAR.W.shortcut = CAR.W.s,
                                     prop.scale = 1 / 5)

fit.affine.constrained <- fit.affine(y = y, z = z, Xz = Xz, Xy = Xy,
                             offset = offset,
                             seed = 0,
                             init = fit.spatial.constrained$final,
                             CAR.D = CAR.D, CAR.W = CAR.W, NB = NB,
                             print.iters = 100, burn.in = burn.in, n.iters = n.iters,
                             CONSTRAIN_SCALES = TRUE,
                             ESTIMATE_RHO = TRUE,
                             SPATIAL = TRUE,
                             USE_PRIOR = TRUE,
                             COND_PRIOR = TRUE,
                             CAR.D.shortcut = CAR.D.s,
                             CAR.W.shortcut = CAR.W.s,
                             prop.scale = 1 / 5)

save(fit.nonspatial,
     fit.spatial.unconstrained,
     fit.spatial.constrained,
     fit.affine.unconstrained,
     fit.affine.constrained,
     file = "fits-semipar.RData")



### ggplot ###

z.coord <- which(stringr::str_detect(colnames(fit.nonspatial$trace), "^beta.z."))

plot.data <- data.frame()

log.acs <- log(fit.nonspatial$trace[keep, "acs.scalar"]) +
  fit.nonspatial$trace[keep, z.coord] %*% t(z.test.spline)

plot.data <- rbind(plot.data,
                   data.frame(Exposure = z.test,
                              Model = "Non-spatial",
                              Estimate = exp(colMeans(log.acs)),
                              Lower = exp(apply(log.acs, 2,
                                                quantile, prob = 0.025)),
                              Upper = exp(apply(log.acs, 2,
                                                quantile, prob = 0.975)))
)

log.acs <- log(fit.spatial.unconstrained$trace[keep, "acs.scalar"]) +
  fit.spatial.unconstrained$trace[keep, z.coord] %*% t(z.test.spline)

plot.data <- rbind(plot.data,
                   data.frame(Exposure = z.test,
                              Model = "Spatial unconstrained",
                              Estimate = exp(colMeans(log.acs)),
                              Lower = exp(apply(log.acs, 2,
                                                quantile, prob = 0.025)),
                              Upper = exp(apply(log.acs, 2,
                                                quantile, prob = 0.975)))
                   )

log.acs <- log(fit.affine.constrained$trace[keep, "acs.scalar"]) +
  fit.affine.constrained$trace[keep, z.coord] %*% t(z.test.spline)

plot.data <- rbind(plot.data,
                   data.frame(Exposure = z.test,
                              Model = "Affine constrained",
                              Estimate = exp(colMeans(log.acs)),
                              Lower = exp(apply(log.acs, 2,
                                                quantile, prob = 0.025)),
                              Upper = exp(apply(log.acs, 2,
                                                quantile, prob = 0.975)))
)

rug.data <- rbind(data.frame(Model = "Non-spatial", PCT_HHNV1MI = data$PCT_HHNV1MI),
                  data.frame(Model = "Spatial unconstrained", PCT_HHNV1MI = data$PCT_HHNV1MI),
                  data.frame(Model = "Affine constrained", PCT_HHNV1MI = data$PCT_HHNV1MI))

rate.scale <- overall.crude.rate * 100000

pdf("result-semipar.pdf", width = 7, height = 3)
ggplot(plot.data %>% filter(Exposure <= 20),
       aes(x = Exposure, y = Estimate * rate.scale, group = Model)) +
  geom_line() +
  geom_ribbon(aes(ymin = Lower * rate.scale,
                  ymax = Upper * rate.scale), alpha = 0.3) +
  facet_grid(cols = vars(Model)) +
  ylab("CVD deaths per 100k") +
  xlab("% HHNV1MI") +
  geom_rug(aes(x = PCT_HHNV1MI, y = NULL),
           data = rug.data %>% filter(PCT_HHNV1MI <= 20),
           alpha = 0.1) +
  theme_light() +
  theme(strip.background = element_rect(fill = "grey85"),
        strip.text = element_text(color = "grey10"))
dev.off()

