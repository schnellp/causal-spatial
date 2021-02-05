source("data-compilation.R")
source("model-fit-linear.R")

library(table1)
library(gridExtra)
library(viridis)

my.render.cont <- function(x) {
  with(stats.apply.rounding(stats.default(x), digits=2), c("",
                                                           "Q1, Median, Q3"=sprintf("%s, %s, %s", Q1, MEDIAN, Q3)))
}

table1( ~ PctUrban + PctWhite + PctBlack + PctHisp + PctHighSchool +
          MedianHHInc + PctPoor + PctFemale + 
          PctMovedIn5 + PctOccupied + MedianHValue +
          PopPerSQM + TotPop +
          SmokeRate, data = data,
        render.continuous = my.render.cont)


############
### maps ###
############


map.counties <- map_data("county")
parishes <- which(map.counties$region == "louisiana")
map.counties$subregion[parishes] <- paste(map.counties$subregion[parishes], "parish")
map.counties$RSR <- paste0(map.counties$region, map.counties$subregion)
map.counties$RSR.comma <- paste0(map.counties$region, ",", map.counties$subregion)
map.counties$RSR.comma <- str_split(map.counties$RSR.comma, " parish", simplify = TRUE)[, 1]

data(county.fips)
county.fips.aug <- county.fips
county.fips.aug$fips <- str_pad(county.fips.aug$fips, 5, pad = "0")
county.fips.aug$polyname <- str_split(county.fips.aug$polyname, ":", simplify = TRUE)[, 1]
county.fips.aug <- county.fips.aug %>% distinct()

map.counties <- left_join(map.counties, county.fips.aug, by = c("RSR.comma" = "polyname"))

data$RSR = paste0(data$State, data$County)
map.counties$FIPS <- map.counties$fips
choro <- merge(map.counties, data, sort = FALSE, by = "FIPS")
choro <- choro[order(choro$order), ]

theme_set(theme_light())



plot.outcome <- ggplot(choro, aes(long, lat)) +
  geom_polygon(aes(group = group, fill = log10(RelativeRisk))) +
  borders(database = "state", regions = STATES, fill = NA,
          colour = "white") +
  scale_fill_viridis(name = "Relative risk",
                     limits = c(-0.5, 0.4),
                     labels = c("0.4", "0.6", "1.0", "1.6", "2.5"),
                     breaks = log10(c(0.4, 0.6, 1.0, 1.6, 2.5))) +
  coord_map("conic", lat = 30) +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "bottom")




plot.exposure <- ggplot(choro, aes(long, lat)) +
  geom_polygon(aes(group = group, fill = log10(PCT_HHNV1MI))) +
  borders(database = "state", regions = STATES, fill = NA,
          colour = "white") +
  scale_fill_viridis(name = "% HHNV1MI",
                     limits = c(-1, 1.51),
                     labels = c("0.1", "0.3", "1.0", "3.2", "10", "32"),
                     breaks = log10(c(0.1, 0.3, 1.0, 3.2, 10, 32))) +
  coord_map("conic", lat = 30) +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "bottom")



pdf("maps.pdf", width = 7, height = 3)
grid.arrange(plot.exposure, plot.outcome, ncol = 2)
dev.off()



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

y <- data$Deaths
offset <- log(data$Expected)
min.nzz <- min(data$PCT_HHNV1MI[data$PCT_HHNV1MI > 0])
z <- log(pmax(data$PCT_HHNV1MI, min.nzz)) - mean(log(pmax(data$PCT_HHNV1MI, min.nzz)))

saved.na.action <- options('na.action')
options(na.action = 'na.pass')
Xy <- model.matrix( ~ PCT_HHNV1MI +
                      SmokeRate +
                      log10(TotPop) + log10(PopPerSQM) +
                      log10(MedianHValue) + log10(MedianHHInc) +
                      PctUrban + PctHighSchool + PctPoor + PctFemale +
                      PctMovedIn5 + PctOccupied + PctWhite + PctBlack + PctHisp,
                    data = data)
colnames(Xy)[2] <- "z"
Xy <- scale(Xy)
Xy[, 1] <- 1

Xz <- model.matrix( ~ 1 +
                      SmokeRate +
                      log10(TotPop) + log10(PopPerSQM) +
                      log10(MedianHValue) + log10(MedianHHInc) +
                      PctUrban + PctHighSchool + PctPoor + PctFemale +
                      PctMovedIn5 + PctOccupied + PctWhite + PctBlack + PctHisp,
                    data = data)
Xz <- scale(Xz)
Xz[, 1] <- 1
options(na.action = saved.na.action$na.action)


z.test <- seq(0, 25, by = 1)

burn.in <- 1000
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
     file = "fits-linear.RData")
                  

result.mcmc <- matrix(NA, nrow = 5, ncol = 3)
colnames(result.mcmc) <- c("mean", "lower", "upper")
rownames(result.mcmc) <- c("nonspatial",
                           "spatial.unconstrained", "spatial.constrained",
                           "affine.unconstrained", "affine.constrained")


result.mcmc["nonspatial", "mean"] <-
  mean(fit.nonspatial$trace[keep, "beta.z"])
result.mcmc["spatial.unconstrained", "mean"] <- 
  mean(fit.spatial.unconstrained$trace[keep, "beta.z"])
result.mcmc["spatial.constrained", "mean"] <-
  mean(fit.spatial.constrained$trace[keep, "beta.z"])
result.mcmc["affine.unconstrained", "mean"] <-
  mean(fit.affine.unconstrained$trace[keep, "beta.z"])
result.mcmc["affine.constrained", "mean"] <-
  mean(fit.affine.constrained$trace[keep, "beta.z"])

result.mcmc["nonspatial", c("lower", "upper")] <-
  quantile(fit.nonspatial$trace[keep, "beta.z"],
           prob = c(0.025, 0.975))
result.mcmc["spatial.unconstrained", c("lower", "upper")] <-
  quantile(fit.spatial.unconstrained$trace[keep, "beta.z"],
           prob = c(0.025, 0.975))
result.mcmc["spatial.constrained", c("lower", "upper")] <-
  quantile(fit.spatial.constrained$trace[keep, "beta.z"],
           prob = c(0.025, 0.975))
result.mcmc["affine.unconstrained", c("lower", "upper")] <-
  quantile(fit.affine.unconstrained$trace[keep, "beta.z"],
           prob = c(0.025, 0.975))
result.mcmc["affine.constrained", c("lower", "upper")] <-
  quantile(fit.affine.constrained$trace[keep, "beta.z"],
           prob = c(0.025, 0.975))

### transformed ###

round(exp(result.mcmc), 3)
mean(exp(fit.nonspatial$trace[keep, "beta.z"]) > 1)
mean(exp(fit.spatial.unconstrained$trace[keep, "beta.z"]) > 1)
mean(exp(fit.spatial.constrained$trace[keep, "beta.z"]) > 1)
mean(exp(fit.affine.unconstrained$trace[keep, "beta.z"]) > 1)
mean(exp(fit.affine.constrained$trace[keep, "beta.z"]) > 1)

d.nonspatial <- density(exp(fit.nonspatial$trace[keep, "beta.z"]))
d.spatial.unconstrained <- density(exp(fit.spatial.unconstrained$trace[keep, "beta.z"]))
d.spatial.constrained <- density(exp(fit.spatial.constrained$trace[keep, "beta.z"]))
d.affine.unconstrained <- density(exp(fit.affine.unconstrained$trace[keep, "beta.z"]))
d.affine.constrained <- density(exp(fit.affine.constrained$trace[keep, "beta.z"]))


### ggplot ###

ppi <- 1 / attributes(Xy)$`scaled:scale`["z"] # how big of a change in standardized PCT_HHNV1MI?

densities <- rbind(
  data.frame(Model = "Non-spatial", beta.z = exp(ppi * fit.nonspatial$trace[keep, "beta.z"])),
  data.frame(Model = "Spatial unconstrained", beta.z = exp(ppi * fit.spatial.unconstrained$trace[keep, "beta.z"])),
  data.frame(Model = "Spatial constrained", beta.z = exp(ppi * fit.spatial.constrained$trace[keep, "beta.z"])),
  data.frame(Model = "Affine unconstrained", beta.z = exp(ppi * fit.affine.unconstrained$trace[keep, "beta.z"])),
  data.frame(Model = "Affine constrained", beta.z = exp(ppi * fit.affine.constrained$trace[keep, "beta.z"]))
)

pdf("densities.pdf", width = 6, height = 4)
ggplot(densities) +
  stat_density(aes(x = beta.z, linetype = Model, color = Model),
               geom = "line", position = "identity", size = 1) +
  xlab("Relative risk") +
  ylab("Density") +
  theme(legend.position = c(0.6, 0.8),
        legend.key.width = grid::unit(3, "lines"),
        legend.box.background = element_rect(),
        axis.text.y = element_blank()) +
  guides(linetype = guide_legend(ncol = 2))
dev.off()




### diagnostics ###

# exposure

# map check of CAR assumption (no more spatial correlation)

# affine constrained
est.ac <- colMeans(fit.affine.constrained$trace[keep, ])
gamma.ac <- est.ac[which(substr(names(est.ac), 1, 6) == "gamma.")]
H.ac <- est.ac["tau.z"] * (CAR.D - est.ac["phi.z"] * CAR.W)
G.ac <- est.ac["tau.u"] * (CAR.D - est.ac["phi.u"] * CAR.W)
Q.ac <- diag.spam(-est.ac["rho"] * sqrt(diag.spam(G.ac) * diag.spam(H.ac)))
H.marginal <- (H.ac - Q.ac %*% solve(G.ac) %*% Q.ac)
z.ac.res <- H.marginal %*% (z - Xz %*% gamma.ac) / sqrt(diag(H.marginal))
plot(z.ac.res ~ fit.affine.constrained$u.mean)

# non-spatial
est.ns <- colMeans(fit.nonspatial$trace[keep, ])
gamma.ns <- est.ac[which(substr(names(est.ns), 1, 6) == "gamma.")]
z.ns.res <- (z - Xz %*% gamma.ns) * sqrt(est.ns["tau.z"])



data.res <- cbind(data, z.ac.res, z.ns.res)

choro <- merge(map.counties, data.res, sort = FALSE, by = "FIPS")
choro <- choro[order(choro$order), ]

plot.exposure.ns.res <- ggplot(choro, aes(long, lat)) +
  geom_polygon(aes(group = group, fill = pnorm(z.ns.res))) +
  borders(database = "state", regions = STATES, fill = NA,
          colour = "white") +
  scale_fill_viridis(name = "Residual log % HHNV1MI",
                     limits = c(0, 1)) +
  coord_map("conic", lat = 30) +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "bottom")

plot.exposure.ac.res <- ggplot(choro, aes(long, lat)) +
  geom_polygon(aes(group = group, fill = pnorm(z.ac.res))) +
  borders(database = "state", regions = STATES, fill = NA,
          colour = "white") +
  scale_fill_viridis(name = "Residual log % HHNV1MI",
                     limits = c(0, 1)) +
  coord_map("conic", lat = 30) +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "bottom")

grid.arrange(plot.exposure.ns.res, plot.exposure.ac.res, ncol = 2)

# linearity check
plot(z.ac.res ~ Xz %*% gamma.ac)


# constant conditional correlation check
e.cond.z.mu <- numeric(length(z))
for (i in 1 : length(z)) {
  j <- NB[[i]]
  e.cond.z.mu[i] <- Xz[i, ] %*% gamma.ac +
    est.ac["phi.z"] / length(j) * sum(z[j] - Xz[j, ] %*% gamma.ac)
}
u.scaled <- fit.affine.constrained$final$u * sqrt(est.ac["tau.u"] / est.ac["tau.z"])
plot((z - e.cond.z.mu - est.ac["rho"] * u.scaled) ~ u.scaled)
summary(lm((z - e.cond.z.mu - est.ac["rho"] * u.scaled) ~ u.scaled))

# outcome

est.ac <- colMeans(fit.affine.constrained$trace[keep, ])
beta.ac <- est.ac[which(substr(names(est.ac), 1, 5) == "beta.")]
lr.ac <- offset + Xy %*% beta.ac + fit.affine.constrained$final$u
y.ac.res <- (y - exp(lr.ac)) / sqrt(exp(lr.ac))

est.ns <- colMeans(fit.nonspatial$trace[keep, ])
beta.ns <- est.ns[which(substr(names(est.ns), 1, 5) == "beta.")]
lr.ns <- offset + Xy %*% beta.ns
y.ns.res <- (y - exp(lr.ns)) / sqrt(exp(lr.ns))


data.res <- cbind(data, y.ac.res, y.ns.res)

choro <- merge(map.counties, data.res, sort = FALSE, by = "FIPS")
choro <- choro[order(choro$order), ]

plot.outcome.ns.res <- ggplot(choro, aes(long, lat)) +
  geom_polygon(aes(group = group, fill = pnorm(y.ns.res))) +
  borders(database = "state", regions = STATES, fill = NA,
          colour = "white") +
  scale_fill_viridis(name = "Residual mortality",
                     # limits = c(-3, 3)) +
                     limits = c(0, 1)) +
  coord_map("conic", lat = 30) +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "bottom")

plot.outcome.ac.res <- ggplot(choro, aes(long, lat)) +
  geom_polygon(aes(group = group, fill = pnorm(y.ac.res))) +
  borders(database = "state", regions = STATES, fill = NA,
          colour = "white") +
  scale_fill_viridis(name = "Residual mortality",
                     # limits = c(-3, 3)) +
                     limits = c(0, 1)) +
  coord_map("conic", lat = 30) +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "bottom")

grid.arrange(plot.outcome.ns.res, plot.outcome.ac.res, ncol = 2)

plot(y.ac.res ~ lr.ac)


# scale constraint check

diff.unconstrained <- fit.affine.unconstrained$trace[keep, "phi.u"] -
  fit.affine.unconstrained$trace[keep, "phi.z"]

mean(diff.unconstrained)
sd(diff.unconstrained)
quantile(diff.unconstrained, prob = c(0.01))
range(diff.unconstrained)

diff.constrained <- fit.affine.constrained$trace[keep, "phi.u"] -
  fit.affine.constrained$trace[keep, "phi.z"]

mean(diff.constrained)
sd(diff.constrained)
quantile(diff.constrained, prob = c(0.01))
range(diff.constrained)

mean(fit.affine.constrained$trace[keep, "rho"])
quantile(fit.affine.unconstrained$trace[keep, "rho"], prob = c(0.025, 0.975))

cor(fit.affine.constrained$trace[keep, "rho"],
    exp(fit.affine.constrained$trace[keep, "beta.z"]))

ggplot(data.frame(fit.affine.constrained$trace[keep, ]), aes(x = rho, y = exp(beta.z))) +
  geom_density_2d() +
  labs(x = expression(rho),
       y = expression(exp(beta[Z])))

