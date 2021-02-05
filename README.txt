These directories contain data sets and code needed to reproduce analyses and simulations in the paper "Mitigating unobserved spatial confounding when estimating the effect of supermarket access on cardiovascular disease deaths."

1. analysis: contains data and analysis code
  - Data files:
    - data-access.csv contains county-level food access data
    - data-census.csv contains county-level aggregated demographics from Census data
    - data-mortality.csv contains county-level CVD mortality counts from the CDC
    - data-combined.csv is the result of linking the food access, demographics, and mortality data, for convenience
    - data-adjacency.csv contains county adjacency information used for constructing spatial structures
  - Analysis files:
    - analysis-linear.R fits the Poisson GLM models to the data set
    - analysis-semipar.R fits the penalized spline extension to the data set
    - model-fit-linear.R defines the MCMC sampler for the Poisson GLM models and is called by analysis-linear.R
    - model-fit-semipar.R defines the MCMC sampler for the penalized spline extension and is called by analysis-semipar.R
    - data-compilation.R loads the individual data sets, constructs (and saves) the combined data set, and constructs county spatial structures; it is called by analysis-linear.R and analysis-semipar.R

2. simulation: contains simulation code
  - linear
    - simulation-linear.R runs the Poisson GLM simulations
    - sim-data-linear.R simulates data according to the Poisson GLM and is called by simulation-linear.R
    - model-fit-linear.R defines the MCMC sampler for the Poisson GLM models
  - semipar: analogous to the "lienar" directory for the semiparametric simulation

