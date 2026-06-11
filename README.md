# bayesdir

`bayesdir` is an experimental R package for Bayesian-style directional data
analysis.

The long-term goal is to support predictive resampling and martingale posterior
methods for circular and spherical models.

## Installation

The package is not on CRAN. Install the development version from GitHub with:

```r
# Install devtools if not already installed
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}

# Install bdynets from GitHub
devtools::install_github("sj156/bayesdir", subdir = "bayesdir", dependencies = TRUE)
```

## Quick example

```r
library(bayesdir)
library(parallel)

#Generate synthetic data.
n <- 100
mu_true <- 2
kappa_true <- 4
nn_cores <- max(detectCores()-1,1)

# Simulate directions on the unit circle for p = 2.
x2 <- generate_vmf_data(n,mu_true, kappa_true)$X

# Draw simple posterior resamples
res_mart <- run_predresamp_vmf(x2,B=500,M=100,n_cores = nn_cores) ## Most effective one. (Hybrid with pooling)
summary(res_mart$samples$kappa_sam)

# Basic plots
plot(density(kappa_sam))
```
If one would like to show the trace plot for each martingale, please refer 5simulation.Rmd.

## Current features

At this stage, `bayesdir` includes simple placeholder tools for:
- working with unit-vector directional data;
- simulating directions from von-mises fisher distribution for both circular angles and cartesian directions;
- fitting a basic placeholder directional model;
- give a trace plot for each martingale.

For more examples, simulation and real data cases, one may see from 5simulation.Rmd and 6realdata_oscar.Rmd.

## Development

This package is developed using a `litr` workflow.

Main source file:

```text
index.Rmd
```

Supporting literate source files:

```text
Rmd/
```

Typical development workflow:

```r
litr::render("index.Rmd")
devtools::document()
devtools::test()
devtools::check()
```

## License

GPL-3

