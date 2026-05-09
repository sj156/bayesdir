# bayesdir

`bayesdir` is an experimental R package for Bayesian-style directional data
analysis.

The long-term goal is to support predictive resampling and martingale posterior
methods for circular and spherical models.

> **Status:** very early development. Current functions are mostly placeholders,
> and the API may change.

## Installation

The package is not on CRAN. Install the development version from GitHub with:

```r
install.packages("remotes")
remotes::install_github("sj156/bayesdir")
```

## Quick example

```r
library(bayesdir)

# Simulate directions on the unit circle
x <- r_unif_sphere(50, p = 2, seed = 1)

# Fit a placeholder directional model
fit <- fit_bayesdir(x)
fit

# Draw simple posterior resamples
post <- posterior_resample(fit, n_draws = 200, seed = 2)
summary(post)

# Basic plots
plot(fit)
plot(post)
```

## Current features

At this stage, `bayesdir` includes simple placeholder tools for:

- working with unit-vector directional data;
- converting circular angles to Cartesian directions;
- simulating uniform directions;
- fitting a basic placeholder directional model;
- simple posterior and predictive resampling.

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
