
## *seroepirecipes*: implementations and tutorials of common models in seroepidemiology <img src="man/figures/logo.png" align="right" width="130"/>

<!-- badges: start -->
<!-- badges: end -->

***seroepirecipes*** implements and links to R packages implementing
commonly used mathematical and statistical models for analyzing
serological data. The package implements a range of methods in R
vignettes, from fitting antibody kinetics models to longitudinal
antibody titer data, estimating the force of infection using
serocatalytic models, and inferring infection histories using
time-since-infection methods. This codebase accompanies a [literature
review](https://osf.io/kqdsn) of analytical methods for
seroepidemiology.

All of the tutorials use either simulated datasets from the
[serosim](https://seroanalytics.github.io/serosim/) R-package or
publicly available datasets.

## DEVELOPER TO DO

- Include Stan models within source code rather than vignette-specific
  scripts.
- Attach example datasets
- Move the jahR antibody kinetics functions here

## Contents

- (WIP) Description of datasets included in ***seroepirecipes***
- (WIP) Common considerations when preparing serological data for
  analysis
- (WIP) [Fitting a hierarchical antibody kinetics model to longitudinal
  antibody
  measurements](https://seroanalytics.github.io/seroepirecipes/articles/antibody_kinetics_models.html)
- (WIP) Using mixture models to classify serostatus
- [Estimating the force of infection using serocatalytic
  models](https://seroanalytics.github.io/seroepirecipes/articles/serocatalytic_model.html)
- (WIP) Inferring infection histories using *serosolver*

Upcoming vignettes:

- (WIP) Estimating the force of infection using antibody acquisition
  models
- (WIP) Infection times using reversible jump MCMC
- (WIP) Time-since-infection machine learning classifier

## Installation

You can install the development version of ***seroepirecipes*** with:

``` r
if(!require("remotes")) install.packages("remotes")
remotes::install_github("seroanalytics/seroepirecipes")
```

## Issues

Please report bugs and requests on the
[issues](https://github.com/seroanalytics/seroepirecipes) link.

## Contribute

\<\< TBC \>\>

We welcome contributions and ideas for additional vignettes and model
implementations. Please follow the [package contributing
guide](https://github.com/seronalytics/seroepirecipes/blob/main/.github/CONTRIBUTING.md).
