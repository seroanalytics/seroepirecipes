
## *seroepirecipes*: implementations and tutorials of common models in seroepidemiology <img src="man/figures/logo.png" align="right" width="130"/>

<!-- badges: start -->
<!-- badges: end -->

***seroepirecipes*** implements and links to R packages implementing
commonly used mathematical and statistical models for analyzing
serological data. We implement a range of methods in R vignettes, from
fitting antibody kinetics models to longitudinal antibody titer data,
estimating the force of infection using serocatalytic models, and
inferring infection histories using time-since-infection methods.

## Contents

- Description of datasets included in ***seroepirecipes***
- Fitting a hierarchical antibody kinetics model to longitudinal
  antibody measurements
- Using mixture models to classify serostatus
- Estimating the force of infection using serocatalytic models
- Infection times using reversible jump MCMC
- Inferring infection histories using *serosolver*

## Installation

You can install the development version of ***seroepirecipes*** from
[GitHub](https://github.com/) with:

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
