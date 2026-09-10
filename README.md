
<!-- README.md is generated from README.Rmd. Please edit that file -->

# geoaddSAE2

<!-- badges: start -->
<!-- badges: end -->

geoaddSAE2 provides functions for fitting area-level Geoadditive Small 
Area Estimation (SAE) models. The Geoadditive SAE model extends the 
classical Fay-Herriot model by incorporating nonlinear covariate effects 
and spatial variation.

## Authors

Amalyah Rizky Khairunnisa Sutarto, Novi Hidayat Pusponegoro

## Maintainer

Amalyah Rizky Khairunnisa Sutarto <222212489@stis.ac.id>

## Installation

You can install the released version of geoaddSAE2 from
[CRAN](https://CRAN.R-project.org) or find my github repository
[Github](https://github.com/amalyahrizky04/geoaddSAE2)

## Example

The following example illustrates the basic use of geosae().


``` r
#Load the library
library(geoaddSAE2)

#Load the dataset
data(simulated_sae)

#Estimate geoadditive SAE
fit <- geosae(
        data = simulated_sae, 
        formula = y ~ x1 + x2, 
        vardir = vardir, 
        nonlinear = "x3", 
        spatial = c("lat", "lon"), 
        bootstrap = TRUE,
        B = 5
      )

#Estimate geoadditive SAE and compare with Fay-Herriot and Spatial Fay-Herriot
fit <- geosae(
        data = simulated_sae, 
        formula = y ~ x1 + x2, 
        vardir = vardir, 
        nonlinear = "x3", 
        spatial = c("lat", "lon"), 
        compare = TRUE,
        B = 5
      ) 
      
print(fit) 
summary(fit)
```

## References

-   Fay, R. E. and Herriot, R. A. (1979). Estimates of income for small 
    places: An application of James-Stein procedures to census data. 
    Journal of the American Statistical Association, 74(366), 269–277.
-   Ruppert, D. (2002). Selecting the Number of Knots for Penalized Splines. 
    Journal of Computational and Graphical Statistics, 11(4), 735–757. 
-   Kammann, E. E., & Wand, M. P. (2003). Geoadditive models. Journal of the 
    Royal Statistical Society. Series C: Applied Statistics, 52(1), 1–18. 
-   Rao, J. N. K. ., & Molina, Isabel. (2015). Small area estimation. 
    John Wiley & Sons, Inc.
-   Wood, S. N. (2025). Generalized additive models. Annual Review of 
    Statistics and Its Application, 12(1), 497-526.
