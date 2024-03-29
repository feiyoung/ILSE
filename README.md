# ILSE

=========================================================================

<!-- badges: start -->

[![](https://www.r-pkg.org/badges/version-ago/ILSE)](https://cran.r-project.org/package=ILSE)
[![](https://cranlogs.r-pkg.org/badges/ILSE?color=orange)](https://cran.r-project.org/package=ILSE)
[![](https://cranlogs.r-pkg.org/badges/grand-total/ILSE?color=orange)](https://cran.r-project.org/package=ILSE)
<!-- badges: end -->

Linear Regression by Iterative Least Square Estimation When Covariates Include Missing Values. In *ILSE* package, we also provide Full Information Maximum Likelihood for Linear Regression *fimlreg* that can handle missing Covariates or missing Response variables. 

Please see our new paper for model details:

[Huazhen Lin, Wei Liu, & Wei Lan (2021). Regression Analysis with individual-specific patterns of missing covariates. Journal of Business & Economic Statistics, 39(1), 179-188.](https://www.tandfonline.com/doi/abs/10.1080/07350015.2019.1635486?needAccess=true&journalCode=ubes20)

# Installation

To install the the packages 'ILSE' from 'Github', firstly, install the 'remotes' package.
```{Rmd}
install.packages("remotes")
remotes::install_github("feiyoung/ILSE")
```
Or install the the packages "ILSE" from 'CRAN'
```{Rmd}
install.packages("ILSE")
```

## Usage
For usage examples and guided walkthroughs, check the `vignettes` directory of the repo. 

* [ILSE for simulated data](https://feiyoung.github.io/ILSE/articles/ILSE.Simu.html)
* [ILSE for a toy real data](https://feiyoung.github.io/ILSE/articles/ILSE.Realdata.html)


# Website of ILSE package

We set up a package website to illustrate the usage of this package. For examples of typical ILSE usage, please see our [Package Website](https://feiyoung.github.io/ILSE/index.html) for a demonstration and overview of the functions included in ILSE.
