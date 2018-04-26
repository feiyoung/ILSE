# ILSE
Corrected Least Square Estimation or Full Information Maximum Likelihood Estimation for Linear Regression When Data Include Missing Values.

## Installation
 If you are on Windows, you can 
* Install MASS 
* install mvnmle by running:
```
  install.packages("MASS")
  install.packages("mvnmle")
```
* Download the package as a zip file and install package from the  zip file in Rgui.exe or Rstudio.

If you are on Linux/Mac, you can 
* Install [Rtools](http://cran.r-project.org/bin/windows/Rtools/)
* Install [devtools](http://cran.r-project.org/web/packages/devtools/index.html) by running 
```
install.packages("devtools")
library(devtools)
install_github('feiyoung/ILSE')
```
* Running the following R code to use it!
```
  library(ILSE)
  example('ilse')
```