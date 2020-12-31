\name{fimlreg}
\alias{fimlreg}
\title{Full Information Maximum Likelihood Linear Regression Whether Covariates Include Missing Values or Response Include Missing Values.}
\description{
  Estimate regression coefficients based on Full Information Maximum Likelihood Estimation, which can couple missing data, including response missing or covariates missing.
}
\usage{
  # S3 method for 'formula'
  fimlreg(formula, data, ...)
}
\arguments{
  \item{formula}{an object of class "formula" (or one that can be coerced to that class): a symbolic description of the model to be fitted. The details of model specification are given under 'Details'.}
  \item{data}{an optional data frame, list or environment (or object coercible by as.data.frame to a data frame) containing the variables in the model. If not found in data, the variables are taken from environment(formula), typically the environment from which clse is called.}
  \item{...}{Optional arguments to be passed to the stats::nlm optimization routine.}
}
\details{
  Note that arguments ... of stats::nlm are the parameters of algorithm, see the details in help file of "nlm". "fimlreg" can cople with any type of missing data.
}
\value{
  Return a list including following components:
  \item{beta}{A named vector of coefficients}
  \item{iterations}{A positive integer, iterative times in total.}
  \item{stop.code}{The stop code returned by nlm.}
}
\author{
  Liu Wei
}
\seealso{
  mvnmle package
}
\examples{
data(nhanes)
## example one: include missing value
fimlreg(age~., data=nhanes, iterlim= 40)

# example two: No missing vlaue
## example two: No missing value
n <- 100
group <- rnorm(n, sd=4)
weight <- 3.2*group + 1.5 + rnorm(n, sd=0.1)
fimllm <- fimlreg(weight~group, data=data.frame(weight=weight, group=group))
fimllm
}
\keyword{Linear Model, Full Information Maximum Likelihood.}
