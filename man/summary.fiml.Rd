\name{summary.fiml}
\alias{summary.fiml}
\title{Summarizing full-information Maximum Likelihood Estimation Regression Whether Covariates Include Missing Values or Response Include Missing Values.}
\description{
  summary method for class "fiml"
}
\usage{
  ## S3 method for class 'fiml'
  summary(object, Nbt=100)
}
\arguments{
  \item{object}{an object of class "fiml".}
  \item{Nbt}{an positive integer, the repeated times of bootstrap to eatimate covariance matrix of regression coefficient.}
}

\value{
  The function summary.fiml computes and returns a named matrix of summary statistics of the fitted linear model given in object by full-information maximum likelihood.
}
\seealso{
  summary.lm
}
\examples{
  # example one
  data(nhanes)
  fimllm <- fimlreg(age~., data=nhanes, iterlim= 40)
  summary(fimllm, Nbt = 50)
}
\keyword{summary}
