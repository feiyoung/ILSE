\name{summary.ilse}
\alias{summary.ilse}
\title{Summarizing Iterative Least Square Estimation Regression When Covariates Include Missing Values}
\description{
  summary method for class "ilse"
}
\usage{
  ## S3 method for class 'ilse'
  summary(object, Nbt=100)
  ##
  Fitted.values(object)
  ##
  Residuals(object)

}
\arguments{
  \item{object}{an object of class "ilse".}
  \item{Nbt}{an positive integer, the repeated times of bootstrap to eatimate covariance matrix of regression coefficient.}
}

\value{
  The function summary.ilse computes and returns a named matrix of summary statistics of the fitted linear model given in object by corrected least square estimation.
  The function Fitted.values return a vector, fitted repsonse vlaues.
  The function Residuals return a vector, residuals.
}
\seealso{
summary.lm
fitted.vlaues
residuals
}
\examples{
# example one
data(nhanes)
NAlm <- ilse(age~., data=nhanes)
summary(NAlm, Nbt=50)
}
\keyword{summary}
