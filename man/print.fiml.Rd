\name{print.fiml}
\alias{print.fiml}
\title{Print the Information of FIML Regression Whether Covariates Include Missing Values or Response Include Missing Values.}
\description{
  print method for class "fiml"..
}
\usage{
  ## S3 method for class 'fiml'
  print(object)
}
\arguments{
  \item{object}{an object of class "fiml".}
}

\value{
  print the basic information of fiml estimation and return a list including

  \item{beta}{A named vector of coefficients}
  \item{iterations}{A positive integer, iterative times in total.}
  \item{stop.code}{The stop code returned by nlm.}
}
\seealso{
  print.lm
}
\examples{

  data(nhanes)
  fimllm <- fimlreg(age~., data=nhanes, iterlim= 40)
  a <- print(fimllm)
  a
}
\keyword{print}
