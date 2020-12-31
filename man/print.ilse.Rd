\name{print.ilse}
\alias{print.ilse}
\title{Print the Information of Iterative Least Square Estimation Regression When Covariates Include Missing Values}
\description{
  print method for class "ilse"..
}
\usage{
  ## S3 method for class 'ilse'
  print(object)
}
\arguments{
  \item{object}{an object of class "ilse".}
}

\value{
  print the basic information of ilse estimation and algorithm and return a list including

\item{beta}{a named vector of coefficients}
\item{Bmat}{a named matrix that summary the estimated beta in every iteration.}
\item{residuals}{the residuals, that is response minus fitted values.}
\item{fitted.values}{the fitted mean values.}
\item{d.fn}{a nonnegative value, vlaue of relative variation rate of objective function value}
\item{d.par}{a nonnegative value, relative variation rate of estimated parametric vector when algorithm stopped.}
\item{K}{a positive integer, iterative times in total.}
}
\seealso{
print.lm
}
\examples{

data(nhanes)
NAlm1 <- ilse(age~., data=nhanes)
a <- print(NAlm1)
a
}
\keyword{print}
