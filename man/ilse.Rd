\name{ilse}
\alias{ilse}
\title{Linear Regression by Iterative Least Square Estimation}
\description{
  Linear regression when covariates include missing values embedding the
    correlation imformation bettwen covariates by Iterative Least Square Estimation.
}
\usage{
  ilse(formula, data,  bw=NULL, intercept=F,k.type=NULL,bw.type='fix.bw',K=NULL,method="Par.cond", max.iter=50, peps=1e-7, feps = 1e-7,infor_output=F)
}
\arguments{
  \item{formula}{an object of class "formula" (or one that can be coerced to that class): a symbolic description of the model to be fitted. The details of model specification are given under 'Details'.}
  \item{data}{an optional data frame, list or environment (or object coercible by as.data.frame to a data frame) containing the variables in the model. If not found in data, the variables are taken from environment(formula), typically the environment from which ilse is called.}
  \item{bw}{a positive value, specify the bandwidth in estimating missing values, default as NULL. When bw=NULL, it is automatically selected by empirical method.}
  \item{intercept}{a optional logical value, indicate whether X or data including intercept term (full-one column), default as 'FALSE'.}
  \item{k.type}{an optional character string, specify the type of kernel used in iterative estimating algorithm and support 'epk', 'biweight', 'triangle', 'guassian', 'triweight', 'tricube', 'cosine', 'uniform' in current version, defualt as 'guassian'.}
  \item{bw.type}{an optional character string, specify the type of kernel bandwidth, only support 'fix.bw'(fixed bandwidth) in current version, may spport variable bandwidth in future version, defualt as 'fix.bw'.}
  \item{K}{an optional positive integer, specify the sample size used in kernel estimates when choose the variable bandwidth.}
  \item{method}{an optional character string, specify the iterative algorithm, only support 'Par.cond' and 'Full.cond' in current version and will support  'SGD' in future version.}
  \item{max.iter}{an optional positive integer, the maximum iterative times, defualt as '50'.}
  \item{peps}{an optional positive value, tolerance vlaue of relative variation rate of estimated parametric vector, defualt as '1e-7'.}
  \item{feps}{an optional positive vlaue, tolerance vlaue of relative variation rate of objective function value, defualt as '1e-7'.}
  \item{infor_output}{an optional logical value, indicate whether output the iterative information, default as 'TRUE'.}
}
\details{
  Models for ilse are specified symbolically. A typical model has the form response ~ terms where response is the (numeric) response vector and terms is a series of terms which specifies a linear predictor for response. A terms specification of the form first + second indicates all the terms in first together with all the terms in second with duplicates removed. A specification of the form first:second indicates the set of terms obtained by taking the interactions of all terms in first with all terms in second. The specification first*second indicates the cross of first and second. This is the same as first + second + first:second.
}
\value{
  ilse returns an object of class "NAlm".

The functions summary and anova are used to obtain and print a summary and analysis of variance table of the results. The generic accessor functions coefficients, effects, fitted.values and residuals extract various useful features of the value returned by lm.

An object of class "NAlm" is a list containing at least the following components:

\item{beta}{a named vector of coefficients}
\item{Bmat}{a named matrix that summary the estimated beta in every iteration.}
\item{d.fn}{a nonnegative value, vlaue of relative variation rate of objective function value}
\item{d.par}{a nonnegative value, relative variation rate of estimated parametric vector when algorithm stopped.}
\item{iterations}{a positive integer, iterative times in total.}
\item{residuals}{the residuals, that is response minus fitted values.}
\item{fitted.values}{the fitted mean values.}
\item{inargs}{a list including all input arguments.}
}
\references{
  Huazhen Lin, Wei Liu, Wei Lan(2019). Regression Analysis with individual-specific
patterns of missing covariates.
}
\author{
Wei Liu
}
\note{
nothing
}
\seealso{
\link{lm}
}
\examples{
## exmaple one: include missing value
data(nhanes)
NAlm1 <- ilse(age~., data=nhanes, method = 'Full.cond')
NAlm1
NAlm2 <- ilse(age~., data=nhanes,bw=1,bw.type = 'fix.bw',method = 'Par.cond', k.type='gaussian', infor_output = T)

## example two: No missing value
n <- 100
group <- rnorm(n, sd=4)
weight <- 3.2*group + 1.5 + rnorm(n, sd=0.1)
NAlm3 <- ilse(weight~group, data=data.frame(weight=weight, group=group), intercept = F)
NAlm3

}
\keyword{Linear Model, Missing Value.}
