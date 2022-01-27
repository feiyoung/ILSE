\name{disProBar}
\alias{disProBar}
\title{Display Process in for loop}
\description{
Display running process in for loop.
}
\usage{
disProBar(k, NN)
}
\arguments{
  \item{k}{a positive integer, current value of loop index.}
  \item{NN}{a positive integer, total times of repetition.}
}
\examples{
  for(i in 1:1000){
  i <- i+1
  disProBar(i, 1000)
  }
}
\value{
  No velue returned.
}
