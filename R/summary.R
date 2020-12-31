summary <- function(x, ...) UseMethod("summary")
summary.ilse <- function(x, Nbt=100){
  Est <- x$beta
  res <- matrix(0, nrow=length(Est), ncol=4)
  Acov <- bootstrap(x, repTimes = Nbt)
  stdErr <- sqrt(diag(Acov))
  Zvalue <- Est / stdErr
  Pvalue <- 2*(1-pnorm(abs(Zvalue)))
  res[,1] <- Est
  res[,2] <- stdErr
  res[,3] <- Zvalue
  res[,4] <- Pvalue
  row.names(res) <- names(Est)
  colnames(res) <- c('Estimate', 'std. Error', 'Z value', 'Pr(>|Z|)')
  res
}

summary.fiml <- function(x, Nbt=100){
  Est <- x$beta
  res <- matrix(0, nrow=length(Est), ncol=4)
  Acov <- bootstrap(x, repTimes = Nbt)
  stdErr <- sqrt(diag(Acov))
  Zvalue <- Est / stdErr
  Pvalue <- 2*(1-pnorm(abs(Zvalue)))
  res[,1] <- Est
  res[,2] <- stdErr
  res[,3] <- Zvalue
  res[,4] <- Pvalue
  row.names(res) <- names(Est)
  colnames(res) <- c('Estimate', 'std. Error', 'Z value', 'Pr(>|Z|)')
  res
}
print <- function(x, ...) UseMethod("print")
print.ilse <- function(x) print(x[1:5])
print.fiml <- function(x) print(x[1:3])
Coef <- function(x) {
  if(!is.element(class(x), c('ilse', 'fiml')))
    stop('x must be class "ilse" or "fiml"!\n')
  return(x$beta)
}
Fitted.values <- function(x){
  if(!is.element(class(x), c('ilse')))
    stop('x must be class "ilse"!\n')
  return(x$fitted.values)
}
Residuals <- function(x){
  if(!is.element(class(x), c('ilse')))
    stop('x must be class "ilse"!\n')
  return(x$residuals)
}
