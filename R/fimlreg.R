fimlreg <- function(x, ...) UseMethod("fimlreg")
fimlreg.numeric <- function(Y, X, ...){
  if(is.null(Y) || is.null(X)) stop('"X","Y" must be given
                                    simutaneously!')
  if (is.null(n <- nrow(X)))
    stop("'X' must be a matrix")
  if (n == 0L)
    stop("0 (non-NA) cases")
  if(!is.matrix(X)) X <- as.matrix(X)
  data <- data.frame(cbind(Y,X))
  resEst <- mlest(data, ...)
  p1 <- length(resEst$muhat)
  muhat <- resEst$muhat
  sigmahat <- resEst$sigmahat
  beta2 <- as.vector(sigmahat[1,2:p1]%*% qr.solve(sigmahat[2:p1,2:p1]))
  beta1 <- muhat[1] - as.vector(sigmahat[1,2:p1]%*% qr.solve(sigmahat[2:p1,2:p1])%*% matrix(muhat[2:p1],nrow=p1-1))
  beta <- c(beta1, beta2)
  names(beta)[1] <- 'intercept'
  names(beta)[2:p1] <- names(data)[2:(p1)]
  return(list(beta=beta, iterations=resEst$iterations, stop.code = resEst$stop.code))
}
fimlreg.formula <- function(formula, data=NULL, ...){
  require(stats)
  if(!is.null(data)){
    data <- model.frame(formula = formula, data = data, na.action=NULL)
  }else{
    XYdat <- model.frame(formula = formula, na.action=NULL)
    p <- ncol(XYdat[[2]])
    XYdat <- cbind(XYdat[[1]], XYdat[[2]])
    colnames(XYdat) <- c("Y", paste0('X', 1:p))
    data <- as.data.frame(XYdat)
    formula <- as.formula('Y~.')
  }
  suppressMessages({resEst <- mlest(data, ...)})
  p1 <- length(resEst$muhat)
  muhat <- resEst$muhat
  sigmahat <- resEst$sigmahat
  beta2 <- as.vector(sigmahat[1,2:p1]%*% qr.solve(sigmahat[2:p1,2:p1]))
  beta1 <- muhat[1] - as.vector(sigmahat[1,2:p1]%*% qr.solve(sigmahat[2:p1,2:p1])%*% matrix(muhat[2:p1],nrow=p1-1))
  beta <- c(beta1, beta2)
  names(beta)[1] <- 'intercept'
  names(beta)[2:p1] <- names(data)[2:(p1)]
  res <- list(beta=beta, iterations=resEst$iterations, stop.code = resEst$stop.code,
              formula=formula, data=data)
  class(res) <- 'fiml'
  return(res)

}
