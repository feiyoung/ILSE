fimlreg <- function(...) UseMethod("fimlreg")
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
  # require(stats)
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


mlest <- function (data, ...)
{
  data <- as.matrix(data)
  sortlist <- mysort(data)
  nvars <- ncol(data)
  nobs <- nrow(data)
  if (nvars > 50)
    stop("mlest cannot handle more than 50 variables.")
  startvals <- getstartvals(data)
  lf <- getclf(data = sortlist$sorted.data, freq = sortlist$freq)
  mle <- nlm(lf, startvals, ...)
  muhat <- mle$estimate[1:nvars]
  del <- make.del(mle$estimate[-(1:nvars)])
  factor <- solve(del, diag(nvars))
  sigmahat <- t(factor) %*% factor
  list(muhat = muhat, sigmahat = sigmahat, value = mle$minimum,
       gradient = mle$gradient, stop.code = mle$code, iterations = mle$iterations)
}

getstartvals <- function (x, eps = 0.001){
  n <- ncol(x)
  startvals <- double(n + n * (n + 1)/2)
  startvals[1:n] <- apply(x, 2, mean, na.rm = TRUE)
  sampmat <- cov(x, use = "p")
  eig <- eigen(sampmat, symmetric = TRUE)
  realvals <- sapply(eig$values, function(y) ifelse(is.complex(y),
                                                    0, y))
  smalleval <- eps * min(realvals[realvals > 0])
  posvals <- pmax(smalleval, realvals)
  mypdmat <- eig$vectors %*% diag(posvals) %*% t(eig$vectors)
  myfact <- chol(mypdmat)
  mydel <- solve(myfact, diag(n))
  signchange <- diag(ifelse(diag(mydel) > 0, 1, -1))
  mydel <- mydel %*% signchange
  startvals[(n + 1):(2 * n)] <- log(diag(mydel))
  for (i in 2:n) {
    startvals[(2 * n + sum(1:(i - 1)) - i + 2):(2 * n + sum(1:(i -
                                                                 1)))] <- mydel[1:(i - 1), i]
  }
  startvals
}

mysort <- function (x) {
  nvars <- ncol(x)
  powers <- as.integer(2^((nvars - 1):0))
  binrep <- ifelse(is.na(x), 0, 1)
  decrep <- binrep %*% powers
  sorted <- x[order(decrep), ]
  decrep <- decrep[order(decrep)]
  list(sorted.data = sorted, freq = as.vector(table(decrep)))
}

getclf <- function (data, freq)
{
  nvars <- ncol(data)
  pars <- double(nvars + nvars * (nvars + 1)/2)
  testdata <- data[cumsum(freq), ]
  presabs <- ifelse(is.na(testdata), 0, 1)
  data <- t(data)
  presabs <- t(presabs)
  dim(presabs) <- NULL
  dim(data) <- NULL
  data <- data[!is.na(data)]
  function(pars) {
    .C("evallf", as.double(data), as.integer(nvars),
       as.integer(freq), as.integer(x = length(freq)), as.integer(presabs),
       as.double(pars), val = double(1), PACKAGE = "ILSE")$val
  }
}

make.del <- function (pars)
{
  k <- floor((-1 + sqrt(1 + 8 * length(pars)))/2)
  mymatrix <- diag(exp(pars[1:k]))
  pars <- pars[-(1:k)]
  if (k > 1) {
    for (i in 2:k) {
      mymatrix[1:(i - 1), i] <- pars[1:(i - 1)]
      pars <- pars[-(1:(i - 1))]
    }
  }
  mymatrix
}
