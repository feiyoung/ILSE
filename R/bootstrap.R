bootstrap <- function(x, ...) UseMethod("bootstrap")

bootstrap.ilse <- function(obj,repTimes=100){
  data <- obj$data
  formula <- obj$formula
  bw <- obj$inargs$bw
  intercept <- obj$inargs$intercept
  k.type<- obj$inargs$k.type
  bw.type <- obj$inargs$bw.type
  K <- obj$inargs$K
  method <- obj$inargs$method
  max.iter <- obj$inargs$max.iter
  peps <- obj$inargs$peps
  feps  <- obj$inargs$ feps
  infor_output <- F
  if(!intercept){
    real_p <- ncol(data)
  }else{
    real_p <- ncol(data)-1
  }
  n <- nrow(data)
  message('===================Start bootstrapping================\n')
  res.par <- matrix(nrow=repTimes, ncol=real_p)
  for(k in 1:repTimes)
  {
    set.seed(k)
    ind <- sample(1:n, n, replace = T)
    data1 <- data[ind, ]
    disProBar(k, repTimes)
    try(coef.par <- ilse(formula, data1, bw,  intercept, k.type,K,
                           bw.type , method, max.iter,
                           peps, feps,infor_output)$beta, silent = T)
    res.par[k, ] <- coef.par
  }
  message('===================Finish bootstrapping================\n')
  return(cov(res.par))
}

bootstrap.fiml <- function(obj, repTimes=100){
  data <- obj$data
  formula <- obj$formula
  n <- nrow(data)
  p <- ncol(data)
  res.par <- matrix(nrow=repTimes, ncol= p)
  message('===================Start bootstrapping================\n')
  for(k in 1:repTimes)
  {
    set.seed(k)
    ind <- sample(1:n, n, replace = T)
    data1 <- data[ind, ]
    disProBar(k, repTimes)
    try(coef.par <- fimlreg(formula, data1)$beta, silent = T)
    res.par[k, ] <- coef.par
  }
  message('===================Finish bootstrapping================\n')
  return(cov(res.par))
}
