ilse <- function(x, ...) UseMethod("ilse")


# methods(generic.function = ilse)
ilse.numeric <- function(Y, X,bw=NULL, intercept=F, k.type=NULL,K,bw.type='fix.bw',method="Par.cond", max.iter=50, peps=1e-5, feps = 1e-7,infor_output=F){
  if(is.null(Y) || is.null(X)) stop('"X","Y" must be given
                                    simutaneously!')
  if (is.null(n <- nrow(X)))
    stop("'X' must be a matrix")
  if (n == 0L)
    stop("0 (non-NA) cases")
  if(!is.null(bw) && !is.numeric(bw)) stop('"bw" must be NULL or a positive scalar!')
  if(!is.matrix(X)) X <- as.matrix(X)
  if(! intercept) {X <- cbind(1,X); colnames(X)[1] <- 'intercept';}
  p <- ncol(X)
  i.com <- which(complete.cases(X))
  Xmat0 <- X
  Xmat0[is.na(X)] <- 0
  Z2  <- Xmat0  # X^tilde
  IDX <- indx.comp(X)

  cc.coef <- as.numeric(lm(Y~.+0, data.frame(X), subset= i.com )$coef) # complete cases' estimate as initial estimate
  if(any(is.na(cc.coef)))  cc.coef[is.na(cc.coef)] <- 1
  if(is.null(bw)) bw  <- n^(-1/3)* sd(Xmat0%*%matrix(cc.coef,p,1))
  bw <- ifelse(bw >1, bw, 1)
  beta <- cc.coef
  Bmat <- beta
  rss <- rss.old <- sqrt(sum((Y-Xmat0%*%matrix(beta,p,1))^2))
  k <- 0
  if(infor_output==T){
    message("iter:", k, "d.fn:",NA, "d.par:", NA, "\n")
    message("par:", as.vector(round(beta,2)),'\n')
  }
  k <- k+1
  NA_ind <- which(apply(is.na(X), 1, sum) !=0)
  while (k < max.iter){
    W <- Xmat0 %*% beta
    #Z1 <- t(vapply(1:n, kern.est, FUN.VALUE=numeric(p), beta = beta, Xmat = X, Y= Y,  IDX = IDX, bw = bw, K=K, bw.type=bw.type))
    #Z2[Xmat0==0] <- Z1[Xmat0==0]
    Z1 <- X
    for(i in NA_ind){
      temp_Z1 <- kern.est(ind=i, beta=beta, Xmat=X, Y=Y,IDX=IDX,bw=bw,k.type=k.type,K=K, bw.type='fix.bw')
      if(is.null(temp_Z1)){return(NULL)}else{
        Z1[i,] <- temp_Z1
      }
    }
    Z2[is.na(X)] <- Z1[is.na(X)]
    if (method =="Full.cond"){
      xtx <- t(Z1) %*% Z2
      xty <- t(Z1)%*%Y
      beta.new <- solve(xtx)%*%xty
    }
    else if (method =="Par.cond"){
      beta.new <- lm(Y~.+0, data.frame(Z2))$coef

    }
    else if (method =="SGD"){
      set.seed(k + 20)
      shf.indx <- sample(n,n)
      beta.new <- as.numeric(lm(Y~.-1, data.frame(X), subset= i.com )$coef)
      for (i in shf.indx){
        message('i=',i, '\n')
        xi <- Z2[i,]
        message('xi=', xi, '\n')
        message('beta.new=', beta.new, '\n')
        ngradi <- 2*(Y[i]-xi%*%beta.new)*xi
        #print(max(abs(ngradi)))
        if (max(abs(ngradi))< 0.1) break
        else beta.new <- beta.new + c(rep(0.01,p/2),rep(0.03,p/2))*ngradi*(Xmat0[i,]!=0)
      }
      print(max(abs(2*t(Z2)%*%(Y-Z2%*%beta.new))))
    }
    rss.new <- sqrt(sum((Y-Xmat0%*%beta.new)^2))

    d.fn <- abs(rss.new - rss.old) / rss.old
    d.par <- max(abs(beta.new-beta)) / max(abs(beta))

    if(infor_output==T){
      message("iter:", k, "d.fn:",d.fn, "d.par:", d.par, "\n")
      message("coefs:", as.vector(round(beta.new,2)),'\n')
    }
    beta <- beta.new
    Bmat <- rbind(Bmat,matrix(beta.new, nrow=1))
    rss.old <- rss.new
    rss <- c(rss,rss.new)
    if ((method != "SGD" & d.par < peps )| (method != "SGD" & d.fn < feps) ){
      break
    }
    else if (method == "SGD" & max(abs(2*t(Z2)%*%(Y-Z2%*%beta.new)))< 0.5){
      break
    }
    k <- k+1
    if(k==max.iter){warning('algorithm may not converge!')}
  }
  beta.new <- as.vector(beta.new)
  names(beta.new) <- colnames(X)
  finalm <- lm(Y~.+0, data.frame(Z2))
  row.names(Bmat) <- paste0('iter', 0:k)
  res <- list(beta=beta.new, Bmat = Bmat, d.fn=d.fn, d.par=d.par, iterations=k,
              residuals = finalm$residuals, fitted.values=finalm$fitted.values,
              inargs=list(bw=bw, intercept=intercept, k.type=k.type, bw.type=bw.type,
                          K=K,method=method, max.iter=max.iter, peps=peps,
                          feps = feps,infor_output=infor_output))
  return(res)
}

ilse.formula <- function(formula, data=NULL, bw=NULL,  intercept=F, k.type=NULL,K=NULL,
                         bw.type='fix.bw', method="Par.cond", max.iter=50, peps=1e-5, feps = 1e-7,infor_output=F){
  if(is.null(formula)) stop('formula must be given!')
  if (!inherits(formula, "formula"))
    stop("method is only for formula objects")

  if(!is.null(formula)){
    require(stats)
    if(!is.null(data)){
      XYdat <- model.frame(formula = formula, data = data, na.action=NULL)
    }else{
      XYdat <- model.frame(formula = formula, na.action=NULL)
      p <- ncol(XYdat[[2]])
      XYdat <- cbind(XYdat[[1]], XYdat[[2]])
      colnames(XYdat) <- c("Y", paste0('X', 1:p))
      data <- as.data.frame(XYdat)
      formula <- as.formula('Y~.')
    }
    np <- dim(XYdat)
    XYdat <- as.matrix(XYdat)
    Y <- XYdat[,1]
    X <- XYdat[,-1]
    n <- length(Y)
    X <- matrix(X, nrow=np[1], ncol=np[2]-1)
    colnames(X) <- colnames(XYdat)[2:np[2]]
  }
  cl <- match.call()
  cl[[1]] <- as.name("clse")
  if(!is.matrix(X)) X <- as.matrix(X)
  res <- ilse.numeric(Y, X,bw=bw, intercept=intercept, k.type=k.type, bw.type=bw.type,
                      K=K,method=method, max.iter=max.iter, peps=peps,
                      feps = feps,infor_output=infor_output)
  res$call <- cl
  res$formula <- formula
  res$data <- data
  class(res) <- 'ilse'
  return(res)
}
