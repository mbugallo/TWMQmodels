
### QUANTILE REGRESSION FITTING ALGORITHM WITH WEIGHTS - IRLS

QRLMweights <- function (x, y, case.weights = rep(1, nrow(x)), var.weights = rep(1, nrow(x)), ...,
                         w = rep(1, nrow(x)), init = "ls", psi = psi.huber, scale.est = c("MAD", 
                         "Huber", "proposal 2"), k2 = 1.345, method = c("M", "MM"), maxit = 20, 
                         acc = 1e-04, test.vec = "resid", q = 0.5, w1)
{
  irls.delta <- function(old, new) sqrt(sum((old - new)^2)/max(1e-20, sum(old^2)))
  irls.rrxwr <- function(x, w, r) {
    w <- sqrt(w)
    max(abs((matrix(r*w,1,length(r)) %*% x)/sqrt(matrix(w,1,length(r)) %*% (x^2))))/sqrt(sum(w*r^2))
  }
  method <- match.arg(method)
  nmx <- deparse(substitute(x))
  if (is.null(dim(x))) {
    x <- as.matrix(x)
    colnames(x) <- nmx
  }
  else x <- as.matrix(x)
  if (is.null(colnames(x))) 
    colnames(x) <- paste("X", seq(ncol(x)), sep = "")
  if (qr(x)$rank < ncol(x)) 
    stop("x is singular: singular fits are not implemented in rlm")
  if (!(any(test.vec == c("resid", "coef", "w", "NULL")) || is.null(test.vec))) 
    stop("invalid testvec")
  if (length(var.weights) != nrow(x)) 
    stop("Length of var.weights must equal number of observations")
  if (any(var.weights < 0)) 
    stop("Negative var.weights value")
  if (length(case.weights) != nrow(x)) 
    stop("Length of case.weights must equal number of observations")
  w <- (w * case.weights)/var.weights
  if (method == "M") {
    scale.est <- match.arg(scale.est)
    if (!is.function(psi)) 
      psi <- get(psi, mode = "function")
    arguments <- list(...)
    if (length(arguments)) {
      pm <- pmatch(names(arguments), names(formals(psi)), nomatch = 0)
      if (any(pm == 0)) 
        warning(paste("some of ... do not match"))
      pm <- names(arguments)[pm > 0]
      formals(psi)[pm] <- unlist(arguments[pm])
    }
    if (is.character(init)) {
      if (init == "ls") 
        temp <- lm.wfit(x, y, w, method = "qr")
      else if (init == "lts") 
        temp <- lqs.default(x, y, intercept = FALSE, nsamp = 200)
      else stop("init method is unknown")
      coef <- temp$coef
      resid <- temp$resid
    }
    else {
      if (is.list(init)) 
        coef <- init$coef
      else coef <- init
      resid <- y - x %*% coef
    }
  }
  else if (method == "MM") {
    scale.est <- "MM"
    temp <- lqs.default(x, y, intercept = FALSE, method = "S", k0 = 1.548)
    coef <- temp$coef
    resid <- temp$resid
    psi <- psi.bisquare
    if (length(arguments <- list(...))) 
      if (match("c", names(arguments), nomatch = FALSE)) {
        c0 <- arguments$c
        if (c0 > 1.548) {
          psi$c <- c0
        }
        else warning("c must be at least 1.548 and has been ignored")
      }
    scale <- temp$scale
  }
  else stop("method is unknown")
  done <- FALSE
  conv <- NULL
  n1 <- nrow(x) - ncol(x)
  if (scale.est != "MM") 
    scale <- mad(resid/sqrt(var.weights), 0)
  theta <- 2 * pnorm(k2) - 1
  gamma <- theta + k2^2 * (1 - theta) - 2 * k2 * dnorm(k2)
  qest <- matrix(0, nrow = ncol(x), ncol = length(q))
  qwt <- matrix(0, nrow = nrow(x), ncol = length(q))
  qfit <- matrix(0, nrow = nrow(x), ncol = length(q))
  qres <- matrix(0, nrow = nrow(x), ncol = length(q))
  for(i in 1:length(q)) {
    for (iiter in 1:maxit) {
      if (!is.null(test.vec)) 
        testpv <- get(test.vec)
      if (scale.est != "MM") {
        if (scale.est == "MAD") 
          scale <- median(abs(resid/sqrt(var.weights)))/0.6745
        else scale <- sqrt(sum(pmin(resid^2/var.weights,(k2*scale)^2))/(n1*gamma))
        if (scale == 0) {
          done <- TRUE
          break
        }
      }
      w <- psi(resid/(scale * sqrt(var.weights))) * case.weights
      ww <- 2 * (1 - q[i]) * w
      ww[resid > 0] <- 2 * q[i] * w[resid > 0]
      w <- ww*diag(w1)
      temp <- lm.wfit(x, y, w, method = "qr")
      coef <- temp$coef
      resid <- temp$residuals
      if (!is.null(test.vec)) 
        convi <- irls.delta(testpv, get(test.vec))
      else convi <- irls.rrxwr(x, wmod, resid)
      conv <- c(conv, convi)
      done <- (convi <= acc)
      if (done) 
        break
    }
    if (!done) 
      warning(paste("rlm failed to converge in", maxit, "steps at q = ", q[i]))
    qest[, i] <- coef
    qwt[, i] <- w
    qfit[, i] <- temp$fitted.values
    qres[,i] <- resid
  }
  list(fitted.values = qfit, residuals = qres, q.values = q, q.weights = qwt, coefficients = qest)
}

# COMPUTING OF THE QUANTILE-ORDERS
"zerovalinter"<-function(y, x)
{
  if(min(y) > 0) {
    xmin <- x[y == min(y)]
    if(length(xmin) > 0)
      xmin <- xmin[length(xmin)]
    xzero <- xmin
  }
  else {
    if(max(y) < 0) {
      xmin <- x[y == max(y)]
      if(length(xmin) > 0)
        xmin <- xmin[1]
      xzero <- xmin
    }
    else {
      y1 <- min(y[y > 0])
      if(length(y1) > 0)
        y1 <- y1[length(y1)]
      y2 <- max(y[y < 0])
      if(length(y2) > 0)
        y2 <- y2[1]
      x1 <- x[y == y1]
      if(length(x1) > 0)
        x1 <- x1[length(x1)]
      x2 <- x[y == y2]
      if(length(x2) > 0)
        x2 <- x2[1]
      xzero <- (x2 * y1 - x1 * y2)/(y1 - y2)
      xmin <- x1
      if(abs(y2) < y1)
        xmin <- x2
    }
  }
  resu <-  xzero
  resu
}

# LINEAR INTERPOLATION FUNCTION
# It assumes that the "zerovalinter" function has been already loaded

"gridfitinter"<-function(y,expectile,Q)
  # Computing of the expectile-order of each observation of y by linear interpolation
{
  nq<-length(Q)
  diff <- y %*% t(as.matrix(rep(1, nq))) - expectile        
  vectordest <- apply(diff, 1, zerovalinter,Q)    
}

# BORROWING STRENGTH FROM TIME: TWMQ WEIGHTS

compute.wt <- function(x.s, y.s, regioncode.s, timecode.s, mod.SAE, x.r){
  m <- length(unique(regioncode.s))
  T <- length(unique(timecode.s))
  p <- dim(x.s)[2]
  n <- dim(x.s)[1]
  
  # Predicted values
  pred.st<-list()
  for (j in 1:m){ pred.st[[j]]<-list()
  for (t in 1:T){ pred.st[[j]][[t]]<-x.s[regioncode.s==j & timecode.s==t,]%*%mod.SAE$coef[,j] }
  }
  
  # Time series analysis
  serie<- ts(aggregate(y.s-unlist(pred.st), by=list(regioncode.s,timecode.s), mean)[,3],frequency=m)
  serie<-serie+abs(min(serie))+1
  
  fit.Arima<-auto.arima(serie, max.p=0, max.q=0, max.Q=0, max.d=0, seasonal=TRUE, stepwise=TRUE, 
                        approximation=FALSE, trace=FALSE) 
  P<-fit.Arima$arma[3]; P
  
  # Definition of time weights
  if (P>=1){
    coefARabs<-as.numeric(abs(fit.Arima$coef[P:1])/sum(abs(fit.Arima$coef[P:1])))
    w<-ww<-list()
    for (t in 1:T){
      if(t <= P){ w[[t]]<-c(coefARabs[1:t],rep(0, T-t))/sum(coefARabs[1:t]) } 
      else {	w[[t]]<-c(rep(0,t-P),coefARabs,rep(0, T-t)) }
      ww[[t]]<-list()
      for (tt in 1:T){ ww[[t]][[tt]]<-rep(w[[t]][tt], each=n.t[tt]) }
      ww[[t]]<-unlist(ww[[t]]) }
    
    order.s<-with(df.s, order(df.s$timecode))
    order.r<-with(df.r, order(df.r$timecode))
    
    x.s<-x.s[order.s, ]
    y.s<-y.s[order.s]
    timecode.s<-timecode.s[order.s]
    regioncode.s<-regioncode.s[order.s]
    
    x.r<-x.r[order.r, ]
    timecode.r<-timecode.r[order.r]
    regioncode.r<-regioncode.r[order.r]
    
    W<-matrix(ncol=n, nrow=n)
    for (i in 1:n){ W[i,]<-ww[[timecode.s[i]]] }
  }
  
  return(list(W = W, P = P))
}

# MSE ESTIMATION FOR PREDICTORS DERIVED FROM TWMQ MODELS

# Valid function for the application to real data
TWMQ.mse<-function(x.s, y.s, x.r, QRLMweights.model, QRLMweights50.model, Ndt.matrix, ndt.matrix, 
                   regioncode.s, timecode.s, regioncode.r, timecode.r, P, grid){
  
  D <- length(unique(regioncode.s))
  T <- length(unique(timecode.s))
  p <- dim(x.s)[2]
  n <- dim(x.s)[1]
  
  nd.t<-rowSums(ndt.matrix)
  
  TWMQ<-TWMQ.BC<-mse1<-mse2<-mse3<-mse4<-mse5<-mse6<-var1<-var2<-
        var3<-var4<-bias<-var.beta.prod2<-opt.var3<-opt.var4<-list()
  k.dt<-i.dt<-matrix(nrow=D,ncol=T)
  
  fun.min1<-fun.min2<-fun.dif<-phi<-list()
  for (i in 1:length(grid)){ fun.min1[[i]]<-fun.min2[[i]]<-fun.dif[[i]]<-phi[[i]]<-matrix(nrow=D, ncol=T) }
  
  for (d in 1:D){
    TWMQ[[d]]<-TWMQ.BC[[d]]<-mse1[[d]]<-mse2[[d]]<-mse3[[d]]<-mse4[[d]]<-
      mse5[[d]]<-mse6[[d]]<-var1[[d]]<-var2[[d]]<-var3[[d]]<-var4[[d]]<-
      bias[[d]]<-var.beta.prod2[[d]]<-opt.var3[[d]]<-opt.var4[[d]]<-list()
    
    for (t in 1:T){
      if ((t == 1)|(P==0)){ incl.t = t } else { incl.t = t-P+1; if (incl.t < 1){ incl.t=1 }}
      
      regioncode.s.t<-regioncode.s[timecode.s %in% incl.t:t]
      timecode.s.t<-timecode.s[timecode.s %in% incl.t:t]
      x.s.t<-x.s[timecode.s %in% incl.t:t, ]
      y.s.t<-y.s[timecode.s %in% incl.t:t]
      
      regioncode.r.t<-regioncode.r[timecode.r %in% incl.t:t]
      timecode.r.t<-timecode.r[timecode.r %in% incl.t:t]
      x.r.t<-x.r[timecode.r %in% incl.t:t, ]
      
      nu.t<-length(timecode.s.t)
      
      wei.out<-(Ndt.matrix[d,t]-ndt.matrix[d,t])/(nd.t[d]-1)
      xdt.s.mean<-colMeans(x.s.t[regioncode.s.t==d & timecode.s.t==t, ])
      xdt.r.sum<-colSums(x.r.t[regioncode.r.t==d & timecode.r.t==t, ])
      xdt.r.mean<-colMeans(x.r.t[regioncode.r.t==d & timecode.r.t==t, ])
      
      weig.dt<-diag(QRLMweights.model[[t]]$q.weights[,d])
      uu.dt<-weig.dt%*%x.s.t%*%solve(t(x.s.t)%*%weig.dt%*%x.s.t)%*%xdt.r.sum
      
      ww.dt<-(uu.dt+(timecode.s.t==t & regioncode.s.t==d))
      lambda.dt<- (uu.dt^2+(regioncode.s.t==d & timecode.s.t==t)*wei.out)
      
      pred.st<-x.s.t[timecode.s.t==t & regioncode.s.t==d,]%*%QRLMweights.model[[t]]$coef[,d]
      pred.rt<-x.r.t[timecode.r.t==t & regioncode.r.t==d,]%*%QRLMweights.model[[t]]$coef[,d]
      res.st<-y.s.t[timecode.s.t==t & regioncode.s.t==d]-pred.st
      s<-fun.MAD(res.st)
      
      # MQ predictor
      TWMQ[[d]][[t]]<-c(ww.dt)%*%y.s.t/(Ndt.matrix[d,t])
      
      # Variance of beta 
      res.d.t<-y.s.t-x.s.t%*%QRLMweights.model[[t]]$coef[,d]
      s.d.t<-fun.MAD(res.d.t)
      
      var.beta<-nu.t^2*s.d.t^2/(nu.t-p)*sum(hub.psi(res.d.t/s.d.t, k=1.345)^2)/
        (sum(der.hub.psi(res.d.t/s.d.t, k=1.345))^2)*solve(t(x.s.t)%*%x.s.t)
      
      suma.var1=suma.var2=suma.var4=suma.bias=0
      
      for (d.k in 1:D){
        index <- (regioncode.s.t==d.k)
        
        # MQ predictor - Variance 1
        suma.var1=suma.var1+sum(lambda.dt[index]*(y.s.t[index]-
                        x.s.t[index,]%*%QRLMweights50.model[[t]]$coef)^2)
        
        suma.var2=suma.var2+sum(lambda.dt[index]*(y.s.t[index]-
                        x.s.t[index,]%*%QRLMweights.model[[t]]$coef[,d.k])^2)
        
        # MQ predictor - Variance 2 and MQ\BC predictor Variance
        suma.var4=suma.var4+sum((y.s.t[index & timecode.s.t==t]-x.s.t[index & timecode.s.t==t,]
                          %*%QRLMweights.model[[t]]$coef[,d.k])^2)
        
        # MQ predictor - Bias
        suma.bias=suma.bias+sum(ww.dt[index]*x.s.t[index,]%*%QRLMweights.model[[t]]$coef[,d.k])
      }
      
      # MQ predictor - Variance 1
      var1[[d]][[t]]<-(1/Ndt.matrix[d,t])^2*suma.var1
      var2[[d]][[t]]<-(1/Ndt.matrix[d,t])^2*suma.var2
      
      # MQ predictor - Variance 2
      var.Nn<-(1-ndt.matrix[d,t]/Ndt.matrix[d,t])^2
      var.beta.prod<-var.Nn*xdt.r.mean%*%var.beta%*%xdt.r.mean
      
      suma.var3=sum(res.st^2)
      opt.var3[[d]][[t]]<-var.Nn/((Ndt.matrix[d,t]-ndt.matrix[d,t])*(nd.t[d]-1))*suma.var3
      var3[[d]][[t]]<-var.beta.prod+opt.var3[[d]][[t]]
      opt.var4[[d]][[t]]<-var.Nn/((Ndt.matrix[d,t]-ndt.matrix[d,t])*(n-D))*suma.var4
      var4[[d]][[t]]<-var.beta.prod+opt.var4[[d]][[t]]
      
      # MQ predictor - Bias
      bias[[d]][[t]]<-1/Ndt.matrix[d,t]*(suma.bias-(sum(pred.st)+sum(pred.rt)))
      
      # MQ predictor - MSE
      mse1[[d]][[t]]<-var1[[d]][[t]]+bias[[d]][[t]]^2
      mse2[[d]][[t]]<-var2[[d]][[t]]+bias[[d]][[t]]^2
      mse3[[d]][[t]]<-var3[[d]][[t]]+bias[[d]][[t]]^2
      mse4[[d]][[t]]<-var4[[d]][[t]]+bias[[d]][[t]]^2
      
      # MQ\RB predictor definitions
      Nn.rob<-(Ndt.matrix[d,t]-ndt.matrix[d,t])/(Ndt.matrix[d,t]*ndt.matrix[d,t])
      var.beta.prod2[[d]][[t]]<-var.Nn*(xdt.r.mean-xdt.s.mean)%*%var.beta%*%(xdt.r.mean-xdt.s.mean)
      
      # Optimal calculation of the robustness parameter
      c.prob<-rep(NA,length(grid))
      for (i in 1:length(grid)){
        fun.min1[[i]][d,t]<-sum(hub.psi(res.st/s, k=grid[i])^2)*(s/ndt.matrix[d,t])^2*var.Nn
        phi[[i]][d,t]<-Nn.rob*sum(s*hub.psi(res.st/s, k=grid[i]))
        fun.min2[[i]][d,t]<-(bias[[d]][[t]] + phi[[i]][d,t])^2
        fun.dif[[i]][d,t]<-fun.min2[[i]][d,t]+fun.min1[[i]][d,t]
        c.prob[i]<-fun.dif[[i]][d,t]
      } 
      k.dt[d,t]<-grid[which.min(abs(c.prob))]
      i<-(which(grid==k.dt[d,t]))
      
      # MQ\RB predictor
      TWMQ.BC[[d]][[t]]<-TWMQ[[d]][[t]] + phi[[i]][d,t]
      
      # MQ\BC predictor variance (MSE) and bias
      mse5[[d]][[t]]<-fun.min2[[i]][d,t]+
        fun.min1[[i]][d,t]+var.beta.prod2[[d]][[t]]+opt.var3[[d]][[t]]
      mse6[[d]][[t]]<-fun.min2[[i]][d,t]+
        fun.min1[[i]][d,t]+var.beta.prod2[[d]][[t]]+opt.var4[[d]][[t]]
      
    }
}    

  return(list(TWMQ=unlist(TWMQ), TWMQ.BC=unlist(TWMQ.BC), mse1=unlist(mse1), 
              mse2=unlist(mse2), mse3=unlist(mse3), mse4=unlist(mse4), 
              bias.sq=unlist(bias), mse5.BC=unlist(mse5), mse6.BC=unlist(mse6), k.dt=k.dt))
}

# General function that calculates the MSE estimator for BTMQ predictors 
# proposed in Section A of the Supplementary Material
TWMQ.mse2<-function(x.s, y.s, x.r, QRLMweights.model, Ndt.matrix, ndt.matrix, 
                    regioncode.s, timecode.s, regioncode.r, timecode.r, P, mqo, c.dt){
  
  D <- length(unique(regioncode.s))
  T <- length(unique(timecode.s))
  p <- dim(x.s)[2]
  n <- dim(x.s)[1]
  
  mse <- bias <- list()
  for (d in 1:D){
    mse[[d]]<-bias[[d]]<-list()
    
    for (t in 1:T){
      if (P==0){ incl.t = t } else { incl.t = t-P+1; if (incl.t < 1){ incl.t=1 }}
      
      regioncode.s.t<-regioncode.s[timecode.s %in% incl.t:t]
      timecode.s.t<-timecode.s[timecode.s %in% incl.t:t]
      x.s.t<-x.s[timecode.s %in% incl.t:t, ]
      y.s.t<-y.s[timecode.s %in% incl.t:t]
      
      regioncode.r.t<-regioncode.r[timecode.r %in% incl.t:t]
      timecode.r.t<-timecode.r[timecode.r %in% incl.t:t]
      x.r.t<-x.r[timecode.r %in% incl.t:t, ]
      
      nu.t<-length(timecode.s.t)
      
      # Variance of beta 
      res.d.t<-y.s.t-x.s.t%*%QRLMweights.model[[t]]$coef[,d]
      s.d.t<-fun.MAD(res.d.t)
      
      var.beta<-nu.t^2*s.d.t^2/(nu.t-p)*sum(hub.psi(res.d.t/s.d.t, k=1.345)^2)/
        (sum(der.hub.psi(res.d.t/s.d.t, k=1.345))^2)*solve(t(x.s.t)%*%x.s.t)
      
      # Atypical and non-atypical subsets
      res.d.t.1<-res.d.t[regioncode.s.t==d & timecode.s.t==t]
      index.d.t.G<-which(abs(res.d.t.1/s.d.t) < c.dt[d,t])
      
      # Estimation of V1
      # TWMQ models
      mod.aux<-QRLMweights(x=x.s.t, y=y.s.t, q=tau, maxit=25, k = 1.345, 
                           w1=W[(n.ts[incl.t]+1):(n.ts[t+1]), (n.ts[incl.t]+1):(n.ts[t+1])])
      qo<-matrix(c(gridfitinter(y.s.t,mod.aux$fitted.values,mod.aux$q.values)), ncol=1)
      qo.s.t<-qo[regioncode.s.t==d & timecode.s.t==t]

      V1.1<-(1-ndt.matrix[d,t]/Ndt.matrix[d,t])^2/ndt.matrix[d,t]^2*var(qo.s.t)*sum(res.d.t.1[index.d.t.G]^2)/mean((qo.s.t[index.d.t.G]-mqo[d])^2)
      V1.1[is.na(V1.1)] <- 0
      V1.2<-(1-ndt.matrix[d,t]/Ndt.matrix[d,t])/(ndt.matrix[d,t]*Ndt.matrix[d,t])*var(qo.s.t)*sum(res.d.t.1^2)/mean((qo.s.t-mqo[d])^2)
      
      # Estimation of V2 and V3
      xdt.s.G.beta<-sum(apply(matrix(x.s.t[regioncode.s.t==d & timecode.s.t==t, ][index.d.t.G,], ncol=p), 1, function(x){ x%*%var.beta%*%x }))
      xdt.r.beta<-sum(apply(x.r.t[regioncode.r.t==d & timecode.r.t==t, ], 1, function(x){ x%*%var.beta%*%x }))
      
      V2<-V3<-(1-ndt.matrix[d,t]/Ndt.matrix[d,t])^2/(ndt.matrix[d,t]^2)*xdt.s.G.beta + 1/Ndt.matrix[d,t]^2*xdt.r.beta
      
      # Estimation of B1
      B1.1<-c.dt[d,t]*sum(sign(res.d.t.1[-index.d.t.G]))
      B1.2<-1/(2*s.d.t)*sum(res.d.t.1[index.d.t.G]^2)
      
      bias[[d]][[t]]<-(1-ndt.matrix[d,t]/Ndt.matrix[d,t])/ndt.matrix[d,t]*(B1.1+B1.2)
      mse[[d]][[t]]<-V1.1+V1.2+V2+V3+bias[[d]][[t]]^2
      
    }  
  }
  return(list(mse=unlist(mse), bias.sq=unlist(bias)))
}  


