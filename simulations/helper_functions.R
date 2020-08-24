# helper functions to calculate class-specific effect estimators (in each class)
IPWestimator <- function(X,Z,Y) {
  ps.fit <- glm(as.formula(paste0("Z~",paste(colnames(X),collapse="+"))), 
                family=binomial("logit"), data=cbind(Z,X))
  ps.hat <- predict(ps.fit,type="response")
  ipw <- (1-Z)/(1-ps.hat) + Z/ps.hat
  gee.fit <- geeglm(Y~Z, weights=ipw, id=1:length(Y),
                    family=gaussian(link = "identity"), corstr="independence")
  tau <- summary(gee.fit)$coef["Z","Estimate"] + 
    c(0,-1,1)*qnorm(0.975)*summary(gee.fit)$coef["Z","Std.err"]
  names(tau) <- c("IPW.gee","ci.low","ci.upp")  
  y1 <- sum(ipw*Z*Y)/sum(ipw*Z)
  y0 <- sum(ipw*(1-Z)*Y)/sum(ipw*(1-Z))
  tau <- c("IPW.raw"=y1-y0,tau)
  return(tau)
}

DeltasClassSpecific <- function(C,X,Z,Y) {
  Cvals <- sort(unique(C))
  sapply(Cvals, function(cl) {
    if (any(c(sum(Z[C==cl])==c(0,sum(C==cl)),sum(Y[C==cl])==c(0,sum(C==cl))))){
      cl_ <- rep(NA,4)
      names(cl_) <- c("IPW.raw","IPW.gee","ci.low","ci.upp")  
      return(cl_)
    } else {
      IPWestimator(X=X[C==cl,,drop=FALSE],Z=Z[C==cl],Y=Y[C==cl])  
    }
  })
}

Expit <- function(x) {
  x_noInf <- pmin(x,-log(.Machine$double.eps)) # avoid Inf in denominator
  exp(x_noInf)/(1+exp(x_noInf))
}

# function to calculate sequence of class-specific effect estimates
# given an (ordered) dataset
DeltasSequence <- function(Data) {
  n <- nrow(Data)
  Xnames <- colnames(Data)[!colnames(Data) %in% c("Z","Y")]
  deltas.l1o <- NULL
  for (i in 1:(n-1)) {
    C.l1o <- rep(2,n)
    C.l1o[1:i] <- 1
    fit.ps <- lapply(1:2, function(cl) {
      # check if PS model can be fitted within each class
      stm <- tryCatch(
        fit <- glm(ps.model,family=binomial("logit"),data=Data[C.l1o==cl,]), 
        error=function(cond) return(NA))
      if (!is.na(stm)) {
        ps.hat <- predict(fit, type="response")
      } else {
        ps.hat <- rep(NA,sum(C.l1o==cl))
      }
      return(ps.hat)
    })
    deltas.l1o[[i]] <- DeltasClassSpecific(
      C=C.l1o,X=Data[,Xnames],Z=Data$Z,Y=Data$Y)
    if (i %% 1000 == 0) {
      cat(i,"/",n,"\n")  
    }
  }
  deltas.l1o <- do.call(rbind,deltas.l1o)
  return(deltas.l1o)
}