# helper functions to calculate class-specific effect estimators (in each class)
IPWestimator <- function(X,Z,Y) {
  tau <- rep(NA,3)
  names(tau) <- c("IPW","ci.low","ci.upp")
  stm <- tryCatch(
    system.time(
      ps.fit <- glm(as.formula(paste0("Z~",paste(colnames(X),collapse="+"))),
                    family=binomial("logit"), data=cbind(Z,X))),
    error=function(cond) return(NA))
  if (any(is.na(stm))) {
    return(tau)
  } else {
    if(!ps.fit$converged || ps.fit$boundary) {
      return(tau)
    }
    ps.hat <- predict(ps.fit,type="response")
    ipw <- (1-Z)/(1-ps.hat) + Z/ps.hat
    y1 <- sum(ipw*Z*Y)/sum(ipw*Z)
    y0 <- sum(ipw*(1-Z)*Y)/sum(ipw*(1-Z))
    tau[1] <- y1-y0
    stm.gee <- tryCatch(
      system.time(
        gee.fit <- geeglm(formula=Y~Z, weights=ipw, id=1:length(Y),
                          family=gaussian(link = "identity"), 
                          corstr="independence")),
      error=function(cond) return(NA))
    if (any(is.na(stm.gee))) {
      return(tau)
    } else {
      tau[2:3] <- summary(gee.fit)$coef["Z","Estimate"] + 
        c(-1,1)*qnorm(0.975)*summary(gee.fit)$coef["Z","Std.err"]
      return(tau)
    }
  }
}

DeltasClassSpecific <- function(C,X,Z,Y) {
  Cvals <- sort(unique(C))
  sapply(Cvals, function(cl) {
    IPWestimator(X=X[C==cl,,drop=FALSE],Z=Z[C==cl],Y=Y[C==cl])
  })
}

Expit <- function(x) {
  x_noInf <- pmin(x,-log(.Machine$double.eps)) # avoid Inf in denominator
  exp(x_noInf)/(1+exp(x_noInf))
}

# function to calculate sequence of class-specific effect estimates
# given a dataset and vector of (posterior) probabilities
DeltasSequence <- function(Data,Predclass,Probs,Xnames) {
  d.traj <- list()
  for (cl in 1:ncol(Probs)) {
    ptm <- proc.time()[3]
    # subset in imputed class
    idx_in <- which(Predclass==cl)
    idx_out <- which(Predclass!=cl)
    # order data according to decreasing probability
    lambda_in <- idx_in[order(Probs[idx_in,cl],na.last=TRUE,decreasing=TRUE)]
    lambda_out <- idx_out[order(Probs[idx_out,cl],na.last=TRUE,decreasing=TRUE)]
    # first add individuals in imputed class, then those from another class
    Data.ordered <- Data[c(lambda_in,lambda_out),]
    # nested sequence of subsets
    d.traj[[cl]] <- t(sapply(1:nrow(Data.ordered), function(idx) {
      IPWestimator(X=Data.ordered[1:idx,Xnames,drop=FALSE],
                   Z=Data.ordered[1:idx,"Z"],
                   Y=Data.ordered[1:idx,"Y"])
    }))
    d.traj[[cl]] <- cbind(d.traj[[cl]],
                          "pr.belong"=Probs[c(lambda_in,lambda_out), cl])
    cat(cl,round((proc.time()[3]-ptm)/60),"mins \n")
  }
  return( d.traj )
}
