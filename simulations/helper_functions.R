# helper functions to calculate class-specific effect estimators (in each class)
AIPWestimator <- function(X,Z,Y,type,ps.form=NULL,out.form=NULL) {
  tau <- rep(NA,3)
  names(tau) <- c("ATE","ci.low","ci.upp")
  # propensity score model ====================================================
  if (is.null(ps.form)) {
    # main effects for covariates
    ps.form <- as.formula(paste0("Z~",paste(colnames(X),collapse="+")))
  }
  if (grepl("logit_Z",type)) {
    # logistic regression model
    stm <- tryCatch(
      system.time(
        ps.fit <- glm(ps.form,family=binomial("logit"), data=cbind(Z,X))),
      error=function(cond) return(NA))
    if (any(is.na(stm))) {
      return(tau)
    } else if(!ps.fit$converged || ps.fit$boundary) {
      return(tau)
    }
  } else if (grepl("CBPS_Z",type)) {
    # CBPS
    stm <- tryCatch(
      system.time(
        ps.fit <- CBPS::CBPS(formula=ps.form, data=cbind(Z,X), ATT=0,
                             iterations=5000)),
      error=function(cond) return(NA))
    if (any(is.na(stm))) {
      return(tau)
    } else if(ps.fit$converged!=0) {
      return(tau)
    }
  }
  ps.hat <- predict(ps.fit,type="response")
  # error if there is an extreme propensity score (up to machine error)
  if(any(pmin(ps.hat,1-ps.hat) < (.Machine$double.eps*1e1))) {
    return(tau)
  }
  
  # effect estimator ==========================================================
  if (grepl("IPWonly",type)) {
    # IPW estimator without any outcome model
    ipw <- (1-Z)/(1-ps.hat) + Z/ps.hat # inverse probability weights
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
      tau[2:3] <- tau[1] + c(-1,1)*qnorm(0.975)*
        summary(gee.fit)$coef["Z","Std.err"]
      return(tau)
    }
  } else if (grepl("_Y",type)) {
    # AIPW estimator with parametric outcome model ============================
    binaryY <- length(sort(unique(Y)))==2
    if (grepl("enet",type)) {
      if (ncol(X)<log2(1/.Machine$double.eps)) {
        # saturated outcome model with all possible interactions
        satu_y <- as.formula(paste0("Y~",paste(
          c("Z",colnames(X)),collapse="*")))
      } else {
        # outcome model with all possible 2-way interactions with treatment
        satu_y <- as.formula(paste0("Y~",paste(
          paste0("Z*",colnames(X),collapse="+"))))
        # satu_y <- as.formula(paste0("Y~",paste(
        #   paste0("Z*",combn(colnames(X),2,FUN=paste, collapse="*")),
        #   collapse="+")))
      }
      modmat <- MatrixModels::model.Matrix(
        object=satu_y,data=cbind(Y,Z,X),sparse=TRUE,drop.unused.levels=TRUE)
      modmat <- modmat[,-1] # remove intercept
      fit_OK <- tryCatch(
        fitY.cv <- cv.glmnet(x=as.matrix(modmat),
                             y=Y,
                             family=ifelse(binaryY,"binomial","gaussian"),
                             type.measure=ifelse(binaryY,"class","deviance"),
                             standardize=FALSE,trace.it=FALSE,nlambda=1000, 
                             alpha=0.5, #elastic net mixing parameter
                             # better to leave main effects unpenalized
                             penalty.factor=grepl(
                               pattern=":",x=colnames(modmat),fixed=TRUE)*1.0),
        error=function(cond) return(NA))
      if(any(is.na(fit_OK))) {
        return(tau)
      } else {
        fitY.coef <- coef(fitY.cv, s="lambda.min")[,1]
        fitY.coef <- fitY.coef[names(fitY.coef) %in% colnames(modmat)]
        # non-zero estimates
        sel.Y <- names(fitY.coef)[abs(fitY.coef)>.Machine$double.eps]
        if (length(sel.Y) >= length(Y)) {
          # more non-zero coefficient estimates than observations
          fitY.coef.ordered <- fitY.coef[order(abs(fitY.coef),decreasing=TRUE)]
          sel.Y <- names(fitY.coef.ordered)[1:(length(Y)-1)]
        } else if (length(sel.Y)==0) {
          sel.Y <- "1"
        }
        # relaxed LASSO: refit with selected predictors
        out.form <- paste0("Y~",paste(sel.Y,collapse="+"))
        rm(fitY.cv,fitY.coef,sel.Y)
      }
    } else if (grepl("lm",type)) {
      # main effects only
      if (is.null(out.form)) {
        out.form <- paste0("Y~",paste(c("Z",colnames(X)),collapse="+"))
      }
    }
    out.form <- as.formula(out.form)
    if (binaryY) {
      # logistic regression model
      stm.out <- tryCatch(
        system.time(
          out.fit <- glm(out.form, family=binomial("logit"), data=cbind(Y,Z,X))),
        error=function(cond) return(NA))
    } else {
      # linear regression model
      stm.out <- tryCatch(
        system.time(
          out.fit <- lm(out.form, data=cbind(Y,Z,X))),
        error=function(cond) return(NA))
    }
    if (any(is.na(stm.out))) {
      return(tau)
    } else if (any(class(out.fit)=="glm")) {
      if(!out.fit$converged | out.fit$boundary) {
        return(tau)
      }
    }
    # for predicting counterfactuals
    mydata.A1 <- mydata.A0 <- as.data.frame(cbind(Z,X))
    mydata.A1[,"Z"] <- 1
    mydata.A0[,"Z"] <- 0
    Y.Z1.hat <- predict(out.fit,newdata=mydata.A1,type="response")
    Y.Z0.hat <- predict(out.fit,newdata=mydata.A0,type="response")
    
    infn.est <- (Z/ps.hat)*(Y-Y.Z1.hat) - ((1-Z)/(1-ps.hat))*(Y-Y.Z0.hat) +
      (Y.Z1.hat-Y.Z0.hat)
    
    # one-step plug-in estimator
    tau[1] <- mean(infn.est)
    # individual influence functions
    infn.est <- infn.est - tau[1]
    tau[2:3] <- tau[1] + 
      c(-1,1)*qnorm(0.975)*sqrt(var(infn.est)/length(infn.est))
    return(tau)
  }
}

# Calculate individual treatment effects using ML-based method
MLestimator <- function(X,Z,Y,ml_type) {
  if (ml_type=="CF") {
    # Causal forest ===========================================================
    # https://github.com/grf-labs/grf
    tau.forest <- grf::causal_forest(X=X, Y=Y, W=Z, tune.parameters = "all")
    # Estimate the conditional average treatment effect on the full sample (CATE)
    tau <- grf::average_treatment_effect(tau.forest, target.sample = "all")
    tau <- tau["estimate"] + c(0,-1,1)*qnorm(0.975)*tau["std.err"]
  } else if (ml_type=="hdm") {
    # Following example in \S4.4 of 
    # https://cran.r-project.org/web/packages/hdm/vignettes/hdm.pdf
    ## all possible higher-order and interaction terms for covariates
    satu_y <- as.formula(paste0("Y~",paste(colnames(X),collapse="*")))
    modmat <- model.matrix(satu_y,data=cbind(Y,X))[,-1]
    doublesel.effect = hdm::rlassoEffect(
      x = modmat, y = Y, d = Z, method = "double selection")
    tau <- doublesel.effect$alpha + c(0,-1,1)*qnorm(0.975)*doublesel.effect$se
  }
  tau <- as.numeric(tau)
  names(tau) <- c("ATE","ci.low","ci.upp")
  return(tau)
}

One_estimator <- function(X,Z,Y,
                          meth.names=c("AIPW.enet_Y.CBPS_Z"),
                          status=FALSE) {
  res <- NULL
  for (meth in meth.names) {
    ptm <- proc.time()[3]
    if (grepl("[.]",meth)) {
      res[[meth]] <- AIPWestimator(X,Z,Y,type=meth)
    } else {
      stm <- tryCatch(
        system.time(
          res[[meth]] <- MLestimator(X,Z,Y,ml_type=meth)),
        error=function(cond) return(NA))
      if(any(is.na(stm))) {
        res[[meth]] <- rep(NA,3)
      }
    }
    if (status==TRUE) {
      cat(meth, "took", proc.time()[3]-ptm, "secs \n")
    }
  }
  return(unlist(res))
}

DeltasClassSpecific <- function(C,X,Z,Y,...) {
  Cvals <- sort(unique(C))
  sapply(Cvals, function(cl) {
    One_estimator(X=X[C==cl,,drop=FALSE],Z=Z[C==cl],Y=Y[C==cl],...)
  })
}

Expit <- function(x) {
  x_noInf <- pmin(x,-log(.Machine$double.eps)) # avoid Inf in denominator
  exp(x_noInf)/(1+exp(x_noInf))
}

# function to calculate sequence of class-specific effect estimates
# given a dataset and vector of (posterior) probabilities
DeltasSequence <- function(Data,Predclass,Probs,Xnames,type,...) {
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
    d.traj[[cl]] <- list()
    pb <- txtProgressBar(min=1, max=nrow(Data.ordered), initial=0, style=3)
    cat("\n")
    for (idx in 1:nrow(Data.ordered)) {
      d.traj[[cl]][[idx]] <- One_estimator(
        X=Data.ordered[1:idx,Xnames,drop=FALSE],
        Z=Data.ordered[1:idx,"Z"],
        Y=Data.ordered[1:idx,"Y"],...)
      setTxtProgressBar(pb,idx);cat("\n")
    }
    close(pb)
    d.traj[[cl]] <- do.call(rbind,d.traj[[cl]])
    d.traj[[cl]] <- cbind(d.traj[[cl]],
                          "pr.belong"=Probs[c(lambda_in,lambda_out), cl])
    cat("Class",cl,"took",round((proc.time()[3]-ptm)/60),"mins \n")
  }
  return( d.traj )
}

# bias-corrected class-specific ATE ===========================================
ATE_biascorrected <- function(C,X,Z,Y,Probs,out_or_ipw="out") {
  Cvals <- sort(unique(C))
  # matrix of classification probabilities
  class.probs <- matrix(NA,nrow=length(Cvals),ncol=length(Cvals))
  for (j in Cvals) {
    for (k in Cvals) {
      class.probs[j,k] <- mean((Probs[,k]*(C==j))/(mean(C==j)))
    }
  }
  # calculate inverse
  class.probs.inv <- chol2inv(chol(class.probs))
  
  if (out_or_ipw=="out") {
    # conditional expectations of observed outcomes: main effects only
    out.form <- as.formula(paste0("Y~",paste(c("Z",colnames(X)),collapse="+")))
    binaryY <- length(sort(unique(Y)))==2
    tau_C <- NULL
    tau_C[["Z0"]] <- tau_C[["Z1"]] <- matrix(NA,nrow=length(Y),ncol=length(Cvals))
    for (cl in Cvals) {
      if (binaryY) {
        # logistic regression model
        stm.out <- tryCatch(
          system.time(
            out.fit <- glm(out.form, family=binomial("logit"), data=cbind(Y,Z,X),
                           subset=C==cl)),
          error=function(cond) return(NA))
      } else {
        # linear regression model
        stm.out <- tryCatch(
          system.time(
            out.fit <- lm(out.form, data=cbind(Y,Z,X),subset=C==cl)),
          error=function(cond) return(NA))
      }
      if (all(!is.na(stm.out))) {
        if (any(class(out.fit)=="glm")) {
          if (out.fit$converged & !out.fit$boundary) {
            # predicted class-specific outcomes for all individuals in the sample
            tau_C[["Z0"]][,cl] <- predict(out.fit,newdata=cbind(Z=0,X),
                                          type="response")
            tau_C[["Z1"]][,cl] <- predict(out.fit,newdata=cbind(Z=1,X),
                                          type="response")
          }
        }
      }
    }
    # sample average of potential outcomes under true class
    tau_Cstar <- lapply(tau_C, function(yhat) 
      rowMeans(tcrossprod(class.probs.inv, yhat)))
    return( Reduce('-',tau_Cstar) )
  } else {
    # IPW estimator without any outcome model
    # PS model: logistic regression model with main effects only
    ps.form <- as.formula(paste0("Z~",paste(colnames(X),collapse="+")))
    tau_Cstar <- rep(NA,length(Cvals))
    for (cl in Cvals) {
      ps.fit <- glm(ps.form, family=binomial("logit"), data=cbind(Z,X),
                    subset=C==cl)
      ps.hat <- predict(ps.fit,newdata=cbind(Z,X),type="response")
      ipw <- (1-Z)/(1-ps.hat) + Z/ps.hat # inverse probability weights
      ipw.misclass <- Probs[,cl]*tcrossprod(
        class.probs.inv[cl,],sapply(Cvals, function(k) (C==k)/(mean(C==k))))[1,]
      ipw <- ipw*ipw.misclass
      y1 <- sum(ipw*Z*Y)/sum(ipw*Z)
      y0 <- sum(ipw*(1-Z)*Y)/sum(ipw*(1-Z))
      tau_Cstar[cl] <- y1-y0
      rm(ipw.misclass,y1,y0)
    }
    return( tau_Cstar )
  }
}

