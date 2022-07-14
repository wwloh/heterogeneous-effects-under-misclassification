# read in preamble ============================================================
rm(list=ls())
source2 <- function(file, start, end, ...) {
  file.lines <- scan(file, what=character(), 
                     skip=start-1, nlines=end-start+1, sep='\n')
  file.lines.collapsed <- paste(file.lines, collapse='\n')
  source(textConnection(file.lines.collapsed), ...)
}
source2("sim-uncertainty-latentclass.R",2,99)
rm(source2)

for (j in 1:M) {
  ptm <- proc.time()[3]
  # check whether mixture model can be fitted
  mixOutOK <- FALSE
  while(!mixOutOK) {
    # observed sample
    obs_idx <- sort(unlist(sapply(1:2, function(cl) 
      sample(x=which(C_true==cl),
             size=as.integer(ifelse(cl==1,Lambda,1-Lambda)*n),
             replace=FALSE)))) # sample without replacement to fix Lambda
    obs_C_true <- C_true[obs_idx]
    ## observed treatment under true class
    Z_true <- sapply(1:2, function(cl) {
      Z_cl <- rep(NA,n)
      ps_true <- (as.matrix(cbind(1,X[obs_idx,])) %*% alpha[,cl])[,1]
      Z_cl[obs_C_true==cl] <- rbinom(sum(obs_C_true==cl),1,
                                     Expit(ps_true[obs_C_true==cl]))
      return(Z_cl)
    })
    Z <- rowSums(Z_true,na.rm=TRUE)
    ## reveal potential outcome under true class for observed treatment
    Y <- unlist(sapply(1:2, function(cl) {
      y_cl <- Y_true[[cl]][obs_idx,][obs_C_true==cl,]
      Z_cl <- Z[obs_C_true==cl]
      Y_cl <- y_cl[,1]
      Y_cl[Z_cl==1] <- y_cl[Z_cl==1,2]
      return(Y_cl)
    }))
    obs_data <- data.frame("i"=1:n,X[obs_idx,],Z,Y)
    row.names(obs_data) <- NULL
    
    obsOK <- c(
      all(table(obs_C_true)/n==table(C_true)/N),
      all(rowSums(is.na(Z_true))==1), # each obs belongs to only one class
      all(table(Z,Y)!=0), # no perfect separation
      all(table(Z,obs_C_true)!=0), # each class has both levels of treatment
      all(table(Y,obs_C_true)!=0) # each class has both levels of outcome
    )
    if (any(!obsOK)) next
    rm(Y,Z,Z_true,obs_idx)
    
    if (meth=="poLCA") {
      # fit latent class analysis
      mixOut <- poLCA(cbind(X1, X2, X3, X4, X5, X6)~1, 
                      data=obs_data, nclass=2,verbose=FALSE)
      fixed_C <- mixOut$predclass
    } else if (meth=="mclust") {
      # model-based clustering: parameterized finite Gaussian mixture models
      mixOut <- densityMclust(data=obs_data[,Xnames],G=2,verbose=FALSE,
                              modelNames="VVV")
      fixed_C <- mixOut$classification
    }
    # check that different criteria for mixture model have been met
    # check variation of X within each class
    var.X.withinclass <- sapply(sort(unique(mixOut$predclass)), function(cl) {
      apply(obs_data[fixed_C==cl,Xnames], 2, function(x) length(unique(x)))
    })
    mixOutOK <- c(all(var.X.withinclass>1),
                  # at least one individual in each class
                  length(unique(fixed_C)==2))
    mixOutOK <- all(mixOutOK)
    if (mixOutOK) {
      # class-specific estimates under estimated class memberships
      est <- DeltasClassSpecific(
        C=fixed_C,X=obs_data[,Xnames],Z=obs_data$Z,Y=obs_data$Y,status=TRUE)
      mixOutOK <- all(!is.na(est))  
    }
  }
    
  res.j <- list()
  # population values
  res.j[["pop"]] <- pop_effs
  # true class memberships
  res.j[["true"]] <- DeltasClassSpecific(
    C=obs_C_true,X=obs_data[,Xnames],Z=obs_data$Z,Y=obs_data$Y,status=TRUE)
  # misclassification
  res.j[["confusion_fixed"]] <- matrix(
    table("true"=obs_C_true,"fixed"=fixed_C)/n,nrow=2)
  res.j[["est"]] <- est
  # switch labels to align with true classes
  if(sum(diag(res.j[["confusion_fixed"]]))<0.5) {
    fixed_C <- 3-fixed_C # switch labels
    res.j[["est"]] <- res.j[["est"]][,2:1]
    res.j[["confusion_fixed"]] <- res.j[["confusion_fixed"]][,2:1]
  }
  rm(est)
  
  # bias-corrected estimates
  if (meth=="poLCA") {
    lambda.hat <- mixOut$posterior
  } else if (meth=="mclust") {
    lambda.hat <- mixOut$z
  }
  if (all(diag(table(apply(lambda.hat, 1, which.max),fixed_C))==0)) {
    # switch labels
    lambda.hat <- lambda.hat[,2:1]
  }
  for (mt in c("out","ipw")) {
    res.j[[paste0("bias_corrected_",mt)]] <- ATE_biascorrected(
      C=fixed_C,X=obs_data[,Xnames],Z=obs_data$Z,Y=obs_data$Y,Probs=lambda.hat,
      out_or_ipw=mt)  
    # use true classes
    res.j[[paste0("true_bias_corrected_",mt)]] <- ATE_biascorrected(
      C=obs_C_true,X=obs_data[,Xnames],Z=obs_data$Z,Y=obs_data$Y,
      Probs=cbind(obs_C_true==1,obs_C_true==2)*1.0,
      out_or_ipw=mt)  
  }
  
  # trajectory based on posterior class membership probabilities
  d.traj <- DeltasSequence(Data=obs_data,
                           Predclass=fixed_C,
                           Probs=lambda.hat,
                           Xnames=Xnames,
                           meth.names=c("AIPW.enet_Y.logit_Z"))
  
  # range of trajectory and whether class-specific population effect was recovered
  res.j[["trajectories"]] <- sapply(1:2, function(cl) {
    res.cl <- list(
      # point estimates
      "range.point"=range(d.traj[[cl]][,1],na.rm=TRUE),
      "captured_between.consecutive_points"=any(sapply(2:n, function(i) {
        i.interval <- sort(d.traj[[cl]][c(i-1,i),1])
        (i.interval[1] <= pop_effs[cl]) && (pop_effs[cl] <= i.interval[2])
        }),na.rm=TRUE),
      # CIs
      "captured_in.some_ci"=any(sapply(1:n, function(i) {
        i.ci <- d.traj[[cl]][i,2:3]
        (i.ci[1] <= pop_effs[cl]) && (pop_effs[cl] <= i.ci[2])
        }),na.rm=TRUE)
    )
    res.cl[["captured_in.range"]] <- 
      (res.cl$range[1] <= pop_effs[cl]) && (pop_effs[cl] <= res.cl$range[2])
    return(unlist(res.cl))
  })

  res.j[["simsetting"]] <- simsettings[seed,]
  res[[j]] <- res.j
  rm(res.j)
  subfolder <- "sim2-trajectories/"
  myfile <- paste0("sim-",seed,".Rdata")
  save(res,file=paste0(subfolder,myfile))
  cat(j,"|",round((proc.time()[3]-ptm)/60),
      "mins  ############################################################## \n")
  # 90 mins per sim
}
q()
