rm(list=ls())
library("geepack")
library("poLCA");library("mclust")
library("mvtnorm"); library("glmnet")
library("grf")
library("CBPS")

source("helper_functions.R")

# simulation settings
simsettings <- expand.grid("meth"=c("poLCA","mclust"),
                           "Lambda"=c(0.4),
                           "direct"=c(0,1))
simsettings

# initialize for parallel cluster jobs
args <- 1
if (!grepl("apple",sessionInfo()[[1]]$platform)) {
  args <- commandArgs(trailingOnly=TRUE) # for CMD BATCH '--args 1'
  simsettings <- simsettings[rep(1:nrow(simsettings),each=200),]
  nrow(simsettings)
}
(seed <- as.integer(args[1]))
rm(args)

# data-generating parameters ==================================================
(meth <- simsettings[seed,"meth"])
(Lambda <- simsettings[seed,"Lambda"])
(direct <- simsettings[seed,"direct"])
M <- 5 # number of simulations
n <- 1000 # sample size
p <- 6 # number of covariates
m <- 100 # number of Monte Carlo or bootstrap draws
N <- n*100 # population size
if (meth=="poLCA") {
  pX1 <- c(0.75, 0.25) # Pr(X=1) for dichotomous covariates in each latent class
} else if (meth=="mclust") {
  mu <- list()
  mu[[1]] <- rep(c(0.75, 0.25),each=p/2)
  mu[[2]] <- -mu[[1]]
  sigma <- list()
  sigma[[1]] <- matrix(0, nrow = p, ncol = p); diag(sigma[[1]]) <- 1
  sigma[[2]] <- sigma[[1]]
}
alpha <- matrix(0,nrow=p+1,ncol=2)
beta <- matrix(0,nrow=p+2,ncol=2)
if (direct==1) {
  confound_idx <- list()
  confound_idx[[1]] <- 1:(p/2)
  confound_idx[[2]] <- (p/2)+(1:(p/2))
  for (cl in 1:2) {
    alpha[1+confound_idx[[cl]],cl] <- beta[2+confound_idx[[cl]],cl] <- 0.7
  }
}
row.names(alpha) <- paste0("a",0:p)
row.names(beta) <- paste0("b",0:(p+1))
beta["b1",] <- c(-2.4,1.8) # treatment coefficient

# generate data ===============================================================
set.seed(9000)
# true subgroup
C_true <- rep(2L,N); C_true[1:as.integer(Lambda*N)] <- 1L
if (meth=="poLCA") {
  # dichotomous indicator variables for measuring latent class membership
  X <- data.frame(do.call(cbind,lapply(1:p, function(u) {
    rbinom(N,1,rowSums(sapply(1:2, function(cl) {
      (C_true==cl)*pX1[cl]
    }))) + 1L # poLCA requires (non-zero) positive integers
  })))
  Xnames <- colnames(X) <- paste0("X",1:p)
} else if (meth=="mclust") {
  # continuous manifest variables for measuring latent subgroup membership
  X_true <- lapply(1:2, function(cl) {
    (C_true==cl)*rmvnorm(N,mu[[cl]],sigma[[cl]])
    })
  X <- data.frame(Reduce('+', X_true))
  Xnames <- colnames(X)
}
# generate potential outcomes
Y_true <- lapply(1:2, function(cl) {
  y_true <- matrix(NA,nrow=N,ncol=2)
  for (z in 0:1) {
    yz_true <- (as.matrix(cbind(1,z,X[C_true==cl,])) %*% beta[,cl])[,1]
    y_true[C_true==cl,z+1] <- rbinom(length(yz_true),1,Expit(yz_true))
  }
  return(y_true)
})
## check potential outcomes generated only under the true class
all.equal(which(rowSums(!is.na(Y_true[[1]]))==2),which(C_true==1))
all.equal(which(rowSums(!is.na(Y_true[[2]]))==2),which(C_true==2))
# true class-specific effects under true classes
(pop_effs <- sapply(Y_true, function(y_true) 
  mean(apply(y_true,1,diff),na.rm=TRUE)))
# population (weighted) average effect
mean(apply(do.call(rbind,Y_true),1,diff),na.rm=TRUE)

# observed samples ============================================================
set.seed(seed)
res <- NULL
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
  
  # parametric bootstrap
  mc.te <- list()
  mc.jit <- list()
  mc.confusion <- list()
  mc_C_resample <- NULL
  resample <- TRUE
  pb <- txtProgressBar(min=0, max=m, initial=0, style=3);cat("\n")
  while(resample) {
    # create simulated data with fitted model’s assumed data-generating process 
    if (meth=="poLCA") {
      mc_sim <- poLCA::poLCA.simdata(N=n,probs=mixOut$probs,P=mixOut$P)
      mc_dat <- mc_sim$dat
      colnames(mc_dat) <- colnames(mixOut$y)
      # fit latent class analysis
      mc_fit <- poLCA(cbind(X1, X2, X3, X4, X5, X6)~1, 
                      data=mc_dat, nclass=2,verbose=FALSE)
      # calculate posterior probabilities for observed data
      # using randomly drawn parameter values from sampling distribution
      mc_dens <- poLCA::poLCA.posterior(lc=mc_fit,y=mixOut$y)
    } else if (meth=="mclust") {
      mc_sim <- mclust::sim(modelName=mixOut$modelName,
                            parameters=mixOut$parameters,n=n)
      mc_dat <- data.frame(mc_sim[,-1])
      colnames(mc_dat) <- colnames(mixOut$data)
      # fit model-based clustering
      mc_fit <- densityMclust(data=mc_dat,G=2,verbose=FALSE,modelNames="VVV")
      # calculate component densities for observed data
      # "z" = conditional probabilities of belonging to each mixture component
      # using randomly drawn parameter values from sampling distribution
      mc_dens <- predict.densityMclust(object=mc_fit,newdata=mixOut$data,what="z")
    }
    
    # optimal Bayes rule (modal assignment)
    mc_C <- as.integer(apply(mc_dens, 1, which.max))
    # check variation of X within each class
    var.X.withinclass <- sapply(sort(unique(mc_C)), function(cl) {
      apply(obs_data[mc_C==cl,Xnames], 2, function(x) 
        length(unique(x)))
    })
    if (all(var.X.withinclass>1) && length(unique(mc_C))==2) {
      # class-specific estimates under resampled class memberships
      mc_est <-DeltasClassSpecific(
        C=mc_C,X=obs_data[,Xnames],Z=obs_data$Z,Y=obs_data$Y)
      # misclassification
      mc_confusion <- matrix(table("true"=obs_C_true,"fixed"=mc_C)/n,nrow=2)
      # switch labels to align with true classes
      if(sum(diag(mc_confusion))<0.5) {
        mc_C <- 3-mc_C # switch labels
        mc_est <- mc_est[,2:1]
        mc_confusion <- mc_confusion[,2:1]
      }
      mc.te <- c(mc.te, list(mc_est))
      mc_C_resample <- cbind(mc_C_resample,mc_C)
      mc.confusion <- c(mc.confusion, list(mc_confusion))
      rm(mc_est,mc_C,mc_confusion)
      
      resample <- (length(mc.te) < m)
      setTxtProgressBar(pb,length(mc.te))
    }
    rm(var.X.withinclass)
  }
  
  res_mc <- lapply(row.names(res.j$est),function(est_) {
    do.call(rbind,lapply(mc.te, "[", est_, ))
  })
  names(res_mc) <- row.names(res.j$est)
  meths_ate <- unlist(strsplit(grep("ATE",names(res_mc),value=TRUE),split=".ATE"))
  for (ci.union_level in c(0,.025)) {
    meths_ate.ci <- NULL
    for (meth_ate in meths_ate) {
      meths_ate.ci[[meth_ate]] <- sapply(1:2, function(cl) {
        c(res.j$est[paste0(meth_ate,".ATE"),cl],
          quantile(res_mc[[paste0(meth_ate,".ci.low")]][,cl],
                   probs=ci.union_level,na.rm=TRUE),
          quantile(res_mc[[paste0(meth_ate,".ci.upp")]][,cl],
                   probs=1-ci.union_level,na.rm=TRUE))
      })
      row.names(meths_ate.ci[[meth_ate]]) <- 
        paste0(meth_ate,c(".ATE",".ci.low",".ci.upp"))
    }
    res.j[[paste0("boot.ci",ifelse(ci.union_level==0,"",".95"))]] <- 
      do.call(rbind,meths_ate.ci)
  }
  res.j[["confusion_boot"]] <- Reduce('+', mc.confusion)/length(mc.confusion)
  rm(obs_C_true,mc_C_resample,mc.confusion)
  
  res.j[["simsetting"]] <- simsettings[seed,]
  res[[j]] <- res.j
  rm(res.j)
  subfolder <- "sim-uncertainty-latentclass/"
  myfile <- paste0("sim-",seed,".Rdata")
  save(res,file=paste0(subfolder,myfile))
  cat(j,"|",round((proc.time()[3]-ptm)/60),
      "mins  ############################################################## \n")
  # 90 mins per sim
}
q()
