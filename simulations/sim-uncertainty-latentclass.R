rm(list=ls())
library("geepack")
library("poLCA")
library("mclust")
library("mvtnorm")

source("helper_functions.R")

# simulation settings
simsettings <- expand.grid("meth"=c("poLCA","mclust"),
                           "Lambda"=c(0.8,0.6))
simsettings

# initialize for parallel cluster jobs
args <- 1
if (!grepl("apple",sessionInfo()[[1]]$platform)) {
  args <- commandArgs(trailingOnly=TRUE) # for CMD BATCH '--args 1'
  simsettings <- simsettings[rep(1:nrow(simsettings),each=10),]
  nrow(simsettings)
}
(seed <- as.integer(args[1]))
rm(args)

# data-generating parameters ==================================================
(meth <- simsettings[seed,"meth"])
(Lambda <- simsettings[seed,"Lambda"])
M <- 100 # number of simulations
n <- 1000 # sample size
m <- 100 # number of Monte Carlo or bootstrap draws
N <- 10^5 # population size
if (meth=="poLCA") {
  pX1 <- c(0.8, 0.2) # Pr(X=1) for dichotomous covariates
} else if (meth=="mclust") {
  mu <- list()
  mu[[1]] <- c(1,1,1)*0.75
  mu[[2]] <- c(1,1,1)*(-0.3)
  sigma <- list()
  sigma[[1]] <- matrix(1/2, nrow = 3, ncol = 3); diag(sigma[[1]]) <- 1
  sigma[[2]] <- matrix(1/4, nrow = 3, ncol = 3); diag(sigma[[2]]) <- 1
}
alpha <- list()
alpha[["a0"]] <- c(0,0)
alpha[["a1"]] <- c(1,0)
alpha[["a2"]] <- -alpha[["a1"]]/2
alpha[["a3"]] <- -alpha[["a1"]]/2
beta <- list()
beta[["b0"]] <- c(0,0)
beta[["b1"]] <- c(-1.4,0.7)
beta[["b2"]] <- c(1,2)
beta[["b3"]] <- -beta[["b2"]]/2
beta[["b4"]] <- -beta[["b2"]]/2

# generate data ===============================================================
set.seed(9000)
C_true <- rbinom(N,1,Lambda) + 1L
if (meth=="poLCA") {
  # dichotomous indicator variables for measuring latent class membership
  X <- data.frame(do.call(cbind,lapply(1:3, function(u) {
    rbinom(N,1,rowSums(sapply(1:2, function(cl) {
      (C_true==cl)*pX1[cl]
    }))) + 1L # poLCA requires (non-zero) positive integers
  })))
  Xnames <- colnames(X) <- paste0("X",1:3)
} else if (meth=="mclust") {
  X_true <- lapply(1:2, function(cl) (C_true==cl)*rmvnorm(N,mu[[cl]],sigma[[cl]]))
  X <- data.frame(X_true[[1]] + X_true[[2]])
  Xnames <- colnames(X)
}
Z_true <- sapply(1:2, function(cl) {
  ps_true <- alpha$a0[cl] + alpha$a1[cl]*X$X1 + alpha$a2[cl]*X$X2 + 
    alpha$a3[cl]*X$X3
  (C_true==cl)*rbinom(N,1,Expit(ps_true))
})
Z <- rowSums(Z_true)
Y_true <- sapply(1:2, function(cl) {
  y_true <- beta$b0[cl] + beta$b1[cl]*Z + beta$b2[cl]*X$X1 + beta$b3[cl]*X$X2 +
    beta$b4[cl]*X$X3
  (C_true==cl)*rgeom(N,Expit(y_true))
})
# higher underlying probability => fewer flips for first success
Y <- as.integer(rowSums(Y_true)==0) 

pop_data <- data.frame("i"=1:N,X,Z,Y)
rm(X,Z,Y,Z_true,Y_true)

## class-specific effects under true classes
pop_effs <- DeltasClassSpecific(C=C_true,X=pop_data[,Xnames],
                                Z=pop_data$Z,Y=pop_data$Y)
pop_effs

# population (weighted) average effect
(ate <- IPWestimator(X=pop_data[,Xnames],Z=pop_data$Z,Y=pop_data$Y))

# observed samples ============================================================
set.seed(seed)
res <- NULL
ptm=proc.time()[3]
for (j in 1:M) {
  # check whether mixture model can be fitted
  mixOutOK <- FALSE
  while(!mixOutOK) {
    # observed sample
    obs_idx <- sort(sample(N,n,replace=FALSE))
    obs_data <- pop_data[obs_idx,]
    obs_C_true <- C_true[obs_idx]
    
    if (meth=="poLCA") {
      # fit latent class analysis
      mixOut <- poLCA(cbind(X1, X2, X3)~1, data=obs_data, nclass=2,verbose=FALSE)
      # class-specific estimates under estimated class memberships
      est <- DeltasClassSpecific(
        C=mixOut$predclass,X=obs_data[,Xnames],Z=obs_data$Z,Y=obs_data$Y)
      mixOutOK <- all(!is.na(est))
    } else if (meth=="mclust") {
      # model-based clustering: parameterized finite Gaussian mixture models
      mixOut <- densityMclust(data=obs_data[,Xnames],G=2,verbose=FALSE)
      # check that different criteria for mixture model have been met
      mixOutOK <- c(
        # at least one individual in each class
        length(unique(mixOut$classification)==2)
      )
      mixOutOK <- all(mixOutOK)
      if (mixOutOK) {
        # class-specific estimates under estimated class memberships
        est <- DeltasClassSpecific(
          C=mixOut$classification,X=obs_data[,Xnames],Z=obs_data$Z,Y=obs_data$Y)
        mixOutOK <- mixOutOK && all(!is.na(est))
      }
    }
  }
  row.names(obs_data) <- NULL
  obs_data[,"i"] <- 1:n  
    
  res.j <- list()
  # true class memberships
  res.j[["true"]] <- DeltasClassSpecific(
    C=obs_C_true,X=obs_data[,Xnames],Z=obs_data$Z,Y=obs_data$Y)
  # switch labels so that class 1 effect < class 2 effect
  res.j[["est"]] <- est[,order(est["IPW.gee",])]
  # misclassification
  if (meth=="poLCA") {
    fixed_C <- mixOut$predclass
  } else if (meth=="mclust") {
    fixed_C <- mixOut$classification
  }
  if (any(order(est["IPW.gee",])!=(1:2))) {
    fixed_C <- 3-fixed_C # switch labels
  }
  res.j[["confusion_fixed"]] <- matrix(
    table("true"=obs_C_true,"fixed"=fixed_C)/n,nrow=2)
  rm(est,fixed_C)
  
  # parametric bootstrap
  mc.te <- list()
  mc.jit <- list()
  mc_C_resample <- NULL
  resample <- TRUE
  while(resample) {
    # create simulated data with fitted model’s assumed data-generating process 
    if (meth=="poLCA") {
      mc_sim <- poLCA::poLCA.simdata(N=n,probs=mixOut$probs,P=mixOut$P)
      mc_dat <- mc_sim$dat
      colnames(mc_dat) <- colnames(mixOut$y)
      # fit latent class analysis
      mc_fit <- poLCA(cbind(X1, X2, X3)~1, data=mc_dat, nclass=2,verbose=FALSE)
      # calculate posterior probabilities for observed data
      # using randomly drawn parameter values from sampling distribution
      mc_dens <- poLCA::poLCA.posterior(lc=mc_fit,y=mixOut$y)
    } else if (meth=="mclust") {
      mc_sim <- mclust::sim(modelName=mixOut$modelName,
                            parameters=mixOut$parameters,n=n)
      mc_dat <- data.frame(mc_sim[,-1])
      colnames(mc_dat) <- colnames(mixOut$data)
      # fit model-based clustering
      mc_fit <- densityMclust(data=mc_dat,G=2,verbose=FALSE)
      # calculate component densities for observed data
      # "z" = conditional probabilities of belonging to each mixture component
      # using randomly drawn parameter values from sampling distribution
      mc_dens <- predict.densityMclust(object=mc_fit,newdata=mixOut$data,what="z")
    }
    
    # optimal Bayes rule
    mc_C <- as.integer(apply(mc_dens, 1, which.max))
    if (length(unique(mc_C))==2) {
      # class-specific estimates under resampled class memberships
      mc_est <-DeltasClassSpecific(
        C=mc_C,X=obs_data[,Xnames],Z=obs_data$Z,Y=obs_data$Y)
      # switch labels so that class 1 effect < class 2 effect
      mc.te <- c(mc.te, list(mc_est[,order(mc_est["IPW.gee",])]))
      if (any(order(mc_est["IPW.gee",])!=(1:2))) {
        mc_C <- 3-mc_C # switch labels
      }
      mc_C_resample <- cbind(mc_C_resample,mc_C)
      rm(mc_est,mc_C)
      
      resample <- (length(mc.te) < m)
    }
  }
  res.j[["boot.ci"]] <- sapply(1:2, function(cl) {
    c("IPW.gee"=mean(unlist(lapply(mc.te, function(te) te["IPW.gee",cl])),na.rm=TRUE),
      "ci.low"=min(unlist(lapply(mc.te, function(te) te["ci.low",cl])),na.rm=TRUE),
      "ci.upp"=max(unlist(lapply(mc.te, function(te) te["ci.upp",cl])),na.rm=TRUE))
  })
  
  res.j[["confusion_boot"]] <- matrix(rowMeans(
    apply(mc_C_resample,2,function(mc_C) {
      table("true"=obs_C_true,"fixed"=mc_C)/n
    })),nrow=2)
  rm(obs_C_true,mc_C_resample)
  
  res.j[["simsetting"]] <- simsettings[seed,]
  res[[j]] <- res.j
  save(res,file=paste0("sim-",meth,"-",seed,".Rdata"))
  cat(j,"|",round((proc.time()[3]-ptm)/60),
      "mins  ############################################################## \n")
}
q()
