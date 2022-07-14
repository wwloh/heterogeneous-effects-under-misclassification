rm(list=ls())
library("Hmisc")
library("geepack")
library("xtable")
library("glmnet")
library("MatrixModels")
library("poLCA")

source("helper_functions.R")

# load processed data
argss <- 1001
if (!grepl("apple",sessionInfo()[[1]]$platform)) {
  argss <- commandArgs(trailingOnly=TRUE) # for CMD BATCH '--args 1'
}
(seed <- as.integer(argss[1]))

m <- 1000 # number of perturbations
if (seed <= m) {
  dataset_name <- "rhc" # RHC example  
} else {
  dataset_name <- "lindner" # example from twang package
}
dataset_name

filename <- paste0("illustration-",dataset_name,"-lca-subgroupATE.Rdata")
load(file=filename)
# reset to correct seed for perturbations (not from observed data)
rm(seed)
(seed <- as.integer(argss[1]))
rm(argss)
set.seed(9000*seed)

# formula for inclusive latent class model
polca.f

ptm <- proc.time()[3]
# calculate perturbed class membership probabilities ==========================
resample <- TRUE
sim_counter <- 0
while(resample) {
  # create simulated data with fitted model’s assumed data-generating process
  mc_sim <- poLCA::poLCA.simdata(N=n,probs=mixLCA$probs,P=mixLCA$P)
  # check levels of each manifest item in simulated data same as observed
  mc_sim.y <- mc_sim$dat[,grep("Y",colnames(mc_sim$dat))]
  colnames(mc_sim.y) <- colnames(mixLCA$y)
  mc_sim.OK <- identical(unlist(lapply(apply(mc_sim.y,2,unique),sort)),
                         unlist(lapply(apply(mixLCA$y,2,unique),sort)))
  if (mc_sim.OK) {
    mc_dat <- mc_sim$dat
    colnames(mc_dat) <- colnames(mixLCA$y)
    # fit latent class analysis
    ## avoid error with randomly-chosen starting values
    ## https://github.com/dlinzer/poLCA/issues/15
    nrep_iter <- 100L
    while(nrep_iter>0) {
      mc_fit_OK <- tryCatch(
        system.time(
          mc_fit <- poLCA(polca.f, data=mc_dat, nclass=length(mixLCA$P),
                          nrep=nrep_iter,verbose=FALSE)),
        error=function(cond) return(NA))
      if(all(!is.na(mc_fit_OK))) {
        # calculate posterior probabilities for observed data
        # using randomly drawn parameter values from sampling distribution
        mc_posterior <- poLCA::poLCA.posterior(lc=mc_fit,
                                               y=mapply(as.numeric,mixLCA$y))
        # optimal Bayes rule
        mc_C <- as.integer(apply(mc_posterior, 1, which.max))
        if (length(unique(mc_C))==length(mixLCA$P) & all(!is.na(mc_C))) {
          # class-specific estimates under resampled class memberships
          if (dataset_name=="rhc") {
            mc_est <- DeltasClassSpecific(
              C=mc_C,X=obs_data[,Xnames],Z=obs_data$Z,Y=obs_data$Y,
              meth.names=c("AIPW.enet_Y.logit_Z"))
          } else if (dataset_name=="lindner") {
            mc_est <- DeltasClassSpecific(
              C=mc_C,X=obs_data[,Xnames],Z=obs_data$Z,Y=obs_data$Y)
          }
          if (all(!is.na(mc_est))) {
            # order classes by increasing lower bound of CI for ATE
            (class_order.mc_est <- order(mc_est[1,],decreasing=FALSE))
            ## predicted class membership (probabilities)
            mc_lambda <- mc_posterior[,class_order.mc_est]
            mc_imputedclass <- as.integer(apply(mc_lambda, 1, which.max))
            resample <- FALSE
            nrep_iter <- 0L
          } else {
            nrep_iter <- nrep_iter-1
            cat(sim_counter, "|", nrep_iter, "============================ \n") 
          }
        }
      } else {
        nrep_iter <- nrep_iter-1
        cat(sim_counter, "|", nrep_iter, "================================ \n")
      }
    }
    if (resample==TRUE) {
      sim_counter <- sim_counter+1
    }
  }
}

mc_est
table(mc_imputedclass)

# trajectory based on posterior class membership probabilities ================
if (dataset_name=="rhc") {
  mc_d.traj <- DeltasSequence(Data=obs_data,
                              Predclass=mc_imputedclass,
                              Probs=mc_lambda,
                              Xnames=Xnames,
                              meth.names=c("AIPW.enet_Y.logit_Z"))
} else if (dataset_name=="lindner") {
  mc_d.traj <- DeltasSequence(Data=obs_data,
                              Predclass=mc_imputedclass,
                              Probs=mc_lambda,
                              Xnames=Xnames)
}

(filename <- gsub("subgroupATE",paste0("subgroupATE-",seed),filename))
save.image(file=paste0("illustrations-all-2-perturb-lca/",filename))
proc.time()[3]-ptm
q()
