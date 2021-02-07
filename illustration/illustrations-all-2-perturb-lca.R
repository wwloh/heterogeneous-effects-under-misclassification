rm(list=ls())
library("Hmisc")
library("geepack")
library("poLCA")
library("Matching")
library("mvtnorm")

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
  dataset_name <- "GerberGreenImai" # Matching example   
}
dataset_name

filename <- paste0("illustration-",dataset_name,"-lca.Rdata")
load(file=filename)
# reset to correct seed for perturbations (not from observed data)
rm(seed,args)
(seed <- as.integer(argss[1]))
rm(argss)
set.seed(9000*seed)

# formula for inclusive latent class model
polca.f

if (dataset_name=="rhc") {
  mixLCA <- mixLCA_k[[5]]
  est <- mixLCA_k_res[[5]][["est"]]
} else if (dataset_name=="GerberGreenImai") {
  mixLCA <- mixLCA_k[[4]]
  est <- mixLCA_k_res[[4]][["est"]]
}

# calculate perturbed class membership probabilities ==========================
resample <- TRUE
sim_counter <- 0
ptm <- proc.time()[3]
while(resample) {
  # create simulated data with fitted model’s assumed data-generating process
  mc_sim <- poLCA::poLCA.simdata(N=n,
                                 probs=mixLCA$probs,
                                 x=as.matrix(mixLCA$x[,-1]),
                                 b=mixLCA$coeff,
                                 P=mixLCA$P)
  # check levels of each manifest item in simulated data same as observed
  mc_sim.y <- mc_sim$dat[,grep("Y",colnames(mc_sim$dat))]
  colnames(mc_sim.y) <- colnames(mixLCA$y)
  if (any(apply(mc_sim.y,2,max) != apply(mixLCA$y,2,max))) {
    resample <- TRUE
    sim_counter <- sim_counter+1
    cat(sim_counter,"\n")
  } else {
    mc_dat <- mc_sim$dat
    colnames(mc_dat) <- c(colnames(mixLCA$y),colnames(mixLCA$x)[-1])
    # fit latent class analysis
    ## https://github.com/dlinzer/poLCA/issues/15
    mc_fit_OK <- NA
    nrep_iter <- 100
    while(any(is.na(mc_fit_OK))) {
      mc_fit_OK <- tryCatch(
        mc_fit <- poLCA(polca.f, data=mc_dat, nclass=length(mixLCA$P),
                        nrep=nrep_iter,verbose=TRUE),
        error=function(cond) return(NA))
      if(any(is.na(mc_fit_OK))) nrep_iter <- max(1,nrep_iter-1)
      cat(nrep_iter, "\n ====================================================")
    }
    # calculate posterior probabilities for observed data
    # using randomly drawn parameter values from sampling distribution
    mc_posterior <- poLCA::poLCA.posterior(lc=mc_fit,
                                           y=mapply(as.numeric,mixLCA$y),
                                           x=as.matrix(mixLCA$x[,-1]))
    # optimal Bayes rule
    mc_C <- as.integer(apply(mc_posterior, 1, which.max))
    if (length(unique(mc_C))==length(mixLCA$P) & all(!is.na(mc_C))) {
      # class-specific estimates under resampled class memberships
      mc_est <- DeltasClassSpecific(
        C=mc_C,X=obs_data[,Xnames],Z=obs_data$Z,Y=obs_data$Y)
      # order classes by increasing lower bound of CI for ATE
      (class_order.mc_est <- order(mc_est["ci.low",],decreasing=FALSE))
      ## predicted class membership (probabilities)
      mc_lambda <- mc_posterior[,class_order.mc_est]
      mc_imputedclass <- as.integer(apply(mc_lambda, 1, which.max))
      resample <- FALSE
    }
  }
}

# trajectory based on posterior class membership probabilities ================
mc_d.traj <- DeltasSequence(Data=obs_data,
                            Predclass=mc_imputedclass,
                            Probs=mc_lambda,
                            Xnames=Xnames)

(filename <- paste0("illustration-",dataset_name,"-perturb-lca-",seed,".Rdata"))
save.image(file=filename)
proc.time()[3]-ptm
q()
