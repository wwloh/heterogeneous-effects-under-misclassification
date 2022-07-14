rm(list=ls())
library("Hmisc")
library("geepack")
library("xtable")
library("glmnet")
library("MatrixModels")
library("poLCA")

source("helper_functions.R")
set.seed(9000)

# load processed data
args <- 1
args <- commandArgs(trailingOnly=TRUE) # for CMD BATCH '--args 1'
(seed <- as.integer(args[1]))
rm(args)
set.seed(9000*seed)

if (seed==1) {
  dataset_name <- "rhc" # RHC example  
} else if (seed==2) {
  dataset_name <- "lindner" # example from twang package
}
dataset_name

## Latent class analysis ======================================================
filename <- paste0("illustration-",dataset_name,"-lca.Rdata")
load(file=filename)

# formula for inclusive latent class model
polca.f

mixLCA_k_res <- lapply(mixLCA_k, function(x) {
  # AIC, BIC
  x.res <- list("AIC"=x$aic,"BIC"=x$bic)
  # proportion in each imputed class
  x.res[["prop"]] <- round(table(x$predclass)/x$Nobs,2)
  # class-specific estimates under estimated class memberships
  if (dataset_name=="rhc") {
    x.res[["est"]] <- DeltasClassSpecific(
      C=x$predclass,X=obs_data[,Xnames],Z=obs_data$Z,Y=obs_data$Y,status=TRUE,
      meth.names=c("AIPW.enet_Y.logit_Z"))
  } else if (dataset_name=="lindner") {
    x.res[["est"]] <- DeltasClassSpecific(
      C=x$predclass,X=obs_data[,Xnames],Z=obs_data$Z,Y=obs_data$Y,status=TRUE)
  }
  return(x.res)
})
lapply(mixLCA_k_res, "[", "prop")
# minimum AIC/BIC (first entry of list is empty)
mixLCA_k_fit <- sapply(paste0(c("A","B"),"IC"), function(ic) {
  c(NA,as.numeric(unlist(lapply(mixLCA_k_res, "[", ic))))
})
apply(mixLCA_k_fit, 2, which.min)

pdf(gsub(".Rdata","-ABIC.pdf",filename),width=6,height=8)
par(mfrow=c(2,1))
for (ic in 1:ncol(mixLCA_k_fit)) {
  plot(mixLCA_k_fit[,ic],
       main=colnames(mixLCA_k_fit)[ic],
       ylab="Criterion",
       xlab="Number of latent classes",
       type="b",xlim=c(2,nrow(mixLCA_k_fit)))
}
dev.off()

# posterior probabilities
lapply(lapply(mixLCA_k, "[[", "posterior"),summary)

# estimates
lapply(mixLCA_k_res, "[", "est")

# select number of latent classes
if (dataset_name=="rhc") {
  mixLCA <- mixLCA_k[[4]]
  est <- mixLCA_k_res[[4]][["est"]]
} else if (dataset_name=="lindner") {
  mixLCA <- mixLCA_k[[2]]
  est <- mixLCA_k_res[[2]][["est"]]
}

# individuals with unit probability of belonging to a particular class
mean(rowSums(mixLCA$posterior > 1-.Machine$double.eps)>0)

# order classes by increasing lower bound of CI for ATE
row.names(est) <- c("ATE","ci.low","ci.high")
(class_order.ATE <- order(est["ci.low",],decreasing=FALSE))
est <- est[,class_order.ATE]
round(est,2)

## predicted class membership (probabilities)
lambda <- mixLCA$posterior[,class_order.ATE]
imputedclass <- as.integer(apply(lambda, 1, which.max))
#### check predicted class identical to hard assignment
sum(diag(table(data.frame(cbind("old"=class_order.ATE[mixLCA$predclass],
                                "new"=imputedclass)))))==n

# proportion in each imputed class
round(table(imputedclass)/mixLCA$Nobs,2)

# summarize characteristics of latent classes =================================
lc_characs <- do.call(rbind,lapply(1:length(mixLCA$probs), function(x) {
  data.frame(
    names(mixLCA$probs)[x],
    t(colSums(t(mixLCA$probs[[x]])*(0:(ncol(mixLCA$probs[[x]])-1)))),
    "pv"=chisq.test(table(imputedclass,obs_data[,names(mixLCA$probs)[x]]))$p.value)
}))
dim(lc_characs)
print.xtable(xtable(lc_characs,digits=3),include.rownames=FALSE)

# trajectory based on posterior class membership probabilities ================
filename <- gsub("lca","lca-subgroupATE",filename)
if (!exists("d.traj")) {
  if (dataset_name=="rhc") {
    d.traj <- DeltasSequence(Data=obs_data,
                             Predclass=imputedclass,
                             Probs=lambda,
                             Xnames=Xnames,
                             meth.names=c("AIPW.enet_Y.logit_Z"))
  } else if (dataset_name=="lindner") {
    d.traj <- DeltasSequence(Data=obs_data,
                             Predclass=imputedclass,
                             Probs=lambda,
                             Xnames=Xnames)
  }
  save.image(file=filename)
}

# bias-corrected estimates (use bootstrap for inference)
ate.bc <- ATE_biascorrected(
  C=imputedclass,X=obs_data[,Xnames],Z=obs_data$Z,Y=obs_data$Y,Probs=lambda)
ate.bc
  
# (weighted) average effect
(ate <- One_estimator(X=obs_data[,Xnames],Z=obs_data$Z,Y=obs_data$Y))

filename <- paste0("illustration-",dataset_name,"-lca-subgroupATE.Rdata")
save.image(file=filename)

q()

# plots of results ############################################################
rm(list=ls())
# dataset_name <- "rhc" # RHC example  
# dataset_name <- "lindner" # example from twang package
filename <- paste0("illustration-",dataset_name,"-lca-subgroupATE.Rdata")
load(file=filename)
# plots =======================================================================
filename.plot <- paste0("plot-",filename)
filename.plot <- gsub("Rdata","pdf",filename.plot)
n_class <- max(imputedclass)
pdf(filename.plot,width=2.5*n_class,height=6)
layout(matrix((1:(n_class*2)), nrow=2, ncol=n_class, byrow = FALSE))
d.traj.std <- lapply(d.traj,function(x) {
  x.std <- x[,1]
  x.std[(x[,"pr.belong"] > 0.999)] <- NA
  return(x.std)
})
x_range <- c(0,1)
for (s in 1:ncol(lambda)) {
  # visualize probabilities of class memberships
  plot(x=(1:n)/n,sort(lambda[,s],decreasing=TRUE,na.last=TRUE),
       main=paste("Class ",s),
       ylab="Posterior prob.",
       xlab="Cumulative proportion",
       ylim=c(0,1),type="l",col=s,lwd=1.5,xlim=x_range)
  abline(v=mean(imputedclass==s),lty=2,col=s)
  y_range <- range(c(est[-1,],range(d.traj.std,na.rm=TRUE)),na.rm=TRUE)
  plot(x=(1:n)/n,d.traj.std[[s]],
       main=paste("Class ",s),
       ylab="Estimate",
       xlab="Cumulative proportion",
       ylim=y_range,type="l",col=s,lwd=1.5,xlim=x_range)
  points(x=mean(imputedclass==s),y=est["ATE",s],cex=1.5,col=s)
  lines(x=rep(mean(imputedclass==s),2), y=est[-1,s],col=s)
  points(x=rep(mean(imputedclass==s),2), y=est[-1,s],pch="-",col=s)
  abline(h=ate[1],col="grey50",lty=3)
  abline(h=ate.bc[s],col=s,lty=2)
}
dev.off()

# Examples of individual subgroups 
if (dataset_name=="rhc") {
  s_i <- c(1,3)
} else if (dataset_name=="lindner") {
  s_i <- 1
}
for (s in s_i) {
  if (dataset_name=="rhc") {
    plot.title <- ifelse(s==1,"Example 2","Example 3")
  } else if (dataset_name=="lindner") {
    plot.title <- "Example 1"
  }
  pdf(gsub("lca",paste0("lca-",gsub(" ","",plot.title)),filename.plot),
      width=4,height=8)
  par(mfrow=c(2,1))
  x_range <- range(which(!is.na(d.traj.std[[s]])))
  x_range[2] <- pmin(sum(lambda[,s]>=0.01),x_range[2])
  x_range <- x_range/n
  # visualize probabilities of class memberships
  plot(x=(1:n)/n,sort(lambda[,s],decreasing=TRUE,na.last=TRUE),
       main=plot.title,
       ylab="Membership probability",
       xlab="Cumulative proportion",
       ylim=c(0,1),type="l",lwd=1.5,xlim=x_range)
  abline(v=mean(imputedclass==s),lty=2)
  y_range <- range(c(est[-1,s],range(d.traj.std[[s]],na.rm=TRUE)),na.rm=TRUE)
  if (min(abs(y_range))==Inf) {
    y_range <- est[-1,s]
  }
  plot(x=(1:n)/n,d.traj.std[[s]],
       main=plot.title,
       ylab="Estimate",
       xlab="Cumulative proportion",
       ylim=y_range,type="l",lwd=1.5,xlim=x_range)
  points(x=mean(imputedclass==s),y=est["ATE",s],cex=1.5)
  lines(x=rep(mean(imputedclass==s),2), y=est[-1,s])
  points(x=rep(mean(imputedclass==s),2), y=est[-1,s],pch="-")
  dev.off()
}

