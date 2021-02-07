rm(list=ls())
library("Hmisc")
library("geepack")
library("poLCA")
library("Matching")
library("mvtnorm")
library("xtable")

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
  dataset_name <- "GerberGreenImai" # Matching example   
}
dataset_name

filename <- paste0("illustration-",dataset_name,"-lca.R")
source(file=filename)
filename <- gsub(".R",".Rdata",filename)

# formula for inclusive latent class model
polca.f

## Latent class analysis ======================================================
# fit each model with different number of classes in turn
if (file.exists(filename)) {
  load(file=filename)
} else {
  mixLCA_k <- list()
  for (k in 2:5) {
    # Estimate multiple times using different initial parameter values
    mixLCA_k[[k]] <- poLCA(polca.f,data=obs_data,nclass=k,nrep=500,verbose=TRUE)
    save.image(file=filename)
  }
}

mixLCA_k_res <- lapply(mixLCA_k, function(x) {
  # AIC, BIC
  x.res <- list("AIC"=x$aic,"BIC"=x$bic)
  # proportion in each imputed class
  x.res[["prop"]] <- round(table(x$predclass)/x$Nobs,2)
  # class-specific estimates under estimated class memberships
  x.res[["est"]] <- DeltasClassSpecific(
    C=x$predclass,X=obs_data[,Xnames],Z=obs_data$Z,Y=obs_data$Y)
  return(x.res)
})
lapply(mixLCA_k_res, "[", "prop")
# minimum AIC/BIC (first entry of list is empty)
which.min(unlist(lapply(mixLCA_k_res, "[", "AIC")))+1
which.min(unlist(lapply(mixLCA_k_res, "[", "BIC")))+1
# estimates
lapply(mixLCA_k_res, "[", "est")

if (dataset_name=="rhc") {
  mixLCA <- mixLCA_k[[5]]
  est <- mixLCA_k_res[[5]][["est"]]
} else if (dataset_name=="GerberGreenImai") {
  mixLCA <- mixLCA_k[[4]]
  est <- mixLCA_k_res[[4]][["est"]]
}

# individuals with unit probability of belonging to a particular class
mean(rowSums(mixLCA$posterior > 1-.Machine$double.eps)>0)

# order classes by increasing lower bound of CI for ATE
(class_order.ATE <- order(est["ci.low",],decreasing=FALSE))
est <- est[,class_order.ATE]
round(est,2)

## predicted class membership (probabilities)
lambda <- mixLCA$posterior[,class_order.ATE]
imputedclass <- as.integer(apply(lambda, 1, which.max))
table(data.frame(cbind("old"=mixLCA$predclass,"new"=imputedclass)))

# proportion in each imputed class
round(table(imputedclass)/mixLCA$Nobs,2)

# summary of latent classes ===================================================
names(mixLCA$probs) <- Lnames$L.names[
  paste0("L.",Lnames$L.idx) %in% names(mixLCA$probs)]
lapply(mixLCA$probs,function(x) {
  x.reord <- x[class_order.ATE,]
  row.names(x.reord) <- paste("class", 1:nrow(x.reord))
  xtable(x.reord,)
})
# multinomial logit coefficient estimates 
# Rows correspond to concomitant variables X. 
# Columns correspond to the second through Rth latent classes
# All logit coefficients are calculated for classes with respect to class 1.
pred_coefs <- data.frame(matrix(NA,nrow=nrow(mixLCA$coeff),
                                ncol=ncol(mixLCA$coeff)+1))
for (i in 1:nrow(mixLCA$coeff)) {
  for (j in 1:ncol(mixLCA$coeff)) {
    pred_coefs[i,j+1] <- paste0(
      format(round(mixLCA$coeff[i,j], 2), nsmall = 2)," (",
      format(round(mixLCA$coeff.se[i,j], 2), nsmall = 2),")")
  }
}
row.names(pred_coefs)[-1] <- Lnames$L.names[
  paste0("L.",Lnames$L.idx) %in% row.names(mixLCA$coeff)[-1]]
xtable(pred_coefs[,class_order.ATE])
# distribution of concomittant predictors across classes
for (i in 2:nrow(mixLCA$coeff)) {
  x_idx <- row.names(mixLCA$coeff)[i]
  x_pred <- mixLCA$x[,x_idx]
  if (length(unique(x_pred))>5) {
    boxplot(x_pred~imputedclass, 
            xlab="Imputed class",
            ylab=Lnames[paste0("L.",Lnames$L.idx)==x_idx,"L.names"])
  } else {
    print(xtable(table(imputedclass,x_pred)/rowSums(table(imputedclass,x_pred))))
  }
}

# trajectory based on posterior class membership probabilities ================
if (!exists("d.traj")) {
  d.traj <- DeltasSequence(Data=obs_data,
                           Predclass=imputedclass,
                           Probs=lambda,
                           Xnames=Xnames)
  save.image(file=filename)
}

# (weighted) average effect
(ate <- IPWestimator(X=obs_data[,Xnames],Z=obs_data$Z,Y=obs_data$Y))

# plots =======================================================================
filename.plot <- paste0("plot-",filename)
filename.plot <- gsub("Rdata","pdf",filename.plot)
n_class <- max(imputedclass)
pdf(filename.plot,width=2*n_class,height=5)
layout(matrix((1:(n_class*2)), nrow=2, ncol=n_class, byrow = FALSE))
d.traj.std <- lapply(d.traj,function(x) {
  x.std <- x[,"IPW"]
  x.std[(x[,"pr.belong"] > 0.99) | (x[,"pr.belong"] < 0.01)] <- NA
  return(x.std)
})
for (s in 1:ncol(lambda)) {
  if(all(is.na(d.traj.std[[s]]))) {
    x_range <- rep(mean(imputedclass==s),2)
  } else {
    x_range <- range(which(!is.na(d.traj.std[[s]])))/n
  }
  # visualize probabilities of class memberships
  plot(x=(1:n)/n,sort(lambda[,s],decreasing=TRUE,na.last=TRUE),
       main=paste("Class ",s),
       ylab="Posterior prob.",
       xlab="Cumu. prop.",
       ylim=c(0,1),type="l",col=s,lwd=1.5,xlim=x_range)
  points(x=mean(imputedclass==s),
         y=min(lambda[((imputedclass==s)*lambda[,s])>0,s]),
         pch=19,col=s,cex=2)
  y_range <- range(c(est[-1,],range(d.traj.std,na.rm=TRUE)),na.rm=TRUE)
  plot(x=(1:n)/n,d.traj.std[[s]],
       main=paste("Class ",s),
       ylab="Estimate",
       xlab="Cumu. prop.",
       ylim=y_range,type="l",col=s,lwd=1.5,xlim=x_range)
  points(x=mean(imputedclass==s),y=est["IPW",s],pch=19,col=s,cex=2)
  lines(x=rep(mean(imputedclass==s),2), y=est[-1,s],col=s)
  points(x=rep(mean(imputedclass==s),2), y=est[-1,s],pch="-",col=s)
  abline(h=ate["IPW"],col="grey50",lty=2)
}
dev.off()

# Examples 
if (dataset_name=="rhc") {
  s_i <- c(5,2)
} else {
  s_i <- 1
}
for (s in s_i) {
  if (dataset_name=="rhc") {
    plot.title <- ifelse(s==2,"Example 3","Example 2")
  } else {
    plot.title <- "Example 1"
  }
  pdf(gsub("lca",paste0("lca-",gsub(" ","",plot.title)),filename.plot),
      width=4,height=8)
  par(mfrow=c(2,1))
  if(all(is.na(d.traj.std[[s]]))) {
    x_range <- rep(mean(imputedclass==s),2)
  } else {
    x_range <- range(which(!is.na(d.traj.std[[s]])))/n
  }
  # visualize probabilities of class memberships
  plot(x=(1:n)/n,sort(lambda[,s],decreasing=TRUE,na.last=TRUE),
       main=plot.title,
       ylab="Membership probability",
       xlab="Cumulative proportion",
       ylim=c(0,1),type="l",col=1,lwd=1.5,xlim=x_range)
  points(x=mean(imputedclass==s),
         y=min(lambda[((imputedclass==s)*lambda[,s])>0,s]),
         pch=19,col=1,cex=2)
  y_range <- range(c(est[-1,s],range(d.traj.std[[s]],na.rm=TRUE)),na.rm=TRUE)
  if (min(abs(y_range))==Inf) {
    y_range <- est[-1,s]
  }
  plot(x=(1:n)/n,d.traj.std[[s]],
       main=plot.title,
       ylab="Estimate",
       xlab="Cumulative proportion",
       ylim=y_range,type="l",col=1,lwd=1.5,xlim=x_range)
  points(x=mean(imputedclass==s),y=est["IPW",s],pch=19,col=1,cex=2)
  lines(x=rep(mean(imputedclass==s),2), y=est[-1,s])
  points(x=rep(mean(imputedclass==s),2), y=est[-1,s],pch="-")
  dev.off()
}

