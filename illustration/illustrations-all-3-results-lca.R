rm(list=ls())
library("combinat")

# load processed data
dataset_name <- "rhc" # RHC example
# dataset_name <- "GerberGreenImai" # Matching example

# observed data ===============================================================
library("geepack")
filename <- paste0("illustration-",dataset_name,"-lca.Rdata")
load(file=filename)
(ate <- IPWestimator(X=obs_data[,Xnames],Z=obs_data$Z,Y=obs_data$Y))

# combine results from parallel jobs for calculating each perturbation ========
deltas.l1o.boot <- NULL
est.boot <- NULL
subfolder <- "illustrations-all-2-perturb-lca/"
myfiles <- list.files(subfolder)
myfiles <- grep(pattern=dataset_name,myfiles,value=TRUE)
for (ll in myfiles) {
  temp <- new.env()
  # load MC objects into temporary workspace
  with(temp, {
    load(paste0(subfolder,ll))
    mc_d.traj.std <- lapply(mc_d.traj,function(x) {
      x.std <- x[,c("ci.low","ci.upp")]
      x.std[(x[,"pr.belong"] > 0.99) | (x[,"pr.belong"] < 0.01),] <- rep(NA,2)
      return(x.std)
    })
    ## all permutations of imputed class size matrix
    all_perms <- permn(class_order.mc_est)
    ## permutation that maximizes similarity to imputed class size
    max_sim <- which.max(unlist(lapply(all_perms, function(class_perm) {
      sum(diag(table("obs"=imputedclass,"mc"=mc_imputedclass)[,class_perm]))
    })))
    chosen_perm <- all_perms[max_sim][[1]]
    class_order.mc_est <- class_order.mc_est[chosen_perm]
    # reorder class estimates
    mc.est <- mc_est[,class_order.mc_est]
    mc.est <- rbind(mc.est,
                    "imputed"=table(mc_imputedclass)[chosen_perm]/
                      length(mc_imputedclass))
    mc_d.traj.std <- mc_d.traj.std[chosen_perm]
  })
  deltas.l1o.boot <- c(deltas.l1o.boot,list(with(temp,mc_d.traj.std)))
  est.boot <- c(est.boot,list(with(temp,mc.est)))
  rm(temp)
  cat(ll,"\n")
}
length(deltas.l1o.boot)

# perturbed CIs
perturbedCIs <- lapply(1:ncol(lambda), function(s) {
  range(apply(do.call(rbind,lapply(est.boot, function(x) x[2:3,s])),
              2,quantile, probs=c(.025,.975),na.rm=TRUE))
})
lapply(perturbedCIs, range)

perturbed_classsize <- lapply(1:ncol(lambda), function(s) {
  matrix(c(
    mean(imputedclass==s),
    quantile(unlist(lapply(est.boot, function(x) 
      x[4,s])),probs=c(.025,.975),na.rm=TRUE)),nrow=1)
})
summ_perturb <- rbind(
  "Proportion"=do.call(cbind,perturbed_classsize)[1,],
  "Fixed"=as.vector(est),
  "Perturbed"=as.vector(rbind(est["IPW",],do.call(cbind,perturbedCIs))))
summ_perturb_formatted <- t(apply(summ_perturb,1,function(x) {
  x_ <- NULL
  for (i in 1:length(x)) {
    x_i <- format(round(x[i], 2), nsmall = 2)
    if (i %% 3 == 1) x_[i] <- x_i
    if (i %% 3 == 2) x_[i] <- paste0("(",x_i,",")
    if (i %% 3 == 0) x_[i] <- paste0(x_i,")")
  }
  return(x_)
}))
library(xtable)
xtable(summ_perturb_formatted)

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
y_range <- range(c(est[-1,],range(d.traj.std,na.rm=TRUE)),na.rm=TRUE)
y_range <- range(c(y_range,range(perturbedCIs)))
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
  plot(x=(1:n)/n,d.traj.std[[s]],
       main=paste("Class ",s),
       ylab="Estimate",
       xlab="Cumu. prop.",
       ylim=y_range,type="l",col=s,lwd=1.5,xlim=x_range)
  points(x=mean(imputedclass==s),y=est["IPW",s],pch=19,col=s,cex=2)
  lines(x=rep(mean(imputedclass==s),2), y=est[-1,s],col=s)
  points(x=rep(mean(imputedclass==s),2), y=est[-1,s],pch="-",col=s)
  abline(h=ate["IPW"],col="grey50",lty=2)
  # perturbed CIs
  abline(h=range(perturbedCIs[[s]]),col=s,lty=3,lwd=2)
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
  pdf(gsub("lca",paste0("lca-perturb",gsub(" ","",plot.title)),filename.plot),
      width=4,height=4)
  xy_range <- do.call(rbind,lapply(est.boot, function(oneboot) {
    c(oneboot["imputed",s],oneboot[c("ci.low","ci.upp"),s])
  }))
  x_range <- range(xy_range[rowSums(is.na(xy_range))==0,1], na.rm=TRUE)
  y_range <- range(xy_range[,-1],na.rm=TRUE)
  plot(x=x_range,y_range,
       main=plot.title,
       ylab="Estimate",
       xlab="Cumu. prop.",
       ylim=y_range,type="n",xlim=x_range)
  # perturbed CIs
  abline(h=range(perturbedCIs[[s]]),lty=3, lwd=2)
  for (bb in 1:length(est.boot)) {
    if (plot.title=="Example 1") {
      xval <- runif(1,est.boot[[bb]]["imputed",s]*.99,
                    est.boot[[bb]]["imputed",s]*1.01)# random jitter  
    } else {
      xval <- est.boot[[bb]]["imputed",s]
    }
    lines(rep(xval,2),est.boot[[bb]][c("ci.low","ci.upp"),s],col="grey",lwd=1)
  }
  # based on imputed class memberships
  points(x=mean(imputedclass==s),y=est["IPW",s],pch=19,col=1,cex=2)
  lines(x=rep(mean(imputedclass==s),2), y=est[-1,s])
  points(x=rep(mean(imputedclass==s),2), y=est[-1,s],pch="-")
  dev.off()
}
