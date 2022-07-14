rm(list=ls())
library("combinat")

# # load processed data
# dataset_name <- "rhc" # RHC example
# dataset_name <- "lindner" # RHC example

# observed data ===============================================================
filename <- paste0("illustration-",dataset_name,"-lca-subgroupATE.Rdata")
load(file=filename)

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
      x.std <- x[,2:3]
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
  "Perturbed"=as.vector(rbind(est[1,],do.call(cbind,perturbedCIs))))
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
pdf(filename.plot,width=2.5*n_class,height=6)
layout(matrix((1:(n_class*2)), nrow=2, ncol=n_class, byrow = FALSE))
d.traj.std <- lapply(d.traj,function(x) {
  x.std <- x[,1]
  x.std[(x[,"pr.belong"] > 0.99 | x[,"pr.belong"] < 0.01)] <- NA
  return(x.std)
})
x_range <- range(lapply(1:ncol(lambda), function(s) {
    range(which(!is.na(d.traj.std[[s]])))/n
  }))
# x_range <- c(0,1)
for (s in 1:ncol(lambda)) {
  # visualize probabilities of class memberships
  plot(x=(1:n)/n,sort(lambda[,s],decreasing=TRUE,na.last=TRUE),
       main=paste("Class ",s),
       ylab="Posterior prob.",
       xlab="Cumulative proportion",
       ylim=c(0,1),type="l",col=s,lwd=1.5,xlim=x_range)
  abline(v=mean(imputedclass==s),lty=2,col=s)
  y_range <- range(c(est[-1,],range(d.traj.std,na.rm=TRUE)))*c(1.5,1)
  plot(x=(1:n)/n,d.traj.std[[s]],
       main=paste("Class ",s),
       ylab="Estimate",
       xlab="Cumulative proportion",
       ylim=y_range,type="l",col=s,lwd=1.5,xlim=x_range)
  points(x=mean(imputedclass==s),y=est["ATE",s],cex=1.5,col=s)
  lines(x=rep(mean(imputedclass==s),2), y=est[-1,s],col=s)
  points(x=rep(mean(imputedclass==s),2), y=est[-1,s],pch="-",col=s)
  abline(h=ate[1],col="grey50")
  # abline(h=ate.bc[s],col=s,lty=2)
  # perturbed CIs
  abline(h=range(perturbedCIs[[s]]),col=s,lty=3,lwd=2)
}
dev.off()

ate
