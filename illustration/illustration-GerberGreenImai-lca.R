rm(list=ls())
library("Matching")
library("poLCA")
library("geepack")
library("mvtnorm")
library("xtable")

# source("../simulations/helper_functions.R") # from GitHub repository
source("helper_functions.R")
set.seed(9000)

# load data ===================================================================
# demo(topic=GerberGreenImai, package = "Matching")
data(GerberGreenImai)
obs_data <- GerberGreenImai
summary(obs_data)
n <- nrow(obs_data)
Xnames <- c("PERSONS","VOTE96.1","NEW","MAJORPTY","AGE","WARD","AGE2")

obs_data <- data.frame("i"=1:n,
                       obs_data[,Xnames],
                       #treatment phone calls
                       "Z"=as.integer(GerberGreenImai$PHN.C1),
                       #outcome, turnout
                       "Y"=as.integer(GerberGreenImai$VOTED98)) 

#replication of Imai's propensity score matching model ========================
ps.model <- as.formula(Z ~ PERSONS + VOTE96.1 + NEW + MAJORPTY + AGE +
                         WARD + PERSONS:VOTE96.1 + PERSONS:NEW + AGE2)
ps.fit <- glm(ps.model, family=binomial("logit"), data=obs_data)
ps.fit.mat <- model.matrix(ps.fit)[,!is.na(coef(ps.fit))]

# (weighted) average effect
ate <- IPWestimator(X=data.frame(ps.fit.mat[,-1]),Z=obs_data$Z,Y=obs_data$Y)
#   IPW.raw   IPW.gee    ci.low    ci.upp 
# 0.1113282 0.1113282 0.0261894 0.1964670

## Latent class analysis ======================================================
#### discretize continuous covariate
summary(ps.fit.mat[,-1])
age.fact <- cut(ps.fit.mat[,"AGE"], 
                quantile(ps.fit.mat[,"AGE"],probs=(0:5)/5),
                ordered_result=TRUE)
levels(age.fact)
age.fact <- as.integer(age.fact)
ps.fit.mat[is.na(age.fact),"AGE"]
age.fact[is.na(age.fact)] <- 1L

ps.fit.mat.dicot <- sapply(colnames(ps.fit.mat)[-1], function(zname) {
  if (zname=="PERSONS") {
    return( as.integer(ps.fit.mat[,zname]) )
  } else if (zname == "AGE") {
    return( age.fact )
  } else if (zname == "AGE2") {
    return( NULL )
  } else if (zname == "PERSONS:VOTE96.1") {
    return( NULL )
  } else if (zname == "PERSONS:NEW") {
    return( NULL )
  } else {
    return( as.integer(ps.fit.mat[,zname]+1) )
  }
})
ps.fit.mat.dicot <- ps.fit.mat.dicot[!unlist(lapply(ps.fit.mat.dicot, is.null))]
ps.fit.mat.dicot <- data.frame(do.call(cbind,ps.fit.mat.dicot))
apply(ps.fit.mat.dicot, 2, table)

# fit basic latent class model
polca.f <- as.formula(
  paste0("cbind(",paste(colnames(ps.fit.mat.dicot),collapse=","),")~1"))
polca.f
# Estimate multiple times using different initial parameter values
mixLCA <- poLCA(polca.f, data=ps.fit.mat.dicot, nclass=2, nrep=20, verbose=TRUE)
lapply(mixLCA$probs,xtable)
## predicted class membership probabilities
table(round(mixLCA$posterior,1))
#     0  0.1  0.2  0.3  0.4  0.5  0.6  0.7  0.8  0.9    1 
# 8056  671  661  576  517  696  517  576  661  671 8056 

mean(rowSums(apply(mixLCA$posterior,2,function(x) x<.Machine$double.eps))==0)
# 0.2657678

round(table(mixLCA$predclass)/mixLCA$Nobs,2)
# 1    2 
# 0.37 0.63

# class-specific estimates under estimated class memberships
est <- DeltasClassSpecific(
  C=mixLCA$predclass,X=data.frame(ps.fit.mat[,-1]),Z=obs_data$Z,Y=obs_data$Y)
est

## predicted class membership probabilities
lambda <- mixLCA$posterior
lambda1 <- lambda[,1]

# visualize probabilities of class memberships
plot(lambda[,1],
     main="Ordered class membership probabilities",
     ylab="Probability",
     xlab="Individuals ordered by their probability of being in class 1",
     ylim=c(0,1),type="n")
for (s in 1:2) {
  lines(lambda[order(lambda[,1],decreasing=TRUE),s],col=s)
}
legend("bottomleft",legend=paste("CLASS = ",1:2),col=1:2,lty=1,cex=.75,bty="n")
rm(ps.fit)

obs_data[,"class"] <- mixLCA$predclass
table("treat"=obs_data$Z,"out"=obs_data$Y,"class"=obs_data$class)

obs_Data <- cbind(data.frame(ps.fit.mat[,-1]),Z=obs_data$Z,Y=obs_data$Y)
head(obs_Data)
rm(obs_data)


# trajectory based on estimated class membership probabilities ================
# ordered by decreasing probability of membership to class 1
ptm <- proc.time()[3]
deltas.l1o <- DeltasSequence(Data=obs_Data[order(lambda1,decreasing=TRUE),])
save(deltas.l1o, file="GerberGreenImai-2class-lca.Rdata")
proc.time()[3]-ptm
# 11298.65/60: 188 mins

# results =====================================================================
load(file="GerberGreenImai-2class-lca.Rdata")
deltas.l1o.all <- deltas.l1o
deltas.l1o <- deltas.l1o[row.names(deltas.l1o)=="IPW.gee",]
# visualize probabilities of class memberships
pdf("plot-illustration-traj.pdf",width=6,height=8)
par(mfrow=c(2,1))
plot(lambda[,1],
     main="Ordered class membership probabilities",
     ylab="Probability",
     xlab="Individuals ordered by their probability of being in class 1",
     ylim=c(0,1),type="n")
for (s in 1:2) {
  lines(lambda[order(lambda[,1],decreasing=TRUE),s],col=s)
}
legend("bottomleft",legend=paste("CLASS = ",1:2),col=1:2,lty=1,cex=.75,bty="n")
abline(v=sum(lambda[,1]>0.5),col="grey50",lwd=0.5)
lambda1.ordered <- lambda1[order(lambda1,decreasing=TRUE)]
my_xlim <- range(which(pmin(lambda1.ordered,1-lambda1.ordered) > 1/n))
my_ylim <- range(deltas.l1o[my_xlim[1]:my_xlim[2],])
plot(deltas.l1o[,1],
     main="Class-specific treatment effects",
     ylab="Treatment effect estimates",
     xlab="Assumed number of individuals in class 1",
     ylim=my_ylim,xlim=my_xlim,
     type="n",col=1)
for (s in 1:2) {
  lines(my_xlim[1]:my_xlim[2],deltas.l1o[my_xlim[1]:my_xlim[2],s],type="s",col=s)  
}
abline(h=ate["IPW.gee"],col="grey50")
abline(h=0,lty=2,col="grey50")
abline(h=est["IPW.gee",1],lty=3,col=1)
abline(h=est["IPW.gee",2],lty=3,col=2)
abline(v=sum(lambda[,1]>0.5),col="grey50",lwd=0.5)
legend("bottomleft",legend=paste("CLASS = ",1:2),col=1:2,lty=1,cex=.75,bty="n")
dev.off()

  
# calculate perturbed class membership probabilities ==========================
resampled_list <- NULL
resample <- TRUE
m <- 1000 # number of perturbations
while(resample) {
  # create simulated data with fitted model’s assumed data-generating process 
  mc_sim <- poLCA::poLCA.simdata(N=n,probs=mixLCA$probs,P=mixLCA$P)
  mc_dat <- mc_sim$dat
  colnames(mc_dat) <- colnames(mixLCA$y)
  # fit latent class analysis
  mc_fit <- poLCA(polca.f, data=mc_dat, nclass=2, verbose=TRUE)
  # calculate posterior probabilities for observed data
  # using randomly drawn parameter values from sampling distribution
  mc_dens <- poLCA::poLCA.posterior(lc=mc_fit,y=mixLCA$y)
  
  # optimal Bayes rule
  mc_C <- as.integer(apply(mc_dens, 1, which.max))
  if (length(unique(mc_C))==2) {
    # class-specific estimates under resampled class memberships
    mc_est <-DeltasClassSpecific(
      C=mc_C,X=obs_Data[,!(colnames(obs_Data) %in% c("Z","Y"))],
      Z=obs_Data$Z,Y=obs_Data$Y)
    # switch labels so that class 1 effect < class 2 effect
    mc_dens <- mc_dens[,order(mc_est["IPW.gee",])]
    
    resampled_list <- c(resampled_list,list(mc_dens))
    resample <- (length(resampled_list) < m)
    cat(length(resampled_list),"| order =",order(mc_est["IPW.gee",]), 
        "################################################################## \n")
  }
}

save(obs_Data,resampled_list,file="GerberGreenImai-2class-perturbed_lca.Rdata")

# parallel jobs to calculate trajetory for each perturbation ==================
# illustration-GerberGreenImai-lca-perturbed.R

# combine results =============================================================
deltas.l1o.boot <- NULL
subfolder <- "lca-pb/"
myfiles <- list.files(subfolder)
myfiles <- grep(pattern=".Rdata",myfiles,value=TRUE)
for (ll in myfiles) {
  load(paste0(subfolder,ll))
  deltas.l1o.boot <- c(deltas.l1o.boot, list(
    deltas.l1o[row.names(deltas.l1o)=="IPW.gee",]))
  rm(deltas.l1o)
}
deltas.l1o.boot <- deltas.l1o.boot[unlist(lapply(deltas.l1o.boot, function(x) 
  max(abs(x),na.rm=TRUE)<1))]
length(deltas.l1o.boot)

deltas.l1o.bootrange <- list()
for (s in 1:2) {
  deltas.l1o.s <- do.call(cbind,lapply(deltas.l1o.boot, function(x) x[,s]))
  deltas.l1o.bootrange[[s]] <- t(apply(deltas.l1o.s, 1, function(x) {
    x.range <- rep(NA,2)
    if (any(!is.na(x))) {
      # x.range <- quantile(x, probs=c(0.005, 0.995), na.rm=TRUE)
      x.range <- range(x, na.rm=TRUE)
    }
    return(x.range)
  }))
}

# results =====================================================================
load(file="GerberGreenImai-2class-lca.Rdata")
deltas.l1o <- deltas.l1o[row.names(deltas.l1o)=="IPW.gee",]
lambda1.ordered <- lambda1[order(lambda1,decreasing=TRUE)]
my_xlim <- range(which(pmin(lambda1.ordered,1-lambda1.ordered) > 1/n))
my_ylim <- range(lapply(deltas.l1o.bootrange, function(x)
  range(x[my_xlim[1]:my_xlim[2],],na.rm=TRUE)))
pdf("plot-illustration-traj-mc.pdf",width=6,height=8)
par(mfrow=c(2,1))
for (s in 1:2) {
  plot(deltas.l1o[,s],
       main=paste0("Class ",s),
       ylab="Treatment effect estimates",
       xlab="Assumed number of individuals in class 1",
       ylim=my_ylim,xlim=my_xlim,type="n",col=1)
  y.poly.l <- deltas.l1o.bootrange[[s]][my_xlim[1]:my_xlim[2],1]
  y.poly.u <- deltas.l1o.bootrange[[s]][my_xlim[1]:my_xlim[2],2]
  x.poly <- c((my_xlim[1]:my_xlim[2])[!is.na(y.poly.l)],
              rev((my_xlim[1]:my_xlim[2])[!is.na(y.poly.u)]))
  y.poly <- c(y.poly.l[!is.na(y.poly.l)],rev(y.poly.u[!is.na(y.poly.u)]))
  polygon(x=x.poly,y=y.poly,col="grey90",border=NA,fillOddEven=TRUE)
  lines(my_xlim[1]:my_xlim[2],deltas.l1o[my_xlim[1]:my_xlim[2],s],type="s",col=s)
  abline(h=ate["IPW.gee"],col="grey50")
  abline(h=0,lty=2,col="grey50")
  abline(h=est["IPW.gee",s],lty=3,col=s)
}
dev.off()
