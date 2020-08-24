# read in preamble ============================================================
rm(list=ls())
source2 <- function(file, start, end, ...) {
  file.lines <- scan(file, what=character(), 
                     skip=start-1, nlines=end-start+1, sep='\n')
  file.lines.collapsed <- paste(file.lines, collapse='\n')
  source(textConnection(file.lines.collapsed), ...)
}
source2("sim-uncertainty-latentclass.R",2,12)

# load saved results ===========================================================
library("xtable")
subfolder <- "sim-results/"
myfiles <- list.files(subfolder)
myfiles <- grep(pattern=".Rdata",myfiles,value=TRUE)
res_all <- vector(mode = "list", length = nrow(simsettings))
res_cm <- NULL

for (ll in myfiles) {
  load(paste0(subfolder,ll))
  # sim setting
  ll_sim <- unique(do.call(rbind,lapply(res, "[[", "simsetting")))
  ll_sim_idx <- which(apply(apply(simsettings, 1, "==", ll_sim),2,all))
  # confusion matrices
  row.names(ll_sim) <- NULL
  cm1 <- cbind("type"="fixed",ll_sim,
               do.call(rbind,lapply(
                 lapply(res, "[[", "confusion_fixed"),as.vector)))
  res_cm <- c(res_cm,list(cm1))
  cm2 <- cbind("type"="boot",ll_sim,
               do.call(rbind,lapply(
                 lapply(res, "[[", "confusion_boot"),as.vector)))
  res_cm <- c(res_cm,list(cm2))
  rm(cm1,cm2)
  res <- lapply(res, function(x) {
    x[c("simsetting","confusion_fixed","confusion_boot")] <- NULL
    x
  })
  res_all[[ll_sim_idx]] <- c(res_all[[ll_sim_idx]],res)
  rm(res,ll_sim,ll_sim_idx)
}

# confusion matrices
res_cm <- do.call(rbind,res_cm)
res_cm <- cbind(res_cm,
                "correct"=rowSums(res_cm[, c(4,7)]),
                "incorrect"=rowSums(res_cm[, c(5,6)]))
res_cm[,c("1","2","3","4")] <- NULL
res_cm_summ <- by(res_cm[,4:5], INDICES=res_cm[,1:3], function(x) 
  round(c("mean"=colMeans(x),"sd"=apply(x,2,sd)),2))
res_cm_summ <- cbind(expand.grid(dimnames(res_cm_summ)),
                     do.call(rbind,res_cm_summ))
rm(res_cm)

# consider each unique simulation setting in turn =============================
for (ss in 1:length(res_all)) {
  res_list <- res_all[[ss]]
  nas <- unlist(lapply(res_list, function(x) any(unlist(lapply(x,is.na)))))
  res <- res_list[!nas]
  rm(res_list)
  length(res)
  
  # calculate population effects
  meth <- as.character(simsettings[ss,"meth"])
  Lambda <- simsettings[ss,"Lambda"]
  source2("sim-uncertainty-latentclass.R",27,92)
  
  # misclassification rates
  res_cm <- res_cm_summ[res_cm_summ$meth==meth & res_cm_summ$Lambda==Lambda,]
  
  # summarize results
  res_summary <- lapply(res, function(onesim) {
    unlist(lapply(onesim, function(onesim_meth) {
      c("est"=onesim_meth["IPW.gee",],
        "bias"=onesim_meth["IPW.gee",]-pop_effs["IPW.gee",],
        "cover"=sapply(1:2, function(cl) {
          as.integer((onesim_meth["ci.low",cl] <= pop_effs["IPW.gee",cl]) &
                       (onesim_meth["ci.upp",cl] >= pop_effs["IPW.gee",cl]))
        }),
        "halfwidth"=sapply(1:2, function(cl) {
          as.numeric((onesim_meth["ci.upp",cl]-onesim_meth["ci.low",cl])/2)
        })
      )
    }))
  })
  res_summary <- do.call(rbind,res_summary)
  colnames(res_summary)
  
  print(simsettings[ss,])
  print(pop_effs["IPW.gee",])
  print(ate)
  print(xtable(t(sapply(c("true","est[.]","boot"), function(mm) {
    mm_av <- colMeans(res_summary[,grep(mm,colnames(res_summary),value=TRUE)])
    mm_se <- apply(res_summary[,grep(mm,colnames(res_summary),value=TRUE)],2,sd)
    if (mm=="true") {
      mm_cm <- rep(NA,2)  
    } else {
      res_cm_mm <- mm
      if (mm=="est[.]") {
        res_cm_mm <- "fixed"
      }
      mm_cm <- res_cm[res_cm$type==res_cm_mm,c("mean.correct","sd.correct")]
    }
    c(mm_av[-c(1:2)],mm_se[grep("[.]est",names(mm_se))],mm_cm)
  }))))
  
  meths <- c("true.est","est.est","boot.ci.est")
  names(meths) <- c("True","Fixed","Perturbed")
  est_density <- apply(res_summary[,sapply(meths, function(mm) 
    grep(mm,colnames(res_summary),value=TRUE))],2,density)
  my_xlim <- lapply(est_density, "[[", "x")
  my_xlim <- sapply(1:2, function(cs) range(my_xlim[grep(cs,names(my_xlim))]))
  my_ylim <- lapply(est_density, "[[", "y")
  my_ylim <- sapply(1:2, function(cs) range(my_ylim[grep(cs,names(my_ylim))]))
  
  filename <- gsub("[.]","_",paste0("sim-",meth,"-",Lambda))
  pdf(paste0(filename,".pdf"),height=9,width=9)
  par(mfrow=c(3,3))
  for (mm in meths) {
    res_bias <- res_summary[,grep(mm,colnames(res_summary),value=TRUE)]
    plot(res_bias,
         main=names(meths)[meths==mm],
         col="grey50",pch=19,cex=.5,
         xlab="Class 1",ylab="Class 2",xlim=my_xlim[,1],ylim=my_xlim[,2])
    abline(v=pop_effs["IPW.gee",1],lty=2,lwd=2)
    abline(h=pop_effs["IPW.gee",2],lty=2,lwd=2)
    for (s in 1:2) {
      plot(density(res_bias[,s]),
           xlim=my_xlim[,s],ylim=my_ylim[,s],
           xlab="Estimated effects",
           main=paste0(names(meths)[meths==mm]," (Class ", s,")"))
      abline(v=pop_effs["IPW.gee",s],lty=2,lwd=2)
    }
  }
  dev.off()
}
