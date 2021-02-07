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
res_cm[,c("2","3")] <- NULL
colnames(res_cm)[4:5] <- paste0("est.class.",1:2)
res_cm$est.class.1 <- res_cm$est.class.1/(1-res_cm$Lambda)
res_cm$est.class.2 <- res_cm$est.class.2/res_cm$Lambda
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
  
  # calculate population effects
  meth <- as.character(simsettings[ss,"meth"])
  Lambda <- simsettings[ss,"Lambda"]
  source2("sim-uncertainty-latentclass.R",27,92)
  
  # correctly classified in each class
  res_cm <- res_cm_summ[res_cm_summ$meth==meth & res_cm_summ$Lambda==Lambda,]
  
  # summarize results
  res_summary <- lapply(res, function(onesim) {
    unlist(lapply(onesim, function(onesim_meth) {
      c("est"=onesim_meth[1,],
        "bias"=onesim_meth[1,]-pop_effs["IPW",],
        "cover"=sapply(1:2, function(cl) {
          as.integer(
            (onesim_meth[grep("ci.low",row.names(onesim_meth)),cl] 
             <= pop_effs["IPW",cl]) &
              (onesim_meth[grep("ci.upp",row.names(onesim_meth)),cl] 
               >= pop_effs["IPW",cl]))
        }),
        "halfwidth"=sapply(1:2, function(cl) {
          as.numeric(
            (onesim_meth[grep("ci.upp",row.names(onesim_meth)),cl]-
               onesim_meth[grep("ci.low",row.names(onesim_meth)),cl])/2)
        })
      )
    }))
  })
  res_summary <- do.call(rbind,res_summary)
  colnames(res_summary)
  
  cat("---------------------------------------------------------------------\n")
  print(simsettings[ss,])
  print(length(res))
  print(pop_effs["IPW",])
  print(ate)
  print(xtable(t(sapply(c("true","est[.]","boot.ci.95"), function(mm) {
    mm_av <- colMeans(res_summary[,grep(mm,colnames(res_summary),value=TRUE)])
    mm_se <- apply(res_summary[,grep(mm,colnames(res_summary),value=TRUE)],2,sd)
    mm_mse <- sqrt(
      (mm_av[grep("bias",names(mm_av),value=TRUE)]^2+
         mm_se[grep("[.]est",names(mm_av),value=TRUE)]^2))
    if (mm=="true") {
      mm_cm <- rep(NA,2)  
    } else {
      if (mm=="est[.]") {
        res_cm_mm <- "fixed"
      } else {
        res_cm_mm <- "boot"
      }
      mm_cm <- res_cm[res_cm$type==res_cm_mm,c("mean.est.class.1",
                                               "mean.est.class.2")]
    }
    # results table
    as.numeric(c(
      mm_av[grep("bias",names(mm_av),value=TRUE)],
      # mm_se[grep("[.]est",names(mm_av),value=TRUE)],
      # mm_mse,
      mm_av[grep("cover",names(mm_av),value=TRUE)],
      mm_av[grep("halfwidth",names(mm_av),value=TRUE)],
      mm_cm))
  }))))
}
