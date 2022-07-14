# read in preamble ============================================================
rm(list=ls())
source2 <- function(file, start, end, ...) {
  file.lines <- scan(file, what=character(), 
                     skip=start-1, nlines=end-start+1, sep='\n')
  file.lines.collapsed <- paste(file.lines, collapse='\n')
  source(textConnection(file.lines.collapsed), ...)
}
source2("sim-uncertainty-latentclass.R",10,14)
rm(source2)

# load saved results ===========================================================
library("xtable")
subfolder <- "sim2-trajectories/"
myfiles <- list.files(subfolder)
myfiles <- grep(pattern=".Rdata",myfiles,value=TRUE)
res_all <- vector(mode = "list", length = nrow(simsettings))
res_cm <- NULL # confusion matrices

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
  rm(cm1)
  res <- lapply(res, function(x) {
    x[c("simsetting","confusion_fixed")] <- NULL
    x
  })
  res_all[[ll_sim_idx]] <- c(res_all[[ll_sim_idx]],res)
  rm(res,ll_sim,ll_sim_idx)
}

# confusion matrices
res_cm <- do.call(rbind,res_cm)
res_cm[,c("2","3")] <- NULL
colnames(res_cm)[ncol(res_cm)-(1:0)] <- paste0("est.class.",1:2)
res_cm$est.class.1 <- 1-abs(res_cm$est.class.1-(1-res_cm$Lambda))
res_cm$est.class.2 <- 1-abs(res_cm$est.class.2-res_cm$Lambda)
res_cm_summ <- by(res_cm[,ncol(res_cm)-(1:0)], INDICES=res_cm[,1:4], function(x) 
  round(c("mean"=colMeans(x),"sd"=apply(x,2,sd)),2))
res_cm_summ <- cbind(expand.grid(dimnames(res_cm_summ)),
                     do.call(rbind,res_cm_summ))
rm(res_cm)

# consider each unique simulation setting in turn =============================
for (ss in 1:length(res_all)) {
  res_list <- res_all[[ss]]
  # nas <- unlist(lapply(res_list, function(x) any(unlist(lapply(x,is.na)))))
  # res <- res_list[!nas]
  # res <- res[!unlist(lapply(res,is.null))]
  # rm(res_list)
  res <- res_list
  
  # population effects
  meth <- as.character(simsettings[ss,"meth"])
  Lambda <- simsettings[ss,"Lambda"]
  direct <- simsettings[ss,"direct"]
  pop_effs <- unique(do.call(rbind,lapply(res, "[[", "pop")))[1,]
  
  # correctly classified in each class
  res_cm <- res_cm_summ[res_cm_summ$meth==meth & res_cm_summ$Lambda==Lambda &
                          res_cm_summ$direct==direct,]
  
  # summarize results
  ## estimators
  est_types <- unlist(strsplit(grep("ATE",row.names(res[[1]]$true),value=TRUE),
                               split=".ATE"))
  ## methods to determine class memberships
  meths_type <- expand.grid("estci"=c("true","est"),
                            "et"=est_types,
                            stringsAsFactors=FALSE)
  
  res <- res[unlist(lapply(res, function(x)
    any(rowSums(sapply(est_types, function(x_) 
      grepl(x_,row.names(x$est))))==1)))]
  
  res_summary <- lapply(res, function(onesim) {
    onesim.res <- NULL
    for (mti in 1:nrow(meths_type)) {
      estci <- meths_type[mti,"estci"]
      et <- meths_type[mti,"et"]
      onesim.estci <- onesim[[estci]]
      onesim_meth_et <- onesim.estci[grep(et,row.names(onesim.estci)),]
      onesim.res[[mti]] <- c(
        "est"=onesim_meth_et[1,],
        "bias"=onesim_meth_et[1,]-pop_effs,
        "cover"=sapply(1:2, function(cl) {
          as.integer(
            (onesim_meth_et[grep("ci.low",row.names(onesim_meth_et)),cl] 
             <= pop_effs[cl]) &
              (onesim_meth_et[grep("ci.upp",row.names(onesim_meth_et)),cl] 
               >= pop_effs[cl]))
        }),
        "halfwidth"=sapply(1:2, function(cl) {
          as.numeric(
            (onesim_meth_et[grep("ci.upp",row.names(onesim_meth_et)),cl]-
               onesim_meth_et[grep("ci.low",row.names(onesim_meth_et)),cl])/2)
        })
      )
    }
    onesim.res <- data.frame(meths_type,do.call(rbind,onesim.res))
    
    onesim.res.bc <- data.frame(
      expand.grid("estci"=paste0("bias_corrected",c("",".true")),
                  "et"=c("outcome","ipw"),
                  stringsAsFactors=FALSE),
      rbind(
        c(onesim$bias_corrected_out,onesim$bias_corrected_out-pop_effs),
        c(onesim$true_bias_corrected_out,onesim$true_bias_corrected_out-pop_effs),
        c(onesim$bias_corrected_ipw,onesim$bias_corrected_ipw-pop_effs),
        c(onesim$true_bias_corrected_ipw,onesim$true_bias_corrected_ipw-pop_effs)),
      matrix(NA,nrow=4,ncol=4)
    )
    colnames(onesim.res.bc) <- colnames(onesim.res)
    onesim.res <- rbind(onesim.res,onesim.res.bc)
    
    # trajectories
    colnames(onesim$trajectories) <- paste0("class",1:ncol(onesim$trajectories))
    onesim.res.traj <- lapply(1:nrow(onesim$trajectories), function(i) {
      onesim$trajectories[i,]
    })
    names(onesim.res.traj) <- row.names(onesim$trajectories)
    onesim.res.traj <- onesim.res.traj[c(3:5,1:2)]
    onesim.res.traj <- unlist(onesim.res.traj)
    
    return(list(onesim.res,onesim.res.traj))
  })
  res_summary.traj <- do.call(rbind,lapply(res_summary,"[[",2L))
  res_summary <- do.call(rbind,lapply(res_summary,"[[",1L))
  
  cat("---------------------------------------------------------------------\n")
  print(simsettings[ss,])
  print(length(res))
  print(pop_effs)
  
  meths_type <- unique(res_summary[,1:2])
  res_tab <- NULL
  for (mti in 1:nrow(meths_type)) {
    res_mti <- res_summary[
      res_summary$estci==meths_type[mti,"estci"] & 
        res_summary$et==meths_type[mti,"et"],]
    
    mm_av <- list(
      colMeans(res_mti[,-c(1:2)],na.rm=TRUE),
      "ese"=apply(res_mti[,paste0("est",1:2)],2,sd,na.rm=TRUE),
      "rmse"=sqrt(colMeans(res_mti[,paste0("bias",1:2)],na.rm=TRUE)^2+
                    apply(res_mti[,paste0("est",1:2)],2,var,na.rm=TRUE)))
    
    if (grepl("bias_corrected",meths_type[mti,"estci"])) {
      mm_cm <- rep(NA,2)  
    } else if (meths_type[mti,"estci"]=="true") {
      mm_cm <- rep(1.0,2)  
    } else {
      if (meths_type[mti,"estci"]=="est") {
        res_cm_mm <- "fixed"
      } else {
        res_cm_mm <- "boot"
      }
      mm_cm <- res_cm[res_cm$type==res_cm_mm,c("mean.est.class.1",
                                               "mean.est.class.2")]
    }
    
    mm_av <- c(mm_av,"correct_class"=list(mm_cm))
    res_tab[[mti]] <- unlist(mm_av)
  }
  res_tab <- do.call(rbind,res_tab)
  
  res_tab <- data.frame(
    "subgroup"=rep(ifelse(simsettings[ss,1]=="poLCA","LCA","GMM"),nrow(res_tab)),
    "Lambda"=rep(simsettings[ss,2],nrow(res_tab)),
    "direct"=rep(ifelse(simsettings[ss,3]==1,"Yes","No"),nrow(res_tab)),
    meths_type,
    res_tab[,as.character(
      sapply(c("bias","cover","ese.est","rmse.bias","correct_class"),
             function(s) paste0(s,1:2)))])
  
  res_tab[,paste0("bias",1:2)] <- res_tab[,paste0("bias",1:2)]*100
  res_tab[,paste0("ese.est",1:2)] <- res_tab[,paste0("ese.est",1:2)]*100
  res_tab[,paste0("rmse.bias",1:2)] <- res_tab[,paste0("rmse.bias",1:2)]*100
  
  print(xtable(res_tab),include.rownames=FALSE)
  
  res_summary.traj <- res_summary.traj[,c(1:7,9,8,10)]
  print(xtable(data.frame(res_tab[1,1:3],
                          t(colMeans(res_summary.traj)[-(5:6)]),
                          t(pop_effs))),
        include.rownames=FALSE)
  
}
