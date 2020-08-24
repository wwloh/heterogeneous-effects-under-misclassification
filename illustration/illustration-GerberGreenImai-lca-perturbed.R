rm(list=ls())
library("geepack")
source("helper_functions.R")
load(file="GerberGreenImai-2class-perturbed_lca.Rdata")

# initialize for parallel cluster jobs
args <- 1
if (!grepl("apple",sessionInfo()[[1]]$platform)) {
  args <- commandArgs(trailingOnly=TRUE) # for CMD BATCH '--args 1'
}
(seed <- as.integer(args[1]))
rm(args)

ptm <- proc.time()[3]
# perturbed class membership probabilities ==================================
lambda <- resampled_list[[seed]]
lambda1 <- lambda[,1]
# ordered by decreasing probability of membership to class 1
deltas.l1o <- DeltasSequence(Data=obs_Data[order(lambda1,decreasing=TRUE),])
save(deltas.l1o, 
     file=paste0("GerberGreenImai-2class-perturbed_lca-",seed,".Rdata"))
proc.time()[3]-ptm
q()
