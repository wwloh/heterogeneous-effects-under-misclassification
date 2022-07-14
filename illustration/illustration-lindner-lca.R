library("poLCA")
# install.packages("twang")
library("twang")
library("data.table")

# data preparation
# observational study
## load data publicly available as part of the twang package
help(lindner)
data(lindner)
dim(lindner) # number of obs
summary(lindner)
sapply(lindner, class)

lindner.dt <- data.table(sapply(lindner, as.numeric))
setkey(lindner.dt)

# covariates
l <- lindner.dt[,!(colnames(lindner.dt) %in% c(
  "lifepres","cardbill","abcix","sixMonthSurvive")),with=FALSE]
summary(l)
l <- as.data.frame(l)
obs_data <- data.frame("i"=1:nrow(l),
                       l,
                       "Z"=as.integer(lindner.dt$abcix),
                       "Y"=as.integer(lindner.dt$sixMonthSurvive))
table(obs_data[,c("Y","Z")])
(n <- nrow(obs_data))
summary(obs_data)

rm(l,lindner,lindner.dt)

## Latent class analysis ======================================================
## explanatory variables of latent class
Xnames_nonindicators <- names(obs_data)[
  unlist(c("i"=1,
           "Z"=which(names(obs_data)=="Z"),
           "Y"=which(names(obs_data)=="Y"),
           sapply(c("height","female"), 
                  grep, x=names(obs_data), fixed=TRUE)))]
Xnames_nonindicators
(Xnames_indicators <- names(obs_data)[!(names(obs_data) %in% Xnames_nonindicators)])

# manifest indicators of latent class
## discretize continuous covariates
Xcont_idx <- names(which(apply(obs_data[,Xnames_indicators], 2, function(x)
  length(unique(x)) > 10)))
length(Xcont_idx)
summary(obs_data[,Xcont_idx])
for (x in Xnames_indicators) {
  if (x %in% Xcont_idx) {
    x.breaks <- unique(quantile(obs_data[,x],probs=(0:5)/5,na.rm=TRUE))
    x.fact <- cut(obs_data[,x],breaks=x.breaks,
                  include.lowest=TRUE,ordered_result=TRUE)
    x.fact <- as.integer(x.fact)
    obs_data[, x] <- x.fact
    rm(x.breaks,x.fact)
  } else {
    obs_data[, x] <- obs_data[, x] + 1L # minimum value of 1
  }
}
summary(obs_data[,Xcont_idx])
rm(x)

# all indicators have at least two unique values
table(apply(obs_data[,Xnames_indicators],2,function(x) length(unique(x))))

Xnames <- names(obs_data)[!(names(obs_data) %in% c("i","Z","Y"))]

# fit basic latent class model
polca.f <- as.formula(paste0("cbind(",paste(Xnames_indicators,collapse=","),")~1"))
polca.f

dataset_name <- "lindner"
filename <- paste0("illustration-",dataset_name,"-lca.Rdata")

# fit each model with different number of classes in turn
mixLCA_k <- list()
for (k in 2:10) {
  # Estimate multiple times using different initial parameter values
  mixLCA_k[[k]] <- poLCA(polca.f,data=obs_data,nclass=k,nrep=500,verbose=TRUE)
  save.image(file=filename)
}
