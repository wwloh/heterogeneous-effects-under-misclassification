# data preparation
# observational study
Hmisc::getHdata(rhc,where="https://hbiostat.org/data/repo/")
dim(rhc) # number of obs
summary(rhc$swang1) # number of treated vs. untreated
rhc <- rhc[,names(rhc)[order(names(rhc))]] # order column names alphabetically
for (cc in 1:ncol(rhc)) {
  cat(names(rhc)[cc], "----------------------------------------------------\n")
  print(summary(rhc[,cc]))
}
rhc$cat2[is.na(rhc$cat2)] <- "Missing" # remove NAs from cat2
## relabelled based on Hirano and Imbens (2002) Table 1
rhc.dt <- data.frame(
  "treat"=rhc$swang1=="RHC",
  "Y"=rhc$death=="Yes",
  "age"=rhc$age,
  "sex"=rhc$sex=="Female",
  # dummy variables for race
  "raceblack"=rhc$race=="black",
  "raceother"=rhc$race=="other",
  "edu"=rhc$edu,
  # dummy variables for income
  "income1"=rhc$income=="$11-$25k",
  "income2"=rhc$income=="$25-$50k",
  "income3"=rhc$income=="> $50k",
  # dummy variables for Medical insurance
  "ins_care"=rhc$ninsclas=="Medicare",
  "ins_pcare"=rhc$ninsclas=="Private & Medicare",
  "ins_caid"=rhc$ninsclas=="Medicaid",
  "ins_no"=rhc$ninsclas=="No insurance",
  "ins_carecaid"=rhc$ninsclas=="Medicare & Medicaid",
  # dummy variables for Primary disease category
  "cat1_copd"=rhc$cat1=="COPD",
  "cat1_mosfsep"=rhc$cat1=="MOSF w/Sepsis",
  "cat1_mosfmal"=rhc$cat1=="MOSF w/Malignancy",
  "cat1_chf"=rhc$cat1=="CHF",
  "cat1_coma"=rhc$cat1=="Coma",
  "cat1_cirr"=rhc$cat1=="Cirrhosis",
  "cat1_lung"=rhc$cat1=="Lung Cancer",
  "cat1_colon"=rhc$cat1=="Colon Cancer",
  # dummy variables for Secondary disease category
  "cat2_mosfsep"=rhc$cat2=="MOSF w/Sepsis",
  "cat2_coma"=rhc$cat2=="Coma",
  "cat2_mosfmal"=rhc$cat2=="MOSF w/Malignancy",
  "cat2_lung"=rhc$cat2=="Lung Cancer",
  "cat2_cirr"=rhc$cat2=="Cirrhosis",
  "cat2_colon"=rhc$cat2=="Colon Cancer", 
  # Categories of admission diagnosis
  "resp"=rhc$resp=="Yes",
  "card"=rhc$card=="Yes",
  "neuro"=rhc$neuro=="Yes",
  "gastr"=rhc$gastr=="Yes",
  "renal"=rhc$renal=="Yes",
  "meta"=rhc$meta=="Yes",
  "hema"=rhc$hema=="Yes",
  "seps"=rhc$seps=="Yes",
  "trauma"=rhc$trauma=="Yes",
  "ortho"=rhc$ortho=="Yes",
  rhc$das2d3pc,
  "dnr1"=rhc$dnr1=="Yes",
  "ca_yes"=rhc$ca=="Yes",
  "ca_meta"=rhc$ca=="Metastatic",
  rhc$surv2md1,
  rhc$aps1,
  rhc$scoma1,
  rhc$wtkilo1,
  rhc$temp1,
  rhc$meanbp1,
  rhc$resp1,
  rhc$hrt1,
  rhc$pafi1,
  rhc$paco21,
  rhc$ph1,
  rhc$wblc1,
  rhc$hema1,
  rhc$sod1,
  rhc$pot1,
  rhc$crea1,
  rhc$bili1,
  rhc$alb1,
  rhc$urin1, # not included in Hirano and Imbens
  # Categories of comorbidities illness:
  rhc$cardiohx,
  rhc$chfhx,
  rhc$dementhx,
  rhc$psychhx,
  rhc$chrpulhx,
  rhc$renalhx,
  rhc$liverhx,
  rhc$gibledhx,
  rhc$malighx,
  rhc$immunhx,
  rhc$transhx,
  rhc$amihx,
  "wt0"=rhc$wtkilo1==0
)
names(rhc.dt) <- sapply(names(rhc.dt), function(x) {
  ifelse(grepl("rhc.",x),strsplit(x,"rhc.")[[1]][2], x)
})
names(rhc.dt)

# remove variable with missing data
(l.rm <- which(colSums(is.na(rhc.dt))>0))
rhc.dt <- rhc.dt[,-l.rm]

# covariates
l <- rhc.dt[,3:ncol(rhc.dt)]
Lnames <- data.frame(cbind("L.idx"=1:ncol(l),"L.names"=names(l)))

# convert any logical covariates to integer
l <- apply(l, 2, function(x) as.numeric(x))
summary(l)
l <- as.data.frame(l)
colnames(l) <- NULL
obs_data <- data.frame("i"=1:nrow(l),
                       "L"=l,
                       "Z"=as.integer(rhc.dt$treat),
                       "Y"=as.integer(rhc.dt$Y))
summary(obs_data)

## Latent class analysis ======================================================
# manifest indicators of latent class
Xnames_indicators <- grep("L",names(obs_data),value=TRUE)[-(1:13)]
Xcont_idx <- names(which(apply(obs_data[,Xnames_indicators], 2, function(x)
  length(unique(x)) > 10)))
length(Xcont_idx)
summary(obs_data[,Xcont_idx])
for (x in Xnames_indicators) {
  if (x %in% Xcont_idx) {
    x.breaks <- unique(quantile(obs_data[,x],probs=(0:5)/5,na.rm=TRUE))
    x.fact <- cut(obs_data[,x],breaks=x.breaks,labels=FALSE,
                  include.lowest=TRUE,right=FALSE,ordered_result=TRUE)
    obs_data[, x] <- x.fact
  } else {
    obs_data[, x] <- obs_data[, x] + 1L # minimum value of 1
  }
  # # new category for missing data
  # obs_data[is.na(obs_data[, x]),x] <- max(obs_data[, x],na.rm=TRUE) + 1L
}
summary(obs_data[,Xcont_idx])

# any indicators have fewer than 2 unique values
any(apply(obs_data[,Xnames_indicators],2,max,na.rm=TRUE)<2)

table(obs_data$Z)
table(obs_data$Y)
(n <- nrow(obs_data))

Xnames <- grep(pattern="L.",x=colnames(obs_data),value=TRUE)
rm(l,rhc.dt)

# fit inclusive latent class model
polca.f <- as.formula(
  paste0("cbind(",paste(Xnames_indicators,collapse=","),")~",
         paste(paste0("L.",1:13),collapse="+")))
