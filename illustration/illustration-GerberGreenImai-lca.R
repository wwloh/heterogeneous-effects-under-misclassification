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
ps.fit.mat <- data.frame(model.matrix(ps.fit)[,-1])

## Latent class analysis ======================================================
summary(ps.fit.mat)
# manifest indicators of latent class
Xnames_indicators <- grep(pattern="AGE",x=colnames(ps.fit.mat),
                          invert=TRUE,value=TRUE)
# no interaction terms
Xnames_indicators <- grep(pattern="PERSONS[.]",x=Xnames_indicators,
                          invert=TRUE,value=TRUE)
summary(ps.fit.mat[,Xnames_indicators])
ps.fit.mat[,-1] <- ps.fit.mat[,-1] + 1L # minimum value of 1

obs_data <- data.frame("i"=obs_data$i,ps.fit.mat,"Z"=obs_data$Z,"Y"=obs_data$Y)
Xnames <- colnames(ps.fit.mat)
rm(ps.fit.mat)

# fit inclusive latent class model
polca.f <- as.formula(
  paste0("cbind(",paste(Xnames_indicators,collapse=","),")~",
         paste(c("AGE","AGE2"),collapse="+")))
