library(dplyr)
library(ggplot2)
library(reshape2)
library(ggpubr)
library(fitdistrplus)

# Setting parameters
dirpath <- file.path("GitHub Data and Codes")
lcl_dat <- read.csv(file.path(dirpath,"LCL_data_040722.csv"))
lcl_dat <- subset(lcl_dat,Chem_index!=0)
lcl_dat$Chem <- factor(lcl_dat$Chem)
lcl_dat$Cell.line <- factor(lcl_dat$Cell.line)
lcl_dat$Chem_index <- factor(lcl_dat$Chem_index)
chems <- levels(lcl_dat$Chem_index)
cells <- levels(lcl_dat$Cell.line)

# Loading GM & GSD data
dirpath <- file.path("GitHub Data and Codes")
gm.gsd.df <- read.csv(file.path(dirpath,"GM_GSD_data.csv"))

# POD for each chemical
POD.med.chem <- gm.gsd.df[,paste("Chem",1:42,"GM",sep=".")]
POD.sens01.chem <- gm.gsd.df[,paste("Chem",1:42,"Sens.01",sep=".")]

# Selecting chemicals with medians>=300
col.300 <- unname(lapply(POD.med.chem,quantile,probs=0.05)<=300)

POD.med.chem <- POD.med.chem[,col.300]
POD.sens01.chem <- POD.sens01.chem[,col.300]

# Converting Mixture GM to original concentration (uM)
POD.med.mix <- gm.gsd.df[,paste("Mix",1:8,"GM",sep=".")]
POD.sens01.mix <- gm.gsd.df[,paste("Mix",1:8,"Sens.01",sep=".")]

POD.med.mix$Mix.1.GM <- POD.med.mix$Mix.1.GM*48.3/100
POD.sens01.mix$Mix.1.Sens.01 <- POD.sens01.mix$Mix.1.Sens.01*48.3/100
POD.med.mix$Mix.2.GM <- POD.med.mix$Mix.2.GM*6236.3/100
POD.sens01.mix$Mix.2.Sens.01 <- POD.sens01.mix$Mix.2.Sens.01*6236.3/100
POD.med.mix$Mix.3.GM <- POD.med.mix$Mix.3.GM*2767.1/100
POD.sens01.mix$Mix.3.Sens.01 <- POD.sens01.mix$Mix.3.Sens.01*2767.1/100
POD.med.mix$Mix.4.GM <- POD.med.mix$Mix.4.GM*21348.4/100
POD.sens01.mix$Mix.4.Sens.01 <- POD.sens01.mix$Mix.4.Sens.01*21348.4/100
POD.med.mix$Mix.5.GM <- POD.med.mix$Mix.5.GM*79.4/100
POD.sens01.mix$Mix.5.Sens.01 <- POD.sens01.mix$Mix.5.Sens.01*79.4/100
POD.med.mix$Mix.6.GM <- POD.med.mix$Mix.6.GM*79.9/100
POD.sens01.mix$Mix.6.Sens.01 <- POD.sens01.mix$Mix.6.Sens.01*79.9/100
POD.med.mix$Mix.7.GM <- POD.med.mix$Mix.7.GM*83.8/100
POD.sens01.mix$Mix.7.Sens.01 <- POD.sens01.mix$Mix.7.Sens.01*83.8/100
POD.med.mix$Mix.8.GM <- POD.med.mix$Mix.8.GM*115.7/100
POD.sens01.mix$Mix.8.Sens.01 <- POD.sens01.mix$Mix.8.Sens.01*115.7/100

# Retrieving fraction of chemicals in each mixture
mix.info <- read.csv("Mix_info.csv",as.is=TRUE)
mix.info <- mix.info[,!(names(mix.info)=="Index")]
mix.info <- as.matrix(mix.info)
mix.info <- mix.info[col.300,]

# Loading PODs for each chemical or mixture
pod.chem.ind.df <- list()
for (chemnum in 1:50){
  fileprefix <- chemnum
  load(paste(fileprefix,"stanfit.Rdata",sep="_"))
  fitparms_df <- as.data.frame(rstan::extract(stan_fit))
  
  temp.chem.ind.df <- fitparms_df[,paste0("ec10.",1:146)]
  pod.chem.ind.df[[chemnum]] <- temp.chem.ind.df
}

# Extracting chemicals' PODs, sampling 4000 if samples > 4000, and selecting chemicals with median >= 300
pod.chem.ind.list <- list()
for (chemnum in 1:42){
  temp.list <- as.data.frame(pod.chem.ind.df[chemnum])
  ifelse(nrow(temp.list)>4000,temp.list<-temp.list[sample(nrow(temp.list),4000),],NA)
  pod.chem.ind.list[[chemnum]]<-temp.list
} 
pod.chem.ind.list <- pod.chem.ind.list[col.300]

# Extracting mixtures' PODs and sampling 4000 if samples > 4000
pod.mix.ind.list <- list()
for (chemnum in 43:50){
  temp.list <- as.data.frame(pod.chem.ind.df[chemnum])
  ifelse(nrow(temp.list)>4000,temp.list<-temp.list[sample(nrow(temp.list),4000),],NA)
  pod.mix.ind.list[[chemnum-42]]<-temp.list
} 

# Extracting PODs for each cell line and Multiplying inverse PODs (32 chems*4000) and fractions (8 mixes*32 chems)
pod.cell.ind.inverse <- list()
for (cellnum in 1:146){
  cellname <- paste("ec10",cellnum,sep=".")
  temp.cell.df <- as.data.frame(lapply(pod.chem.ind.list.300, function(x) x[,cellname]))
  temp.cell.df.inverse <- as.matrix(1/temp.cell.df)
  temp.pod.df <- temp.cell.df.inverse %*% mix.info
  pod.cell.ind.inverse[[cellnum]]<-temp.pod.df
}

# Censored normal distribution of log-transformed CA Mixture PODs
POD.mix1.CA.Indiv.inverse <- as.data.frame(lapply(pod.cell.ind.inverse, function(x) x[,"Mix.1"]))
colnames(POD.mix1.CA.Indiv.inverse) <- unique(lcl_dat$Cell.line)

POD.mix1.CA.Indiv <- 1/POD.mix1.CA.Indiv.inverse
cens.mix1.df <- as.data.frame(matrix(NA,nrow=4000,ncol=2))

for (i in 1:4000){
  temp.censor <- as.data.frame(log(t(POD.mix1.CA.Indiv[i,])))
  colnames(temp.censor) <- "value"
  temp.censor$right <- ifelse(temp.censor$value<log(100),temp.censor$right<-temp.censor$value,NA)
  temp.censor$left <- ifelse(temp.censor$value>=log(100),temp.censor$left<-log(100),temp.censor$left<-temp.censor$value)
  temp.censor <- temp.censor[,c("left","right")]
  fln <- fitdistcens(temp.censor,"norm")
  fln.gm <- exp(unname(fln$estimate)[1])
  fln.gsd <- exp(unname(fln$estimate)[2])
  cens.mix1.df[i,1] <- fln.gm
  cens.mix1.df[i,2] <- fln.gm * (fln.gsd)^-2.326  # 1st percentile
}
colnames(cens.mix1.df) <- c("Med","Perc.01")

POD.mix2.CA.Indiv.inverse <- as.data.frame(lapply(pod.cell.ind.inverse, function(x) x[,"Mix.2"]))
colnames(POD.mix2.CA.Indiv.inverse) <- unique(lcl_dat$Cell.line)

POD.mix2.CA.Indiv <- 1/POD.mix2.CA.Indiv.inverse
cens.mix2.df <- as.data.frame(matrix(NA,nrow=4000,ncol=2))

for (i in 1:4000){
  temp.censor <- as.data.frame(log(t(POD.mix2.CA.Indiv[i,])))
  colnames(temp.censor) <- "value"
  temp.censor$right <- ifelse(temp.censor$value<log(100),temp.censor$right<-temp.censor$value,NA)
  temp.censor$left <- ifelse(temp.censor$value>=log(100),temp.censor$left<-log(100),temp.censor$left<-temp.censor$value)
  temp.censor <- temp.censor[,c("left","right")]
  fln <- fitdistcens(temp.censor,"norm")
  fln.gm <- exp(unname(fln$estimate)[1])
  fln.gsd <- exp(unname(fln$estimate)[2])
  cens.mix2.df[i,1] <- fln.gm
  cens.mix2.df[i,2] <- fln.gm * (fln.gsd)^-2.326  # 1st percentile
}
colnames(cens.mix2.df) <- c("Med","Perc.01")

POD.mix3.CA.Indiv.inverse <- as.data.frame(lapply(pod.cell.ind.inverse, function(x) x[,"Mix.3"]))
colnames(POD.mix3.CA.Indiv.inverse) <- unique(lcl_dat$Cell.line)

POD.mix3.CA.Indiv <- 1/POD.mix3.CA.Indiv.inverse
cens.mix3.df <- as.data.frame(matrix(NA,nrow=4000,ncol=2))

# Dividing into 2 groups, because sample no. 3851 failed running
for (i in 1:3850){
  temp.censor <- as.data.frame(log(t(POD.mix3.CA.Indiv[i,])))
  colnames(temp.censor) <- "value"
  temp.censor$right <- ifelse(temp.censor$value<log(100),temp.censor$right<-temp.censor$value,NA)
  temp.censor$left <- ifelse(temp.censor$value>=log(100),temp.censor$left<-log(100),temp.censor$left<-temp.censor$value)
  temp.censor <- temp.censor[,c("left","right")]
  fln <- fitdistcens(temp.censor,"norm")
  fln.gm <- exp(unname(fln$estimate)[1])
  fln.gsd <- exp(unname(fln$estimate)[2])
  cens.mix3.df[i,1] <- fln.gm
  cens.mix3.df[i,2] <- fln.gm * (fln.gsd)^-2.326  # 1st percentile
}

for (i in 3852:4000){
  temp.censor <- as.data.frame(log(t(POD.mix3.CA.Indiv[i,])))
  colnames(temp.censor) <- "value"
  temp.censor$right <- ifelse(temp.censor$value<log(100),temp.censor$right<-temp.censor$value,NA)
  temp.censor$left <- ifelse(temp.censor$value>=log(100),temp.censor$left<-log(100),temp.censor$left<-temp.censor$value)
  temp.censor <- temp.censor[,c("left","right")]
  fln <- fitdistcens(temp.censor,"norm")
  fln.gm <- exp(unname(fln$estimate)[1])
  fln.gsd <- exp(unname(fln$estimate)[2])
  cens.mix3.df[i,1] <- fln.gm
  cens.mix3.df[i,2] <- fln.gm * (fln.gsd)^-2.326  # 1st percentile
}
colnames(cens.mix3.df) <- c("Med","Perc.01")

# Applying median of the other data to empty values at sample no. 3851
cens.mix3.df[3851,] <- apply(cens.mix3.df,2,median)
rownames(cens.mix3.df) <- 1:4000

POD.mix4.CA.Indiv.inverse <- as.data.frame(lapply(pod.cell.ind.inverse, function(x) x[,"Mix.4"]))
colnames(POD.mix4.CA.Indiv.inverse) <- unique(lcl_dat$Cell.line)

POD.mix4.CA.Indiv <- 1/POD.mix4.CA.Indiv.inverse
cens.mix4.df <- as.data.frame(matrix(NA,nrow=4000,ncol=2))

for (i in 1:4000){
  temp.censor <- as.data.frame(log(t(POD.mix4.CA.Indiv[i,])))
  colnames(temp.censor) <- "value"
  temp.censor$right <- ifelse(temp.censor$value<log(100),temp.censor$right<-temp.censor$value,NA)
  temp.censor$left <- ifelse(temp.censor$value>=log(100),temp.censor$left<-log(100),temp.censor$left<-temp.censor$value)
  temp.censor <- temp.censor[,c("left","right")]
  fln <- fitdistcens(temp.censor,"norm")
  fln.gm <- exp(unname(fln$estimate)[1])
  fln.gsd <- exp(unname(fln$estimate)[2])
  cens.mix4.df[i,1] <- fln.gm
  cens.mix4.df[i,2] <- fln.gm * (fln.gsd)^-2.326  # 1st percentile
}
colnames(cens.mix4.df) <- c("Med","Perc.01")

POD.mix5.CA.Indiv.inverse <- as.data.frame(lapply(pod.cell.ind.inverse, function(x) x[,"Mix.5"]))
colnames(POD.mix5.CA.Indiv.inverse) <- unique(lcl_dat$Cell.line)

POD.mix5.CA.Indiv <- 1/POD.mix5.CA.Indiv.inverse
cens.mix5.df <- as.data.frame(matrix(NA,nrow=4000,ncol=2))

for (i in 1:4000){
  temp.censor <- as.data.frame(log(t(POD.mix5.CA.Indiv[i,])))
  colnames(temp.censor) <- "value"
  temp.censor$right <- ifelse(temp.censor$value<log(100),temp.censor$right<-temp.censor$value,NA)
  temp.censor$left <- ifelse(temp.censor$value>=log(100),temp.censor$left<-log(100),temp.censor$left<-temp.censor$value)
  temp.censor <- temp.censor[,c("left","right")]
  fln <- fitdistcens(temp.censor,"norm")
  fln.gm <- exp(unname(fln$estimate)[1])
  fln.gsd <- exp(unname(fln$estimate)[2])
  cens.mix5.df[i,1] <- fln.gm
  cens.mix5.df[i,2] <- fln.gm * (fln.gsd)^-2.326  # 1st percentile
}
colnames(cens.mix5.df) <- c("Med","Perc.01")

POD.mix6.CA.Indiv.inverse <- as.data.frame(lapply(pod.cell.ind.inverse, function(x) x[,"Mix.6"]))
colnames(POD.mix6.CA.Indiv.inverse) <- unique(lcl_dat$Cell.line)

POD.mix6.CA.Indiv <- 1/POD.mix6.CA.Indiv.inverse
cens.mix6.df <- as.data.frame(matrix(NA,nrow=4000,ncol=2))

for (i in 1:4000){
  temp.censor <- as.data.frame(log(t(POD.mix6.CA.Indiv[i,])))
  colnames(temp.censor) <- "value"
  temp.censor$right <- ifelse(temp.censor$value<log(100),temp.censor$right<-temp.censor$value,NA)
  temp.censor$left <- ifelse(temp.censor$value>=log(100),temp.censor$left<-log(100),temp.censor$left<-temp.censor$value)
  temp.censor <- temp.censor[,c("left","right")]
  fln <- fitdistcens(temp.censor,"norm")
  fln.gm <- exp(unname(fln$estimate)[1])
  fln.gsd <- exp(unname(fln$estimate)[2])
  cens.mix6.df[i,1] <- fln.gm
  cens.mix6.df[i,2] <- fln.gm * (fln.gsd)^-2.326  # 1st percentile
}
colnames(cens.mix6.df) <- c("Med","Perc.01")

POD.mix7.CA.Indiv.inverse <- as.data.frame(lapply(pod.cell.ind.inverse, function(x) x[,"Mix.7"]))
colnames(POD.mix7.CA.Indiv.inverse) <- unique(lcl_dat$Cell.line)

POD.mix7.CA.Indiv <- 1/POD.mix7.CA.Indiv.inverse
cens.mix7.df <- as.data.frame(matrix(NA,nrow=4000,ncol=2))

for (i in 1:4000){
  temp.censor <- as.data.frame(log(t(POD.mix7.CA.Indiv[i,])))
  colnames(temp.censor) <- "value"
  temp.censor$right <- ifelse(temp.censor$value<log(100),temp.censor$right<-temp.censor$value,NA)
  temp.censor$left <- ifelse(temp.censor$value>=log(100),temp.censor$left<-log(100),temp.censor$left<-temp.censor$value)
  temp.censor <- temp.censor[,c("left","right")]
  fln <- fitdistcens(temp.censor,"norm")
  fln.gm <- exp(unname(fln$estimate)[1])
  fln.gsd <- exp(unname(fln$estimate)[2])
  cens.mix7.df[i,1] <- fln.gm
  cens.mix7.df[i,2] <- fln.gm * (fln.gsd)^-2.326  # 1st percentile
}
colnames(cens.mix7.df) <- c("Med","Perc.01")

POD.mix8.CA.Indiv.inverse <- as.data.frame(lapply(pod.cell.ind.inverse, function(x) x[,"Mix.8"]))
colnames(POD.mix8.CA.Indiv.inverse) <- unique(lcl_dat$Cell.line)

POD.mix8.CA.Indiv <- 1/POD.mix8.CA.Indiv.inverse
cens.mix8.df <- as.data.frame(matrix(NA,nrow=4000,ncol=2))

for (i in 1:4000){
  temp.censor <- as.data.frame(log(t(POD.mix8.CA.Indiv[i,])))
  colnames(temp.censor) <- "value"
  temp.censor$right <- ifelse(temp.censor$value<log(100),temp.censor$right<-temp.censor$value,NA)
  temp.censor$left <- ifelse(temp.censor$value>=log(100),temp.censor$left<-log(100),temp.censor$left<-temp.censor$value)
  temp.censor <- temp.censor[,c("left","right")]
  fln <- fitdistcens(temp.censor,"norm")
  fln.gm <- exp(unname(fln$estimate)[1])
  fln.gsd <- exp(unname(fln$estimate)[2])
  cens.mix8.df[i,1] <- fln.gm
  cens.mix8.df[i,2] <- fln.gm * (fln.gsd)^-2.326  # 1st percentile
}
colnames(cens.mix8.df) <- c("Med","Perc.01")

# Extracting CA for median and sensitive 1st percentile populations
POD.med.CA.Indiv <- as.data.frame(lapply(list(cens.mix1.df,cens.mix2.df,cens.mix3.df,cens.mix4.df,
                                         cens.mix5.df,cens.mix6.df,cens.mix7.df,cens.mix8.df), function(x) x$Med))
colnames(POD.med.CA.Indiv)<- c(paste0("Mix.",1:8))
POD.sens01.CA.Indiv <- as.data.frame(lapply(list(cens.mix1.df,cens.mix2.df,cens.mix3.df,cens.mix4.df,
                                            cens.mix5.df,cens.mix6.df,cens.mix7.df,cens.mix8.df), function(x) x$Perc.01))
colnames(POD.sens01.CA.Indiv)<- c(paste0("Mix.",1:8))

# write.csv(POD.med.CA.Indiv,file.path(dirpath,"POD_CAIndiv_Median.csv"),row.names=FALSE)
# write.csv(POD.sens01.CA.Indiv,file.path(dirpath,"POD_CAIndiv_Sens01.csv"),row.names=FALSE)

# Deriving Additivity Index
AI.CA.Indiv.med <- POD.med.mix/POD.med.CA.Indiv
AI.CA.Indiv.sens01 <- POD.sens01.mix/POD.sens01.CA.Indiv
colnames(AI.CA.Indiv.med) <- c(paste0("Mix.",1:8))
colnames(AI.CA.Indiv.sens01) <- c(paste0("Mix.",1:8))

# write.csv(AI.CA.Indiv.med,file.path(dirpath,"AI_CAIndiv_Median.csv"),row.names=FALSE)
# write.csv(AI.CA.Indiv.sens01,file.path(dirpath,"AI_CAIndiv_Sens01.csv"),row.names=FALSE)