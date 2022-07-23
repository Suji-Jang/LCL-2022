library(dplyr)
library(ggplot2)
library(reshape2)
library(ggpubr)
library(lognorm)

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

# Creating GM only dataframe
gm.only.df <- gm.gsd.df[,paste("Chem",1:42,"GM",sep=".")]
gm.only.df <- gm.only.df[,col.300]

# Creating Sigma only dataframe
sigma.only.df <- log(gm.gsd.df[,paste("Chem",1:42,"GSD",sep=".")])
sigma.only.df <- sigma.only.df[,col.300]
colnames(sigma.only.df) <- gsub("GSD","sigma",colnames(sigma.only.df))

# Sampling 4000 if sample size > 4000 for both GM and Sigma
ifelse(nrow(gm.only.df)>4000,gm.only.df<-gm.only.df[sample(nrow(gm.only.df),4000),],NA)
ifelse(nrow(sigma.only.df)>4000,sigma.only.df<-sigma.only.df[sample(nrow(sigma.only.df),4000),],NA)

# Assigning empty dataframes for Mu data
mu.mix1.df <- as.data.frame(matrix(NA,nrow=4000,ncol=32))
mu.mix2.df <- as.data.frame(matrix(NA,nrow=4000,ncol=32))
mu.mix3.df <- as.data.frame(matrix(NA,nrow=4000,ncol=32))
mu.mix4.df <- as.data.frame(matrix(NA,nrow=4000,ncol=32))
mu.mix5.df <- as.data.frame(matrix(NA,nrow=4000,ncol=32))
mu.mix6.df <- as.data.frame(matrix(NA,nrow=4000,ncol=32))
mu.mix7.df <- as.data.frame(matrix(NA,nrow=4000,ncol=32))
mu.mix8.df <- as.data.frame(matrix(NA,nrow=4000,ncol=32))

# Converting GM to mu
for(chemnum in 1:32){
  frac <- mix.info[chemnum,'Mix.1']
  temp.chem.df <- log(frac/(gm.only.df[,chemnum]))
  mu.mix1.df[,chemnum] <- temp.chem.df
}
colnames(mu.mix1.df) <- gsub("GM","mu",colnames(gm.only.df))

for(chemnum in 1:32){
  frac <- mix.info[chemnum,'Mix.2']
  temp.chem.df <- log(frac/(gm.only.df[,chemnum]))
  mu.mix2.df[,chemnum] <- temp.chem.df
}
colnames(mu.mix2.df) <- gsub("GM","mu",colnames(gm.only.df))

for(chemnum in 1:32){
  frac <- mix.info[chemnum,'Mix.3']
  temp.chem.df <- log(frac/(gm.only.df[,chemnum]))
  mu.mix3.df[,chemnum] <- temp.chem.df
}
colnames(mu.mix3.df) <- gsub("GM","mu",colnames(gm.only.df))

for(chemnum in 1:32){
  frac <- mix.info[chemnum,'Mix.4']
  temp.chem.df <- log(frac/(gm.only.df[,chemnum]))
  mu.mix4.df[,chemnum] <- temp.chem.df
}
colnames(mu.mix4.df) <- gsub("GM","mu",colnames(gm.only.df))

for(chemnum in 1:32){
  frac <- mix.info[chemnum,'Mix.5']
  temp.chem.df <- log(frac/(gm.only.df[,chemnum]))
  mu.mix5.df[,chemnum] <- temp.chem.df
}
colnames(mu.mix5.df) <- gsub("GM","mu",colnames(gm.only.df))

for(chemnum in 1:32){
  frac <- mix.info[chemnum,'Mix.6']
  temp.chem.df <- log(frac/(gm.only.df[,chemnum]))
  mu.mix6.df[,chemnum] <- temp.chem.df
}
colnames(mu.mix6.df) <- gsub("GM","mu",colnames(gm.only.df))

for(chemnum in 1:32){
  frac <- mix.info[chemnum,'Mix.7']
  temp.chem.df <- log(frac/(gm.only.df[,chemnum]))
  mu.mix7.df[,chemnum] <- temp.chem.df
}
colnames(mu.mix7.df) <- gsub("GM","mu",colnames(gm.only.df))

for(chemnum in 1:32){
  frac <- mix.info[chemnum,'Mix.8']
  temp.chem.df <- log(frac/(gm.only.df[,chemnum]))
  mu.mix8.df[,chemnum] <- temp.chem.df
}
colnames(mu.mix8.df) <- gsub("GM","mu",colnames(gm.only.df))

# Adding sigma values into mu dataframe
mu.sig.mix1.df <- cbind(mu.mix1.df,sigma.only.df)
mu.sig.mix2.df <- cbind(mu.mix2.df,sigma.only.df)
mu.sig.mix3.df <- cbind(mu.mix3.df,sigma.only.df)
mu.sig.mix4.df <- cbind(mu.mix4.df,sigma.only.df)
mu.sig.mix5.df <- cbind(mu.mix5.df,sigma.only.df)
mu.sig.mix6.df <- cbind(mu.mix6.df,sigma.only.df)
mu.sig.mix7.df <- cbind(mu.mix7.df,sigma.only.df)
mu.sig.mix8.df <- cbind(mu.mix8.df,sigma.only.df)

# Assigning empty columns for Mu and Sigma
mu.sig.mix1.df$CA.Mu <- NA
mu.sig.mix1.df$CA.Sigma <- NA
mu.sig.mix2.df$CA.Mu <- NA
mu.sig.mix2.df$CA.Sigma <- NA
mu.sig.mix3.df$CA.Mu <- NA
mu.sig.mix3.df$CA.Sigma <- NA
mu.sig.mix4.df$CA.Mu <- NA
mu.sig.mix4.df$CA.Sigma <- NA
mu.sig.mix5.df$CA.Mu <- NA
mu.sig.mix5.df$CA.Sigma <- NA
mu.sig.mix6.df$CA.Mu <- NA
mu.sig.mix6.df$CA.Sigma <- NA
mu.sig.mix7.df$CA.Mu <- NA
mu.sig.mix7.df$CA.Sigma <- NA
mu.sig.mix8.df$CA.Mu <- NA
mu.sig.mix8.df$CA.Sigma <- NA

# Sum of Lognormal distributions
for(i in 1:4000){
  mu.vec <- as.vector(t(mu.sig.mix1.df[i,gsub("GM","mu",colnames(gm.only.df))]))
  sig.vec <- as.vector(t(mu.sig.mix1.df[i,colnames(sigma.only.df)]))
  
  mu.sig.mix1.df[i,"CA.Mu"] <- unname(estimateSumLognormal(mu.vec,sig.vec)[1])
  mu.sig.mix1.df[i,"CA.Sigma"] <- unname(estimateSumLognormal(mu.vec,sig.vec)[2])
}

for(i in 1:4000){
  mu.vec <- as.vector(t(mu.sig.mix2.df[i,gsub("GM","mu",colnames(gm.only.df))]))
  sig.vec <- as.vector(t(mu.sig.mix2.df[i,colnames(sigma.only.df)]))
  
  mu.sig.mix2.df[i,"CA.Mu"] <- unname(estimateSumLognormal(mu.vec,sig.vec)[1])
  mu.sig.mix2.df[i,"CA.Sigma"] <- unname(estimateSumLognormal(mu.vec,sig.vec)[2])
}

for(i in 1:4000){
  mu.vec <- as.vector(t(mu.sig.mix3.df[i,gsub("GM","mu",colnames(gm.only.df))]))
  sig.vec <- as.vector(t(mu.sig.mix3.df[i,colnames(sigma.only.df)]))
  
  mu.sig.mix3.df[i,"CA.Mu"] <- unname(estimateSumLognormal(mu.vec,sig.vec)[1])
  mu.sig.mix3.df[i,"CA.Sigma"] <- unname(estimateSumLognormal(mu.vec,sig.vec)[2])
}

for(i in 1:4000){
  mu.vec <- as.vector(t(mu.sig.mix4.df[i,gsub("GM","mu",colnames(gm.only.df))]))
  sig.vec <- as.vector(t(mu.sig.mix4.df[i,colnames(sigma.only.df)]))
  
  mu.sig.mix4.df[i,"CA.Mu"] <- unname(estimateSumLognormal(mu.vec,sig.vec)[1])
  mu.sig.mix4.df[i,"CA.Sigma"] <- unname(estimateSumLognormal(mu.vec,sig.vec)[2])
}

for(i in 1:4000){
  mu.vec <- as.vector(t(mu.sig.mix5.df[i,gsub("GM","mu",colnames(gm.only.df))]))
  sig.vec <- as.vector(t(mu.sig.mix5.df[i,colnames(sigma.only.df)]))
  
  mu.sig.mix5.df[i,"CA.Mu"] <- unname(estimateSumLognormal(mu.vec,sig.vec)[1])
  mu.sig.mix5.df[i,"CA.Sigma"] <- unname(estimateSumLognormal(mu.vec,sig.vec)[2])
}

for(i in 1:4000){
  mu.vec <- as.vector(t(mu.sig.mix6.df[i,gsub("GM","mu",colnames(gm.only.df))]))
  sig.vec <- as.vector(t(mu.sig.mix6.df[i,colnames(sigma.only.df)]))
  
  mu.sig.mix6.df[i,"CA.Mu"] <- unname(estimateSumLognormal(mu.vec,sig.vec)[1])
  mu.sig.mix6.df[i,"CA.Sigma"] <- unname(estimateSumLognormal(mu.vec,sig.vec)[2])
}

for(i in 1:4000){
  mu.vec <- as.vector(t(mu.sig.mix7.df[i,gsub("GM","mu",colnames(gm.only.df))]))
  sig.vec <- as.vector(t(mu.sig.mix7.df[i,colnames(sigma.only.df)]))
  
  mu.sig.mix7.df[i,"CA.Mu"] <- unname(estimateSumLognormal(mu.vec,sig.vec)[1])
  mu.sig.mix7.df[i,"CA.Sigma"] <- unname(estimateSumLognormal(mu.vec,sig.vec)[2])
}

for(i in 1:4000){
  mu.vec <- as.vector(t(mu.sig.mix8.df[i,gsub("GM","mu",colnames(gm.only.df))]))
  sig.vec <- as.vector(t(mu.sig.mix8.df[i,colnames(sigma.only.df)]))
  
  mu.sig.mix8.df[i,"CA.Mu"] <- unname(estimateSumLognormal(mu.vec,sig.vec)[1])
  mu.sig.mix8.df[i,"CA.Sigma"] <- unname(estimateSumLognormal(mu.vec,sig.vec)[2])
}

# Converting CA-mu&sigma to CA-GM&GSD
mu.sig.mix1.df$CA.GM <- 1/exp(mu.sig.mix1.df$CA.Mu)
mu.sig.mix1.df$CA.GSD <- exp(mu.sig.mix1.df$CA.Sigma)
mu.sig.mix2.df$CA.GM <- 1/exp(mu.sig.mix2.df$CA.Mu)
mu.sig.mix2.df$CA.GSD <- exp(mu.sig.mix2.df$CA.Sigma)
mu.sig.mix3.df$CA.GM <- 1/exp(mu.sig.mix3.df$CA.Mu)
mu.sig.mix3.df$CA.GSD <- exp(mu.sig.mix3.df$CA.Sigma)
mu.sig.mix4.df$CA.GM <- 1/exp(mu.sig.mix4.df$CA.Mu)
mu.sig.mix4.df$CA.GSD <- exp(mu.sig.mix4.df$CA.Sigma)
mu.sig.mix5.df$CA.GM <- 1/exp(mu.sig.mix5.df$CA.Mu)
mu.sig.mix5.df$CA.GSD <- exp(mu.sig.mix5.df$CA.Sigma)
mu.sig.mix6.df$CA.GM <- 1/exp(mu.sig.mix6.df$CA.Mu)
mu.sig.mix6.df$CA.GSD <- exp(mu.sig.mix6.df$CA.Sigma)
mu.sig.mix7.df$CA.GM <- 1/exp(mu.sig.mix7.df$CA.Mu)
mu.sig.mix7.df$CA.GSD <- exp(mu.sig.mix7.df$CA.Sigma)
mu.sig.mix8.df$CA.GM <- 1/exp(mu.sig.mix8.df$CA.Mu)
mu.sig.mix8.df$CA.GSD <- exp(mu.sig.mix8.df$CA.Sigma)

# CA - Sensitive 1st percentile data
mu.sig.mix1.df$CA.Sens01 <- mu.sig.mix1.df$CA.GM*(mu.sig.mix1.df$CA.GSD)^(-2.326)
mu.sig.mix2.df$CA.Sens01 <- mu.sig.mix2.df$CA.GM*(mu.sig.mix2.df$CA.GSD)^(-2.326)
mu.sig.mix3.df$CA.Sens01 <- mu.sig.mix3.df$CA.GM*(mu.sig.mix3.df$CA.GSD)^(-2.326)
mu.sig.mix4.df$CA.Sens01 <- mu.sig.mix4.df$CA.GM*(mu.sig.mix4.df$CA.GSD)^(-2.326)
mu.sig.mix5.df$CA.Sens01 <- mu.sig.mix5.df$CA.GM*(mu.sig.mix5.df$CA.GSD)^(-2.326)
mu.sig.mix6.df$CA.Sens01 <- mu.sig.mix6.df$CA.GM*(mu.sig.mix6.df$CA.GSD)^(-2.326)
mu.sig.mix7.df$CA.Sens01 <- mu.sig.mix7.df$CA.GM*(mu.sig.mix7.df$CA.GSD)^(-2.326)
mu.sig.mix8.df$CA.Sens01 <- mu.sig.mix8.df$CA.GM*(mu.sig.mix8.df$CA.GSD)^(-2.326)

# Extracting CA for median and sensitive 1st percentile populations
POD.med.CA.LNSum <- as.data.frame(lapply(list(mu.sig.mix1.df,mu.sig.mix2.df,mu.sig.mix3.df,mu.sig.mix4.df,
                                         mu.sig.mix5.df,mu.sig.mix6.df,mu.sig.mix7.df,mu.sig.mix8.df), 
                                    function(x) x$CA.GM))
colnames(POD.med.CA.LNSum) <- paste("Mix",1:8,sep=".")

POD.sens01.CA.LNSum <- as.data.frame(lapply(list(mu.sig.mix1.df,mu.sig.mix2.df,mu.sig.mix3.df,mu.sig.mix4.df,
                                            mu.sig.mix5.df,mu.sig.mix6.df,mu.sig.mix7.df,mu.sig.mix8.df), 
                                       function(x) x$CA.Sens01))
colnames(POD.sens01.CA.LNSum) <- paste("Mix",1:8,sep=".")

# write.csv(POD.med.CA.LNSum,file.path(dirpath,"POD_CALNSum_Median.csv"),row.names=FALSE)
# write.csv(POD.sens01.CA.LNSum,file.path(dirpath,"POD_CALNSum_Sens01.csv"),row.names=FALSE)

# Deriving Additivity Index
AI.CA.LNSum.med <- POD.med.mix/POD.med.CA.LNSum
AI.CA.LNSum.sens01 <- POD.sens01.mix/POD.sens01.CA.LNSum
colnames(AI.CA.LNSum.med) <- c(paste0("Mix.",1:8))
colnames(AI.CA.LNSum.sens01) <- c(paste0("Mix.",1:8))

# write.csv(AI.CA.LNSum.med,file.path(dirpath,"AI_CALNSum_Median.csv"),row.names=FALSE)
# write.csv(AI.CA.LNSum.sens01,file.path(dirpath,"AI_CALNSum_Sens01.csv"),row.names=FALSE)
