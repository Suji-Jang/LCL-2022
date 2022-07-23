library(ggplot2)
library(reshape2)
library(viridis)
library(dplyr)
library(ggpubr)

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

# Converting mixture PODs to original concentration (uM)
pod.mix.ind.list[[1]] <- pod.mix.ind.list[[1]]*48.3/100
pod.mix.ind.list[[2]] <- pod.mix.ind.list[[2]]*6236.3/100
pod.mix.ind.list[[3]] <- pod.mix.ind.list[[3]]*2767.1/100
pod.mix.ind.list[[4]] <- pod.mix.ind.list[[4]]*21348.4/100
pod.mix.ind.list[[5]] <- pod.mix.ind.list[[5]]*79.4/100
pod.mix.ind.list[[6]] <- pod.mix.ind.list[[6]]*79.9/100
pod.mix.ind.list[[7]] <- pod.mix.ind.list[[7]]*83.8/100
pod.mix.ind.list[[8]] <- pod.mix.ind.list[[8]]*115.7/100

# Deriving median of individual
ind.chem.med.df <- as.data.frame(matrix(NA,nrow=146,ncol=32))
for (chemnum in 1:32){
  temp.ind.med <- apply(as.data.frame(pod.chem.ind.list[chemnum]),2,median)
  ind.chem.med.df[,chemnum] <- temp.ind.med
}
ind.chem.med.df <- t(ind.chem.med.df)

ind.mix.med.df <- as.data.frame(matrix(NA,nrow=146,ncol=8))
for (mixnum in 1:8){
  temp.ind.med <- apply(as.data.frame(pod.mix.ind.list[mixnum]),2,median)
  ind.mix.med.df[,mixnum] <- temp.ind.med
}

# CA for individual
ca.ind.chem.med.df <- as.data.frame(matrix(NA,nrow=146,ncol=8))
for (cellnum in 1:146){
  temp.ind.ca <- t(mix.info[,"Mix.1"]) %*% as.matrix(1/(ind.chem.med.df[,cellnum]))
  ca.ind.chem.med.df[cellnum,1] <- temp.ind.ca
  temp.ind.ca <- t(mix.info[,"Mix.2"]) %*% as.matrix(1/(ind.chem.med.df[,cellnum]))
  ca.ind.chem.med.df[cellnum,2] <- temp.ind.ca
  temp.ind.ca <- t(mix.info[,"Mix.3"]) %*% as.matrix(1/(ind.chem.med.df[,cellnum]))
  ca.ind.chem.med.df[cellnum,3] <- temp.ind.ca
  temp.ind.ca <- t(mix.info[,"Mix.4"]) %*% as.matrix(1/(ind.chem.med.df[,cellnum]))
  ca.ind.chem.med.df[cellnum,4] <- temp.ind.ca
  temp.ind.ca <- t(mix.info[,"Mix.5"]) %*% as.matrix(1/(ind.chem.med.df[,cellnum]))
  ca.ind.chem.med.df[cellnum,5] <- temp.ind.ca
  temp.ind.ca <- t(mix.info[,"Mix.6"]) %*% as.matrix(1/(ind.chem.med.df[,cellnum]))
  ca.ind.chem.med.df[cellnum,6] <- temp.ind.ca
  temp.ind.ca <- t(mix.info[,"Mix.7"]) %*% as.matrix(1/(ind.chem.med.df[,cellnum]))
  ca.ind.chem.med.df[cellnum,7] <- temp.ind.ca
  temp.ind.ca <- t(mix.info[,"Mix.8"]) %*% as.matrix(1/(ind.chem.med.df[,cellnum]))
  ca.ind.chem.med.df[cellnum,8] <- temp.ind.ca
}
colnames(ca.ind.chem.med.df) <- paste0("Mix.",1:8)

# Converting inverse CA POD to CA POD
ca.ind.chem.med.df <- 1/ca.ind.chem.med.df

# Merging CA EC10 and Measured EC10
ca.mea.plot.df <- as.data.frame(matrix(NA,nrow=0,ncol=3))

for (mixnum in 1:8){
  mix.name <- paste0("Mix.",mixnum)
  temp.df <- data.frame(Mix=rep(mix.name,146),CA.EC10.med=ca.ind.chem.med.df[,mixnum],Mea.EC10.med=ind.mix.med.df[,mixnum])
  ca.mea.plot.df <- rbind(ca.mea.plot.df,temp.df)
}
# Changing mixture names
ca.mea.plot.df <- ca.mea.plot.df %>%
  mutate(Mix = case_when(
    Mix=="Mix.1" ~ "AC50-L",
    Mix=="Mix.2" ~ "AC50-H",
    Mix=="Mix.3" ~ "POD-L",
    Mix=="Mix.4" ~ "POD-H",
    Mix=="Mix.5" ~ "Expo-L",
    Mix=="Mix.6" ~ "Expo-H",
    Mix=="Mix.7" ~ "RfD-L",
    Mix=="Mix.8" ~ "RfD-H",
  ))

ca.mea.plot.df$Mix <- factor(ca.mea.plot.df$Mix, levels=c("AC50-L","AC50-H","POD-L","POD-H","Expo-L","Expo-H","RfD-L","RfD-H"))

# Figure 2A - Comparing PODs for measured and CA with scatterplot
options(scipen=999)
min(ca.mea.plot.df[,c(2,3)])
max(ca.mea.plot.df[,c(2,3)])
p1 <- ggplot(ca.mea.plot.df,aes(x=Mea.EC10.med,y=CA.EC10.med,color=Mix)) + geom_point() + 
  scale_x_log10(labels = function(x) ifelse(x == 0, "0", x),limits=c(0.2,2400), breaks=10^seq(0,3)) + 
  scale_y_log10(labels = function(x) ifelse(x == 0, "0", x),limits=c(0.2,2400), breaks=10^seq(0,3)) +
  xlab("Measured Individual Mixture EC10 (uM)") + ylab("CA Individual Mixture EC10 (uM)") +
  geom_abline(linetype="dashed") + 
  geom_abline(linetype="dotted",colour="gray",size=0.8,aes(slope=log10(10),intercept=-1)) + 
  geom_abline(linetype="dotted",colour="gray",size=0.8,aes(slope=log10(10),intercept=1)) +
  scale_color_viridis(option="turbo",discrete=TRUE) + theme_classic() + 
  theme(legend.title=element_blank(),legend.position="none")

# Deriving Additivity Index (AI)
ai.ind.chem.med.df <- ind.mix.med.df/ca.ind.chem.med.df
colnames(ai.ind.chem.med.df) <- paste0("Mix.",1:8)

ai.ind.med.melted <- melt(ai.ind.chem.med.df)

ai.ind.med.melted <- ai.ind.med.melted %>%
  mutate(variable = case_when(
    variable=="Mix.1" ~ "AC50-L",
    variable=="Mix.2" ~ "AC50-H",
    variable=="Mix.3" ~ "POD-L",
    variable=="Mix.4" ~ "POD-H",
    variable=="Mix.5" ~ "Expo-L",
    variable=="Mix.6" ~ "Expo-H",
    variable=="Mix.7" ~ "RfD-L",
    variable=="Mix.8" ~ "RfD-H",
  ))

ai.ind.med.melted$variable <- factor(ai.ind.med.melted$variable, 
                                     levels=c("RfD-H","RfD-L","Expo-H","Expo-L","POD-H","POD-L","AC50-H","AC50-L"))

# Figure 2B - Plotting boxplot for AI of each mixture
min(ai.ind.med.melted$value)
max(ai.ind.med.melted$value)
p2 <- ggplot(ai.ind.med.melted,aes(x=value,y=variable,fill=variable)) + geom_boxplot(outlier.shape = NA) +
  scale_x_log10(limits=c(5e-3,5e2),breaks=10^seq(-4,4),labels = function(x) ifelse(x == 0, "0", x)) + 
  theme_classic() + geom_vline(xintercept=1,linetype="dashed") +
  geom_vline(xintercept=0.1,linetype="dotted",colour="gray") +
  geom_vline(xintercept=10,linetype="dotted",colour="gray") +
  xlab("LAI = (Measured EC10) / (CA EC10)") +
  theme(axis.title.y=element_blank(),legend.title=element_blank()) +
  scale_fill_viridis(option="turbo",discrete=TRUE,direction=-1,
                     breaks=c("AC50-L","AC50-H","POD-L","POD-H","Expo-L","Expo-H","RfD-L","RfD-H"))

p1.2 <- ggarrange(p1,p2,widths=c(1,1.2))
p1.2

# Figure 3C - Heatmap for mixtures * individual cell lines
library(gplots)
library(RColorBrewer)
library(colorspace)

coolwarm_hcl <- colorspace::diverging_hcl(100, h = c(250, 10), c = 100, l = c(37, 88), power = c(0.7, 1.7))
warmcool_hcl <- rev(coolwarm_hcl)
ai.ind.chem.med.df.trans <- as.data.frame(t(ai.ind.chem.med.df))
heatmap.df <- as.matrix(log10(ai.ind.chem.med.df.trans))
colnames(heatmap.df) <- cells
rownames(heatmap.df) <- c("AC50-L","AC50-H","POD-L","POD-H","Expo-L","Expo-H","RfD-L","RfD-H")
pdf("Figure 2C - Heatmap.pdf",height=5,width=20)
heatmap.2(heatmap.df,trace="n",cexRow=1,scale="none",lhei=c(2,8.5),lwid=c(0.5,8.5),col= warmcool_hcl,
          key.title=NA,key.xlab=NA,key.ylab=NA)(100)
dev.off()

# ANOVA and eta-squared
library(lsr)

melted.heatmap.df <- melt(as.data.frame(heatmap.df))
melted.heatmap.df$Mix <- rep(c("AC50-L","AC50-H","POD-L","POD-H","Expo-L","Expo-H","RfD-L","RfD-H"),146)

res <- lm(value ~ variable+Mix,data=melted.heatmap.df)
aovres <- aov(res)
print(summary(res))
print(summary(aovres)) 
print(etaSquared(aovres))