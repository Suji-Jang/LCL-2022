library(dplyr)
library(ggplot2)
library(reshape2)
library(ggpubr)

# Loading GM & GSD data
dirpath <- file.path("GitHub Data and Codes")
gm.gsd.df <- read.csv(file.path(dirpath,"GM_GSD_data.csv"))

# POD for each chemical
POD.med.chem <- gm.gsd.df[,paste("Chem",1:42,"GM",sep=".")]
POD.sens01.chem <- gm.gsd.df[,paste("Chem",1:42,"Sens.01",sep=".")]
col.300 <- unname(lapply(POD.med.chem,quantile,probs=0.05)<=300)

# Selecting chemicals with medians>=300
POD.med.chem <- POD.med.chem[,col.300]
POD.sens01.chem <- POD.sens01.chem[,col.300]

# write.csv(POD.med.chem,"EC10_Median_chemical.csv",row.names=FALSE)
# write.csv(POD.sens01.chem,"EC10_Sens01_chemical.csv",row.names=FALSE)

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
mix.info <- read.csv(file.path(dirpath,"Mix_info.csv"))
mix.info <- mix.info[,!(names(mix.info)=="Index")]
mix.info <- as.matrix(mix.info)
mix.info <- mix.info[col.300,]

# Inverse POD: 1/POD
POD.med.chem.inverse <- as.matrix(1/POD.med.chem)
POD.sens01.chem.inverse <- as.matrix(1/POD.sens01.chem)

# CA
mix.info <- as.matrix(mix.info)
POD.med.CA.Default.inverse <- POD.med.chem.inverse %*% mix.info
POD.sens01.CA.Default.inverse <- POD.sens01.chem.inverse %*% mix.info

# Inverse CA to get CA-Default POD
POD.med.CA.Default <- as.data.frame(1/POD.med.CA.Default.inverse)
POD.sens01.CA.Default <- as.data.frame(1/POD.sens01.CA.Default.inverse)

# write.csv(POD.med.CA.Default,file.path(dirpath,"POD_CADefault_Median.csv"),row.names=FALSE)
# write.csv(POD.sens01.CA.Default,file.path(dirpath,"POD_CADefault_Sens01.csv"),row.names=FALSE)

# Calculating AI = Measured POD / CA POD
AI.CA.Default.med <- POD.med.mix/POD.med.CA.Default
AI.CA.Default.sens01 <- POD.sens01.mix/POD.sens01.CA.Default
colnames(AI.CA.Default.med) <- c(paste0("Mix.",1:8))
colnames(AI.CA.Default.sens01) <- c(paste0("Mix.",1:8))

# write.csv(AI.CA.Default.med,file.path(dirpath,"AI_CADefault_Median.csv"),row.names=FALSE)
# write.csv(AI.CA.Default.sens01,file.path(dirpath,"AI_CADefault_Sens01.csv"),row.names=FALSE)