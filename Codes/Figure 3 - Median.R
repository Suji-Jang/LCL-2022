library(reshape2)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(viridis)

# Loading POD and AI data
dirpath <- file.path("GitHub Data and Codes")
POD.med.CA.Indiv <- read.csv(file.path(dirpath,"POD_CAIndiv_Median.csv"))
POD.med.CA.LNSum <- read.csv(file.path(dirpath,"POD_CALNSum_Median.csv"))
POD.med.CA.Default <- read.csv(file.path(dirpath,"POD_CADefault_Median.csv"))

AI.med.CA.Indiv <- read.csv(file.path(dirpath,"AI_CAIndiv_Median.csv"))
AI.med.CA.LNSum <- read.csv(file.path(dirpath,"AI_CALNSum_Median.csv"))
AI.med.CA.Default <- read.csv(file.path(dirpath,"AI_CADefault_Median.csv"))

# Loading GM & GSD data
gm.gsd.df <- read.csv(file.path(dirpath,"GM_GSD_data.csv"))
POD.med.mix <- gm.gsd.df[,paste("Mix",1:8,"GM",sep=".")]

# Converting Mixture GM to original concentration (uM)
POD.med.mix$Mix.1.GM <- POD.med.mix$Mix.1.GM*48.3/100
POD.med.mix$Mix.2.GM <- POD.med.mix$Mix.2.GM*6236.3/100
POD.med.mix$Mix.3.GM <- POD.med.mix$Mix.3.GM*2767.1/100
POD.med.mix$Mix.4.GM <- POD.med.mix$Mix.4.GM*21348.4/100
POD.med.mix$Mix.5.GM <- POD.med.mix$Mix.5.GM*79.4/100
POD.med.mix$Mix.6.GM <- POD.med.mix$Mix.6.GM*79.9/100
POD.med.mix$Mix.7.GM <- POD.med.mix$Mix.7.GM*83.8/100
POD.med.mix$Mix.8.GM <- POD.med.mix$Mix.8.GM*115.7/100

colnames(POD.med.mix) <- paste0("Mix.",1:8)
measured.med.melted <- melt(POD.med.mix)

# Merging POD dataframes for measured and CA
CA.Indiv.med.melted <- melt(POD.med.CA.Indiv)
CA.Indiv.med.plot.df <- cbind(measured.med.melted,CA.Indiv.med.melted)
CA.Indiv.med.plot.df <- CA.Indiv.med.plot.df[,-3]
CA.Indiv.med.plot.df$CA <- "CA-Indiv"
colnames(CA.Indiv.med.plot.df) <- c("Mix","Measured.EC10","CA.EC10","CA")

CA.LNSum.med.melted <- melt(POD.med.CA.LNSum)
CA.LNSum.med.plot.df <- cbind(measured.med.melted,CA.LNSum.med.melted)
CA.LNSum.med.plot.df <- CA.LNSum.med.plot.df[,-3]
CA.LNSum.med.plot.df$CA <- "CA-LNSum"
colnames(CA.LNSum.med.plot.df) <- c("Mix","Measured.EC10","CA.EC10","CA")

CA.Default.med.melted <- melt(POD.med.CA.Default)
CA.Default.med.plot.df <- cbind(measured.med.melted,CA.Default.med.melted)
CA.Default.med.plot.df <- CA.Default.med.plot.df[,-3]
CA.Default.med.plot.df$CA <- "CA-Default"
colnames(CA.Default.med.plot.df) <- c("Mix","Measured.EC10","CA.EC10","CA")

scat.plot.df <- rbind(CA.Indiv.med.plot.df,CA.LNSum.med.plot.df,CA.Default.med.plot.df)

# Changing mixture names
scat.plot.df <- scat.plot.df %>%
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

scat.plot.df$Mix <- factor(scat.plot.df$Mix, levels=c("AC50-L","AC50-H","POD-L","POD-H","Expo-L","Expo-H","RfD-L","RfD-H"))
scat.plot.df$CA <- factor(scat.plot.df$CA, levels=c("CA-Indiv","CA-LNSum","CA-Default"))

# Figure 3A - Comparing PODs for measured and CA with scatterplot
min(scat.plot.df[,c(2,3)],na.rm=TRUE)
max(scat.plot.df[,c(2,3)],na.rm=TRUE)

p1 <- ggplot(scat.plot.df,aes(x=Measured.EC10,y=CA.EC10,color=Mix)) + geom_point() +
  scale_x_log10(labels = function(x) ifelse(x == 0, "0", x),limits=c(3,3500), breaks=10^seq(1,3)) + 
  scale_y_log10(labels = function(x) ifelse(x == 0, "0", x),limits=c(3,3500), breaks=10^seq(1,3)) +
  xlab("Measured Population Median Mixture EC10 (uM)") + ylab("CA Population Median Mixture EC10 (uM)") +
  theme_classic() + geom_abline(linetype="dashed") + facet_wrap(CA~.,nrow=3)+
  geom_abline(linetype="dotted",colour="gray",size=0.8,aes(slope=log10(10),intercept=-1)) + 
  geom_abline(linetype="dotted",colour="gray",size=0.8,aes(slope=log10(10),intercept=1)) +
  scale_color_viridis(option="turbo",discrete=TRUE) + theme(legend.title=element_blank(),legend.position="none")

# Deriving Additivity Index (AI)
AI.med.CA.Indiv.melted <- melt(AI.med.CA.Indiv)
AI.med.CA.Indiv.melted$group <- "CA-Indiv"
AI.med.CA.LNSum.melted <- melt(AI.med.CA.LNSum)
AI.med.CA.LNSum.melted$group <- "CA-LNSum"
AI.med.CA.Default.melted <- melt(AI.med.CA.Default)
AI.med.CA.Default.melted$group <- "CA-Default"

AI.box.plot.df <- rbind(AI.med.CA.Indiv.melted,AI.med.CA.LNSum.melted,AI.med.CA.Default.melted)

# Changing mixture names
AI.box.plot.df <- AI.box.plot.df %>%
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

AI.box.plot.df$variable <- factor(AI.box.plot.df$variable, 
                                  levels=c("RfD-H","RfD-L","Expo-H","Expo-L","POD-H","POD-L","AC50-H","AC50-L"))
AI.box.plot.df$group <- factor(AI.box.plot.df$group,levels=c("CA-Indiv","CA-LNSum","CA-Default"))

# Figure 3B - Plotting boxplot for AI of each mixture
min(AI.box.plot.df$value,na.rm=TRUE)
max(AI.box.plot.df$value,na.rm=TRUE)

p2 <- ggplot(AI.box.plot.df,aes(x=value,y=variable,fill=variable)) + geom_boxplot(outlier.shape = NA) +
  scale_x_log10(limits=c(0.05,20),breaks=c(0.1,1,10),labels = function(x) ifelse(x == 0, "0", x)) +
  geom_vline(xintercept=0.1,linetype="dotted",colour="gray") +
  geom_vline(xintercept=10,linetype="dotted",colour="gray") + xlab("LAI = (Measured EC10) / (CA EC10)") +
  theme_classic() + facet_wrap(group~.,ncol=1) + geom_vline(xintercept = 1,linetype="dashed") + 
  theme(axis.title.y=element_blank(),legend.title=element_blank()) + 
  scale_fill_viridis(option="turbo",discrete=TRUE,direction=-1,
                     breaks=c("AC50-L","AC50-H","POD-L","POD-H","Expo-L","Expo-H","RfD-L","RfD-H"))

p1.2 <- ggarrange(p1,p2,widths=c(1,1.3))
p1.2
