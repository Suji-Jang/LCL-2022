library(reshape2)
library(ggplot2)
library(ggpubr)

# Loading POD and AI data
dirpath <- file.path("GitHub Data and Codes")
POD.med.CA.Indiv <- read.csv(file.path(dirpath,"POD_CAIndiv_Median.csv"))
POD.med.CA.LNSum <- read.csv(file.path(dirpath,"POD_CALNSum_Median.csv"))
POD.med.CA.Default <- read.csv(file.path(dirpath,"POD_CADefault_Median.csv"))
POD.sens01.CA.Indiv <- read.csv(file.path(dirpath,"POD_CAIndiv_Sens01.csv"))
POD.sens01.CA.LNSum <- read.csv(file.path(dirpath,"POD_CALNSum_Sens01.csv"))
POD.sens01.CA.Default <- read.csv(file.path(dirpath,"POD_CADefault_Sens01.csv"))

# Loading GM & GSD data
gm.gsd.df <- read.csv(file.path(dirpath,"GM_GSD_data.csv"))
POD.med.mix <- gm.gsd.df[,paste("Mix",1:8,"GM",sep=".")]
POD.sens01.mix <- gm.gsd.df[,paste("Mix",1:8,"Sens.01",sep=".")]

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

POD.sens01.mix$Mix.1.Sens.01 <- POD.sens01.mix$Mix.1.Sens.01*48.3/100
POD.sens01.mix$Mix.2.Sens.01 <- POD.sens01.mix$Mix.2.Sens.01*6236.3/100
POD.sens01.mix$Mix.3.Sens.01 <- POD.sens01.mix$Mix.3.Sens.01*2767.1/100
POD.sens01.mix$Mix.4.Sens.01 <- POD.sens01.mix$Mix.4.Sens.01*21348.4/100
POD.sens01.mix$Mix.5.Sens.01 <- POD.sens01.mix$Mix.5.Sens.01*79.4/100
POD.sens01.mix$Mix.6.Sens.01 <- POD.sens01.mix$Mix.6.Sens.01*79.9/100
POD.sens01.mix$Mix.7.Sens.01 <- POD.sens01.mix$Mix.7.Sens.01*83.8/100
POD.sens01.mix$Mix.8.Sens.01 <- POD.sens01.mix$Mix.8.Sens.01*115.7/100

colnames(POD.sens01.mix) <- paste0("Mix.",1:8)

# Deriving TDVF01
TDVF01.fitted.df <- POD.med.mix/POD.sens01.mix
colnames(TDVF01.fitted.df) <- paste0("Mix.",1:8)
TDVF01.CA1.df <- POD.med.CA.Indiv/POD.sens01.CA.Indiv
TDVF01.CA2.df <- POD.med.CA.LNSum/POD.sens01.CA.LNSum
TDVF01.CA3.df <- POD.med.CA.Default/POD.sens01.CA.Default

# Deriving AI
AI.TDVF.CA1.df <- TDVF01.CA1.df/TDVF01.fitted.df
AI.TDVF.CA2.df <- TDVF01.CA2.df/TDVF01.fitted.df
AI.TDVF.CA3.df <- TDVF01.CA3.df/TDVF01.fitted.df

TDVF01.fitted.melted <- melt(TDVF01.fitted.df)
TDVF01.CA1.melted <- melt(TDVF01.CA1.df)
TDVF01.CA2.melted <- melt(TDVF01.CA2.df)
TDVF01.CA3.melted <- melt(TDVF01.CA3.df)

TDVF01.fitted.melted$group <- "Measured"
TDVF01.CA1.melted$group <- "CA-Indiv"
TDVF01.CA2.melted$group <- "CA-LNSum"
TDVF01.CA3.melted$group <- "CA-Default"

plot.TDVF.df <- rbind(TDVF01.fitted.melted,TDVF01.CA1.melted,TDVF01.CA2.melted,TDVF01.CA3.melted)
plot.TDVF.df$group <- factor(plot.TDVF.df$group, levels=c("CA-Default","CA-LNSum","CA-Indiv","Measured"))

# Changing mixture names
plot.TDVF.df <- plot.TDVF.df %>%
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
plot.TDVF.df$variable <- factor(plot.TDVF.df$variable, levels=c("AC50-L","POD-L","Expo-L","RfD-L","AC50-H","POD-H","Expo-H","RfD-H"))

# Color setting
turbo.palette <- viridis::turbo(n=8)
turbo.palette.reorderd <- c("#30123BFF","#1BD0D5FF","#D2E935FF","#DB3A07FF","#4777EFFF","#62FC6BFF","#FE9B2DFF","#7A0403FF")

# Figure 5 - Plotting boxplot for TDVF of measured and each CA method
ggplot(plot.TDVF.df,aes(x=value,y=group,fill=variable)) + geom_boxplot(outlier.shape = NA) +
  scale_x_log10() + theme_classic() + facet_wrap(variable~.,ncol=4) + 
  xlab("TDVF01 = (EC10 Median) / (EC10 Sensitive 1st percentile)") +
  geom_vline(xintercept=3.16,linetype="dashed") + geom_vline(xintercept=1) +
  theme(axis.title.y=element_blank(),legend.position="none") +
  scale_fill_manual(values=turbo.palette.reorderd)