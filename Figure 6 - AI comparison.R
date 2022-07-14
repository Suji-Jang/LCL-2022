# Nanhung data
load("C:/Users/lover/OneDrive - Texas A&M University/Projects/LCL/ec10_pred.rda")
nh.mix1 <- as.data.frame(ec10.list$`AC50-L`)
nh.mix2 <- as.data.frame(ec10.list$`AC50-H`)
nh.mix3 <- as.data.frame(ec10.list$`POD-L`)
nh.mix4 <- as.data.frame(ec10.list$`POD-H`)
nh.mix5 <- as.data.frame(ec10.list$`Expo-L`)
nh.mix6 <- as.data.frame(ec10.list$`Expo-H`)
nh.mix7 <- as.data.frame(ec10.list$`RFD-L`)
nh.mix8 <- as.data.frame(ec10.list$`RFD-H`)

cf.df <- rbind(subset(nh.mix1,Method=="CF"),subset(nh.mix2,Method=="CF"),subset(nh.mix3,Method=="CF"),
               subset(nh.mix4,Method=="CF"),subset(nh.mix5,Method=="CF"),subset(nh.mix6,Method=="CF"),
               subset(nh.mix7,Method=="CF"),subset(nh.mix8,Method=="CF"))
ca.df <- rbind(subset(nh.mix1,Method=="CA"),subset(nh.mix2,Method=="CA"),subset(nh.mix3,Method=="CA"),
               subset(nh.mix4,Method=="CA"),subset(nh.mix5,Method=="CA"),subset(nh.mix6,Method=="CA"),
               subset(nh.mix7,Method=="CA"),subset(nh.mix8,Method=="CA"))


cf.cardio.df <- subset(cf.df,Celltype=="iCell Cardiomyocytes")
cf.endo.df <- subset(cf.df,Celltype=="iCell Endothelial cells")
cf.hepato.df <- subset(cf.df,Celltype=="iCell Hepatocytes")
cf.huvec.df <- subset(cf.df,Celltype=="HUVECs")
cf.neuron.df <- subset(cf.df,Celltype=="iCell Neurons")

ca.cardio.df <- subset(ca.df,Celltype=="iCell Cardiomyocytes")
ca.endo.df <- subset(ca.df,Celltype=="iCell Endothelial cells")
ca.hepato.df <- subset(ca.df,Celltype=="iCell Hepatocytes")
ca.huvec.df <- subset(ca.df,Celltype=="HUVECs")
ca.neuron.df <- subset(ca.df,Celltype=="iCell Neurons")

cf.cardio.cell.df <- subset(cf.cardio.df,Phenotype=="Cell Number")
cf.cardio.other.df <- subset(cf.cardio.df,Phenotype!="Cell Number")
cf.endo.cell.df <- subset(cf.endo.df,Phenotype=="Cell Number")
cf.endo.other.df <- subset(cf.endo.df,Phenotype!="Cell Number")
cf.hepato.cell.df <- subset(cf.hepato.df,Phenotype=="Cell Number")
cf.hepato.other.df <- subset(cf.hepato.df,Phenotype!="Cell Number")
cf.huvec.cell.df <- subset(cf.huvec.df,Phenotype=="Cell Number")
cf.huvec.other.df <- subset(cf.huvec.df,Phenotype!="Cell Number")
cf.neuron.cell.df <- subset(cf.neuron.df,Phenotype=="Cell Number")
cf.neuron.other.df <- subset(cf.neuron.df,Phenotype!="Cell Number")

ca.cardio.cell.df <- subset(ca.cardio.df,Phenotype=="Cell Number")
ca.cardio.other.df <- subset(ca.cardio.df,Phenotype!="Cell Number")
ca.endo.cell.df <- subset(ca.endo.df,Phenotype=="Cell Number")
ca.endo.other.df <- subset(ca.endo.df,Phenotype!="Cell Number")
ca.hepato.cell.df <- subset(ca.hepato.df,Phenotype=="Cell Number")
ca.hepato.other.df <- subset(ca.hepato.df,Phenotype!="Cell Number")
ca.huvec.cell.df <- subset(ca.huvec.df,Phenotype=="Cell Number")
ca.huvec.other.df <- subset(ca.huvec.df,Phenotype!="Cell Number")
ca.neuron.cell.df <- subset(ca.neuron.df,Phenotype=="Cell Number")
ca.neuron.other.df <- subset(ca.neuron.df,Phenotype!="Cell Number")

ai.cardio.cell <- as.data.frame(cf.cardio.cell.df$EC10/ca.cardio.cell.df$EC10)
ai.cardio.other <- as.data.frame(cf.cardio.other.df$EC10/ca.cardio.other.df$EC10)
ai.endo.cell <- as.data.frame(cf.endo.cell.df$EC10/ca.endo.cell.df$EC10)
ai.endo.other <- as.data.frame(cf.endo.other.df$EC10/ca.endo.other.df$EC10)
ai.hepato.cell <- as.data.frame(cf.hepato.cell.df$EC10/ca.hepato.cell.df$EC10)
ai.hepato.other <- as.data.frame(cf.hepato.other.df$EC10/ca.hepato.other.df$EC10)
ai.huvec.cell <- as.data.frame(cf.huvec.cell.df$EC10/ca.huvec.cell.df$EC10)
ai.huvec.other <- as.data.frame(cf.huvec.other.df$EC10/ca.huvec.other.df$EC10)
ai.neuron.cell <- as.data.frame(cf.neuron.cell.df$EC10/ca.neuron.cell.df$EC10)
ai.neuron.other <- as.data.frame(cf.neuron.other.df$EC10/ca.neuron.other.df$EC10)

ai.cardio.cell$group <- "Cardio-Cell#"
ai.cardio.other$group <- "Cardio-Others"
ai.endo.cell$group <- "Endo-Cell#"
ai.endo.other$group <- "Endo-Others"
ai.hepato.cell$group <- "Hepato-Cell#"
ai.hepato.other$group <- "Hepato-Others"
ai.huvec.cell$group <- "HUVEC-Cell#"
ai.huvec.other$group <- "HUVEC-Others"
ai.neuron.cell$group <- "Neuron-Cell#"
ai.neuron.other$group <- "Neuron-Others"

colnames(ai.cardio.cell) <- c("value","group")
colnames(ai.cardio.other) <- c("value","group")
colnames(ai.endo.cell) <- c("value","group")
colnames(ai.endo.other) <- c("value","group")
colnames(ai.hepato.cell) <- c("value","group")
colnames(ai.hepato.other) <- c("value","group")
colnames(ai.huvec.cell) <- c("value","group")
colnames(ai.huvec.other) <- c("value","group")
colnames(ai.neuron.cell) <- c("value","group")
colnames(ai.neuron.other) <- c("value","group")

# Data setting
library(reshape2)
library(ggplot2)
setwd("C:/Users/lover/OneDrive - Texas A&M University/Projects/LCL")

AI1.med <- read.csv("AI_CA1_Median.csv")
AI2.med <- read.csv("AI_CA2_Median.csv")
AI3.med <- read.csv("AI_CA3_Median.csv")
AI1.sens01 <- read.csv("AI_CA1_Sens01.csv")
AI2.sens01 <- read.csv("AI_CA2_Sens01.csv")
AI3.sens01 <- read.csv("AI_CA3_Sens01.csv")
AI1.tdvf <- read.csv("AI_TDVF01_CA1.csv")
AI2.tdvf <- read.csv("AI_TDVF01_CA2.csv")
AI3.tdvf <- read.csv("AI_TDVF01_CA3.csv")

AI1.med.melted <- melt(AI1.med)
AI2.med.melted <- melt(AI2.med)
AI3.med.melted <- melt(AI3.med)
AI1.med.melted$group <- "LCL-PopMed-CA-Indiv"
AI2.med.melted$group <- "LCL-PopMed-CA-LNSum"
AI3.med.melted$group <- "LCL-PopMed-CA-Default"

AI1.sens01.melted <- melt(AI1.sens01)
AI2.sens01.melted <- melt(AI2.sens01)
AI3.sens01.melted <- melt(AI3.sens01)
AI1.sens01.melted$group <- "LCL-Sens01-CA-Indiv"
AI2.sens01.melted$group <- "LCL-Sens01-CA-LNSum"
AI3.sens01.melted$group <- "LCL-Sens01-CA-Default"

AI1.tdvf.melted <- melt(AI1.tdvf)
AI2.tdvf.melted <- melt(AI2.tdvf)
AI3.tdvf.melted <- melt(AI3.tdvf)
AI1.tdvf.melted$group <- "LCL-TDVF01-CA-Indiv"
AI2.tdvf.melted$group <- "LCL-TDVF01-CA-LNSum"
AI3.tdvf.melted$group <- "LCL-TDVF01-CA-Default"

#boxplot.df <- rbind(ai.cardio.cell,ai.cardio.other,ai.endo.cell,ai.endo.other,ai.hepato.cell,
#                    ai.hepato.other,ai.huvec.cell,ai.huvec.other,ai.neuron.cell,ai.neuron.other,
#                    AI1.med.melted,AI2.med.melted,AI3.med.melted,
#                    AI1.sens01.melted,AI2.sens01.melted,AI3.sens01.melted,
#                    AI1.tdvf.melted,AI2.tdvf.melted,AI3.tdvf.melted)
#boxplot.df$group <- factor(boxplot.df$group,levels=c("LCL-TDVF01-CA-Default","LCL-TDVF01-CA-LNSum","LCL-TDVF01-CA-Indiv",
                                               "LCL-Sens01-CA-Default","LCL-Sens01-CA-LNSum","LCL-Sens01-CA-Indiv",
                                               "LCL-PopMed-CA-Default","LCL-PopMed-CA-LNSum","LCL-PopMed-CA-Indiv",
                                               "Neuron-Others","Neuron-Cell#","HUVEC-Others","HUVEC-Cell#",
                                               "Hepato-Others","Hepato-Cell#","Endo-Others","Endo-Cell#","Cardio-Others","Cardio-Cell#"))

boxplot.df <- rbind(ai.cardio.cell,ai.endo.cell,ai.hepato.cell,ai.huvec.cell,ai.neuron.cell,
                    AI1.med.melted,AI1.sens01.melted,AI1.tdvf.melted)
boxplot.df$group <- factor(boxplot.df$group,levels=c("LCL-TDVF01-CA-Indiv","LCL-Sens01-CA-Indiv","LCL-PopMed-CA-Indiv",
                                                     "Neuron-Cell#","HUVEC-Cell#","Hepato-Cell#","Endo-Cell#","Cardio-Cell#"))


min(boxplot.df$value,na.rm=TRUE)
max(boxplot.df$value,na.rm=TRUE)

ggplot(boxplot.df,aes(x=value,y=group)) + geom_boxplot(outlier.shape = NA) +
  geom_jitter(alpha=0.1) +
  scale_x_log10(limits=c(1.5e-5,1.5e5),breaks=10^seq(-4,4),labels = function(x) ifelse(x == 0, "0", x)) + #coord_cartesian() +
  geom_vline(xintercept=0.1,linetype="dotted",colour="gray") +
  geom_vline(xintercept=10,linetype="dotted",colour="gray") + xlab("AI") +
  theme_classic() + geom_vline(xintercept = 1,linetype="dashed") +
  theme(axis.title.y=element_blank())


# Table
mix.uf.table <- as.data.frame(t(unname(quantile(AI1.med.melted$value,probs=c(0.05,0.5,0.95),na.rm=TRUE))))
mix.uf.table <- cbind(mix.uf.table,as.data.frame(t(unname(quantile(AI2.med.melted$value,probs=c(0.05,0.5,0.95),na.rm=TRUE)))),
                      as.data.frame(t(unname(quantile(AI3.med.melted$value,probs=c(0.05,0.5,0.95),na.rm=TRUE)))))
temp.values <- cbind(as.data.frame(t(unname(quantile(AI1.sens01.melted$value,probs=c(0.05,0.5,0.95),na.rm=TRUE)))),
                     as.data.frame(t(unname(quantile(AI2.sens01.melted$value,probs=c(0.05,0.5,0.95),na.rm=TRUE)))),
                     as.data.frame(t(unname(quantile(AI3.sens01.melted$value,probs=c(0.05,0.5,0.95),na.rm=TRUE)))))
mix.uf.table <- rbind(mix.uf.table,temp.values)
temp.values <- cbind(as.data.frame(t(unname(quantile(AI1.tdvf.melted$value,probs=c(0.05,0.5,0.95),na.rm=TRUE)))),
                     as.data.frame(t(unname(quantile(AI2.tdvf.melted$value,probs=c(0.05,0.5,0.95),na.rm=TRUE)))),
                     as.data.frame(t(unname(quantile(AI3.tdvf.melted$value,probs=c(0.05,0.5,0.95),na.rm=TRUE)))))
mix.uf.table <- rbind(mix.uf.table,temp.values)
temp.values <- cbind(as.data.frame(t(unname(quantile(ai.cardio.cell$value,probs=c(0.05,0.5,0.95),na.rm=TRUE)))),
                     as.data.frame(t(unname(quantile(ai.cardio.other$value,probs=c(0.05,0.5,0.95),na.rm=TRUE)))),NA,NA,NA)
colnames(temp.values) <- c("V1","V2","V3","V1","V2","V3","V1","V2","V3")
mix.uf.table <- rbind(mix.uf.table,temp.values)
temp.values <- cbind(as.data.frame(t(unname(quantile(ai.endo.cell$value,probs=c(0.05,0.5,0.95),na.rm=TRUE)))),
                     as.data.frame(t(unname(quantile(ai.endo.other$value,probs=c(0.05,0.5,0.95),na.rm=TRUE)))),NA,NA,NA)
colnames(temp.values) <- c("V1","V2","V3","V1","V2","V3","V1","V2","V3")
mix.uf.table <- rbind(mix.uf.table,temp.values)
temp.values <- cbind(as.data.frame(t(unname(quantile(ai.hepato.cell$value,probs=c(0.05,0.5,0.95),na.rm=TRUE)))),
                     as.data.frame(t(unname(quantile(ai.hepato.other$value,probs=c(0.05,0.5,0.95),na.rm=TRUE)))),NA,NA,NA)
colnames(temp.values) <- c("V1","V2","V3","V1","V2","V3","V1","V2","V3")
mix.uf.table <- rbind(mix.uf.table,temp.values)
temp.values <- cbind(as.data.frame(t(unname(quantile(ai.huvec.cell$value,probs=c(0.05,0.5,0.95),na.rm=TRUE)))),
                     as.data.frame(t(unname(quantile(ai.huvec.other$value,probs=c(0.05,0.5,0.95),na.rm=TRUE)))),NA,NA,NA)
colnames(temp.values) <- c("V1","V2","V3","V1","V2","V3","V1","V2","V3")
mix.uf.table <- rbind(mix.uf.table,temp.values)
temp.values <- cbind(as.data.frame(t(unname(quantile(ai.neuron.cell$value,probs=c(0.05,0.5,0.95),na.rm=TRUE)))),
                     as.data.frame(t(unname(quantile(ai.neuron.other$value,probs=c(0.05,0.5,0.95),na.rm=TRUE)))),NA,NA,NA)
colnames(temp.values) <- c("V1","V2","V3","V1","V2","V3","V1","V2","V3")
mix.uf.table <- rbind(mix.uf.table,temp.values)


# Reproducibility to paper
cf.df <- subset(nh.mix1,Method=="CF")
ca.df <- subset(nh.mix1,Method=="CA")

cf.cardio.df <- subset(cf.df,Celltype=="iCell Cardiomyocytes")
cf.endo.df <- subset(cf.df,Celltype=="iCell Endothelial cells")
cf.hepato.df <- subset(cf.df,Celltype=="iCell Hepatocytes")
cf.huvec.df <- subset(cf.df,Celltype=="HUVECs")
cf.neuron.df <- subset(cf.df,Celltype=="iCell Neurons")

ca.cardio.df <- subset(ca.df,Celltype=="iCell Cardiomyocytes")
ca.endo.df <- subset(ca.df,Celltype=="iCell Endothelial cells")
ca.hepato.df <- subset(ca.df,Celltype=="iCell Hepatocytes")
ca.huvec.df <- subset(ca.df,Celltype=="HUVECs")
ca.neuron.df <- subset(ca.df,Celltype=="iCell Neurons")

cf.cardio.cell.df <- subset(cf.cardio.df,Phenotype=="Cell Number")
cf.endo.cell.df <- subset(cf.endo.df,Phenotype=="Cell Number")
cf.hepato.cell.df <- subset(cf.hepato.df,Phenotype=="Cell Number")
cf.huvec.cell.df <- subset(cf.huvec.df,Phenotype=="Cell Number")
cf.neuron.cell.df <- subset(cf.neuron.df,Phenotype=="Cell Number")

ca.cardio.cell.df <- subset(ca.cardio.df,Phenotype=="Cell Number")
ca.endo.cell.df <- subset(ca.endo.df,Phenotype=="Cell Number")
ca.hepato.cell.df <- subset(ca.hepato.df,Phenotype=="Cell Number")
ca.huvec.cell.df <- subset(ca.huvec.df,Phenotype=="Cell Number")
ca.neuron.cell.df <- subset(ca.neuron.df,Phenotype=="Cell Number")

newplot.df <- subset(nh.mix1,Method!="IA")
library(dplyr)
newplot.df <- newplot.df %>%
  mutate(Celltype = case_when(
    Celltype=="iCell Cardiomyocytes" ~ "Cardio",
    Celltype=="iCell Endothelial cells" ~ "Endo",
    Celltype=="iCell Hepatocytes" ~ "Hepato",
    Celltype=="HUVECs" ~ "HUVEC",
    Celltype=="iCell Neurons" ~ "Neuron",
  ))

newplot.df$Cell.Pheno <- paste(newplot.df$Celltype,newplot.df$Phenotype)
bymedian <- with(newplot.df, reorder(Cell.Pheno, EC10, median))

ggplot(newplot.df,aes(x=EC10,y=bymedian,fill=Method)) + geom_boxplot(outlier.shape = NA) + scale_x_log10() +
  geom_vline(xintercept=0.1,linetype="dotted",colour="gray") + coord_cartesian() +
  geom_vline(xintercept=10,linetype="dotted",colour="gray") + xlab("EC10") +
  theme_classic() + geom_vline(xintercept = 1,linetype="dashed") +
  theme(axis.title.y=element_blank())
