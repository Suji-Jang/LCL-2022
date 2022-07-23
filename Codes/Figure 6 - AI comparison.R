library(scales)
library(viridis)
# Mixture names
mixnames <- c("AC50-L","AC50-H","POD-L","POD-H","Expo-L","Expo-H","RfD-L","RfD-H")

# Nanhung data
dirpath <- file.path("CA Manuscript","2022 July - R files and Figures")
load(file.path(dirpath,"ec10_pred.rda"))
nh.mix1 <- as.data.frame(ec10.list$`AC50-L`)
nh.mix1$Mixture <- "AC50-L"
nh.mix2 <- as.data.frame(ec10.list$`AC50-H`)
nh.mix2$Mixture <- "AC50-H"
nh.mix3 <- as.data.frame(ec10.list$`POD-L`)
nh.mix3$Mixture <- "POD-L"
nh.mix4 <- as.data.frame(ec10.list$`POD-H`)
nh.mix4$Mixture <- "POD-H"
nh.mix5 <- as.data.frame(ec10.list$`Expo-L`)
nh.mix5$Mixture <- "Expo-L"
nh.mix6 <- as.data.frame(ec10.list$`Expo-H`)
nh.mix6$Mixture <- "Expo-H"
nh.mix7 <- as.data.frame(ec10.list$`RFD-L`)
nh.mix7$Mixture <- "RfD-L"
nh.mix8 <- as.data.frame(ec10.list$`RFD-H`)
nh.mix8$Mixture <- "RfD-H"

# Curve fitting = measured EC10
cf.df <- rbind(subset(nh.mix1,Method=="CF"),
               subset(nh.mix2,Method=="CF"),
               subset(nh.mix3,Method=="CF"),
               subset(nh.mix4,Method=="CF"),
               subset(nh.mix5,Method=="CF"),
               subset(nh.mix6,Method=="CF"),
               subset(nh.mix7,Method=="CF"),
               subset(nh.mix8,Method=="CF"))
# Concentration addition EC10
ca.df <- rbind(subset(nh.mix1,Method=="CA"),
               subset(nh.mix2,Method=="CA"),
               subset(nh.mix3,Method=="CA"),
               subset(nh.mix4,Method=="CA"),
               subset(nh.mix5,Method=="CA"),
               subset(nh.mix6,Method=="CA"),
               subset(nh.mix7,Method=="CA"),
               subset(nh.mix8,Method=="CA"))

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

ai.cardio.cell <- cbind(cf.cardio.cell.df[,-(1:2)],data.frame(AI=cf.cardio.cell.df$EC10/ca.cardio.cell.df$EC10))
ai.cardio.other <- cbind(cf.cardio.other.df[,-(1:2)],data.frame(AI=cf.cardio.other.df$EC10/ca.cardio.other.df$EC10))
ai.endo.cell <- cbind(cf.endo.cell.df[,-(1:2)],data.frame(AI=cf.endo.cell.df$EC10/ca.endo.cell.df$EC10))
ai.endo.other <- cbind(cf.endo.other.df[,-(1:2)],data.frame(AI=cf.endo.other.df$EC10/ca.endo.other.df$EC10))
ai.hepato.cell <- cbind(cf.hepato.cell.df[,-(1:2)],data.frame(AI=cf.hepato.cell.df$EC10/ca.hepato.cell.df$EC10))
ai.hepato.other <- cbind(cf.hepato.other.df[,-(1:2)],data.frame(AI=cf.hepato.other.df$EC10/ca.hepato.other.df$EC10))
ai.huvec.cell <- cbind(cf.huvec.cell.df[,-(1:2)],data.frame(AI=cf.huvec.cell.df$EC10/ca.huvec.cell.df$EC10))
ai.huvec.other <- cbind(cf.huvec.other.df[,-(1:2)],data.frame(AI=cf.huvec.other.df$EC10/ca.huvec.other.df$EC10))
ai.neuron.cell <- cbind(cf.neuron.cell.df[,-(1:2)],data.frame(AI=cf.neuron.cell.df$EC10/ca.neuron.cell.df$EC10))
ai.neuron.other <- cbind(cf.neuron.other.df[,-(1:2)],data.frame(AI=cf.neuron.other.df$EC10/ca.neuron.other.df$EC10))

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

# Data setting
library(reshape2)
library(ggplot2)
# setwd("C:/Users/lover/OneDrive - Texas A&M University/Projects/LCL")

AI1.med <- read.csv(file.path(dirpath,"AI_CA1_Median.csv"))
AI2.med <- read.csv(file.path(dirpath,"AI_CA2_Median.csv"))
AI3.med <- read.csv(file.path(dirpath,"AI_CA3_Median.csv"))
AI1.sens01 <- read.csv(file.path(dirpath,"AI_CA1_Sens01.csv"))
AI2.sens01 <- read.csv(file.path(dirpath,"AI_CA2_Sens01.csv"))
AI3.sens01 <- read.csv(file.path(dirpath,"AI_CA3_Sens01.csv"))
AI1.tdvf <- read.csv(file.path(dirpath,"AI_TDVF01_CA1.csv"))
AI2.tdvf <- read.csv(file.path(dirpath,"AI_TDVF01_CA2.csv"))
AI3.tdvf <- read.csv(file.path(dirpath,"AI_TDVF01_CA3.csv"))
names(AI1.med) <- mixnames
names(AI2.med) <- mixnames
names(AI3.med) <- mixnames
names(AI1.sens01) <- mixnames
names(AI2.sens01) <- mixnames
names(AI3.sens01) <- mixnames
names(AI1.tdvf) <- mixnames
names(AI2.tdvf) <- mixnames
names(AI3.tdvf) <- mixnames

AI1.med.melted <- melt(AI1.med,value.name = "AI")
AI2.med.melted <- melt(AI2.med,value.name = "AI")
AI3.med.melted <- melt(AI3.med,value.name = "AI")
AI1.med.melted$group <- "LCL-PopMed-CA-Indiv"
AI2.med.melted$group <- "LCL-PopMed-CA-LNSum"
AI3.med.melted$group <- "LCL-PopMed-CA-Default"
names(AI1.med.melted)[1]<-"Mixture"
names(AI2.med.melted)[1]<-"Mixture"
names(AI3.med.melted)[1]<-"Mixture"

AI1.sens01.melted <- melt(AI1.sens01,value.name = "AI")
AI2.sens01.melted <- melt(AI2.sens01,value.name = "AI")
AI3.sens01.melted <- melt(AI3.sens01,value.name = "AI")
AI1.sens01.melted$group <- "LCL-Sens01-CA-Indiv"
AI2.sens01.melted$group <- "LCL-Sens01-CA-LNSum"
AI3.sens01.melted$group <- "LCL-Sens01-CA-Default"
names(AI1.sens01.melted)[1]<-"Mixture"
names(AI2.sens01.melted)[1]<-"Mixture"
names(AI3.sens01.melted)[1]<-"Mixture"

AI1.tdvf.melted <- melt(AI1.tdvf,value.name = "AI")
AI2.tdvf.melted <- melt(AI2.tdvf,value.name = "AI")
AI3.tdvf.melted <- melt(AI3.tdvf,value.name = "AI")
AI1.tdvf.melted$group <- "LCL-TDVF01-CA-Indiv"
AI2.tdvf.melted$group <- "LCL-TDVF01-CA-LNSum"
AI3.tdvf.melted$group <- "LCL-TDVF01-CA-Default"
names(AI1.tdvf.melted)[1]<-"Mixture"
names(AI2.tdvf.melted)[1]<-"Mixture"
names(AI3.tdvf.melted)[1]<-"Mixture"

# boxplot.df <- rbind(ai.cardio.cell[,-(1:2)],
#                     ai.cardio.other[,-(1:2)],
#                     ai.endo.cell[,-(1:2)],
#                     ai.endo.other[,-(1:2)],
#                     ai.hepato.cell[,-(1:2)],
#                     ai.hepato.other[,-(1:2)],
#                     ai.huvec.cell[,-(1:2)],
#                     ai.huvec.other[,-(1:2)],
#                     ai.neuron.cell[,-(1:2)],
#                     ai.neuron.other[,-(1:2)],
#                     AI1.med.melted,AI2.med.melted,AI3.med.melted,
#                     AI1.sens01.melted,AI2.sens01.melted,AI3.sens01.melted,
#                     AI1.tdvf.melted,AI2.tdvf.melted,AI3.tdvf.melted)
# boxplot.df$Mixture <- factor(boxplot.df$Mixture,
#                              levels=mixnames)
# boxplot.df$group <- factor(boxplot.df$group,levels=c("LCL-TDVF01-CA-Default","LCL-TDVF01-CA-LNSum","LCL-TDVF01-CA-Indiv",
#                                                "LCL-Sens01-CA-Default","LCL-Sens01-CA-LNSum","LCL-Sens01-CA-Indiv",
#                                                "LCL-PopMed-CA-Default","LCL-PopMed-CA-LNSum","LCL-PopMed-CA-Indiv",
#                                                "Neuron-Others","Neuron-Cell#","HUVEC-Others","HUVEC-Cell#",
#                                                "Hepato-Others","Hepato-Cell#","Endo-Others","Endo-Cell#","Cardio-Others","Cardio-Cell#"))


##

quantiles_95CI <- function(x) {
  r <- quantile(x, probs=c(0.025, 0.25, 0.5, 0.75, 0.975))
  names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
  r
}

boxplotall.df <- rbind(ai.cardio.cell[,c(3:5,1,2)],
                    ai.cardio.other[,c(3:5,1,2)],
                    ai.endo.cell[,c(3:5,1,2)],
                    ai.endo.other[,c(3:5,1,2)],
                    ai.hepato.cell[,c(3:5,1,2)],
                    ai.hepato.other[,c(3:5,1,2)],
                    ai.huvec.cell[,c(3:5,1,2)],
                    ai.huvec.other[,c(3:5,1,2)],
                    ai.neuron.cell[,c(3:5,1,2)],
                    ai.neuron.other[,c(3:5,1,2)],
                    data.frame(AI1.med.melted,Phenotype="Viability (Median)",Celltype="LCL (CA-Indiv)"),
                    data.frame(AI2.med.melted,Phenotype="Viability (Median)",Celltype="LCL (CA-LNSum)"),
                    data.frame(AI3.med.melted,Phenotype="Viability (Median)",Celltype="LCL (CA-Default)"),
                    data.frame(AI1.sens01.melted,Phenotype="Viability (Sens01)",Celltype="LCL (CA-Indiv)"),
                    data.frame(AI2.sens01.melted,Phenotype="Viability (Sens01)",Celltype="LCL (CA-LNSum)"),
                    data.frame(AI3.sens01.melted,Phenotype="Viability (Sens01)",Celltype="LCL (CA-Default)"),
                    data.frame(AI1.tdvf.melted,Phenotype="Viability (TDVF01)",Celltype="LCL (CA-Indiv)"),
                    data.frame(AI2.tdvf.melted,Phenotype="Viability (TDVF01)",Celltype="LCL (CA-LNSum)"),
                    data.frame(AI3.tdvf.melted,Phenotype="Viability (TDVF01)",Celltype="LCL (CA-Default)"))
boxplotall.df$Mixture <- factor(boxplotall.df$Mixture,
                             levels=mixnames)
boxplotall.df$group <- paste(boxplotall.df$Celltype,boxplotall.df$Phenotype)
boxplotall.df$group <- factor(boxplotall.df$group,
                              levels=rev(unique(boxplotall.df$group)))
set.seed(3.14159)
sampindxall <- sample.int(nrow(boxplotall.df),
                       size=nrow(boxplotall.df)*0.1) # 10% sample

supfigAIall<-ggplot(boxplotall.df,aes(x=AI,y=group)) + 
  geom_jitter(data=boxplotall.df[sampindxall,],
              aes(color=Mixture,y=group),height=0.1,size=0.1,alpha=0.5) +
  stat_summary(fun.data = quantiles_95CI, geom="boxplot",fill=NA)+
  scale_x_log10(breaks = 10^seq(-4,4),
                labels = trans_format("log10", math_format(10^.x)))+
  coord_cartesian(xlim=c(1e-4,1e4))+
  geom_vline(xintercept=0.1,linetype="dotted",colour="gray") +
  geom_vline(xintercept=10,linetype="dotted",colour="gray") + 
  xlab("LAI = (measured EC10 or TDVF01)/(CA EC10 or TDVF01)") +
  theme_classic() + geom_vline(xintercept = 1,linetype="dashed") +
  scale_color_viridis(option="turbo",discrete=TRUE)+
  guides(color = guide_legend(override.aes = list(size=1,alpha=1)))+
  theme(axis.title.y=element_blank())
print(supfigAIall)
ggsave(file.path(dirpath,"SuppFigureAIall.pdf"),supfigAIall,
       dpi=600,height=8,width=6,scale=1.1)

### ANOVA 2-way
library(lsr)
boxplotall.df$log10AI <- log10(boxplotall.df$AI)
res <- lm(log10AI~Mixture+group,data=subset(boxplotall.df,!(Celltype %in% c("LCL (CA-LNSum)","LCL (CA-Default)"))))
aovres <- aov(res)
print(summary(res))
print(summary(aovres))
print(etaSquared(aovres))


boxplot.cytotox.df <- subset(boxplotall.df,
                             Phenotype == "Cell Number" | Celltype == "LCL (CA-Indiv)")
set.seed(3.14159)
sampindx <- sample.int(nrow(boxplot.cytotox.df),
                       size=nrow(boxplot.cytotox.df)*0.1) # 10% sample

fig6<-ggplot(boxplot.cytotox.df,aes(x=AI,y=group)) + 
  geom_jitter(data=boxplot.cytotox.df[sampindx,],
              aes(color=Mixture,y=group),height=0.1,size=0.1,alpha=0.5) +
  stat_summary(fun.data = quantiles_95CI, geom="boxplot",fill=NA)+
  scale_x_log10(breaks = 10^seq(-4,4),
                labels = trans_format("log10", math_format(10^.x)))+
  coord_cartesian(xlim=c(1e-4,1e4))+
  geom_vline(xintercept=0.1,linetype="dotted",colour="gray") +
  geom_vline(xintercept=10,linetype="dotted",colour="gray") + 
  xlab("LAI = (measured EC10 or TDVF01)/(CA EC10 or TDVF01)") +
  theme_classic() + geom_vline(xintercept = 1,linetype="dashed") +
  scale_color_viridis(option="turbo",discrete=TRUE)+
  guides(color = guide_legend(override.aes = list(size=1,alpha=1)))+
  theme(axis.title.y=element_blank())
print(fig6)
ggsave(file.path(dirpath,"Figure 6.pdf"),fig6,
       dpi=600,height=4,width=6,scale=1.1)


# #
# 
# boxplot.cytotox.df <- rbind(ai.cardio.cell[,-(1:2)],
#                             ai.endo.cell[,-(1:2)],
#                             ai.hepato.cell[,-(1:2)],
#                             ai.huvec.cell[,-(1:2)],
#                             ai.neuron.cell[,-(1:2)],
#                             AI1.med.melted,
#                             AI1.sens01.melted,
#                             AI1.tdvf.melted)
# boxplot.cytotox.df$Mixture <- factor(boxplot.cytotox.df$Mixture,
#                                      levels=mixnames)
# boxplot.cytotox.df$group <- factor(boxplot.cytotox.df$group,levels=c("LCL-TDVF01-CA-Indiv","LCL-Sens01-CA-Indiv","LCL-PopMed-CA-Indiv",
#                                                                      "Neuron-Cell#","HUVEC-Cell#","Hepato-Cell#","Endo-Cell#","Cardio-Cell#"))
# 
# set.seed(3.14159)
# sampindx <- sample.int(nrow(boxplot.cytotox.df),
#                        size=nrow(boxplot.cytotox.df)*0.1) # 10% sample
# 
# fig6<-ggplot(boxplot.cytotox.df,aes(x=AI,y=group)) + 
#   geom_jitter(data=boxplot.cytotox.df[sampindx,],
#               aes(color=Mixture,y=group),height=0.1,size=0.1,alpha=0.5) +
#   stat_summary(fun.data = quantiles_95CI, geom="boxplot",fill=NA)+
#   scale_x_log10(breaks = 10^seq(-4,4),
#                 labels = trans_format("log10", math_format(10^.x)))+
#   coord_cartesian(xlim=c(1e-4,1e4))+
#   geom_vline(xintercept=0.1,linetype="dotted",colour="gray") +
#   geom_vline(xintercept=10,linetype="dotted",colour="gray") + 
#   xlab("LAI = (measured EC10 or TDVF01)/(CA EC10 or TDVF01)") +
#   theme_classic() + geom_vline(xintercept = 1,linetype="dashed") +
#   scale_color_viridis(option="turbo",discrete=TRUE)+
#   guides(color = guide_legend(override.aes = list(size=1,alpha=1)))+
#   theme(axis.title.y=element_blank())
# print(fig6)
# ggsave(file.path(dirpath,"Figure6.pdf"),fig6,dpi=600,height=4,width=6)
# 

# Table
mix.ai.cyto.table <- 
aggregate(log10(AI) ~ group,
          data=subset(boxplotall.df,
                      Phenotype == "Cell Number" | Celltype == "LCL (CA-Indiv)"
                      | Celltype == "LCL (CA-LNSum)" | Celltype == "LCL (CA-Default)"
                      ),FUN=quantile,prob=c(0.025,0.25,0.5,0.75,0.975))
write.csv(mix.ai.cyto.table,file.path(dirpath,"Mixture Log10 AI Cytotox.csv"))

mix.ai.noncyto.table <- 
  aggregate(log10(AI) ~ Celltype,data=
              subset(boxplotall.df,
                     Phenotype != "Cell Number" & Celltype != "LCL (CA-Indiv)"
                     & Celltype != "LCL (CA-LNSum)" & Celltype != "LCL (CA-Default)"
                     ),FUN=quantile,prob=c(0.025,0.25,0.5,0.75,0.975))
write.csv(mix.ai.noncyto.table,file.path(dirpath,"Mixture Log10 AI Non-Cytotox.csv"))


mix.ai.all.table <- 
  aggregate(log10(AI) ~ group,data=boxplotall.df,FUN=quantile,prob=c(0.025,0.25,0.5,0.75,0.975))
write.csv(mix.ai.all.table,file.path(dirpath,"Mixture Log10 AI All Phenotypes.csv"))

# 
# mix.uf.table <- as.data.frame(t(unname(quantile(AI1.med.melted$AI,probs=c(0.05,0.5,0.95),na.rm=TRUE))))
# mix.uf.table <- cbind(mix.uf.table,as.data.frame(t(unname(quantile(AI2.med.melted$AI,probs=c(0.05,0.5,0.95),na.rm=TRUE)))),
#                       as.data.frame(t(unname(quantile(AI3.med.melted$AI,probs=c(0.05,0.5,0.95),na.rm=TRUE)))))
# temp.AIs <- cbind(as.data.frame(t(unname(quantile(AI1.sens01.melted$AI,probs=c(0.05,0.5,0.95),na.rm=TRUE)))),
#                      as.data.frame(t(unname(quantile(AI2.sens01.melted$AI,probs=c(0.05,0.5,0.95),na.rm=TRUE)))),
#                      as.data.frame(t(unname(quantile(AI3.sens01.melted$AI,probs=c(0.05,0.5,0.95),na.rm=TRUE)))))
# mix.uf.table <- rbind(mix.uf.table,temp.AIs)
# temp.AIs <- cbind(as.data.frame(t(unname(quantile(AI1.tdvf.melted$AI,probs=c(0.05,0.5,0.95),na.rm=TRUE)))),
#                      as.data.frame(t(unname(quantile(AI2.tdvf.melted$AI,probs=c(0.05,0.5,0.95),na.rm=TRUE)))),
#                      as.data.frame(t(unname(quantile(AI3.tdvf.melted$AI,probs=c(0.05,0.5,0.95),na.rm=TRUE)))))
# mix.uf.table <- rbind(mix.uf.table,temp.AIs)
# temp.AIs <- cbind(as.data.frame(t(unname(quantile(ai.cardio.cell$AI,probs=c(0.05,0.5,0.95),na.rm=TRUE)))),
#                      as.data.frame(t(unname(quantile(ai.cardio.other$AI,probs=c(0.05,0.5,0.95),na.rm=TRUE)))),NA,NA,NA)
# colnames(temp.AIs) <- c("V1","V2","V3","V1","V2","V3","V1","V2","V3")
# mix.uf.table <- rbind(mix.uf.table,temp.AIs)
# temp.AIs <- cbind(as.data.frame(t(unname(quantile(ai.endo.cell$AI,probs=c(0.05,0.5,0.95),na.rm=TRUE)))),
#                      as.data.frame(t(unname(quantile(ai.endo.other$AI,probs=c(0.05,0.5,0.95),na.rm=TRUE)))),NA,NA,NA)
# colnames(temp.AIs) <- c("V1","V2","V3","V1","V2","V3","V1","V2","V3")
# mix.uf.table <- rbind(mix.uf.table,temp.AIs)
# temp.AIs <- cbind(as.data.frame(t(unname(quantile(ai.hepato.cell$AI,probs=c(0.05,0.5,0.95),na.rm=TRUE)))),
#                      as.data.frame(t(unname(quantile(ai.hepato.other$AI,probs=c(0.05,0.5,0.95),na.rm=TRUE)))),NA,NA,NA)
# colnames(temp.AIs) <- c("V1","V2","V3","V1","V2","V3","V1","V2","V3")
# mix.uf.table <- rbind(mix.uf.table,temp.AIs)
# temp.AIs <- cbind(as.data.frame(t(unname(quantile(ai.huvec.cell$AI,probs=c(0.05,0.5,0.95),na.rm=TRUE)))),
#                      as.data.frame(t(unname(quantile(ai.huvec.other$AI,probs=c(0.05,0.5,0.95),na.rm=TRUE)))),NA,NA,NA)
# colnames(temp.AIs) <- c("V1","V2","V3","V1","V2","V3","V1","V2","V3")
# mix.uf.table <- rbind(mix.uf.table,temp.AIs)
# temp.AIs <- cbind(as.data.frame(t(unname(quantile(ai.neuron.cell$AI,probs=c(0.05,0.5,0.95),na.rm=TRUE)))),
#                      as.data.frame(t(unname(quantile(ai.neuron.other$AI,probs=c(0.05,0.5,0.95),na.rm=TRUE)))),NA,NA,NA)
# colnames(temp.AIs) <- c("V1","V2","V3","V1","V2","V3","V1","V2","V3")
# mix.uf.table <- rbind(mix.uf.table,temp.AIs)
# 

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
