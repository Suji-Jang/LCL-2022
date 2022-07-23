##### Parameters setting #####
dirpath <- file.path("GitHub Data and Codes")
lcl_dat <- read.csv(file.path(dirpath,"LCL_data_040722.csv"))
lcl_dat <- subset(lcl_dat,Chem_index!=0)
lcl_dat$Chem <- factor(lcl_dat$Chem)
lcl_dat$Cell.line <- factor(lcl_dat$Cell.line)
lcl_dat$Chem_index <- factor(lcl_dat$Chem_index)
chems <- levels(lcl_dat$Chem_index)
cells <- levels(lcl_dat$Cell.line)
scale_factor <- as.integer(100)

##### GM and GSD data collecting #####
gm.gsd.df <- as.data.frame(matrix(NA,nrow=4000,ncol=0))
pod.chem.ind.df <- list()
for (chemnum in 1:50){
  fileprefix <- chemnum
  load(file.path(dirpath,paste(fileprefix,"stanfit.Rdata",sep="_")))
  fitparms_df <- as.data.frame(rstan::extract(stan_fit))
  
  temp.gm.gsd.df <- exp(fitparms_df[,c("m_x10","sd_x10")])
  temp.gm.gsd.df$sens.05 <- (temp.gm.gsd.df$m_x10)*(temp.gm.gsd.df$sd_x10)^(-1.645)  # Z(5%) =-1.645
  temp.gm.gsd.df$sens.01 <- (temp.gm.gsd.df$m_x10)*(temp.gm.gsd.df$sd_x10)^(-2.326)  # Z(1%) =-2.326
  colnames(temp.gm.gsd.df) <- c(paste0("Chem.",chemnum,".GM"),paste0("Chem.",chemnum,".GSD"),
                                paste0("Chem.",chemnum,".Sens.05"),paste0("Chem.",chemnum,".Sens.01"))
  gm.gsd.df <- cbind(gm.gsd.df,temp.gm.gsd.df)
  gm.gsd.df <- as.data.frame(gm.gsd.df)
  
  temp.chem.ind.df <- fitparms_df[,paste0("ec10.",1:146)]
  pod.chem.ind.df[[chemnum]] <- temp.chem.ind.df
}
colnames(gm.gsd.df) = gsub("Chem.43","Mix.1",colnames(gm.gsd.df))
colnames(gm.gsd.df) = gsub("Chem.44","Mix.2",colnames(gm.gsd.df))
colnames(gm.gsd.df) = gsub("Chem.45","Mix.3",colnames(gm.gsd.df))
colnames(gm.gsd.df) = gsub("Chem.46","Mix.4",colnames(gm.gsd.df))
colnames(gm.gsd.df) = gsub("Chem.47","Mix.5",colnames(gm.gsd.df))
colnames(gm.gsd.df) = gsub("Chem.48","Mix.6",colnames(gm.gsd.df))
colnames(gm.gsd.df) = gsub("Chem.49","Mix.7",colnames(gm.gsd.df))
colnames(gm.gsd.df) = gsub("Chem.50","Mix.8",colnames(gm.gsd.df))

ifelse(nrow(gm.gsd.df)>4000,gm.gsd.df<-gm.gsd.df[sample(nrow(gm.gsd.df),4000),],NA)

write.csv(gm.gsd.df,file.path(dirpath,"GM_GSD_data.csv"),row.names=FALSE)
