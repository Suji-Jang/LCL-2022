remove.packages(c("StanHeaders", "rstan"))
install.packages("rstan", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))

library(tidyr)
library(dplyr)
library(readr)
library(stringi)
library(ggplot2)
library(parallel)
library(rstan)
library(MASS)
library(Hmisc)
library(lattice)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

# Preparing data for stan fitting
dirpath <- file.path("GitHub Data and Codes")
lcl_dat <- read.csv(file.path(dirpath,"LCL_data_040722.csv"))
lcl_dat <- subset(lcl_dat,Chem_index!=0)  # Removing DMSO data
lcl_dat$Chem <- factor(lcl_dat$Chem)
lcl_dat$Cell.line <- factor(lcl_dat$Cell.line)
lcl_dat$Chem_index <- factor(lcl_dat$Chem_index)
chems <- levels(lcl_dat$Chem_index)
controls <- NA
cells <- levels(lcl_dat$Cell.line)
scale_factor <- as.integer(100)

quants <- c(0.01,0.025,0.05,0.1,0.25,0.5,0.75,0.9,0.95,0.975,0.99)
xplot<-10^(seq(-20,25)/10)
scale_factor <- 100

for (chemnum in chems) {
  onechemdat <- subset(lcl_dat,(Chem_index == chemnum))
  lclstan_dat <- list(
    scale_factor = scale_factor,
    Ni = length(unique(onechemdat$Cell.line)),
    Nj = nrow(onechemdat),
    x = onechemdat$conc,
    ys = onechemdat$value/scale_factor,
    cell = as.numeric(onechemdat$Cell.line),
    quants = quants,
    Nquants = length(quants)
  )
 fname <- paste(chems[chemnum],"stan_dat.R",sep="_")
  with(lclstan_dat, {
    stan_rdump(names(lclstan_dat),file=fname)
  })
}

# rstan fitting
seed = 314159
iter = 2000
for (chemnum in chems) {
  fileprefix <- chemnum
  onechemdat <- subset(lcl_dat,(Chem_index == chemnum))
  lclstan_dat <- list(
    scale_factor = scale_factor,
    Ni = length(unique(onechemdat$Cell.line)),
    Nj = nrow(onechemdat),
    x = onechemdat$conc,
    ys = onechemdat$value/scale_factor,
    cell = as.numeric(onechemdat$Cell.line),
    quants = quants,
    Nquants = length(quants)
  )
  
  source(paste(fileprefix,"stan_dat.R",sep="_"))
  time.start<-proc.time()
  stan_fit <- stan(file="conc_resp_zero_ec10.stan",
                    data = lclstan_dat,
                    seed=seed,iter=iter, 
                    chains = 4,
                    sample_file=paste(fileprefix,"samples",sep="_"))
   time.end<-proc.time()
   print(paste("Working on",fileprefix))
   print(time.end-time.start)
   save(stan_fit,file=paste(fileprefix,"stanfit.RData",sep="_"))
   rhat.list <- stan_rhat(stan_fit)
   pdf(file=paste0(fileprefix,"-Rhat.pdf"))
   plot(rhat.list)
   dev.off()
}
  
# Convergence check and doubling iterations if rhat is >= 1.2
seed = 314159
chem.rhat <- c()

for (chemnum in chems) {
  fileprefix <- chemnum
  load(paste(fileprefix,"stanfit.Rdata",sep="_"))
  rhat.list <- stan_rhat(stan_fit)
  rhat.dat <-  rhat.list$data$stat
  Max.Rhat <- max(rhat.dat,na.rm=TRUE)
  iter=4000
  onechemdat <- subset(lcl_dat,(Chem_index == chemnum))
  lclstan_dat <- list(
    scale_factor = scale_factor,
    Ni = length(unique(onechemdat$Cell.line)),
    Nj = nrow(onechemdat),
    x = onechemdat$conc,
    ys = onechemdat$value/scale_factor,
    cell = as.numeric(onechemdat$Cell.line),
    quants = quants,
    Nquants = length(quants)
  )
  while (Max.Rhat >= 1.2 && iter < 32000){
    source(paste(fileprefix,"stan_dat.R",sep="_"))
    time.start<-proc.time()
    
    stan_fit <- stan(file="conc_resp_zero_ec10.stan",
                     data = lclstan_dat,
                     seed=seed,iter=iter, 
                     chains = 4,
                     sample_file=paste(fileprefix,"samples",sep="_"));
    
    time.end<-proc.time()
    print(paste("Working on",fileprefix))
    print(time.end-time.start)
    save(stan_fit,file=paste(fileprefix,"stanfit2.RData",sep="_"))
    rhat.list <- stan_rhat(stan_fit)
    rhat.dat <-  rhat.list$data$stat
    Max.Rhat <- max(rhat.dat,na.rm=TRUE)
    chem.rhat[chemnum] <- Max.Rhat
    iter <- iter*2
  }
  if(iter == 32000 | Max.Rhat >= 1.2){
    cat(chemnum,"\n", "MODELING FAILED", "\n")
  } else {
    cat(chemnum,"\n", "MODELING SUCCESSFUL", "\n")
  }
} 
seed = 314159
iter = 500

for (chemnum in rhat.redone) {
  fileprefix <- chemnum
  load(paste(fileprefix,"stanfit2.Rdata",sep="_"))
  rhat.list <- stan_rhat(stan_fit)
  pdf(file=paste0(fileprefix,"-Rhat_redone.pdf"))
  plot(rhat.list)
  dev.off()
}
