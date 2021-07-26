########################################################
### Calculate hypervolumes of individual trees in FAB ##
### Anna K Schweiger, March 2019, updates July 2021 #####

library(hypervolume)
library(FactoMineR)
library(tidyverse)

### read data
spec <- read.csv("./FAB_ind_spec.csv")

### PC axes
spec_pca <- PCA(spec[,grepl("^X",colnames(spec))],scale.unit = F, graph = F)

#### extract PC's
pcs <- as.data.frame(spec_pca$ind$coord)
out <- cbind(spec[,!(grepl("^X",colnames(spec)))],pcs) 
# write.csv(out, "./FAB_ind_specPCs.csv", row.names = F)

### Calculate hypervolumes #############
sp_list <- as.character(unique(out$IndID))

volumes <- list()
for (a in 1:length(sp_list)){
  hv_comm = new("HypervolumeList")
  hv_comm@HVList = vector(mode="list",length=2)
  set.seed(1984)
  for (b in 1:100){ ### run 100 repeats
    dat = out[out$IndID==sp_list[a],c(5:7)] ### 3 PC axes
    hv_comm@HVList[[b]] <- hypervolume_gaussian(dat, name=as.character(sp_list[a]),
                                                samples.per.point = 1000)
  }
  vols <- matrix(0,nrow = 100, ncol=2)
  for (i in 1:length(hv_comm@HVList)){
    vols[i,1] <- hv_comm@HVList[[i]]@Name
    vols[i,2] <- hv_comm@HVList[[i]]@Volume
  }
  vols <- as.data.frame(vols)
  colnames(vols) <- c("Plot","Volume")
  vols$Volume <- as.numeric(as.character(vols$Volume))
  
  volumes[[a]] <- vols
}

volcombi <- do.call(rbind, volumes)
# write.csv(volcombi, "./FAB_ind_volcombi.csv", row.names = F)

############################################
### Tree growth ~ individual spectral space
vols <- volcombi %>% group_by(Plot) %>% 
  summarize(vol_m=mean(Volume), vol_sd =sd(Volume), 
            vol_m_log=mean(log(Volume)), vol_sd_log=sd(log(Volume))) %>%
  rename(IndID=Plot)
  
datx <- read.csv("./FAB_ind_growth.csv")
dat <- vols %>% merge(datx,by="IndID", all.x=T)

### Height 2016 
plot(H_16 ~ vol_m_log, data=dat)
mod <- lm(H_16~vol_m_log, data=dat)
summary(mod)

### Height 2017
plot(H_17 ~ vol_m_log, data=dat)
mod <- lm(H_17~vol_m_log, data=dat)
summary(mod)

### Diameter 2016
plot(D_16 ~ vol_m_log, data=dat)
mod <- lm(D_16~vol_m_log, data=dat)
summary(mod)

### Diameter 2017
plot(D_17 ~ vol_m_log, data=dat)
mod <- lm(D_17~vol_m_log, data=dat)
summary(mod)

### END ################

