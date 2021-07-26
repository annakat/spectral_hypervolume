##########################################################
### Calculate hypervolumes of plant communities in FAB ###
### Anna K Schweiger, March 2019, updates July 2021 ######

library(hypervolume)
library(FactoMineR)
library(tidyverse)
spec <- read.csv("./FAB_ind_spec.csv")

### 50 random selections of 9 measurements per plot
rndid <- list()
for (i in 1:50){
  rndid [[i]]<- with(spec, ave(1:nrow(spec), Plot, FUN=function(x) {sample.int(length(x))}))
}
# saveRDS(rndid, "./R_output/FAB_rndid.rds")

pc_out <- list()
volcombi <- list()

for (i in 1:50){
  insample <- rndid[[i]] <=9
  sub <- spec[insample,]

  ### PCA
  spec_pca <- PCA(sub[,grepl("X", colnames(sub))],scale.unit = F, graph = F) 
  pcs <- as.data.frame(spec_pca$ind$coord) ### extract PC's (scores)
  pc_out[[i]] <- cbind(sub[,!(grepl("X", colnames(sub)))],pcs) 
  
  ### Hypervolume calculations
  hypdata <- pc_out[[i]]
  plot_list <- as.character(unique(hypdata$Plot))
  volumes <- list()
  
  for (a in 1:length(plot_list)){ 
    hv_comm = new("HypervolumeList")
    hv_comm@HVList = vector(mode="list",length=2)
    set.seed(1984)
    for (b in 1:25){ 
      dat = hypdata[hypdata$Plot==plot_list[a],c("Dim.1", "Dim.2","Dim.3")] ### 3 PC axes
      hv_comm@HVList[[b]] <- hypervolume_gaussian(dat, name=as.character(plot_list[a]),
                                                  samples.per.point = 1000)
    }
    vols <- matrix(0,nrow = 25, ncol=2)
    for (j in 1:length(hv_comm@HVList)){
      vols[j,1] <- hv_comm@HVList[[j]]@Name
      vols[j,2] <- hv_comm@HVList[[j]]@Volume
    }
    vols <- as.data.frame(vols)
    colnames(vols) <- c("plot","volume")
    vols$volume <- as.numeric(as.character(vols$volume))
    volumes[[a]] <- vols
  }
  volcombi[[i]] <- do.call(rbind, volumes)
}

# saveRDS(rndid, "./R_output/rndid.rds")
# saveRDS(pc_out, "./R_output/FAB_pc_out.rds")
# saveRDS(volcombi, "./R_output/FAB_volcombi.rds")

### compile output
vols <- volcombi %>% 
  lapply(.%>% group_by(plot))%>%
  lapply(.%>%summarize(vol_m=mean(volume), vol_sd=sd(volume), 
         vol_m_log=mean(log(volume)), vol_sd_log=sd(log(volume))))

volss <- do.call(rbind,vols)

volsss <- volss %>% group_by(plot) %>%
  summarize(vol_mm=mean(vol_m), vol_ss=sd(vol_sd),
            vol_mm_log=mean(vol_m_log), vol_ss_log=mean(vol_sd_log))%>%
  rename(Plot=plot)

summary(volsss$vol_mm_log)
summary(volsss$vol_ss_log)


### community productivity ~ spectral space
dat <- read.csv("./FAB_comm_productivity.csv")
volx <- merge(volsss, dat, by="Plot")

### NBE, CE, SE does not make sense for monocultures
vols1 <- volx[-which(volx$des.S==1),]

### net biodiversity effect
plot(NBE_dST_sqrt ~ vol_mm_log, data=vols1)
mod <- lm(NBE_dST_sqrt ~ vol_mm_log, data=vols1)
summary(mod)

### complementarity effect 
plot(CE_dST_sqrt ~ vol_mm_log, data=vols1)
mod <- lm(CE_dST_sqrt ~ vol_mm_log, data=vols1)
summary(mod)

### selection effect
plot(SE_dST_sqrt ~ vol_mm_log, data=vols1)
mod <- lm(SE_dST_sqrt ~ vol_mm_log, data=vols1)
summary(mod)

### species richness
plot(des.S ~ vol_mm_log, data=volx)
mod <- lm(des.S ~ vol_mm_log, data=volx)
summary(mod)

#### END ############
