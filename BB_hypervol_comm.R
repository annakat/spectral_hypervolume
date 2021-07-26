##########################################################
### Calculate hypervolumes of plant communities in FAB ###
### Anna K Schweiger, March 2019, updates July 2021 ######

library(hypervolume)
library(FactoMineR)
library(tidyverse)

spec <- read.csv("./BB_ind_spec.csv")

### 50 random selections of 4 subplots per plot and 3 individuals per subplot
sp <- split(spec, list(spec$Plot))
sam <- list()

set.seed(1840)
for(i in 1:50){
  ### 4 random subplots per plot
  samples <- lapply(sp, function(x) x[x$Subplot %in% sample(unique(x$Subplot),4,replace = F),])
  out <- do.call(rbind, samples)
  row.names(out) <- NULL
  ### 3 random individuals per subplot 
  ssamples <- out %>% group_by(sub_ID) %>% sample_n(size=3)
  sam[[i]] <- as.data.frame(ssamples)
}

# saveRDS(sam, "./R_output/BB_sam.rds")
  
pc_out <- list()
volcombi <- list()

for (i in 1:50){
  sub <- sam[[i]]
  
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

# saveRDS(pc_out, "./R_output/BB_pc_out.rds")
# saveRDS(volcombi, "./R_output/BB_volcombi.rds")

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
dat <- read.csv("./BB_comm_productivity.csv")
volx <- merge(volsss, dat, by="Plot")

### biomass
plot(total_biom~vol_mm_log, data=volx)
mod <- lm(total_biom ~ vol_mm_log, data=volx)
summary(mod)

### species richness
plot(Sp ~vol_mm_log, data=volx)
mod <- lm(Sp ~vol_mm_log, data=volx)
summary(mod)

##### END #######


