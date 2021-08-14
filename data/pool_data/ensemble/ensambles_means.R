library(raster)
library(dplyr)

#rfs <-readRDS("Z:\\karger\\ch_wc_cru\\output\\results\\stacks\\regional\\mammals\\prob.stkch_ts_sp_rf.RDS")
#gams<-readRDS("Z:\\karger\\ch_wc_cru\\output\\results\\stacks\\regional\\mammals\\prob.stkch_ts_sp_gam.RDS")
#glms<-readRDS("Z:\\karger\\ch_wc_cru\\output\\results\\stacks\\regional\\mammals\\prob.stkch_ts_sp_glm.RDS")

setwd("C:/Users/topoa/OneDrive/capII/ensemble")

rfs <- readRDS ("prob.stkch_ts_sp_rf.RDS")
gams <- readRDS ("prob.stkch_ts_sp_gam.RDS")
glms <- readRDS ("prob.stkch_ts_sp_glm.RDS")

spec_all<-Reduce(intersect, list(names(rfs),names(glms),names(gams)))
rfs1 <-subset(rfs,spec_all)
gams1<-subset(gams,spec_all)
glms1<-subset(glms,spec_all)

ens_mean<-stack()
#ens_sd  <-stack()
for (n in 1:length(names(rfs1)))
{
  ens<-stack(glms1[[n]],gams1[[n]],rfs1[[n]])
  
  mm1<-calc(ens,mean)
  #sd1<-calc(ens,sd)
  ens_mean<-stack(ens_mean,mm1)
  #ens_sd<-stack(ens_sd,sd1)
  print(n)
}

speclist<-read.csv("Z:\\karger\\ch_wc_cru\\output\\results\\metadata\\species_names\\mammals.csv",sep=";")
speclist$ID_ch_wc_cru
joinspe<-as.data.frame(names(glms1))
colnames(joinspe)<-c("ID_ch_wc_cru")
joinspe$ID_ch_wc_cru<-as.numeric(gsub("id_","",joinspe$ID_ch_wc_cru))
specnames<-merge(joinspe,speclist,by.x="ID_ch_wc_cru",by.y="ID_ch_wc_cru")

names(ens_mean)<-specnames[,3]
#names(ens_sd)<-specnames

