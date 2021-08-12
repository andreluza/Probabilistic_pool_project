
## load packages
source ("R/Packages.R")
## load functions
source ("R/Functions.R")

#################################################

#               LOCAL COMMUNITY DATA

#################################################

# -----------------------------------
# load data from Luza et al. (2019)
data<-read.csv (here ("data","community_data","AppendixS1- Small_mammal_data.csv"), h=T, sep=";")
data<-data [,1:40]

# interactions between factors to obtain sites
data$site_sub<- as.factor (interaction(data$REFERENCE,data$REGION, data$SITE))
data$hab_frag<-as.factor (interaction(data$HABITAT,data$FRAGMENTED))
data$hab_cut_edge<- as.factor (interaction(data$HABITAT,data$CLEAR_CUT,data$FOREST_EDGE, data$GRASSLAND_FOREST_EDGE))
data$hab_cut_edge<-droplevels (data$hab_cut_edge)

# removing sites without details on trap type and sampling effort
data<-subset (data, data$TRAP_TYPE != "NA", drop=T) ### not removed from database
data<- subset (data, data$EFFORT_PER_HABITAT !="NA") ### not removed from database

## filtering out  non-small mammal orders
data<-subset (data, subset = data$ORDER != "Xenarthra" & 
                data$ORDER != "Monotremata" & 
                data$ORDER != "Macroscelidea" & 
                data$ORDER != "Pholidota")

## filtering data by trap type
data<-data [which (data$BOX.LIKE ==1 | 
                     data$SNAP.LIKE ==1 | 
                     data$WIRE_MESHED == 1 | 
                     data$PITFALL == 1),]

# total number of sites
# length(unique (data$site_sub))
# and of studies 
# length(unique (data$REFERENCE))

## adjusting habitat levels
levels (data$hab_cut_edge) [levels (data$hab_cut_edge)== "edge.0.0.1"]<- "edge.0.0.0"

## list of non-identified spp
list_names <-unique(data$SPECIES) [grep ("\\_sp.", unique(data$SPECIES))] [-c(1,4,22,23,24,25,29)]
### removing non-identified species
data <- data [which (data$SPECIES %in% list_names ==F),]
# check if some one was not detected by grep 
unique (data$SPECIES)[grep ("sp",unique (data$SPECIES))]
list_names<-c(list_names,c("Apodemus_sp.","Neacomys_sp."))
### removing non-identified species, again
data <- data [which (data$SPECIES %in% list_names ==F),]
# no one remains
unique (data$SPECIES)[grep ("sp",unique (data$SPECIES))]

# -----------------------------------
# load PREDICTS' data
data_pred<-read.csv (here ("data","community_data",'hudson_data.csv'),h=T,sep=",")

# filtering by trapping type
data_pred2<-data_pred [data_pred$Sampling_method == "live_traps" | 
                        data_pred$Sampling_method == "baited_sherman_traps_along_transects" | 
                         data_pred$Sampling_method == "pit-fall_traps_with_drift_fences" | 
                         data_pred$Sampling_method == "multiple" 	| 
                         data_pred$Sampling_method == "baited_traps" | 
                         data_pred$Sampling_method == "various"   | 
                         data_pred$Sampling_method == "small_mammal_trapping" | 
                         data_pred$Sampling_method == "mixed_traps",]

# interaction to have sites
data_pred2<-cbind (data_pred2, hab_cut_edge= as.factor (interaction (data_pred2$Habitat,data_pred2$Clearcut,data_pred2$forest_edge, data_pred2$GF_edge)))
data_pred2$hab_cut_edge<-droplevels (data_pred2$hab_cut_edge)
data_pred2$site_sub<-as.factor (interaction (data_pred2$Reference,data_pred2$Site_name))
# number of sites
# length(unique (data_pred2$site_sub))
# and studies
# length(unique (data_pred2$Reference))

# Let's bind the two data sets
## select that interesting columns from Luza et al.  
luza_data <- data [, c(41, 17, 30,29, 31, 32, 43,39)]
#  colnames ( luza_data)

# matching PREDICTS and Luza
data_pred2<-data_pred2 [data_pred2$Presence >0,]

## occurrence matrix
predicts_data <- data_pred2 [, c(74,28,42,43,69,65,73,71)]
colnames(predicts_data)<-colnames(luza_data)

## rbind these datasets
matching<-rbind (luza_data, predicts_data)

## filtering data considering three values of sampling effort
matching_100<-subset (matching, subset= EFFORT_PER_HABITAT >= 100) # minimum of 100 trap-nights
matching_500<-subset (matching, subset= EFFORT_PER_HABITAT >= 500) # minimum of 500 trap-nights
matching_1000<-subset (matching, subset= EFFORT_PER_HABITAT >= 1000)  # minimum of 1000 trap-nights

# occurrence matrices for each sampling effort value
# 100 trap-nights
data_100_traps<-cast(matching_100, 
                     site_sub+hab_cut_edge+EFFORT_PER_HABITAT+LONG+LAT~SPECIES, 
                     fun.aggregate="sum", value="PRESENCE")
data_100_traps_var<-data_100_traps [,-c(1:5)]
data_100_traps_var [data_100_traps_var >= 2] <- 1
colnames (data_100_traps_var) <- gsub ("_",".", colnames (data_100_traps_var))

## 500 trap-nights
data_500_traps<-cast(matching_500, site_sub+hab_cut_edge+EFFORT_PER_HABITAT+LONG+LAT~SPECIES, fun.aggregate="sum", value="PRESENCE")
data_500_traps_var<-data_500_traps [,-c(1:5)]
data_500_traps_var [data_500_traps_var >= 2] <- 1
colnames (data_500_traps_var) <- gsub ("_",".", colnames (data_500_traps_var))

## 1000 trap-nights
data_1000_traps <-cast(matching_1000, site_sub+hab_cut_edge+EFFORT_PER_HABITAT+LONG+LAT~SPECIES, fun.aggregate="sum", value="PRESENCE")
data_1000_traps_var <-data_1000_traps [,-c(1:5)]
data_1000_traps_var [data_1000_traps_var >= 2] <- 1
colnames (data_1000_traps_var) <- gsub ("_",".", colnames (data_1000_traps_var))

# save these occurrence matrices
save (data_100_traps,
      data_500_traps,
      data_1000_traps,
      data_100_traps_var,
      data_500_traps_var,
      data_1000_traps_var,
      file=here ("data","community_data","community_matrices_effort.RData"))

# end