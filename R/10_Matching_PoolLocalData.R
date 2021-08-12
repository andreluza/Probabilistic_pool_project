
## load packages
source ("R/Packages.R")
## load functions
source ("R/Functions.R")

# ---------------------------------------------
#               LOCAL POOL DATA
# ---------------------------------------------

load(here ("data","pool_data","probabilistic_pool.RData"))

# ---------------------------------------------
#               LOCAL COMMUNITY DATA
# ---------------------------------------------

load(here ("data","community_data","community_matrices_effort.RData"))

# list of occurrence matrices
# select sites with more than five species

list_matrices_comp <- lapply(list(data_100_traps_var, 
                                  data_500_traps_var,
                                  data_1000_traps_var), function (i) {
                	a <- i[which(rowSums (i) >=5),]
                	b <- a [,which(colSums(a) >0)]
                    ; # return
                	(b)
})

# how much each order contribute to data
# tab_freq <- matching_1000 [which(gsub ("_",".",matching_1000$SPECIES) %in% colnames(list_matrices_comp[[3]])),]
# tab_freq <- (cast (tab_freq,
#       SPECIES ~ORDER))[,-1]
# contribution of each order to the complete dataset
# (colSums (tab_freq>0))/nrow(tab_freq) *100

# ----------------------
# get coordinates of local sites

coords_100 <- data.frame(LAT=data_100_traps$LAT, LONG=data_100_traps$LONG)
coordinates(coords_100)<- ~ LONG + LAT
proj4string(coords_100) <- "+proj=longlat +datum=WGS84 +no_defs"

# extract pool probabilities for each site
# 100 trap nights
extract_100_traps <-lapply (probabilistic_pool, function (i) # for each pool 
                          extract (i, coords_100, # extract
                                   method='simple',
                                   fun=mean, small=T,df=T,na.rm=T) [,-1]) 

## 500 traps
#coords_500 <- data.frame(LAT=data_pred3$LAT, LONG=data_pred3$LONG)
#coordinates(coords_500)<- ~ LONG + LAT
#proj4string(coords_500)<- proj4string (ENV4)
### extrair
#traps_500<-lapply (smlista2, function (i) 
#	extract (i, coords_500, method='simple',fun=mean, small=T,df=T,na.rm=T) [,-1]) ## it works!!

## 1000 traps
# coords_1000<- data.frame(LAT=data_pred1000$LAT, LONG=data_pred1000$LONG)
# coordinates(coords_1000)<- ~ LONG + LAT
# proj4string(coords_1000)<- proj4string (ENV4)
# extrair 
# traps_1000 <-lapply (smlista2, function (i) 
#	extract (i, coords_1000, method='simple',fun=mean, small=T,df=T,na.rm=T) [,-1]) ## it works!!
## 
# save (traps_100,traps_500, traps_1000,file= "dados_extract.RData")

# fazer um mapa sobrepondo os pontos amostrados com 500 e 1000 traps (pontos vazados) com pontos com 100 traps (pontos solidos)
##
require(maps)
par(mfrow=c(1,1),family="serif")
png("mapa_test.png",width=18, height=15, units="cm", res=600, family="serif")
map("world", interior=F, col="gray70", fill=F, type="l", plot=T, lwd=1)
points(x= coordinates(coords_100)[,1], coordinates(coords_100)[,2], pch=1,col=as.numeric(data_pred100$hab_cut_edge)+7,
	cex=0.9)
points(x= coordinates(coords_500)[,1], coordinates(coords_500)[,2], pch=3,col=as.numeric(data_pred3$hab_cut_edge)+7,
	cex=0.9)
points(x= coordinates(coords_1000)[,1], coordinates(coords_1000)[,2], pch=19,col=as.numeric(data_pred1000$hab_cut_edge)+7,
	cex=0.9)
dev.off()
#text (-143.6677,10.74326, expression(paste(Delta, ""^i, Psi, "=", " ", ""^i, Psi[DxE],"-", Psi[OBS])))
#legend (-168.6677,5.74326, legend=c("1","3","6","9"), horiz=F,cex=0.8,bty = "n", bg=NULL, pch=19, col="black",
#	pt.cex= c(1*0.3,3*0.3,6*0.3, 9*0.3))


## carregar o pool extraido de cada sitio

setwd("C:/Users/topoa/OneDrive/capII")
load ("dados_extract.RData")

### colocar as especies que nao tem no pool mas tem nas comunidades
# site_list_100<-lapply (traps_100, function (i) 
#	cbind (data_pred100A [which (colnames (data_pred100A) %notin% colnames (traps_100[[1]]))], i))
# site_list<-lapply (traps_500, function (i) 
#	cbind (data_pred4 [which (colnames (data_pred4) %notin% colnames (traps_500[[1]]))], i))
# site_list_1000<-lapply (traps_1000, function (i) 
#	cbind (data_pred5 [which (colnames (data_pred5) %notin% colnames (traps_1000[[1]]))], i))

####
# save (site_list_100,site_list,site_list_1000, file="pool_probabilistico_add_especies.RData")

### abrir o RData com o pool probabilistico com especies adicionadas
load ("pool_probabilistico_add_especies.RData")

## ranks

list_spp <-  gsub ("\\."," ",colnames (site_list_1000[[1]]))
x <-  classification(list_spp,db = 'ncbi')

# rm not found
x <-x[which(unlist(lapply (x, length)) > 1)]

table(unlist(lapply (x, function (i) 
  i[which(i$rank == "order"),"name"]
)))/length(x)

### abrir os dados de comunidades (matriz de deteccao e nao deteccao das spp. em um dado sitio

load ("matriz_comunidades.RData")

######
## SEM NECESSIDADE DE UM MODELO NULO PARA RIQUEZA

## 100 traps 
#local_100 <- rowSums (data_pred100A)
#diff_riqueza_100 <- lapply (site_list_100, function (i) 
#	rowSums (i) - local_100)## diferenca entre o observado e esperado

## 500 traps
#local <- rowSums (data_pred4)
#diff_riqueza_500 <- lapply (site_list, function (i) 
#	rowSums (i) - local)## diferenca entre o observado e esperado

### 1000 traps
#local_1000 <- rowSums (data_pred5)
#diff_riqueza_1000 <- lapply (site_list_1000, function (i) 
#	rowSums (i) - local_1000)## diferenca entre o observado e esperado

# save (diff_riqueza_100,diff_riqueza_500,diff_riqueza_1000,
#	file="RESULTADO_diferenca_riqueza.RData")
####################################### 

#### FUNCTIONAL
#library("hypervolume")
#https://benjaminblonder.org/hypervolume_faq.html
#setwd ("C:/Users/topoa/OneDrive/capII/trait_data")
#penone<-read.csv("Penone_et_al_2016_mammal_trait_data_imputed.csv",h=T,sep=",")
### ELTONTRAITS
#elton<-read.table ("elton1.0_my_mod.csv",h=T, sep=",")
#setwd("C:/Users/User/Documents/AndreLuza/capII/data")

#penone1 <- penone [order (penone$IUCN.binomial, decreasing=F),]
#penone1$IUCN.binomial<-droplevels (penone1$IUCN.binomial)
#penone1 <-subset (penone1, subset = penone1$Body.mass.g <= 5000)
#elton1 <- elton [order (elton$binomial, decreasing=F),]
#elton1$binomial<-droplevels (elton1$binomial)
#elton1 <-subset (elton1, subset = elton1$BodyMass_Value <= 5000)

#penone2<-penone1 [which (penone1$IUCN.binomial %in% elton1$binomial),]
#penone2$IUCN.binomial<-droplevels (penone2$IUCN.binomial)
#elton2 <- elton1 [which (elton1$binomial %in% penone2$IUCN.binomial),]
#elton2$binomial<-droplevels (elton2$binomial)

#pairs(data.frame (elton2$ForStrat_Value,elton2$Activity_Nocturnal,
#	elton2$Diet_Inv, elton2$Diet_PlantO, elton2$Diet_Fruit, elton2$Diet_Seed, penone2$LitPerYear,
#	penone2$LitSz), 
#	panel=panel.smooth, lwd=2, pch=19)

#range(
#	cor(data.frame (elton2$ForStrat_Value,elton2$Activity_Nocturnal,
#		elton2$Diet_Inv, elton2$Diet_PlantO, elton2$Diet_Fruit, elton2$Diet_Seed,penone2$LitPerYear,penone2$LitSz)[,-1]))

#traits_complete <- cbind (penone2, elton2)

#traits <- subset (traits_complete, 
#	gsub ("_",".",traits_complete$IUCN.binomial) %in% colnames(site_list[[1]]))
#traits$ForStrat_Value<-as.numeric (traits$ForStrat_Value)
#rownames (traits)<-traits$binomial

#traits <- traits [, c(3,4,6,9,10,13)]

#cor(traits)
#apply (apply (traits, 2, duplicated),2,table)

#### procedimento de imputing dos dados faltantes
#### pegar as especies faltantes e estimar

#sp_para_imputar <- colnames(site_list[[1]]) [which(colnames(site_list[[1]]) %in% gsub ("_",".",traits_complete$IUCN.binomial) == F)]
#sp_para_imputar %in%  gsub ("_",".",traits_complete$IUCN.binomial)

## media das sp do mesmo genero
## botar os atributos de Dasyurus_albopunctatus em maculatus

## "Dasycercus.blythi"
# d_blythi <- traits_complete [grep("Dasycercus", traits_complete$IUCN.binomial),c(3,4,6,9,10,13)]
# d_blythi<- colMeans (d_blythi)

##  "Dasyurus.maculatus" 
# d_maculatus <- traits_complete [grep("Dasyurus", traits_complete$IUCN.binomial),c(3,4,6,9,10,13)]
# d_maculatus <- colMeans (d_maculatus)

## Fukomys nao ha
## https://www.researchgate.net/publication/226224782_African_Mole-rats_Bathyergidae_A_Complex_Radiation_in_Tropical_Soils/figures?lo=1
# f_boca <- traits_complete [grep("Cryptomys", traits_complete$IUCN.binomial),c(3,4,6,9,10,13)]
# f_boca <- colMeans (f_boca)

## hystrix sumatrae em Hystrix.crassispinis
# h_criss <- traits_complete [grep("Hystrix", traits_complete$IUCN.binomial),c(3,4,6,9,10,13)]
# h_criss <-  colMeans (h_criss)

## flaviventris em marmota caligata
# m_caligata <- traits_complete [grep("Marmota", traits_complete$IUCN.binomial),c(3,4,6,9,10,13)]
# m_caligata <-  colMeans (m_caligata)

## maxomys rajah em m taju
## https://academic.oup.com/jmammal/article/94/6/1412/907529
# m_taju <- traits_complete [grep("Maxomys", traits_complete$IUCN.binomial),c(3,4,6,9,10,13)]
# m_taju <- colMeans(m_taju)

## microtus  no lugar de hyperboleus
## parece ocorrer numa regiao similar
# m_hyper <- traits_complete [grep("Microtus", traits_complete$IUCN.binomial),c(3,4,6,9,10,13)]
# m_hyper <-  colMeans(m_hyper)

## monodelphis kunsi no lugar de peruviana
## de acordo com proximidade filogenetic
# m_peru <- traits_complete [grep("Monodelphis", traits_complete$IUCN.binomial),c(3,4,6,9,10,13)]
# m_peru <- colMeans(m_peru)

## nanospalax - nao ha
## https://www.researchgate.net/publication/272830161_How_can_scientific_researches_change_conservation_priorities_A_review_of_decade-long_research_on_blind_mole-rats_Rodentia_Spalacinae_in_the_Carpathian_Basin/figures?lo=1
# n_ehre <- traits_complete [grep("Spalax", traits_complete$IUCN.binomial),c(3,4,6,9,10,13)]
# n_ehre <- colMeans(n_ehre)

##ochotona
# o_mant <- traits_complete [grep("Ochotona", traits_complete$IUCN.binomial),c(3,4,6,9,10,13)]
# o_mant <- colMeans(o_mant)

## pronolagus
# p_saun <- traits_complete [grep("Pronolagus", traits_complete$IUCN.binomial),c(3,4,6,9,10,13)]
# p_saun <- colMeans(p_saun)

# Spilocuscus.maculatus - nao ha
## https://academic.oup.com/jmammal/article/85/5/825/858786
## https://www.researchgate.net/publication/302954153_Phylogenetic_Relationship_of_Cuscuses_Marsupialia_Phalangeridae_from_Papua_and_Maluku_Based_on_Mitochondrial_Sequences_of_NADH_Dehydrogenase_Sub-unit_1_Gene
# s_macu <- traits_complete [grep("Phalanger", traits_complete$IUCN.binomial),c(3,4,6,9,10,13)]
# s_macu <- colMeans(s_macu )

# Tamiops
# t_macc <- traits_complete [grep("Tamiops", traits_complete$IUCN.binomial),c(3,4,6,9,10,13)]
# t_macc <- colMeans(t_macc)

# Thylamys
# t_pulc <- traits_complete [grep("Thylamys", traits_complete$IUCN.binomial),c(3,4,6,9,10,13)]
# t_pulc <- colMeans(t_pulc)

#imputar_estas<-rbind(d_blythi, d_maculatus,f_boca,f_boca,f_boca,f_boca,
#				h_criss,m_caligata,m_taju,
#				m_hyper,m_peru,n_ehre,n_ehre,o_mant,p_saun,s_macu,t_macc,t_pulc)

#rownames(imputar_estas) <- c(
#			"Dasycercus.blythi",
#			"Dasyurus.maculatus" ,
#			"Fukomys.bocagei",
#			"Fukomys.damarensis",
#			"Fukomys.mechowi",
#			"Fukomys.ochraceocinereus",
#			"Hystrix.crassispinis",
#			"Marmota.caligata",
#			"Maxomys.tajuddinii",
#			"Microtus.hyperboreus",
#			"Monodelphis.peruviana",
#			"Nannospalax.ehrenbergi",
#			"Nannospalax.xanthodon",
#			"Ochotona.mantchurica",
#			"Pronolagus.saundersiae" ,
#			"Spilocuscus.maculatus", 
#			"Tamiops.macclellandii",
 # 			"Thylamys.pulchellus")

#traits_imput <- rbind(traits,imputar_estas)
#traits <- decostand(traits_imput, "standardize")
#rownames(traits) <- gsub ("_",".",rownames(traits))
#colnames(site_list[[1]]) [which(colnames(site_list[[1]]) %in% gsub ("_",".",rownames(traits)) == F)]

#### 
# save (traits, file="imputed_trait_media.RData")
setwd("C:/Users/topoa/OneDrive/capII")
load ("imputed_trait_media.RData")

### functional diversity of the pool ( FD for each cell)

### HERE'S YOUR SUGGESTION
niterations <- 50
LPSP_composition <- values (smlista2[[1]])
LPSP_composition_high0 <- LPSP_composition[which(round (rowSums (LPSP_composition,na.rm=T))== 5),]

pool_sample <- replicate (niterations , 
	lapply (seq(1,dim(LPSP_composition_high0)[1]), function (cell) ## for each cell
	sample (LPSP_composition_high0 [cell,], size= round(sum(LPSP_composition_high0 [cell,],na.rm=T)), 
	prob= LPSP_composition_high0 [cell,])), simplify = "array" )

## pool_sample [1,1] ## [species specific pool 1, sample 1]

traits_pool <- traits[which(rownames (traits) %in% colnames (LPSP_composition_high0)),]

## took the traits of each species in the pool
#cl <- makeCluster(6) ## number of cores = generally ncores -1
#clusterEvalQ(cl, library(hypervolume))
# export your data and function
#clusterExport(cl, c("pool_sample","traits"))

#traits_pool <- parLapply (cl,seq(1,nrow(pool_sample)), function (pool)
traits_pool <- lapply (seq(1,nrow(pool_sample)), function (cell)
			lapply (seq(1,ncol(pool_sample)), function (iteration)
			traits [which (rownames(traits) %in% names(pool_sample [cell,iteration][[1]])),]))
#stopCluster(cl)

# observed pool hypervolume 
## 10 ITERATIONS TOOK 2.5 DAYS TO RUN -USING 3 OF 4 CORES OF MY PC
cl <- makeCluster(5) ## number of cores = generally ncores -1
clusterEvalQ(cl, library(hypervolume))
# export your data and function
clusterExport(cl, c("traits_pool"))

hypervolume_pool_obs <- parLapply(cl,traits_pool, function (cell)# for each sample
	lapply(cell, function (iteration)# for each community ##   
      	hypervolume_gaussian (iteration,
		kde.bandwidth=estimate_bandwidth(iteration,method="silverman"))@Volume))

stopCluster(cl)

## mean pool hypervolume (averaged across n interations)

mean_hypervolume_pool_5sp <-lapply (hypervolume_pool_obs, function (cell) # for each cell
	# calculate the mean pool hypervolume for each community (averaged across iterations)
	rowMeans (matrix (unlist(lapply (cell, unlist)),ncol=niterations,byrow=T))) 

save(hypervolume_pool_obs ,mean_hypervolume_pool_5sp,file="mean_hypervolume_pool_5sp.RData")

## abrir os resultdos do hypervolume do pool

load("mean_hypervolume_pool_3_10.RData")
load("mean_hypervolume_pool_first.RData")
load("mean_hypervolume_pool_scd.RData")
load("mean_hypervolume_pool_mais_5.RData")
load("mean_hypervolume_pool_five.RData")
load("mean_hypervolume_pool_5sp.RData")

FD_POOL <- rowMeans (cbind(
	unlist(mean_hypervolume_pool_first),
	unlist(mean_hypervolume_pool_scd),
	unlist(mean_hypervolume_pool),
	unlist (mean_hypervolume_pool_mais_5,
	unlist(mean_hypervolume_pool_five))))

all_cells <- data.frame (
	m5 = round (rowSums (LPSP_composition,na.rm=T))>5)
all_cells$igual5 <- round (rowSums (LPSP_composition,na.rm=T))== 5

all_cells$hyper <- ifelse (all_cells$m5 == FALSE,
		0,TRUE)

all_cells$hyper[which(round (rowSums (LPSP_composition,na.rm=T))>5)] <- FD_POOL
all_cells$hyper[which(round (rowSums (LPSP_composition,na.rm=T))== 5)] <- unlist(mean_hypervolume_pool_5sp)


## conferir
which(all_cells$hyper > 0) == which (round (rowSums (LPSP_composition,na.rm=T))>5)

cp_pool <- sum(smlista2[[1]])
values (cp_pool ) <- all_cells$hyper 

writeRaster(cp_pool,filename= "FD_pool_final", bylayer=TRUE, format="GTiff")

#################################
### Analises para comunidades

########################################
########## FRIC, TODAS AS ESPECIES #####
########################################

########## null model
# pool = matrix with pool probabilities
# local = a vector with local richness (rowSums of the local matrix)

nsamples <- 100

sampling_prevalence_NM <- function (pool, local,nsamples){## 
	sampling <- lapply (as.list(seq(1, nrow (pool))), function (k) ## over all sites
	lapply (as.list(rep(1,nsamples)), function (i) ## replicate 1 x niter
	replicate(i,sample (x=pool[k,],size=local[k],prob=pool[k,], replace=F)))) ## do the prevalence sampling
	# obtain only the sampled probabilities
	result_list <- vector ("list")
	result_list$probabilities <- lapply (sampling, function (i) do.call (cbind, i))
	# obtain the random list of sampled species
	result_list$species <- lapply (sampling, function (k) do.call (cbind, lapply (k, function (i) dimnames(i)[[2]])))
	# return results
	return (result_list)
}


#########################
# NULL MODEL AND HYPERVOLUME ANALYSIS
##########################

#100 traps

local_100 <- rowSums (data_pred100A)
nsamples<- 1000 ## a few iterations to test the function

### run the null model to sort species from the pool (the same number of specis as observed in the local communities)

require(parallel)
cl <- makeCluster(4) ## number of cores = generally ncores -1
clusterEvalQ(cl, library(phytools))
clusterEvalQ(cl, library(phylobase))
clusterEvalQ(cl, library(dplyr))

# export your data and function
clusterExport(cl, c("local_100","nsamples", "site_list_100", "sampling_prevalence_NM"))

nm_teste <- parLapply(cl,site_list_100, function (i) 
	sampling_prevalence_NM (i,local=local_100,nsamples=nsamples))## sample with the size of the obs community

stopCluster(cl)

#### took the traits of the sampled species
#traits
list_local_traits <- lapply (as.list(seq (1, length (nm_teste))), function (j)## for each pool
	lapply (as.list (seq (1,length (local_100))), function (k) ### site
	lapply (as.list (seq(1,nsamples)), function (i) ## and sample
	traits[which (rownames(traits) %in% nm_teste[[j]]$species[[k]] [,i]),]))) ## took the trait

### removing sites with less than 4 species
list_local_traits <- lapply (as.list(seq (1, length (nm_teste))), function (j)
	list_local_traits [[j]] [which (local_100 >= 4)])
	
## NULL FUNCTIONAL COMPLETENESS

library("hypervolume")

### calculate null hypervolume for that random communities
#cl <- makeCluster(4) ## number of cores = generally ncores -1
#clusterEvalQ(cl, library(hypervolume))
#export your data and function
#clusterExport(cl, c("nsamples","nm_teste", "list_local_traits"))

#hypervolume_local_null <-parLapply (cl, as.list(seq (1, length (nm_teste))), function (j)# each pool
#	lapply (as.list(seq(1,length (list_local_traits [[1]]))), function (k)# site
#	unlist (lapply (as.list (seq(1,nsamples)), function (i) ## sample
#	hypervolume_gaussian(list_local_traits [[j]][[k]][[i]],
#	kde.bandwidth=estimate_bandwidth(list_local_traits [[j]][[k]][[i]],method="silverman"))@Volume)))) # calculate hyp and extract the object "Volume"
#stopCluster (cl)

teste <- lapply (list_local_traits, function (i) 
		lapply (i, function (k)
			lapply (k, function (j) 
				hypervolume_gaussian (j,
					kde.bandwidth=estimate_bandwidth(j,
					method="silverman"))@Volume)))

#save(hypervolume_local_null , file="hypervolume_local_null.RData")

# OBSERVED FUNCTIONAL COMPLETENESS

### data for observed communities
local_obs <- data.matrix(data_pred100A [which (local_100 >= 4),order (colnames(data_pred100A), decreasing=F)]) # removing sites with less than 4 sp
local_obs <- local_obs [,-which (colSums (local_obs)==0)] # removing absent species
list_local_obs <- apply (local_obs, 1, list) # it's faster to work with lists

## remove absent species
list_local_obs <- lapply (as.list (seq(1,length (list_local_obs))), function (i) 
	list_local_obs [[i]] [[1]] [which (list_local_obs [[i]] [[1]] == 1)])

### took the traits of species occurring in local communities
list_local_traits_obs <- lapply (as.list (seq(1,length (list_local_obs))), function (i)
	traits[which (rownames(traits) %in% names (list_local_obs [[i]])),])

## calculate observed hypervolume

cl <- makeCluster(6) ## number of cores = generally ncores -1
clusterEvalQ(cl, library(hypervolume))
#export your data and function
clusterExport(cl, c("list_local_traits_obs"))

hypervolume_local_obs<- parLapply (cl,list_local_traits_obs, function (i)
	hypervolume_gaussian(i,	kde.bandwidth=estimate_bandwidth(i,method="silverman"))@Volume) # calculate and extract the object "Volume"

stopCluster (cl)

cl <- makeCluster(6) ## number of cores = generally ncores -1
clusterEvalQ(cl, library(hypervolume))
#export your data and function
clusterExport(cl, c("list_local_traits_obs"))

bandwidths_local <- parLapply (cl,list_local_traits_obs, function (i)
	hypervolume_gaussian(i,	
	kde.bandwidth=estimate_bandwidth(i,method="silverman"))@Parameters) # calculate bandwidths

stopCluster (cl)

### values of bandwidth
bandwidth_kde_local<-lapply (bandwidths_local, function (i)i$kde.bandwidth)
colMeans (do.call (rbind.data.frame, bandwidth_kde_local))
apply (do.call (rbind.data.frame, bandwidth_kde_local),2,sd)
sd_count_local<-lapply (bandwidths_local, function (i)i$sd.count)
## samples per point
sample_per_point<-lapply (bandwidths_local, function (i)i$samples.per.point)
mean (unlist(sample_per_point))
sd(unlist(sample_per_point))

## calculate hypervolume for the pool
pool_obs <- lapply (site_list_100, function (i) i [which (local_100 >=4),order(colnames(i), decreasing=F)])
traits_pool<-traits[which(rownames(traits) %in% colnames (pool_obs[[1]])),]
pool_obs <- lapply (pool_obs, function (i) i[,which (colnames(i) %in% rownames(traits_pool))])
pool_obs <- lapply (pool_obs, function (i) apply (i, 1, list))

### HERE'S THE THRESHOLD WE DISCUSS
#pool_obs <- lapply (seq (1, length (pool_obs)), function (j)
#	lapply (seq(1, length (pool_obs[[1]])), function (k)
#	pool_obs [[j]][[k]][[1]][which (pool_obs [[j]][[k]][[1]] >0.1)])) ## 
#
### HERE'S YOUR SUGGESTION
niterations <- 100
pool_sample <- replicate (niterations , lapply (seq(1,length (pool_obs)), function (pool) ## for each pool
	lapply (seq(1,length (pool_obs[[1]])), function (comm) ## for each community
	sample (pool_obs[[pool]][[comm]][[1]], size=round(sum(pool_obs[[pool]][[comm]][[1]])), 
	prob=(pool_obs[[pool]][[comm]][[1]])))), simplify = "array" )

## pool_sample [1,1] ## [species specific pool 1, sample 1]

## took the traits of each species in the pool
#cl <- makeCluster(6) ## number of cores = generally ncores -1
#clusterEvalQ(cl, library(hypervolume))
# export your data and function
#clusterExport(cl, c("pool_sample","traits"))

#traits_pool <- parLapply (cl,seq(1,nrow(pool_sample)), function (pool)
traits_pool <- lapply (seq(1,nrow(pool_sample)), function (pool)
			lapply (seq(1,ncol(pool_sample)), function (sample)
			lapply (seq(1,length (pool_sample[1,1][[1]])), function (comm)
				traits [which (rownames(traits) %in% names(pool_sample [pool,sample][[1]][[comm]])),])))
#stopCluster(cl)

# observed pool hypervolume 
## 10 ITERATIONS TOOK 2.5 DAYS TO RUN -USING 3 OF 4 CORES OF MY PC
#cl <- makeCluster(6) ## number of cores = generally ncores -1
#clusterEvalQ(cl, library(hypervolume))
# export your data and function
#clusterExport(cl, c("traits_pool"))

#hypervolume_pool_obs <- parLapply(cl, as.list(seq (1, length (traits_pool))), function (pool)# for #each pool
hypervolume_pool_obs <- lapply(traits_pool, function (pool)# for each pool
	lapply (pool, function (sample)# for each sample
	lapply (sample, function (comm)# for each community ##   
      	hypervolume_gaussian (comm,
		kde.bandwidth=estimate_bandwidth(comm,method="silverman"))@Volume)))

#stopCluster(cl)

## mean pool hypervolume (averaged across n interations)

mean_hypervolume_pool<-lapply (hypervolume_pool_obs, function (pool) # for each pool
	# calculate the mean pool hypervolume for each community (averaged across iterations)
	rowMeans (matrix (unlist(lapply (pool, unlist)),ncol=niterations,byrow=T))) 

save (hypervolume_pool_obs, file="hypervolumePOOL_100traps_observed.RData")
save (mean_hypervolume_pool, file="mean_hypervolume_pool_100traps_observed.RData")

# bandwidths used to calculate the pool hypervolume
cl <- makeCluster(6) ## number of cores = generally ncores -1
clusterEvalQ(cl, library(hypervolume))
# export your data and function
clusterExport(cl, c("traits_pool"))

bandwidths_pool_obs <- parLapply(cl, as.list(seq (1, length (traits_pool))), function (pool)# for each pool
	lapply (seq(1, length (traits_pool[[1]])), function (sample)# for each sample
	lapply (seq(1, 2), function (comm)# for each community ## length (traits_pool[[1]][[1]])
      	hypervolume_gaussian (traits_pool[[pool]][[sample]][[comm]],
			kde.bandwidth=estimate_bandwidth(traits_pool[[pool]][[sample]][[comm]],method="silverman"))@Parameters)))

stopCluster(cl)

# Other parameters important to report
bandwidth_kde_pool<-lapply (bandwidths_pool_obs, function (pool)
	lapply (seq (1, length (bandwidths_pool_obs [[1]])), function (sample) 
	lapply (seq (1, length (bandwidths_pool_obs [[1]][[1]])), function (comm) 
	pool[[sample]][[comm]]$kde.bandwidth)))

matrix_bandwidth <- lapply (bandwidth_kde_pool, function (pool) # for each pool
	# calculate the mean pool hypervolume for each community (averaged across iterations)
	matrix (unlist(lapply (pool, unlist)),ncol=niterations,byrow=T))

matrix_bandwidth <- lapply (lapply (bandwidth_kde_pool [[1]], unlist), matrix,ncol=6,nrow=niterations,byrow=T)

# mean bandwidths over samples
lapply(matrix_bandwidth,colMeans) 
lapply(matrix_bandwidth,function (i) apply(i,2,sd)) 

## SD - nice value would be 3
sd_count_pool <- lapply (bandwidths_pool_obs, function (pool)
	lapply (seq (1, length (bandwidths_pool_obs [[1]])), function (sample) 
	lapply (seq (1, length (bandwidths_pool_obs [[1]][[1]])), function (comm) 
	pool[[sample]][[comm]]$sd.count)))

## samples (shoots) per point
sample_per_point_pool<- lapply (bandwidths_pool_obs, function (pool)
	lapply (seq (1, length (bandwidths_pool_obs [[1]])), function (sample) 
	lapply (seq (1, length (bandwidths_pool_obs [[1]][[1]])), function (comm) 
	pool[[sample]][[comm]]$samples.per.point)))

## mean and sd per pool
lapply(sample_per_point_pool, function (i) mean(unlist(i)))
lapply(sample_per_point_pool, function (i) sd(unlist(i)))

#save (bandwidths_pool, bandwidths_local, file="bandwidth_100.RData")


#######################
######    500 traps  ##
#######################

local_500 <- rowSums (data_pred4) # local richness (sample size)

##########
# run null model
require(parallel)
cl <- makeCluster(4) ## number of cores = generally ncores -1
clusterEvalQ(cl, library(phytools))
clusterEvalQ(cl, library(phylobase))
clusterEvalQ(cl, library(dplyr))
# export your data and function
clusterExport(cl, c("local_500","nsamples", "site_list", "sampling_prevalence_NM"))

nm_teste <- parLapply(cl,site_list, function (i) sampling_prevalence_NM (i,local=local_500,nsamples=nsamples))## sample with the size of the obs community

stopCluster(cl)

#### took the traits of the sampled species
#traits
list_local_traits <- lapply (as.list(seq (1, length (nm_teste))), function (j)# for each pool
	lapply (as.list (seq (1,length (local_500))), function (k) ### for each community
	lapply (as.list (seq(1,nsamples)), function (i) # for each iteration
	traits[which (rownames(traits) %in% nm_teste[[j]]$species[[k]] [,i]),]))) ## took the traits

### removing the sites with less than 4 species
list_local_traits <- lapply (as.list(seq (1, length (nm_teste))), function (j)
	list_local_traits [[j]] [which (local_500 >= 4)])
	
### apply to all communities and iterations
#require(parallel)
#cl <- makeCluster(6) ## number of cores = generally ncores -1
#clusterEvalQ(cl, library(hypervolume))
#export your data and function
#clusterExport(cl, c("nsamples","nm_teste", "list_local_traits"))

hypervolume_local_null_500traps <- lapply (list_local_traits, function (i) 
		lapply (i, function (k)
			lapply (k, function (j) 
				hypervolume_gaussian (j,
					kde.bandwidth=estimate_bandwidth(j,
					method="silverman"))@Volume)))

#stopCluster (cl)

save(hypervolume_local_null_500traps , file="hypervolume_local_null_500traps.RData")

# OBSERVED
### data of local species composition
local_obs <- data.matrix(data_pred4 [which (local_500 >= 4),order (colnames(data_pred4), decreasing=F)])
local_obs <- local_obs [,-which (colSums (local_obs)==0)]
list_local_obs <- apply (local_obs, 1, list)

### removing absent species
list_local_obs <- lapply (as.list (seq(1,length (list_local_obs))), function (i) 
	list_local_obs [[i]] [[1]] [which (list_local_obs [[i]] [[1]] == 1)])

### took the traits of the species of each communities
list_local_traits_obs <- lapply (as.list (seq(1,length (list_local_obs))), function (i)
	traits[which (rownames(traits) %in% names (list_local_obs [[i]])),])

## calculate local community hypervolume
cl <- makeCluster(6) ## number of cores = generally ncores -1
clusterEvalQ(cl, library(hypervolume))
#export your data and function
clusterExport(cl, c("list_local_traits_obs"))

hypervolume_local_obs<- parLapply (cl,list_local_traits_obs, function (i)
	hypervolume_gaussian(i,	kde.bandwidth=estimate_bandwidth(i,method="silverman"))@Volume)

stopCluster (cl)

## bandwidth
cl <- makeCluster(6) ## number of cores = generally ncores -1
clusterEvalQ(cl, library(hypervolume))
#export your data and function
clusterExport(cl, c("list_local_traits_obs"))

bandwidths_local_500 <- parLapply (cl,list_local_traits_obs, function (i)
	hypervolume_gaussian(i,	
	kde.bandwidth=estimate_bandwidth(i,method="silverman"))@Parameters)

stopCluster (cl)

### bandwidth values (important to report)
bandwidth_kde_local<-lapply (bandwidths_local_500, function (i)i$kde.bandwidth)
colMeans (do.call (rbind.data.frame, bandwidth_kde_local))
apply (do.call (rbind.data.frame, bandwidth_kde_local),2,sd)
sd_count_local<-lapply (bandwidths_local_500, function (i)i$sd.count)

sample_per_point<-lapply (bandwidths_local_500, function (i)i$samples.per.point)
mean (unlist(sample_per_point))
sd(unlist(sample_per_point))

## pool composition data
pool_obs <- lapply (site_list, function (i) i [which (local_500 >=4),order(colnames(i), decreasing=F)])
traits_pool<-traits[which(rownames(traits) %in% colnames (pool_obs[[1]])),]
pool_obs <- lapply (pool_obs, function (i) i[,which (colnames(i) %in% rownames(traits_pool))])
pool_obs <- lapply (pool_obs, function (i) apply (i, 1, list))

### HERE'S THE THRESHOLD WE DISCUSS
#pool_obs <- lapply (seq (1, length (pool_obs)), function (j)
#	lapply (seq(1, length (pool_obs[[1]])), function (k)
#	pool_obs [[j]][[k]][[1]][which (pool_obs [[j]][[k]][[1]] >0.1)])) ## 
#
### HERE'S YOUR SUGGESTION
#niterations <- 2
pool_sample <- replicate (niterations , lapply (seq(1,length (pool_obs)), function (pool) ## for each pool
	lapply (seq(1,length (pool_obs[[1]])), function (comm) ## for each community
	sample (pool_obs[[pool]][[comm]][[1]], size=round(sum(pool_obs[[pool]][[comm]][[1]])), 
	prob=(pool_obs[[pool]][[comm]][[1]])))), simplify = "array" )

## pool_sample [1,1] ## [species specific pool 1, sample 1]

## took the traits of each species in the pool
cl <- makeCluster(6) ## number of cores = generally ncores -1
clusterEvalQ(cl, library(hypervolume))
# export your data and function
clusterExport(cl, c("pool_sample","traits"))

traits_pool <- parLapply (cl,seq(1,nrow(pool_sample)), function (pool)
			lapply (seq(1,ncol(pool_sample)), function (sample)
			lapply (seq(1,length (pool_sample[1,1][[1]])), function (comm)
				traits [which (rownames(traits) %in% names(pool_sample [pool,sample][[1]][[comm]])),])))

stopCluster(cl)

# observed pool hypervolume 
## 10 ITERATIONS TOOK 2.5 DAYS TO RUN -USING 3 OF 4 CORES OF MY PC
cl <- makeCluster(6) ## number of cores = generally ncores -1
clusterEvalQ(cl, library(hypervolume))
# export your data and function
clusterExport(cl, c("traits_pool"))

hypervolume_pool_obs <- parLapply(cl, as.list(seq (1, length (traits_pool))), function (pool)# for each pool
	lapply (seq(1, length (traits_pool[[1]])), function (sample)# for each sample
	lapply (seq(1,length (traits_pool[[1]][[1]])), function (comm)# for each community ##   
      	hypervolume_gaussian (traits_pool[[pool]][[sample]][[comm]],
		kde.bandwidth=estimate_bandwidth(traits_pool[[pool]][[sample]][[comm]],method="silverman"))@Volume)))

stopCluster(cl)

## mean pool hypervolume (averaged across n interations)

mean_hypervolume_pool<-lapply (hypervolume_pool_obs, function (pool) # for each pool
	# calculate the mean pool hypervolume for each community (averaged across iterations)
	rowMeans (matrix (unlist(lapply (pool, unlist)),ncol=niterations,byrow=T))) 

#save (hypervolume_local_obs, hypervolume_pool_obs, file="hypervolume_500traps_observed.RData")
save (hypervolume_pool_obs, file="hypervolumePOOL_500traps_observed.RData")
save (mean_hypervolume_pool, file="mean_hypervolume_pool_500traps_observed.RData")


# bandwidths we used to calculate the pool hypervolume
cl <- makeCluster(6) ## number of cores = generally ncores -1
clusterEvalQ(cl, library(hypervolume))
# export your data and function
clusterExport(cl, c("traits_pool"))

bandwidths_pool_obs <- parLapply(cl, as.list(seq (1, length (traits_pool))), function (pool)# for each pool
	lapply (seq(1, length (traits_pool[[1]])), function (sample)# for each sample
	lapply (seq(1, 2), function (comm)# for each community ## length (traits_pool[[1]][[1]])
      	hypervolume_gaussian (traits_pool[[pool]][[sample]][[comm]],
			kde.bandwidth=estimate_bandwidth(traits_pool[[pool]][[sample]][[comm]],method="silverman"))@Parameters)))

stopCluster(cl)

# Other parameters important to report
bandwidth_kde_pool<-lapply (bandwidths_pool_obs, function (pool)
	lapply (seq (1, length (bandwidths_pool_obs [[1]])), function (sample) 
	lapply (seq (1, length (bandwidths_pool_obs [[1]][[1]])), function (comm) 
	pool[[sample]][[comm]]$kde.bandwidth)))

matrix_bandwidth <- lapply (bandwidth_kde_pool, function (pool) # for each pool
	# calculate the mean pool hypervolume for each community (averaged across iterations)
	matrix (unlist(lapply (pool, unlist)),ncol=niterations,byrow=T))

matrix_bandwidth <- lapply (lapply (bandwidth_kde_pool [[1]], unlist), matrix,ncol=6,nrow=niterations,byrow=T)

# mean bandwidths over samples
lapply(matrix_bandwidth,colMeans) 
lapply(matrix_bandwidth,function (i) apply(i,2,sd)) 

## SD - nice value would be 3
sd_count_pool <- lapply (bandwidths_pool_obs, function (pool)
	lapply (seq (1, length (bandwidths_pool_obs [[1]])), function (sample) 
	lapply (seq (1, length (bandwidths_pool_obs [[1]][[1]])), function (comm) 
	pool[[sample]][[comm]]$sd.count)))

## samples (shoots) per point
sample_per_point_pool<- lapply (bandwidths_pool_obs, function (pool)
	lapply (seq (1, length (bandwidths_pool_obs [[1]])), function (sample) 
	lapply (seq (1, length (bandwidths_pool_obs [[1]][[1]])), function (comm) 
	pool[[sample]][[comm]]$samples.per.point)))

## mean and sd per pool
lapply(sample_per_point_pool, function (i) mean(unlist(i)))
lapply(sample_per_point_pool, function (i) sd(unlist(i)))

#save (bandwidths_pool_obs, bandwidths_local_500, file="bandwidth_500.RData")

######################
####  1,000 traps  ###
######################

local_1000 <- rowSums (data_pred5) ## local richness

########### run null model
require(parallel)
cl <- makeCluster(6) ## number of cores = generally ncores -1
clusterEvalQ(cl, library(phytools))
clusterEvalQ(cl, library(phylobase))
clusterEvalQ(cl, library(dplyr))
# export your data and function
clusterExport(cl, c("local_1000","nsamples", "site_list_1000", "sampling_prevalence_NM"))

nm_teste <- parLapply(cl,site_list_1000, function (i) sampling_prevalence_NM (i,local=local_1000,nsamples=nsamples))## sample with the size of the obs community

stopCluster(cl)

#### took the traits of the sampled species
#traits
list_local_traits <- lapply (as.list(seq (1, length (nm_teste))), function (j)
	lapply (as.list (seq (1,length (local_1000))), function (k) ### para cada sitio
	lapply (as.list (seq(1,nsamples)), function (i)
	traits[which (rownames(traits) %in% nm_teste[[j]]$species[[k]] [,i]),])))

### removing sites with less than 4 sp
list_local_traits <- lapply (as.list(seq (1, length (nm_teste))), function (j)
	list_local_traits [[j]] [which (local_1000 >= 4)])
	
### apply hypervolume function to all sites and iterations
#require(parallel)
cl <- makeCluster(6) ## number of cores = generally ncores -1
clusterEvalQ(cl, library(hypervolume))
#export your data and function
clusterExport(cl, c("nsamples","nm_teste", "list_local_traits"))

hypervolume_local_null <-parLapply (cl, as.list(seq (1, length (nm_teste))), function (j)# para cada pool
	lapply (as.list(seq(1,length (list_local_traits [[1]]))), function (k)# para cada comunidade
	unlist (lapply (as.list (seq(1,nsamples)), function (i) ## para cada simula??o
	hypervolume_gaussian(list_local_traits [[j]][[k]][[i]],
	kde.bandwidth=estimate_bandwidth(list_local_traits [[j]][[k]][[i]],method="silverman"))@Volume))))

stopCluster (cl)

#save(hypervolume_local_null , file="hypervolume_local_null_1000traps.RData")

# OBSERVED
### observed local composition data
local_obs <- data.matrix(data_pred5 [which (local_1000 >= 4),order (colnames(data_pred5), decreasing=F)])
local_obs <- local_obs [,-which (colSums (local_obs)==0)]
list_local_obs <- apply (local_obs, 1, list)

### removing absent species
list_local_obs <- lapply (as.list (seq(1,length (list_local_obs))), function (i) 
	list_local_obs [[i]] [[1]] [which (list_local_obs [[i]] [[1]] == 1)])

lapply (lapply (list_local_traits_obs, function (i) apply (i, 2, duplicated)), function (k) apply (k,2,table))

### took the traits of species present in the communities
list_local_traits_obs <- lapply (as.list (seq(1,length (list_local_obs))), function (i)
	traits[which (rownames(traits) %in% names (list_local_obs [[i]])),])

## estimate local community hypervolume
cl <- makeCluster(6) ## number of cores = generally ncores -1
clusterEvalQ(cl, library(hypervolume))
#export your data and function
clusterExport(cl, c("list_local_traits_obs"))

hypervolume_local_obs<- parLapply (cl,list_local_traits_obs, function (i)
	hypervolume_gaussian(i,	kde.bandwidth=estimate_bandwidth(i,method="silverman"))@Volume)
stopCluster (cl)

## bandwidths
cl <- makeCluster(6) ## number of cores = generally ncores -1
clusterEvalQ(cl, library(hypervolume))
#export your data and function
clusterExport(cl, c("list_local_traits_obs"))

bandwidth_local_1000<- parLapply (cl,list_local_traits_obs, function (i)
	hypervolume_gaussian(i,	kde.bandwidth=estimate_bandwidth(i,method="silverman"))@Parameters)
stopCluster (cl)

### bandwidth values 
bandwidth_kde_local<-lapply (bandwidth_local_1000, function (i)i$kde.bandwidth)
colMeans (do.call (rbind.data.frame, bandwidth_kde_local))
apply (do.call (rbind.data.frame, bandwidth_kde_local),2,sd)
sd_count_local<-lapply (bandwidth_local_1000, function (i)i$sd.count)

sample_per_point<-lapply (bandwidth_local_1000, function (i)i$samples.per.point)
mean (unlist(sample_per_point))
sd(unlist(sample_per_point))

## pool composition data
pool_obs <- lapply (site_list_1000, function (i) i [which (local_1000 >=4),order(colnames(i), decreasing=F)])
traits_pool<-traits[which(rownames(traits) %in% colnames (pool_obs[[1]])),]
pool_obs <- lapply (pool_obs, function (i) i[,which (colnames(i) %in% rownames(traits_pool))])
pool_obs <- lapply (pool_obs, function (i) apply (i, 1, list))

### HERE'S THE THRESHOLD WE DISCUSS
#pool_obs <- lapply (seq (1, length (pool_obs)), function (j)
#	lapply (seq(1, length (pool_obs[[1]])), function (k)
#	pool_obs [[j]][[k]][[1]][which (pool_obs [[j]][[k]][[1]] >0.1)])) ## 
#
### HERE'S YOUR SUGGESTION
pool_sample <- replicate (niterations , lapply (seq(1,length (pool_obs)), function (pool) ## for each pool
	lapply (seq(1,length (pool_obs[[1]])), function (comm) ## for each community
	sample (pool_obs[[pool]][[comm]][[1]], size=round(sum(pool_obs[[pool]][[comm]][[1]])), 
	prob=(pool_obs[[pool]][[comm]][[1]])))), simplify = "array" )

## pool_sample [1,1] ## [species specific pool 1, sample 1]

## took the traits of each species in the pool

cl <- makeCluster(6) ## number of cores = generally ncores -1
clusterEvalQ(cl, library(hypervolume))
# export your data and function
clusterExport(cl, c("pool_sample","traits"))

traits_pool <- parLapply (cl,seq(1,nrow(pool_sample)), function (pool)
			lapply (seq(1,ncol(pool_sample)), function (sample)
			lapply (seq(1,length (pool_sample[1,1][[1]])), function (comm)
				traits [which (rownames(traits) %in% names(pool_sample [pool,sample][[1]][[comm]])),])))
stopCluster(cl)

# observed pool hypervolume 
## 10 ITERATIONS TOOK 2.5 DAYS TO RUN -USING 3 OF 4 CORES OF MY PC
cl <- makeCluster(6) ## number of cores = generally ncores -1
clusterEvalQ(cl, library(hypervolume))
# export your data and function
clusterExport(cl, c("traits_pool"))

hypervolume_pool_obs <- parLapply(cl, as.list(seq (1, length (traits_pool))), function (pool)# for each pool
	lapply (seq(1, length (traits_pool[[1]])), function (sample)# for each sample
	lapply (seq(1,length (traits_pool[[1]][[1]])), function (comm)# for each community ##   
      	hypervolume_gaussian (traits_pool[[pool]][[sample]][[comm]],
		kde.bandwidth=estimate_bandwidth(traits_pool[[pool]][[sample]][[comm]],method="silverman"))@Volume)))

stopCluster(cl)

## mean pool hypervolume (averaged across n interations)

mean_hypervolume_pool<-lapply (hypervolume_pool_obs, function (pool) # for each pool
	# calculate the mean pool hypervolume for each community (averaged across iterations)
	rowMeans (matrix (unlist(lapply (pool, unlist)),ncol=niterations,byrow=T))) 

#save (hypervolume_local_obs, hypervolume_pool_obs, file="hypervolume_1000traps_observed.RData")
save (hypervolume_pool_obs, file="hypervolumePOOL_1000traps_observed.RData")
save (mean_hypervolume_pool, file="mean_hypervolume_pool_1000traps_observed.RData")



# bandwidths used to calculate the pool hypervolume
cl <- makeCluster(6) ## number of cores = generally ncores -1
clusterEvalQ(cl, library(hypervolume))
# export your data and function
clusterExport(cl, c("traits_pool"))

bandwidths_pool_1000 <- parLapply(cl, as.list(seq (1, length (traits_pool))), function (pool)# for each pool
	lapply (seq(1, length (traits_pool[[1]])), function (sample)# for each sample
	lapply (seq(1, 2), function (comm)# for each community ## length (traits_pool[[1]][[1]])
      	hypervolume_gaussian (traits_pool[[pool]][[sample]][[comm]],
			kde.bandwidth=estimate_bandwidth(traits_pool[[pool]][[sample]][[comm]],method="silverman"))@Parameters)))

stopCluster(cl)

# Other parameters important to report
bandwidth_kde_pool<-lapply (bandwidths_pool_obs, function (pool)
	lapply (seq (1, length (bandwidths_pool_obs [[1]])), function (sample) 
	lapply (seq (1, length (bandwidths_pool_obs [[1]][[1]])), function (comm) 
	pool[[sample]][[comm]]$kde.bandwidth)))

matrix_bandwidth <- lapply (bandwidth_kde_pool, function (pool) # for each pool
	# calculate the mean pool hypervolume for each community (averaged across iterations)
	matrix (unlist(lapply (pool, unlist)),ncol=niterations,byrow=T))

matrix_bandwidth <- lapply (lapply (bandwidth_kde_pool [[1]], unlist), matrix,ncol=6,nrow=niterations,byrow=T)

# mean bandwidths over samples
lapply(matrix_bandwidth,colMeans) 
lapply(matrix_bandwidth,function (i) apply(i,2,sd)) 

## SD - nice value would be 3
sd_count_pool <- lapply (bandwidths_pool_obs, function (pool)
	lapply (seq (1, length (bandwidths_pool_obs [[1]])), function (sample) 
	lapply (seq (1, length (bandwidths_pool_obs [[1]][[1]])), function (comm) 
	pool[[sample]][[comm]]$sd.count)))

## samples (shoots) per point
sample_per_point_pool<- lapply (bandwidths_pool_obs, function (pool)
	lapply (seq (1, length (bandwidths_pool_obs [[1]])), function (sample) 
	lapply (seq (1, length (bandwidths_pool_obs [[1]][[1]])), function (comm) 
	pool[[sample]][[comm]]$samples.per.point)))

## mean and sd per pool
lapply(sample_per_point_pool, function (i) mean(unlist(i)))
lapply(sample_per_point_pool, function (i) sd(unlist(i)))

#save (bandwidth_pool_1000, bandwidth_local_1000, file="bandwidth_1000.RData")
































local_100 <- rowSums (data_pred100A)
nsamples<- 2 ## a few iterations to test the function

### run the null model to sort species from the pool (the same number of specis as observed in the local communities)

require(parallel)
cl <- makeCluster(6) ## number of cores = generally ncores -1
clusterEvalQ(cl, library(phytools))
clusterEvalQ(cl, library(phylobase))
clusterEvalQ(cl, library(dplyr))

# export your data and function
clusterExport(cl, c("local_100","nsamples", "site_list_100", "sampling_prevalence_NM"))

nm_teste <- parLapply(cl,site_list_100, function (i) sampling_prevalence_NM (i,local=local_100,nsamples=nsamples))## sample with the size of the obs community

stopCluster(cl)

#### took the traits of the sampled species
#traits
list_local_traits <- lapply (as.list(seq (1, length (nm_teste))), function (j)## for each pool
	lapply (as.list (seq (1,length (local_100))), function (k) ### site
	lapply (as.list (seq(1,nsamples)), function (i) ## and sample
	traits[which (rownames(traits) %in% nm_teste[[j]]$species[[k]] [,i]),]))) ## took the trait

### removing sites with less than 4 species
list_local_traits <- lapply (as.list(seq (1, length (nm_teste))), function (j)
	list_local_traits [[j]] [which (local_100 >= 4)])
	
## NULL FUNCTIONAL COMPLETENESS
### calculate null hypervolume for that random communities
cl <- makeCluster(6) ## number of cores = generally ncores -1
clusterEvalQ(cl, library(hypervolume))
#export your data and function
clusterExport(cl, c("nsamples","nm_teste", "list_local_traits"))

hypervolume_local_null <-parLapply (cl, as.list(seq (1, length (nm_teste))), function (j)# each pool
	lapply (as.list(seq(1,length (list_local_traits [[1]]))), function (k)# site
	unlist (lapply (as.list (seq(1,nsamples)), function (i) ## sample
	hypervolume_gaussian(list_local_traits [[j]][[k]][[i]],
	kde.bandwidth=estimate_bandwidth(list_local_traits [[j]][[k]][[i]],method="silverman"))@Volume)))) # calculate hyp and extract the object "Volume"

stopCluster (cl)

#save(hypervolume_local_null , file="hypervolume_local_null.RData")

# OBSERVED FUNCTIONAL COMPLETENESS

### data for observed communities
local_obs <- data.matrix(data_pred100A [which (local_100 >= 4),order (colnames(data_pred100A), decreasing=F)]) # removing sites with less than 4 sp
local_obs <- local_obs [,-which (colSums (local_obs)==0)] # removing absent species
list_local_obs <- apply (local_obs, 1, list) # it's faster to work with lists

## remove absent species
list_local_obs <- lapply (as.list (seq(1,length (list_local_obs))), function (i) 
	list_local_obs [[i]] [[1]] [which (list_local_obs [[i]] [[1]] == 1)])

### took the traits of species occurring in local communities
list_local_traits_obs <- lapply (as.list (seq(1,length (list_local_obs))), function (i)
	traits[which (rownames(traits) %in% names (list_local_obs [[i]])),])

## calculate observed hypervolume

cl <- makeCluster(6) ## number of cores = generally ncores -1
clusterEvalQ(cl, library(hypervolume))
#export your data and function
clusterExport(cl, c("list_local_traits_obs"))

hypervolume_local_obs<- parLapply (cl,list_local_traits_obs, function (i)
	hypervolume_gaussian(i,	kde.bandwidth=estimate_bandwidth(i,method="silverman"))@Volume) # calculate and extract the object "Volume"

stopCluster (cl)

cl <- makeCluster(6) ## number of cores = generally ncores -1
clusterEvalQ(cl, library(hypervolume))
#export your data and function
clusterExport(cl, c("list_local_traits_obs"))

bandwidths_local <- parLapply (cl,list_local_traits_obs, function (i)
	hypervolume_gaussian(i,	
	kde.bandwidth=estimate_bandwidth(i,method="silverman"))@Parameters) # calculate bandwidths

stopCluster (cl)

### values of bandwidth
bandwidth_kde_local<-lapply (bandwidths_local, function (i)i$kde.bandwidth)
colMeans (do.call (rbind.data.frame, bandwidth_kde_local))
apply (do.call (rbind.data.frame, bandwidth_kde_local),2,sd)
sd_count_local<-lapply (bandwidths_local, function (i)i$sd.count)
## samples per point
sample_per_point<-lapply (bandwidths_local, function (i)i$samples.per.point)
mean (unlist(sample_per_point))
sd(unlist(sample_per_point))

## calculate hypervolume for the pool
pool_obs <- lapply (site_list_100, function (i) i [which (local_100 >=4),order(colnames(i), decreasing=F)])
traits_pool<-traits[which(rownames(traits) %in% colnames (pool_obs[[1]])),]
pool_obs <- lapply (pool_obs, function (i) i[,which (colnames(i) %in% rownames(traits_pool))])
pool_obs <- lapply (pool_obs, function (i) apply (i, 1, list))

### HERE'S THE THRESHOLD WE DISCUSS
#pool_obs <- lapply (seq (1, length (pool_obs)), function (j)
#	lapply (seq(1, length (pool_obs[[1]])), function (k)
#	pool_obs [[j]][[k]][[1]][which (pool_obs [[j]][[k]][[1]] >0.1)])) ## 
#
### HERE'S YOUR SUGGESTION
niterations <- 50
pool_sample <- replicate (niterations , lapply (seq(1,length (pool_obs)), function (pool) ## for each pool
	lapply (seq(1,length (pool_obs[[1]])), function (comm) ## for each community
	sample (pool_obs[[pool]][[comm]][[1]], size=round(sum(pool_obs[[pool]][[comm]][[1]])), 
	prob=(pool_obs[[pool]][[comm]][[1]])))), simplify = "array" )

## pool_sample [1,1] ## [species specific pool 1, sample 1]

## took the traits of each species in the pool
cl <- makeCluster(6) ## number of cores = generally ncores -1
clusterEvalQ(cl, library(hypervolume))
# export your data and function
clusterExport(cl, c("pool_sample","traits"))

traits_pool <- parLapply (cl,seq(1,nrow(pool_sample)), function (pool)
			lapply (seq(1,ncol(pool_sample)), function (sample)
			lapply (seq(1,length (pool_sample[1,1][[1]])), function (comm)
				traits [which (rownames(traits) %in% names(pool_sample [pool,sample][[1]][[comm]])),])))
stopCluster(cl)

# observed pool hypervolume 
## 10 ITERATIONS TOOK 2.5 DAYS TO RUN -USING 3 OF 4 CORES OF MY PC
cl <- makeCluster(6) ## number of cores = generally ncores -1
clusterEvalQ(cl, library(hypervolume))
# export your data and function
clusterExport(cl, c("traits_pool"))

hypervolume_pool_obs <- parLapply(cl, as.list(seq (1, length (traits_pool))), function (pool)# for each pool
	lapply (seq(1, length (traits_pool[[1]])), function (sample)# for each sample
	lapply (seq(1,length (traits_pool[[1]][[1]])), function (comm)# for each community ##   
      	hypervolume_gaussian (traits_pool[[pool]][[sample]][[comm]],
		kde.bandwidth=estimate_bandwidth(traits_pool[[pool]][[sample]][[comm]],method="silverman"))@Volume)))

stopCluster(cl)

## mean pool hypervolume (averaged across n interations)

mean_hypervolume_pool<-lapply (hypervolume_pool_obs, function (pool) # for each pool
	# calculate the mean pool hypervolume for each community (averaged across iterations)
	rowMeans (matrix (unlist(lapply (pool, unlist)),ncol=niterations,byrow=T))) 

save (hypervolume_pool_obs, file="hypervolumePOOL_100traps_observed.RData")
save (mean_hypervolume_pool, file="mean_hypervolume_pool_100traps_observed.RData")

# bandwidths used to calculate the pool hypervolume
cl <- makeCluster(6) ## number of cores = generally ncores -1
clusterEvalQ(cl, library(hypervolume))
# export your data and function
clusterExport(cl, c("traits_pool"))

bandwidths_pool_obs <- parLapply(cl, as.list(seq (1, length (traits_pool))), function (pool)# for each pool
	lapply (seq(1, length (traits_pool[[1]])), function (sample)# for each sample
	lapply (seq(1, 2), function (comm)# for each community ## length (traits_pool[[1]][[1]])
      	hypervolume_gaussian (traits_pool[[pool]][[sample]][[comm]],
			kde.bandwidth=estimate_bandwidth(traits_pool[[pool]][[sample]][[comm]],method="silverman"))@Parameters)))

stopCluster(cl)

# Other parameters important to report
bandwidth_kde_pool<-lapply (bandwidths_pool_obs, function (pool)
	lapply (seq (1, length (bandwidths_pool_obs [[1]])), function (sample) 
	lapply (seq (1, length (bandwidths_pool_obs [[1]][[1]])), function (comm) 
	pool[[sample]][[comm]]$kde.bandwidth)))

matrix_bandwidth <- lapply (bandwidth_kde_pool, function (pool) # for each pool
	# calculate the mean pool hypervolume for each community (averaged across iterations)
	matrix (unlist(lapply (pool, unlist)),ncol=niterations,byrow=T))

matrix_bandwidth <- lapply (lapply (bandwidth_kde_pool [[1]], unlist), matrix,ncol=6,nrow=niterations,byrow=T)

# mean bandwidths over samples
lapply(matrix_bandwidth,colMeans) 
lapply(matrix_bandwidth,function (i) apply(i,2,sd)) 

## SD - nice value would be 3
sd_count_pool <- lapply (bandwidths_pool_obs, function (pool)
	lapply (seq (1, length (bandwidths_pool_obs [[1]])), function (sample) 
	lapply (seq (1, length (bandwidths_pool_obs [[1]][[1]])), function (comm) 
	pool[[sample]][[comm]]$sd.count)))

## samples (shoots) per point
sample_per_point_pool<- lapply (bandwidths_pool_obs, function (pool)
	lapply (seq (1, length (bandwidths_pool_obs [[1]])), function (sample) 
	lapply (seq (1, length (bandwidths_pool_obs [[1]][[1]])), function (comm) 
	pool[[sample]][[comm]]$samples.per.point)))

## mean and sd per pool
lapply(sample_per_point_pool, function (i) mean(unlist(i)))
lapply(sample_per_point_pool, function (i) sd(unlist(i)))

#save (bandwidths_pool, bandwidths_local, file="bandwidth_100.RData")


#######################
######    500 traps  ##
#######################

local_500 <- rowSums (data_pred4) # local richness (sample size)

##########
# run null model
require(parallel)
cl <- makeCluster(6) ## number of cores = generally ncores -1
clusterEvalQ(cl, library(phytools))
clusterEvalQ(cl, library(phylobase))
clusterEvalQ(cl, library(dplyr))
# export your data and function
clusterExport(cl, c("local_500","nsamples", "site_list", "sampling_prevalence_NM"))

nm_teste <- parLapply(cl,site_list, function (i) sampling_prevalence_NM (i,local=local_500,nsamples=nsamples))## sample with the size of the obs community

stopCluster(cl)

#### took the traits of the sampled species
#traits
list_local_traits <- lapply (as.list(seq (1, length (nm_teste))), function (j)# for each pool
	lapply (as.list (seq (1,length (local_500))), function (k) ### for each community
	lapply (as.list (seq(1,nsamples)), function (i) # for each iteration
	traits[which (rownames(traits) %in% nm_teste[[j]]$species[[k]] [,i]),]))) ## took the traits

### removing the sites with less than 4 species
list_local_traits <- lapply (as.list(seq (1, length (nm_teste))), function (j)
	list_local_traits [[j]] [which (local_500 >= 4)])
	
### apply to all communities and iterations
#require(parallel)
cl <- makeCluster(6) ## number of cores = generally ncores -1
clusterEvalQ(cl, library(hypervolume))
#export your data and function
clusterExport(cl, c("nsamples","nm_teste", "list_local_traits"))

hypervolume_local_null <-parLapply (cl, as.list(seq (1, length (nm_teste))), function (j)# para cada pool
	lapply (as.list(seq(1,length (list_local_traits [[1]]))), function (k)# para cada comunidade
	unlist (lapply (as.list (seq(1,nsamples)), function (i) ## para cada simula??o
	hypervolume_gaussian(list_local_traits [[j]][[k]][[i]],
	kde.bandwidth=estimate_bandwidth(list_local_traits [[j]][[k]][[i]],method="silverman"))@Volume))))

stopCluster (cl)

#save(hypervolume_local_null , file="hypervolume_local_null_500traps.RData")

# OBSERVED
### data of local species composition
local_obs <- data.matrix(data_pred4 [which (local_500 >= 4),order (colnames(data_pred4), decreasing=F)])
local_obs <- local_obs [,-which (colSums (local_obs)==0)]
list_local_obs <- apply (local_obs, 1, list)

### removing absent species
list_local_obs <- lapply (as.list (seq(1,length (list_local_obs))), function (i) 
	list_local_obs [[i]] [[1]] [which (list_local_obs [[i]] [[1]] == 1)])

### took the traits of the species of each communities
list_local_traits_obs <- lapply (as.list (seq(1,length (list_local_obs))), function (i)
	traits[which (rownames(traits) %in% names (list_local_obs [[i]])),])

## calculate local community hypervolume
cl <- makeCluster(6) ## number of cores = generally ncores -1
clusterEvalQ(cl, library(hypervolume))
#export your data and function
clusterExport(cl, c("list_local_traits_obs"))

hypervolume_local_obs<- parLapply (cl,list_local_traits_obs, function (i)
	hypervolume_gaussian(i,	kde.bandwidth=estimate_bandwidth(i,method="silverman"))@Volume)

stopCluster (cl)

## bandwidth
cl <- makeCluster(6) ## number of cores = generally ncores -1
clusterEvalQ(cl, library(hypervolume))
#export your data and function
clusterExport(cl, c("list_local_traits_obs"))

bandwidths_local_500 <- parLapply (cl,list_local_traits_obs, function (i)
	hypervolume_gaussian(i,	
	kde.bandwidth=estimate_bandwidth(i,method="silverman"))@Parameters)

stopCluster (cl)

### bandwidth values (important to report)
bandwidth_kde_local<-lapply (bandwidths_local_500, function (i)i$kde.bandwidth)
colMeans (do.call (rbind.data.frame, bandwidth_kde_local))
apply (do.call (rbind.data.frame, bandwidth_kde_local),2,sd)
sd_count_local<-lapply (bandwidths_local_500, function (i)i$sd.count)

sample_per_point<-lapply (bandwidths_local_500, function (i)i$samples.per.point)
mean (unlist(sample_per_point))
sd(unlist(sample_per_point))

## pool composition data
pool_obs <- lapply (site_list, function (i) i [which (local_500 >=4),order(colnames(i), decreasing=F)])
traits_pool<-traits[which(rownames(traits) %in% colnames (pool_obs[[1]])),]
pool_obs <- lapply (pool_obs, function (i) i[,which (colnames(i) %in% rownames(traits_pool))])
pool_obs <- lapply (pool_obs, function (i) apply (i, 1, list))

### HERE'S THE THRESHOLD WE DISCUSS
#pool_obs <- lapply (seq (1, length (pool_obs)), function (j)
#	lapply (seq(1, length (pool_obs[[1]])), function (k)
#	pool_obs [[j]][[k]][[1]][which (pool_obs [[j]][[k]][[1]] >0.1)])) ## 
#
### HERE'S YOUR SUGGESTION
#niterations <- 2
pool_sample <- replicate (niterations , lapply (seq(1,length (pool_obs)), function (pool) ## for each pool
	lapply (seq(1,length (pool_obs[[1]])), function (comm) ## for each community
	sample (pool_obs[[pool]][[comm]][[1]], size=round(sum(pool_obs[[pool]][[comm]][[1]])), 
	prob=(pool_obs[[pool]][[comm]][[1]])))), simplify = "array" )

## pool_sample [1,1] ## [species specific pool 1, sample 1]

## took the traits of each species in the pool
cl <- makeCluster(6) ## number of cores = generally ncores -1
clusterEvalQ(cl, library(hypervolume))
# export your data and function
clusterExport(cl, c("pool_sample","traits"))

traits_pool <- parLapply (cl,seq(1,nrow(pool_sample)), function (pool)
			lapply (seq(1,ncol(pool_sample)), function (sample)
			lapply (seq(1,length (pool_sample[1,1][[1]])), function (comm)
				traits [which (rownames(traits) %in% names(pool_sample [pool,sample][[1]][[comm]])),])))

stopCluster(cl)

# observed pool hypervolume 
## 10 ITERATIONS TOOK 2.5 DAYS TO RUN -USING 3 OF 4 CORES OF MY PC
cl <- makeCluster(6) ## number of cores = generally ncores -1
clusterEvalQ(cl, library(hypervolume))
# export your data and function
clusterExport(cl, c("traits_pool"))

hypervolume_pool_obs <- parLapply(cl, as.list(seq (1, length (traits_pool))), function (pool)# for each pool
	lapply (seq(1, length (traits_pool[[1]])), function (sample)# for each sample
	lapply (seq(1,length (traits_pool[[1]][[1]])), function (comm)# for each community ##   
      	hypervolume_gaussian (traits_pool[[pool]][[sample]][[comm]],
		kde.bandwidth=estimate_bandwidth(traits_pool[[pool]][[sample]][[comm]],method="silverman"))@Volume)))

stopCluster(cl)

## mean pool hypervolume (averaged across n interations)

mean_hypervolume_pool<-lapply (hypervolume_pool_obs, function (pool) # for each pool
	# calculate the mean pool hypervolume for each community (averaged across iterations)
	rowMeans (matrix (unlist(lapply (pool, unlist)),ncol=niterations,byrow=T))) 

#save (hypervolume_local_obs, hypervolume_pool_obs, file="hypervolume_500traps_observed.RData")
save (hypervolume_pool_obs, file="hypervolumePOOL_500traps_observed.RData")
save (mean_hypervolume_pool, file="mean_hypervolume_pool_500traps_observed.RData")


# bandwidths we used to calculate the pool hypervolume
cl <- makeCluster(6) ## number of cores = generally ncores -1
clusterEvalQ(cl, library(hypervolume))
# export your data and function
clusterExport(cl, c("traits_pool"))

bandwidths_pool_obs <- parLapply(cl, as.list(seq (1, length (traits_pool))), function (pool)# for each pool
	lapply (seq(1, length (traits_pool[[1]])), function (sample)# for each sample
	lapply (seq(1, 2), function (comm)# for each community ## length (traits_pool[[1]][[1]])
      	hypervolume_gaussian (traits_pool[[pool]][[sample]][[comm]],
			kde.bandwidth=estimate_bandwidth(traits_pool[[pool]][[sample]][[comm]],method="silverman"))@Parameters)))

stopCluster(cl)

# Other parameters important to report
bandwidth_kde_pool<-lapply (bandwidths_pool_obs, function (pool)
	lapply (seq (1, length (bandwidths_pool_obs [[1]])), function (sample) 
	lapply (seq (1, length (bandwidths_pool_obs [[1]][[1]])), function (comm) 
	pool[[sample]][[comm]]$kde.bandwidth)))

matrix_bandwidth <- lapply (bandwidth_kde_pool, function (pool) # for each pool
	# calculate the mean pool hypervolume for each community (averaged across iterations)
	matrix (unlist(lapply (pool, unlist)),ncol=niterations,byrow=T))

matrix_bandwidth <- lapply (lapply (bandwidth_kde_pool [[1]], unlist), matrix,ncol=6,nrow=niterations,byrow=T)

# mean bandwidths over samples
lapply(matrix_bandwidth,colMeans) 
lapply(matrix_bandwidth,function (i) apply(i,2,sd)) 

## SD - nice value would be 3
sd_count_pool <- lapply (bandwidths_pool_obs, function (pool)
	lapply (seq (1, length (bandwidths_pool_obs [[1]])), function (sample) 
	lapply (seq (1, length (bandwidths_pool_obs [[1]][[1]])), function (comm) 
	pool[[sample]][[comm]]$sd.count)))

## samples (shoots) per point
sample_per_point_pool<- lapply (bandwidths_pool_obs, function (pool)
	lapply (seq (1, length (bandwidths_pool_obs [[1]])), function (sample) 
	lapply (seq (1, length (bandwidths_pool_obs [[1]][[1]])), function (comm) 
	pool[[sample]][[comm]]$samples.per.point)))

## mean and sd per pool
lapply(sample_per_point_pool, function (i) mean(unlist(i)))
lapply(sample_per_point_pool, function (i) sd(unlist(i)))

#save (bandwidths_pool_obs, bandwidths_local_500, file="bandwidth_500.RData")

######################
####  1,000 traps  ###
######################

local_1000 <- rowSums (data_pred5) ## local richness

########### run null model
require(parallel)
cl <- makeCluster(6) ## number of cores = generally ncores -1
clusterEvalQ(cl, library(phytools))
clusterEvalQ(cl, library(phylobase))
clusterEvalQ(cl, library(dplyr))
# export your data and function
clusterExport(cl, c("local_1000","nsamples", "site_list_1000", "sampling_prevalence_NM"))

nm_teste <- parLapply(cl,site_list_1000, function (i) sampling_prevalence_NM (i,local=local_1000,nsamples=nsamples))## sample with the size of the obs community

stopCluster(cl)

#### took the traits of the sampled species
#traits
list_local_traits <- lapply (as.list(seq (1, length (nm_teste))), function (j)
	lapply (as.list (seq (1,length (local_1000))), function (k) ### para cada sitio
	lapply (as.list (seq(1,nsamples)), function (i)
	traits[which (rownames(traits) %in% nm_teste[[j]]$species[[k]] [,i]),])))

### removing sites with less than 4 sp
list_local_traits <- lapply (as.list(seq (1, length (nm_teste))), function (j)
	list_local_traits [[j]] [which (local_1000 >= 4)])
	
### apply hypervolume function to all sites and iterations
#require(parallel)
cl <- makeCluster(6) ## number of cores = generally ncores -1
clusterEvalQ(cl, library(hypervolume))
#export your data and function
clusterExport(cl, c("nsamples","nm_teste", "list_local_traits"))

hypervolume_local_null <-parLapply (cl, as.list(seq (1, length (nm_teste))), function (j)# para cada pool
	lapply (as.list(seq(1,length (list_local_traits [[1]]))), function (k)# para cada comunidade
	unlist (lapply (as.list (seq(1,nsamples)), function (i) ## para cada simula??o
	hypervolume_gaussian(list_local_traits [[j]][[k]][[i]],
	kde.bandwidth=estimate_bandwidth(list_local_traits [[j]][[k]][[i]],method="silverman"))@Volume))))

stopCluster (cl)

#save(hypervolume_local_null , file="hypervolume_local_null_1000traps.RData")

# OBSERVED
### observed local composition data
local_obs <- data.matrix(data_pred5 [which (local_1000 >= 4),order (colnames(data_pred5), decreasing=F)])
local_obs <- local_obs [,-which (colSums (local_obs)==0)]
list_local_obs <- apply (local_obs, 1, list)

### removing absent species
list_local_obs <- lapply (as.list (seq(1,length (list_local_obs))), function (i) 
	list_local_obs [[i]] [[1]] [which (list_local_obs [[i]] [[1]] == 1)])

lapply (lapply (list_local_traits_obs, function (i) apply (i, 2, duplicated)), function (k) apply (k,2,table))

### took the traits of species present in the communities
list_local_traits_obs <- lapply (as.list (seq(1,length (list_local_obs))), function (i)
	traits[which (rownames(traits) %in% names (list_local_obs [[i]])),])

## estimate local community hypervolume
cl <- makeCluster(6) ## number of cores = generally ncores -1
clusterEvalQ(cl, library(hypervolume))
#export your data and function
clusterExport(cl, c("list_local_traits_obs"))

hypervolume_local_obs<- parLapply (cl,list_local_traits_obs, function (i)
	hypervolume_gaussian(i,	kde.bandwidth=estimate_bandwidth(i,method="silverman"))@Volume)
stopCluster (cl)

## bandwidths
cl <- makeCluster(6) ## number of cores = generally ncores -1
clusterEvalQ(cl, library(hypervolume))
#export your data and function
clusterExport(cl, c("list_local_traits_obs"))

bandwidth_local_1000<- parLapply (cl,list_local_traits_obs, function (i)
	hypervolume_gaussian(i,	kde.bandwidth=estimate_bandwidth(i,method="silverman"))@Parameters)
stopCluster (cl)

### bandwidth values 
bandwidth_kde_local<-lapply (bandwidth_local_1000, function (i)i$kde.bandwidth)
colMeans (do.call (rbind.data.frame, bandwidth_kde_local))
apply (do.call (rbind.data.frame, bandwidth_kde_local),2,sd)
sd_count_local<-lapply (bandwidth_local_1000, function (i)i$sd.count)

sample_per_point<-lapply (bandwidth_local_1000, function (i)i$samples.per.point)
mean (unlist(sample_per_point))
sd(unlist(sample_per_point))

## pool composition data
pool_obs <- lapply (site_list_1000, function (i) i [which (local_1000 >=4),order(colnames(i), decreasing=F)])
traits_pool<-traits[which(rownames(traits) %in% colnames (pool_obs[[1]])),]
pool_obs <- lapply (pool_obs, function (i) i[,which (colnames(i) %in% rownames(traits_pool))])
pool_obs <- lapply (pool_obs, function (i) apply (i, 1, list))

### HERE'S THE THRESHOLD WE DISCUSS
#pool_obs <- lapply (seq (1, length (pool_obs)), function (j)
#	lapply (seq(1, length (pool_obs[[1]])), function (k)
#	pool_obs [[j]][[k]][[1]][which (pool_obs [[j]][[k]][[1]] >0.1)])) ## 
#
### HERE'S YOUR SUGGESTION
pool_sample <- replicate (niterations , lapply (seq(1,length (pool_obs)), function (pool) ## for each pool
	lapply (seq(1,length (pool_obs[[1]])), function (comm) ## for each community
	sample (pool_obs[[pool]][[comm]][[1]], size=round(sum(pool_obs[[pool]][[comm]][[1]])), 
	prob=(pool_obs[[pool]][[comm]][[1]])))), simplify = "array" )

## pool_sample [1,1] ## [species specific pool 1, sample 1]

## took the traits of each species in the pool

cl <- makeCluster(6) ## number of cores = generally ncores -1
clusterEvalQ(cl, library(hypervolume))
# export your data and function
clusterExport(cl, c("pool_sample","traits"))

traits_pool <- parLapply (cl,seq(1,nrow(pool_sample)), function (pool)
			lapply (seq(1,ncol(pool_sample)), function (sample)
			lapply (seq(1,length (pool_sample[1,1][[1]])), function (comm)
				traits [which (rownames(traits) %in% names(pool_sample [pool,sample][[1]][[comm]])),])))
stopCluster(cl)

# observed pool hypervolume 
## 10 ITERATIONS TOOK 2.5 DAYS TO RUN -USING 3 OF 4 CORES OF MY PC
cl <- makeCluster(6) ## number of cores = generally ncores -1
clusterEvalQ(cl, library(hypervolume))
# export your data and function
clusterExport(cl, c("traits_pool"))

hypervolume_pool_obs <- parLapply(cl, as.list(seq (1, length (traits_pool))), function (pool)# for each pool
	lapply (seq(1, length (traits_pool[[1]])), function (sample)# for each sample
	lapply (seq(1,length (traits_pool[[1]][[1]])), function (comm)# for each community ##   
      	hypervolume_gaussian (traits_pool[[pool]][[sample]][[comm]],
		kde.bandwidth=estimate_bandwidth(traits_pool[[pool]][[sample]][[comm]],method="silverman"))@Volume)))

stopCluster(cl)

## mean pool hypervolume (averaged across n interations)

mean_hypervolume_pool<-lapply (hypervolume_pool_obs, function (pool) # for each pool
	# calculate the mean pool hypervolume for each community (averaged across iterations)
	rowMeans (matrix (unlist(lapply (pool, unlist)),ncol=niterations,byrow=T))) 

#save (hypervolume_local_obs, hypervolume_pool_obs, file="hypervolume_1000traps_observed.RData")
save (hypervolume_pool_obs, file="hypervolumePOOL_1000traps_observed.RData")
save (mean_hypervolume_pool, file="mean_hypervolume_pool_1000traps_observed.RData")



# bandwidths used to calculate the pool hypervolume
cl <- makeCluster(6) ## number of cores = generally ncores -1
clusterEvalQ(cl, library(hypervolume))
# export your data and function
clusterExport(cl, c("traits_pool"))

bandwidths_pool_1000 <- parLapply(cl, as.list(seq (1, length (traits_pool))), function (pool)# for each pool
	lapply (seq(1, length (traits_pool[[1]])), function (sample)# for each sample
	lapply (seq(1, 2), function (comm)# for each community ## length (traits_pool[[1]][[1]])
      	hypervolume_gaussian (traits_pool[[pool]][[sample]][[comm]],
			kde.bandwidth=estimate_bandwidth(traits_pool[[pool]][[sample]][[comm]],method="silverman"))@Parameters)))

stopCluster(cl)

# Other parameters important to report
bandwidth_kde_pool<-lapply (bandwidths_pool_obs, function (pool)
	lapply (seq (1, length (bandwidths_pool_obs [[1]])), function (sample) 
	lapply (seq (1, length (bandwidths_pool_obs [[1]][[1]])), function (comm) 
	pool[[sample]][[comm]]$kde.bandwidth)))

matrix_bandwidth <- lapply (bandwidth_kde_pool, function (pool) # for each pool
	# calculate the mean pool hypervolume for each community (averaged across iterations)
	matrix (unlist(lapply (pool, unlist)),ncol=niterations,byrow=T))

matrix_bandwidth <- lapply (lapply (bandwidth_kde_pool [[1]], unlist), matrix,ncol=6,nrow=niterations,byrow=T)

# mean bandwidths over samples
lapply(matrix_bandwidth,colMeans) 
lapply(matrix_bandwidth,function (i) apply(i,2,sd)) 

## SD - nice value would be 3
sd_count_pool <- lapply (bandwidths_pool_obs, function (pool)
	lapply (seq (1, length (bandwidths_pool_obs [[1]])), function (sample) 
	lapply (seq (1, length (bandwidths_pool_obs [[1]][[1]])), function (comm) 
	pool[[sample]][[comm]]$sd.count)))

## samples (shoots) per point
sample_per_point_pool<- lapply (bandwidths_pool_obs, function (pool)
	lapply (seq (1, length (bandwidths_pool_obs [[1]])), function (sample) 
	lapply (seq (1, length (bandwidths_pool_obs [[1]][[1]])), function (comm) 
	pool[[sample]][[comm]]$samples.per.point)))

## mean and sd per pool
lapply(sample_per_point_pool, function (i) mean(unlist(i)))
lapply(sample_per_point_pool, function (i) sd(unlist(i)))

#save (bandwidth_pool_1000, bandwidth_local_1000, file="bandwidth_1000.RData")






























































########### rodar o modelo nulo que sorteia especies com base na prevalencia no pool
require(parallel)
cl <- makeCluster(6) ## number of cores = generally ncores -1
clusterEvalQ(cl, library(phytools))
clusterEvalQ(cl, library(phylobase))
clusterEvalQ(cl, library(dplyr))
# export your data and function
clusterExport(cl, c("local_100","nsamples", "site_list_100", "sampling_prevalence_NM"))

nm_teste <- parLapply(cl,site_list_100, function (i) 
	sampling_prevalence_NM (i,local=local_100,nsamples=nsamples))
	## sample with the size of the obs community
#nm_teste <- lapply(site_list, function (i) sampling_prevalence_NM (i,local=local,nsamples=nsamples))## sample with the size of the obs community

stopCluster(cl)

#### pegar os atributos das especies sorteadas
#traits
list_local_traits <- lapply (as.list(seq (1, length (nm_teste))), function (j)
	lapply (as.list (seq (1,length (local_100))), function (k) ### para cada sitio
	lapply (as.list (seq(1,nsamples)), function (i)
	traits[which (rownames(traits) %in% nm_teste[[j]]$species[[k]] [,i]),])))

### removendo sitios sem sp ou com uma especie
list_local_traits <- lapply (as.list(seq (1, length (nm_teste))), function (j)
	list_local_traits [[j]] [which (local_100 >= 4)])
	
### aplicar pra todas os sitios e simulacoes
#require(parallel)
cl <- makeCluster(6) ## number of cores = generally ncores -1
clusterEvalQ(cl, library(hypervolume))
#export your data and function
clusterExport(cl, c("nsamples","nm_teste", "list_local_traits"))

hypervolume_local_null <-parLapply (cl, as.list(seq (1, length (nm_teste))), function (j)# para cada pool
	lapply (as.list(seq(1,length (list_local_traits [[1]]))), function (k)# para cada comunidade
	unlist (lapply (as.list (seq(1,nsamples)), function (i) ## para cada simula??o
	hypervolume_gaussian(list_local_traits [[j]][[k]][[i]],
	kde.bandwidth=estimate_bandwidth(list_local_traits [[j]][[k]][[i]],method="silverman"))@Volume))))

stopCluster (cl)

save(hypervolume_local_null , file="hypervolume_local_null.RData")

############# traits for species of the pool = NAO PRECISA, O POOL JA ? UMA EXPECTATIVA NULA
#list_pool_traits <- lapply (as.list(seq (1, length (nm_teste))), function (j) ## para cada pool
#	lapply (as.list (seq (1,length (local_100))), function (k) ### para cada sitio
#	#lapply (as.list (seq(1,nsamples)), function (i) ## para cada simulacao
#	traits[which (rownames (traits) %in% colnames (site_list_100[[j]][k,][which (site_list_100[[j]][k,] > 0.1)])),]))#)
#
# pegar dados das comunidades com 3 ou mais especeis
#list_pool_traits_null <- lapply (as.list(seq (1, length (nm_teste))), function (j)
#	list_pool_traits [[j]] [which (local_100 >= 4)])
#
### aplicar pra todas os sitios e simulacoes
#require(parallel)
#cl <- makeCluster(6) ## number of cores = generally ncores -1
#clusterEvalQ(cl, library(hypervolume))
#export your data and function
#clusterExport(cl, c("nsamples","nm_teste", "list_pool_traits"))
#hypervolume_pool_null <-parLapply (cl, as.list(seq (1, length (nm_teste))), function (j)
#	lapply (as.list(seq(1,length (list_pool_traits [[1]]))), function (k)
#	lapply (as.list (seq(1,nsamples)), function (i) ## para cada simulacao
#	hypervolume_gaussian(list_pool_traits [[j]][[k]][[i]],
#	kde.bandwidth=estimate_bandwidth(list_pool_traits [[j]][[k]][[i]],method="silverman"))@Volume)))
#
#stopCluster (cl)
#save (hypervolume_pool_null, file= "hypervolume_pool_100traps_NULL.RData")

# OBSERVADOS
### dados de composicao observados das comunidades
local_obs <- data.matrix(data_pred100A [which (local_100 >= 4),order (colnames(data_pred100A), decreasing=F)])
local_obs <- local_obs [,-which (colSums (local_obs)==0)]
list_local_obs <- apply (local_obs, 1, list)
### removendo especies ausentes
list_local_obs <- lapply (as.list (seq(1,length (list_local_obs))), function (i) 
	list_local_obs [[i]] [[1]] [which (list_local_obs [[i]] [[1]] == 1)])

lapply (lapply (list_local_obs, function (i) apply (i, 2, duplicated)), function (k) apply (k,2,table))

### pegar os atributos das especies de cada comunidade
list_local_traits_obs <- lapply (as.list (seq(1,length (list_local_obs))), function (i)
	traits[which (rownames(traits) %in% names (list_local_obs [[i]])),])

## calcular o volume observado
cl <- makeCluster(6) ## number of cores = generally ncores -1
clusterEvalQ(cl, library(hypervolume))
#export your data and function
clusterExport(cl, c("list_local_traits_obs"))

hypervolume_local_obs<- parLapply (cl,list_local_traits_obs, function (i)
	hypervolume_gaussian(i,	kde.bandwidth=estimate_bandwidth(i,method="silverman"))@Volume)
stopCluster (cl)

## composicao do pool
pool_obs <- lapply (site_list_100, function (i) i [which (local_100 >=4),order(colnames(i), decreasing=F)])
traits_pool<-traits[which(rownames(traits) %in% colnames (pool_obs[[1]])),]
pool_obs <- lapply (pool_obs, function (i) i[,which (colnames(i) %in% rownames(traits_pool))])
pool_obs <- lapply (pool_obs, function (i) apply (i, 1, list))

pool_obs <- lapply (seq (1, length (pool_obs)), function (j)
	lapply (seq(1, length (pool_obs[[1]])), function (k)
	pool_obs [[j]][[k]][[1]][which (pool_obs [[j]][[k]][[1]] >0.1)]))

## pegar os atributos das especies compondo o pool
traits_pool <- lapply (seq (1, length (pool_obs)), function (j)
	lapply (seq(1, length (pool_obs[[1]])), function (k)
	traits_pool [which (rownames(traits_pool) %in% names (pool_obs[[j]][[k]])),]))

# calcular o volume observado do pool
cl <- makeCluster(6) ## number of cores = generally ncores -1
clusterEvalQ(cl, library(hypervolume))
# export your data and function
clusterExport(cl, c("pool_obs","traits_pool"))

hypervolume_pool_obs <- parLapply(cl, as.list(seq (1, length (pool_obs))), function (j)#
	lapply (seq(1, length (pool_obs[[1]])), function (k)#
	hypervolume_gaussian (traits_pool[[j]][[k]],kde.bandwidth=estimate_bandwidth(traits_pool[[j]][[k]],method="silverman"))@Volume))

stopCluster(cl)

save (hypervolume_local_obs, hypervolume_pool_obs, file="hypervolume_100traps_observed.RData")


#500 traps

local_500 <- rowSums (data_pred4)
nsamples<-100

########### rodar o modelo nulo que sorteia especies com base na prevalencia no pool
require(parallel)
cl <- makeCluster(6) ## number of cores = generally ncores -1
clusterEvalQ(cl, library(phytools))
clusterEvalQ(cl, library(phylobase))
clusterEvalQ(cl, library(dplyr))
# export your data and function
clusterExport(cl, c("local_500","nsamples", "site_list", "sampling_prevalence_NM"))

nm_teste <- parLapply(cl,site_list, function (i) sampling_prevalence_NM (i,local=local_500,nsamples=nsamples))## sample with the size of the obs community

stopCluster(cl)

#### pegar os atributos das especies sorteadas
#traits
list_local_traits <- lapply (as.list(seq (1, length (nm_teste))), function (j)
	lapply (as.list (seq (1,length (local_500))), function (k) ### para cada sitio
	lapply (as.list (seq(1,nsamples)), function (i)
	traits[which (rownames(traits) %in% nm_teste[[j]]$species[[k]] [,i]),])))

### removendo sitios sem sp ou com uma especie
list_local_traits <- lapply (as.list(seq (1, length (nm_teste))), function (j)
	list_local_traits [[j]] [which (local_500 >= 4)])
	
### aplicar pra todas os sitios e simulacoes
#require(parallel)
cl <- makeCluster(6) ## number of cores = generally ncores -1
clusterEvalQ(cl, library(hypervolume))
#export your data and function
clusterExport(cl, c("nsamples","nm_teste", "list_local_traits"))

hypervolume_local_null <-parLapply (cl, as.list(seq (1, length (nm_teste))), function (j)# para cada pool
	lapply (as.list(seq(1,length (list_local_traits [[1]]))), function (k)# para cada comunidade
	unlist (lapply (as.list (seq(1,nsamples)), function (i) ## para cada simula??o
	hypervolume_gaussian(list_local_traits [[j]][[k]][[i]],
	kde.bandwidth=estimate_bandwidth(list_local_traits [[j]][[k]][[i]],method="silverman"))@Volume))))

stopCluster (cl)

save(hypervolume_local_null , file="hypervolume_local_null_500traps.RData")

# OBSERVADOS
### dados de composicao observados das comunidades
local_obs <- data.matrix(data_pred4 [which (local_500 >= 4),order (colnames(data_pred4), decreasing=F)])
local_obs <- local_obs [,-which (colSums (local_obs)==0)]
list_local_obs <- apply (local_obs, 1, list)
### removendo especies ausentes
list_local_obs <- lapply (as.list (seq(1,length (list_local_obs))), function (i) 
	list_local_obs [[i]] [[1]] [which (list_local_obs [[i]] [[1]] == 1)])

lapply (lapply (list_local_traits_obs, function (i) apply (i, 2, duplicated)), function (k) apply (k,2,table))

### pegar os atributos das especies de cada comunidade
list_local_traits_obs <- lapply (as.list (seq(1,length (list_local_obs))), function (i)
	traits[which (rownames(traits) %in% names (list_local_obs [[i]])),])

## calcular o volume observado
cl <- makeCluster(6) ## number of cores = generally ncores -1
clusterEvalQ(cl, library(hypervolume))
#export your data and function
clusterExport(cl, c("list_local_traits_obs"))

hypervolume_local_obs<- parLapply (cl,list_local_traits_obs, function (i)
	hypervolume_gaussian(i,	kde.bandwidth=estimate_bandwidth(i,method="silverman"))@Volume)
stopCluster (cl)

## composicao do pool
pool_obs <- lapply (site_list, function (i) i [which (local_500 >=4),order(colnames(i), decreasing=F)])
traits_pool<-traits[which(rownames(traits) %in% colnames (pool_obs[[1]])),]
pool_obs <- lapply (pool_obs, function (i) i[,which (colnames(i) %in% rownames(traits_pool))])
pool_obs <- lapply (pool_obs, function (i) apply (i, 1, list))

pool_obs <- lapply (seq (1, length (pool_obs)), function (j)
	lapply (seq(1, length (pool_obs[[1]])), function (k)
	pool_obs [[j]][[k]][[1]][which (pool_obs [[j]][[k]][[1]] >0.1)]))

## pegar os atributos das especies compondo o pool
traits_pool <- lapply (seq (1, length (pool_obs)), function (j)
	lapply (seq(1, length (pool_obs[[1]])), function (k)
	traits_pool [which (rownames(traits_pool) %in% names (pool_obs[[j]][[k]])),]))

# calcular o volume observado do pool
cl <- makeCluster(6) ## number of cores = generally ncores -1
clusterEvalQ(cl, library(hypervolume))
# export your data and function
clusterExport(cl, c("pool_obs","traits_pool"))

hypervolume_pool_obs <- parLapply(cl, as.list(seq (1, length (pool_obs))), function (j)#
	lapply (seq(1, length (pool_obs[[1]])), function (k)#
	hypervolume_gaussian (traits_pool[[j]][[k]],kde.bandwidth=estimate_bandwidth(traits_pool[[j]][[k]],method="silverman"))@Volume))

stopCluster(cl)

save (hypervolume_local_obs, hypervolume_pool_obs, file="hypervolume_500traps_observed.RData")


#################### 1,000 traps

local_1000 <- rowSums (data_pred5)
nsamples<-100

########### rodar o modelo nulo que sorteia especies com base na prevalencia no pool
require(parallel)
cl <- makeCluster(6) ## number of cores = generally ncores -1
clusterEvalQ(cl, library(phytools))
clusterEvalQ(cl, library(phylobase))
clusterEvalQ(cl, library(dplyr))
# export your data and function
clusterExport(cl, c("local_1000","nsamples", "site_list_1000", "sampling_prevalence_NM"))

nm_teste <- parLapply(cl,site_list_1000, function (i) sampling_prevalence_NM (i,local=local_1000,nsamples=nsamples))## sample with the size of the obs community

stopCluster(cl)

#### pegar os atributos das especies sorteadas
#traits
list_local_traits <- lapply (as.list(seq (1, length (nm_teste))), function (j)
	lapply (as.list (seq (1,length (local_1000))), function (k) ### para cada sitio
	lapply (as.list (seq(1,nsamples)), function (i)
	traits[which (rownames(traits) %in% nm_teste[[j]]$species[[k]] [,i]),])))

### removendo sitios sem sp ou com uma especie
list_local_traits <- lapply (as.list(seq (1, length (nm_teste))), function (j)
	list_local_traits [[j]] [which (local_1000 >= 4)])
	
### aplicar pra todas os sitios e simulacoes
#require(parallel)
cl <- makeCluster(6) ## number of cores = generally ncores -1
clusterEvalQ(cl, library(hypervolume))
#export your data and function
clusterExport(cl, c("nsamples","nm_teste", "list_local_traits"))

hypervolume_local_null <-parLapply (cl, as.list(seq (1, length (nm_teste))), function (j)# para cada pool
	lapply (as.list(seq(1,length (list_local_traits [[1]]))), function (k)# para cada comunidade
	unlist (lapply (as.list (seq(1,nsamples)), function (i) ## para cada simula??o
	hypervolume_gaussian(list_local_traits [[j]][[k]][[i]],
	kde.bandwidth=estimate_bandwidth(list_local_traits [[j]][[k]][[i]],method="silverman"))@Volume))))

stopCluster (cl)

save(hypervolume_local_null , file="hypervolume_local_null_1000traps.RData")

# OBSERVADOS
### dados de composicao observados das comunidades
local_obs <- data.matrix(data_pred5 [which (local_1000 >= 4),order (colnames(data_pred5), decreasing=F)])
local_obs <- local_obs [,-which (colSums (local_obs)==0)]
list_local_obs <- apply (local_obs, 1, list)
### removendo especies ausentes
list_local_obs <- lapply (as.list (seq(1,length (list_local_obs))), function (i) 
	list_local_obs [[i]] [[1]] [which (list_local_obs [[i]] [[1]] == 1)])

lapply (lapply (list_local_traits_obs, function (i) apply (i, 2, duplicated)), function (k) apply (k,2,table))

### pegar os atributos das especies de cada comunidade
list_local_traits_obs <- lapply (as.list (seq(1,length (list_local_obs))), function (i)
	traits[which (rownames(traits) %in% names (list_local_obs [[i]])),])

## calcular o volume observado
cl <- makeCluster(6) ## number of cores = generally ncores -1
clusterEvalQ(cl, library(hypervolume))
#export your data and function
clusterExport(cl, c("list_local_traits_obs"))

hypervolume_local_obs<- parLapply (cl,list_local_traits_obs, function (i)
	hypervolume_gaussian(i,	kde.bandwidth=estimate_bandwidth(i,method="silverman"))@Volume)
stopCluster (cl)

## composicao do pool
pool_obs <- lapply (site_list_1000, function (i) i [which (local_1000 >=4),order(colnames(i), decreasing=F)])
traits_pool<-traits[which(rownames(traits) %in% colnames (pool_obs[[1]])),]
pool_obs <- lapply (pool_obs, function (i) i[,which (colnames(i) %in% rownames(traits_pool))])
pool_obs <- lapply (pool_obs, function (i) apply (i, 1, list))

pool_obs <- lapply (seq (1, length (pool_obs)), function (j)
	lapply (seq(1, length (pool_obs[[1]])), function (k)
	pool_obs [[j]][[k]][[1]][which (pool_obs [[j]][[k]][[1]] >0.1)]))

## pegar os atributos das especies compondo o pool
traits_pool <- lapply (seq (1, length (pool_obs)), function (j)
	lapply (seq(1, length (pool_obs[[1]])), function (k)
	traits_pool [which (rownames(traits_pool) %in% names (pool_obs[[j]][[k]])),]))

# calcular o volume observado do pool
cl <- makeCluster(6) ## number of cores = generally ncores -1
clusterEvalQ(cl, library(hypervolume))
# export your data and function
clusterExport(cl, c("pool_obs","traits_pool"))

hypervolume_pool_obs <- parLapply(cl, as.list(seq (1, length (pool_obs))), function (j)#
	lapply (seq(1, length (pool_obs[[1]])), function (k)#
	hypervolume_gaussian (traits_pool[[j]][[k]],kde.bandwidth=estimate_bandwidth(traits_pool[[j]][[k]],method="silverman"))@Volume))

stopCluster(cl)

save (hypervolume_local_obs, hypervolume_pool_obs, file="hypervolume_1000traps_observed.RData")








































































## 500 traps

local <- rowSums (data_pred4)
nsamples<-2

require(parallel)
cl <- makeCluster(6) ## number of cores = generally ncores -1
clusterEvalQ(cl, library(phytools))
clusterEvalQ(cl, library(phylobase))
clusterEvalQ(cl, library(dplyr))
# export your data and function
clusterExport(cl, c("local","nsamples", "site_list", "sampling_prevalence_NM"))

nm_teste <- parLapply(cl,site_list, function (i) sampling_prevalence_NM (i,local=local,nsamples=nsamples))## sample with the size of the obs community
#nm_teste <- lapply(site_list, function (i) sampling_prevalence_NM (i,local=local,nsamples=nsamples))## sample with the size of the obs community

stopCluster(cl)

#### list of null communities 
list_null_communities <- lapply (as.list(seq (1, length (nm_teste))), function (j)
	lapply (as.list (seq (1,length (local))), function (k) ### para cada sitio
	lapply (as.list (seq(1,nsamples)), function (i)
	data.frame (null_sp=nm_teste[[j]]$species[[k]][,i], pres=as.numeric(rep(1, length(nm_teste[[j]]$species[[k]][,i])))))))

list_null_communities <- lapply (as.list(seq (1, length (nm_teste))), function (j)
	lapply (as.list (seq (1,length (local))), function (k) ### para cada sitio
	lapply (as.list (seq(1,nsamples)), function (i)
	as.matrix(t(list_null_communities[[j]][[k]][[i]]))))) ## fazer uma matriz

list_null_communities <- lapply (as.list(seq (1, length (nm_teste))), function (j)
	lapply (as.list (seq (1,length (local))), function (k) ### para cada sitio
	lapply (as.list (seq(1,nsamples)), function (i)
	matrix (list_null_communities[[j]][[k]][[i]],byrow=F,nrow=2,dimnames=list(NULL, list_null_communities[[j]][[k]][[i]][1,])))))

list_null_communities <- lapply (as.list(seq (1, length (nm_teste))), function (j)
	lapply (as.list (seq (1,length (local))), function (k) ### para cada sitio
	lapply (as.list (seq(1,nsamples)), function (i)
	list_null_communities [[j]][[k]][[i]][-1,])))

list_null_communities <- lapply (as.list(seq (1, length (nm_teste))), function (j)
	lapply (as.list (seq (1,length (local))), function (k) ### para cada sitio
	lapply (as.list (seq(1,nsamples)), function (i)
	as.matrix(t(list_null_communities[[j]][[k]][[i]])))))

list_null_communities <- lapply (as.list(seq (1, length (nm_teste))), function (j)
	lapply (as.list (seq (1,length (local))), function (k) ### para cada sitio
	lapply (as.list (seq(1,nsamples)), function (i)
	list_null_communities[[j]][[k]][[i]] [,order (as.character(colnames(list_null_communities[[j]][[k]][[i]])), decreasing=F)])))

list_null_communities <- lapply (as.list(seq (1, length (nm_teste))), function (j)
	lapply (as.list (seq (1,length (local))), function (k) ### para cada sitio
	lapply (as.list (seq(1,nsamples)), function (i)
	as.matrix(t(list_null_communities[[j]][[k]][[i]])))))

list_null_communities <- lapply (as.list(seq (1, length (nm_teste))), function (j)
	lapply (as.list (seq (1,length (local))), function (k) ### para cada sitio
	lapply (as.list (seq(1,nsamples)), function (i)
	{class(list_null_communities[[j]][[k]][[i]])<- "numeric";list_null_communities[[j]][[k]][[i]]})))

#traits
list_local_traits <- lapply (as.list(seq (1, length (nm_teste))), function (j)
	lapply (as.list (seq (1,length (local))), function (k) ### para cada sitio
	lapply (as.list (seq(1,nsamples)), function (i)
	traits[which (rownames(traits) %in% colnames(list_null_communities[[j]][[k]][[i]])),])))

### removendo sitios sem sp ou com uma especie
list_null_communities <- lapply (as.list(seq (1, length (nm_teste))), function (j)
	list_null_communities[[j]] [which (local >= 4)])

list_local_traits <- lapply (as.list(seq (1, length (nm_teste))), function (j)
	list_local_traits [[j]] [which (local >= 3)])
	
list_null_communities <- lapply (as.list(seq (1, length (nm_teste))), function (j)
	lapply (as.list (seq (1,length (local[which (local >= 3)]))), function (k) ### para cada sitio
	lapply (as.list (seq(1,nsamples)), function (i)
	list_null_communities [[j]][[k]][[i]][,which (colnames(list_null_communities [[j]][[k]][[i]]) %in% rownames(list_local_traits [[j]][[k]][[i]]))])))

### conferindo se tudo esta com as mesmas dimensoes
unlist(lapply (as.list(seq (1, length (nm_teste))), function (j)
	lapply (as.list (seq (1,length (local[which (local >= 3)]))), function (k) ### para cada sitio
	lapply (as.list (seq(1,nsamples)), function (i)
	dim(list_local_traits[[j]][[k]][[i]])[1]  == length (list_null_communities[[j]][[k]][[i]])))))

####
############# traits for species of the pool
list_pool_traits <- lapply (as.list(seq (1, length (nm_teste))), function (j)
	lapply (as.list (seq (1,length (local))), function (k) ### para cada sitio
	#lapply (as.list (seq(1,nsamples)), function (i)
	traits[which (rownames (traits) %in% colnames (site_list[[j]][k,][which (site_list[[j]][k,] > 0.1)])),]))

### species in the sites matching with species with traits
list_pool_sites <- lapply (as.list(seq (1, length (nm_teste))), function (j)
	lapply (as.list (seq (1,length (local))), function (k) ### para cada sitio
	site_list[[j]][k,][which (site_list[[j]][k,] > 0.1)]))

list_pool_sites <- lapply (as.list(seq (1, length (nm_teste))), function (j)
	lapply (as.list (seq (1,length (local))), function (k) ### para cada sitio
	list_pool_sites [[j]][[k]][,which (colnames(list_pool_sites [[j]][[k]]) %in% rownames (list_pool_traits [[j]][[k]]))]))

list_pool_sites <- lapply (as.list(seq (1, length (nm_teste))), function (j)
	lapply (as.list (seq (1,length (local))), function (k) ### para cada sitio
	list_pool_sites [[j]][[k]] [,order (colnames (list_pool_sites [[j]][[k]][which (list_pool_sites [[j]][[k]] > 0.1)]),decreasing=F)]))

### removendo sitios sem sp ou com uma especie
list_pool_sites <- lapply (as.list(seq (1, length (nm_teste))), function (j)
	list_pool_sites[[j]] [which (local >= 3)])

list_pool_traits <- lapply (as.list(seq (1, length (nm_teste))), function (j)
	list_pool_traits[[j]] [which (local >= 3)])

## conferindo se o numero de especies ? similar
unlist(lapply (list_pool_traits[[1]],function (i) dim(i)[1]))==unlist(lapply (list_pool_sites[[1]],function (i) dim(i)[2]))

### aplicar pra todas os sitios e simulacoes
require(parallel)
cl <- makeCluster(6) ## number of cores = generally ncores -1
clusterEvalQ(cl, library(FD))
# export your data and function
clusterExport(cl, c("local","nsamples","nm_teste", "traits", "list_pool_traits","list_pool_sites","list_null_communities","list_local_traits"))

resultados_500 <- parLapply(cl, as.list(seq (1, length (nm_teste))), function (j)
	lapply (as.list (seq (1,length (local [which (local >= 3)]))), function (k) ### para cada sitio
	unlist (lapply (as.list (seq(1,nsamples)), function (i) ## para cada simula??o
	dbFD (x=list_pool_traits[[j]][[k]], a=list_pool_sites [[j]][[k]],stand.x=T,calc.FRic = TRUE,calc.FGR=F,calc.CWM=F,calc.FDiv=F,corr="cailliez",
	w.abun=F, m = "max",messages=T)$FRic - dbFD (x= list_local_traits[[j]][[k]][[i]],
	a=list_null_communities[[j]][[k]][[i]],stand.x=T,calc.FRic = TRUE,calc.FGR=F,calc.CWM=F,calc.FDiv=F, w.abun=F,m = "max",corr="cailliez",messages=T)$FRic))))

stopCluster (cl)
save(resultados_500, file="resultados_rich_500_func.RData")
#resultados_500<- readRDS ("resultados_rich_500_func.rda")

# OBSERVADOS
local_obs <- data.matrix(data_pred4 [which (local >= 4),order (colnames(data_pred4), decreasing=F)])
local_obs <- local_obs [,-which (colSums (local_obs)==0)]
list_local_obs <- apply (local_obs, 1, list)

### FRic observado na escala local
FD_local<-dbFD(traits[which (rownames(traits) %in% colnames (local_obs)),c(33,36,21,27,29)], 
	local_obs, stand.x=T,calc.FRic = TRUE,calc.FGR=F,calc.CWM=F,calc.FDiv=F,corr="cailliez",
	w.abun=F, m = "max",messages=T)$FRic

pool_obs <- lapply (site_list, function (i) i [which (local >=4),order(colnames(i), decreasing=F)])
traits_pool<-traits[which(rownames(traits) %in% colnames (pool_obs[[1]])),c(33,36,21,27,29)]
pool_obs <- lapply (pool_obs, function (i) i[,which (colnames(i) %in% rownames(traits_pool))])
pool_obs <- lapply (pool_obs, function (i) apply (i, 1, list))

pool_obs <- lapply (seq (1, length (pool_obs)), function (j)
	lapply (seq(1, length (pool_obs[[1]])), function (k)
	pool_obs [[j]][[k]][[1]][which (pool_obs [[j]][[k]][[1]] >0.1)]))

traits_pool <- lapply (seq (1, length (pool_obs)), function (j)
	lapply (seq(1, length (pool_obs[[1]])), function (k)
	traits_pool [which (rownames(traits_pool) %in% names (pool_obs[[j]][[k]])),]))

#require(parallel)
cl <- makeCluster(6) ## number of cores = generally ncores -1
clusterEvalQ(cl, library(FD))
# export your data and function
clusterExport(cl, c("pool_obs","traits_pool"))
 
### FRic observada para o pool
FD_pool<- parLapply (cl, seq (1, length (pool_obs)), function (j)
	lapply (seq(1, length (pool_obs[[1]])), function (k)
	dbFD (traits_pool[[j]][[k]],pool_obs[[j]][[k]], stand.x=T,calc.FRic = TRUE,calc.FGR=F,calc.CWM=F,calc.FDiv=F,corr="cailliez",
	w.abun=F, m = "max",messages=T)$FRic))

stopCluster(cl)

save (FD_pool, FD_local, file= "FRIC_obs_500traps_rich.RData")

### 1000 traps

local_1000 <- rowSums (data_pred5)
nsamples<-100

require(parallel)
cl <- makeCluster(6) ## number of cores = generally ncores -1
clusterEvalQ(cl, library(phytools))
clusterEvalQ(cl, library(phylobase))
clusterEvalQ(cl, library(dplyr))
# export your data and function
clusterExport(cl, c("local_1000","nsamples", "site_list_1000", "sampling_prevalence_NM"))

nm_teste <- parLapply(cl,site_list_1000, function (i) sampling_prevalence_NM (i,local=local_1000,nsamples=nsamples))## sample with the size of the obs community
#nm_teste <- lapply(site_list, function (i) sampling_prevalence_NM (i,local=local,nsamples=nsamples))## sample with the size of the obs community

stopCluster(cl)

#### list of null communities 
list_null_communities <- lapply (as.list(seq (1, length (nm_teste))), function (j)
	lapply (as.list (seq (1,length (local_1000))), function (k) ### para cada sitio
	lapply (as.list (seq(1,nsamples)), function (i)
	data.frame (null_sp=nm_teste[[j]]$species[[k]][,i], pres=as.numeric(rep(1, length(nm_teste[[j]]$species[[k]][,i])))))))

list_null_communities <- lapply (as.list(seq (1, length (nm_teste))), function (j)
	lapply (as.list (seq (1,length (local_1000))), function (k) ### para cada sitio
	lapply (as.list (seq(1,nsamples)), function (i)
	as.matrix(t(list_null_communities[[j]][[k]][[i]]))))) ## fazer uma matriz

list_null_communities <- lapply (as.list(seq (1, length (nm_teste))), function (j)
	lapply (as.list (seq (1,length (local_1000))), function (k) ### para cada sitio
	lapply (as.list (seq(1,nsamples)), function (i)
	matrix (list_null_communities[[j]][[k]][[i]],byrow=F,nrow=2,dimnames=list(NULL, list_null_communities[[j]][[k]][[i]][1,])))))

list_null_communities <- lapply (as.list(seq (1, length (nm_teste))), function (j)
	lapply (as.list (seq (1,length (local_1000))), function (k) ### para cada sitio
	lapply (as.list (seq(1,nsamples)), function (i)
	list_null_communities [[j]][[k]][[i]][-1,])))

list_null_communities <- lapply (as.list(seq (1, length (nm_teste))), function (j)
	lapply (as.list (seq (1,length (local_1000))), function (k) ### para cada sitio
	lapply (as.list (seq(1,nsamples)), function (i)
	as.matrix(t(list_null_communities[[j]][[k]][[i]])))))

list_null_communities <- lapply (as.list(seq (1, length (nm_teste))), function (j)
	lapply (as.list (seq (1,length (local_1000))), function (k) ### para cada sitio
	lapply (as.list (seq(1,nsamples)), function (i)
	list_null_communities[[j]][[k]][[i]] [,order (as.character(colnames(list_null_communities[[j]][[k]][[i]])), decreasing=F)])))

list_null_communities <- lapply (as.list(seq (1, length (nm_teste))), function (j)
	lapply (as.list (seq (1,length (local_1000))), function (k) ### para cada sitio
	lapply (as.list (seq(1,nsamples)), function (i)
	as.matrix(t(list_null_communities[[j]][[k]][[i]])))))

list_null_communities <- lapply (as.list(seq (1, length (nm_teste))), function (j)
	lapply (as.list (seq (1,length (local_1000))), function (k) ### para cada sitio
	lapply (as.list (seq(1,nsamples)), function (i)
	{class(list_null_communities[[j]][[k]][[i]])<- "numeric";list_null_communities[[j]][[k]][[i]]})))

#traits
list_local_traits <- lapply (as.list(seq (1, length (nm_teste))), function (j)
	lapply (as.list (seq (1,length (local_1000))), function (k) ### para cada sitio
	lapply (as.list (seq(1,nsamples)), function (i)
	traits[which (rownames(traits) %in% colnames(list_null_communities[[j]][[k]][[i]])),c(33,36,21,27,29)])))

### removendo sitios sem sp ou com uma especie
list_null_communities <- lapply (as.list(seq (1, length (nm_teste))), function (j)
	list_null_communities[[j]] [which (local_1000 >= 4)])

list_local_traits <- lapply (as.list(seq (1, length (nm_teste))), function (j)
	list_local_traits [[j]] [which (local_1000 >= 4)])
	
dist_local<-lapply (as.list(seq (1, length (nm_teste))), function (j)
	lapply (as.list (seq (1,length (local_1000[which (local_1000 >= 4)]))), function (k) ### para cada sitio
	lapply (as.list (seq(1,nsamples)), function (i)
	vegdist (list_local_traits[[j]][[k]][[i]],"gower"))))

list_null_communities <- lapply (as.list(seq (1, length (nm_teste))), function (j)
	lapply (as.list (seq (1,length (local_1000[which (local_1000 >= 4)]))), function (k) ### para cada sitio
	lapply (as.list (seq(1,nsamples)), function (i)
	list_null_communities [[j]][[k]][[i]][,which (colnames(list_null_communities [[j]][[k]][[i]]) %in% rownames(list_local_traits [[j]][[k]][[i]]))])))


### conferindo se tudo esta com as mesmas dimensoes
unlist(lapply (as.list(seq (1, length (nm_teste))), function (j)
	lapply (as.list (seq (1,length (local_1000[which (local_1000 >= 4)]))), function (k) ### para cada sitio
	lapply (as.list (seq(1,nsamples)), function (i)
	dim(list_local_traits[[j]][[k]][[i]])[1]  == length (list_null_communities[[j]][[k]][[i]])))))

####

############# traits for species of the pool
list_pool_traits <- lapply (as.list(seq (1, length (nm_teste))), function (j)
	lapply (as.list (seq (1,length (local_1000))), function (k) ### para cada sitio
	#lapply (as.list (seq(1,nsamples)), function (i)
	traits[which (rownames (traits) %in% colnames (site_list_1000[[j]][k,][which (site_list_1000[[j]][k,] > 0.1)])),c(33,36,21,27,29)]))


### species in the sites matching with species with traits
list_pool_sites <- lapply (as.list(seq (1, length (nm_teste))), function (j)
	lapply (as.list (seq (1,length (local_1000))), function (k) ### para cada sitio
	site_list_1000[[j]][k,][which (site_list_1000[[j]][k,] > 0.1)]))

list_pool_sites <- lapply (as.list(seq (1, length (nm_teste))), function (j)
	lapply (as.list (seq (1,length (local_1000))), function (k) ### para cada sitio
	list_pool_sites [[j]][[k]][,which (colnames(list_pool_sites [[j]][[k]]) %in% rownames (list_pool_traits [[j]][[k]]))]))

list_pool_sites <- lapply (as.list(seq (1, length (nm_teste))), function (j)
	lapply (as.list (seq (1,length (local_1000))), function (k) ### para cada sitio
	list_pool_sites [[j]][[k]] [,order (colnames (list_pool_sites [[j]][[k]][which (list_pool_sites [[j]][[k]] > 0.1)]),decreasing=F)]))

### removendo sitios sem sp ou com uma especie
list_pool_sites <- lapply (as.list(seq (1, length (nm_teste))), function (j)
	list_pool_sites[[j]] [which (local_1000 >= 4)])

list_pool_traits <- lapply (as.list(seq (1, length (nm_teste))), function (j)
	list_pool_traits[[j]] [which (local_1000 >= 4)])

dist_pool<-lapply (as.list(seq (1, length (nm_teste))), function (j)
	lapply (as.list (seq (1,length (local_1000[which (local_1000 >= 4)]))), function (k) ### para cada sitio
	#lapply (as.list (seq(1,nsamples)), function (i)
	vegdist (list_pool_traits[[j]][[k]],"gower")))

## conferindo se o numero de especies ? similar
unlist(lapply (list_pool_traits[[1]],function (i) dim(i)[1]))==unlist(lapply (list_pool_sites[[1]],function (i) dim(i)[2]))

### aplicar pra todas os sitios e simulacoes
#require(parallel)
cl <- makeCluster(6) ## number of cores = generally ncores -1
clusterEvalQ(cl, library(FD))
#export your data and function
clusterExport(cl, c("local_1000","nsamples","nm_teste", "traits","list_pool_sites","list_null_communities","list_pool_traits","list_local_traits"))

resultados_1000_local <- parLapply (cl, as.list(seq (1, length (nm_teste))), function (j)
	lapply (as.list (seq (1,length (local_1000 [which (local_1000 >= 4)]))), function (k) ### para cada sitio
	unlist (lapply (as.list (seq(1,nsamples)), function (i) ## para cada simula??o
	dbFD (x= list_local_traits[[j]][[k]][[i]],
	a=list_null_communities[[j]][[k]][[i]],stand.x=T,calc.FRic = TRUE,calc.FGR=F,calc.CWM=F,calc.FDiv=F, w.abun=F,m = "max",corr="cailliez",messages=T)$FRic))))
	
stopCluster (cl)

cl <- makeCluster(6) ## number of cores = generally ncores -1
clusterEvalQ(cl, library(FD))
#export your data and function
clusterExport(cl, c("local_1000","nsamples","nm_teste", "traits","list_pool_sites","list_null_communities","list_pool_traits","list_local_traits"))

resultados_1000_regional<-parLapply (cl, as.list(seq (1, length (nm_teste))), function (j)
	lapply (as.list (seq (1,length (local_1000 [which (local_1000 >= 4)]))), function (k) ### para cada sitio
	dbFD (x=list_pool_traits[[j]][[k]], a=list_pool_sites [[j]][[k]],stand.x=T,calc.FRic = TRUE,calc.FGR=F,calc.CWM=F,calc.FDiv=F, w.abun=F,m = "max",corr="cailliez",messages=T)$FRic))#
		
stopCluster (cl)

save (resultados_1000_local, resultados_1000_regional, file= "FRIC_1000traps_NULL.RData")

# OBSERVADOS
local_obs <- data.matrix(data_pred5 [which (local_1000 >= 4),order (colnames(data_pred5), decreasing=F)])
local_obs <- local_obs [,-which (colSums (local_obs)==0)]
list_local_obs <- apply (local_obs, 1, list)

FD_local_1000<-dbFD(traits[which (rownames(traits) %in% colnames (local_obs)),c(33,36,21,27,29)], 
	local_obs, stand.x=T,calc.FRic = TRUE,calc.FGR=F,calc.CWM=F,calc.FDiv=F,corr="cailliez",
	w.abun=F, m = "max",messages=T)$FRic

pool_obs <- lapply (site_list_1000, function (i) i [which (local_1000 >=4),order(colnames(i), decreasing=F)])
traits_pool<-traits[which(rownames(traits) %in% colnames (pool_obs[[1]])),c(33,36,21,27,29)]
pool_obs <- lapply (pool_obs, function (i) i[,which (colnames(i) %in% rownames(traits_pool))])
pool_obs <- lapply (pool_obs, function (i) apply (i, 1, list))

pool_obs <- lapply (seq (1, length (pool_obs)), function (j)
	lapply (seq(1, length (pool_obs[[1]])), function (k)
	pool_obs [[j]][[k]][[1]][which (pool_obs [[j]][[k]][[1]] >0.1)]))

traits_pool <- lapply (seq (1, length (pool_obs)), function (j)
	lapply (seq(1, length (pool_obs[[1]])), function (k)
	traits_pool [which (rownames(traits_pool) %in% names (pool_obs[[j]][[k]])),]))

cl <- makeCluster(6) ## number of cores = generally ncores -1
clusterEvalQ(cl, library(FD))
# export your data and function
clusterExport(cl, c("pool_obs","traits_pool"))

FD_pool_1000<- parLapply(cl, seq (1, length (pool_obs)), function (j)
	lapply (seq(1, length (pool_obs[[1]])), function (k)
	dbFD (traits_pool[[j]][[k]],pool_obs[[j]] [[k]], stand.x=T,calc.FRic = TRUE,calc.FGR=F,calc.CWM=F,calc.FDiv=F,corr="cailliez",
	w.abun=F, m = "max",messages=T)$FRic))

stopCluster(cl)

save (FD_local_1000, FD_pool_1000, file="FRIC_1000traps_observed.RData")


###############################################
####### FUNCTIONAL COMPOSITION  ###############
###############################################

## selecionar florestais
for_sp_composition_pool<-lapply  (site_list, function (i) i [colnames (i) %in% iucn_hab [iucn_hab$forest >= 1,]$sp])
for_sp_composition_local <- data_pred4 [, which (colnames (data_pred4) %in% iucn_hab [iucn_hab$forest >= 1,]$sp)]

## funcao modelo nulo
nsamples <- 1
local <- rowSums (for_sp_composition_local)

require(parallel)
cl <- makeCluster(6) ## number of cores = generally ncores -1
clusterEvalQ(cl, library(phytools))
clusterEvalQ(cl, library(phylobase))
clusterEvalQ(cl, library(dplyr))
# export your data and function
clusterExport(cl, c("local","nsamples", "for_sp_composition_pool", "sampling_prevalence_NM"))

nm_teste <- parLapply(cl,site_list, function (i) sampling_prevalence_NM (i,local=local,nsamples=nsamples))## sample with the size of the obs community
#nm_teste <- lapply(site_list, function (i) sampling_prevalence_NM (i,local=local,nsamples=nsamples))## sample with the size of the obs community

stopCluster(cl)

#### list of null communities 
list_null_communities <- lapply (as.list(seq (1, length (nm_teste))), function (j)
	lapply (as.list (seq (1,length (local))), function (k) ### para cada sitio
	lapply (as.list (seq(1,nsamples)), function (i)
	data.frame (null_sp=nm_teste[[j]]$species[[k]][,i], pres=as.numeric(rep(1, length(nm_teste[[j]]$species[[k]][,i])))))))

list_null_communities <- lapply (as.list(seq (1, length (nm_teste))), function (j)
	lapply (as.list (seq (1,length (local))), function (k) ### para cada sitio
	lapply (as.list (seq(1,nsamples)), function (i)
	as.matrix(t(list_null_communities[[j]][[k]][[i]]))))) ## fazer uma matriz

list_null_communities <- lapply (as.list(seq (1, length (nm_teste))), function (j)
	lapply (as.list (seq (1,length (local))), function (k) ### para cada sitio
	lapply (as.list (seq(1,nsamples)), function (i)
	matrix (list_null_communities[[j]][[k]][[i]],byrow=F,nrow=2,dimnames=list(NULL, list_null_communities[[j]][[k]][[i]][1,])))))

list_null_communities <- lapply (as.list(seq (1, length (nm_teste))), function (j)
	lapply (as.list (seq (1,length (local))), function (k) ### para cada sitio
	lapply (as.list (seq(1,nsamples)), function (i)
	list_null_communities [[j]][[k]][[i]][-1,])))

list_null_communities <- lapply (as.list(seq (1, length (nm_teste))), function (j)
	lapply (as.list (seq (1,length (local))), function (k) ### para cada sitio
	lapply (as.list (seq(1,nsamples)), function (i)
	as.matrix(t(list_null_communities[[j]][[k]][[i]])))))

list_null_communities <- lapply (as.list(seq (1, length (nm_teste))), function (j)
	lapply (as.list (seq (1,length (local))), function (k) ### para cada sitio
	lapply (as.list (seq(1,nsamples)), function (i)
	list_null_communities[[j]][[k]][[i]] [,order (as.character(colnames(list_null_communities[[j]][[k]][[i]])), decreasing=F)])))

list_null_communities <- lapply (as.list(seq (1, length (nm_teste))), function (j)
	lapply (as.list (seq (1,length (local))), function (k) ### para cada sitio
	lapply (as.list (seq(1,nsamples)), function (i)
	as.matrix(t(list_null_communities[[j]][[k]][[i]])))))

list_null_communities <- lapply (as.list(seq (1, length (nm_teste))), function (j)
	lapply (as.list (seq (1,length (local))), function (k) ### para cada sitio
	lapply (as.list (seq(1,nsamples)), function (i)
	{class(list_null_communities[[j]][[k]][[i]])<- "numeric";list_null_communities[[j]][[k]][[i]]})))

#traits
list_local_traits <- lapply (as.list(seq (1, length (nm_teste))), function (j)
	lapply (as.list (seq (1,length (local))), function (k) ### para cada sitio
	lapply (as.list (seq(1,nsamples)), function (i)
	traits[which (rownames(traits) %in% colnames(list_null_communities[[j]][[k]][[i]])),c(33,36,21,27,29)])))

### removendo sitios sem sp ou com uma especie
list_null_communities <- lapply (as.list(seq (1, length (nm_teste))), function (j)
	list_null_communities[[j]] [which (local >= 4)])

list_local_traits <- lapply (as.list(seq (1, length (nm_teste))), function (j)
	list_local_traits [[j]] [which (local >= 4)])
	
list_null_communities <- lapply (as.list(seq (1, length (nm_teste))), function (j)
	lapply (as.list (seq (1,length (local[which (local >= 4)]))), function (k) ### para cada sitio
	lapply (as.list (seq(1,nsamples)), function (i)
	list_null_communities [[j]][[k]][[i]][,which (colnames(list_null_communities [[j]][[k]][[i]]) %in% rownames(list_local_traits [[j]][[k]][[i]]))])))

### conferindo se tudo esta com as mesmas dimensoes
unlist(lapply (as.list(seq (1, length (nm_teste))), function (j)
	lapply (as.list (seq (1,length (local[which (local >= 4)]))), function (k) ### para cada sitio
	lapply (as.list (seq(1,nsamples)), function (i)
	dim(list_local_traits[[j]][[k]][[i]])[1]  == length (list_null_communities[[j]][[k]][[i]])))))

####
############# traits for species of the pool
list_pool_traits <- lapply (as.list(seq (1, length (nm_teste))), function (j)
	lapply (as.list (seq (1,length (local))), function (k) ### para cada sitio
	#lapply (as.list (seq(1,nsamples)), function (i)
	traits[which (rownames (traits) %in% colnames (site_list[[j]][k,][which (site_list[[j]][k,] > 0.1)])),c(33,36,21,27,29)]))

### species in the sites matching with species with traits
list_pool_sites <- lapply (as.list(seq (1, length (nm_teste))), function (j)
	lapply (as.list (seq (1,length (local))), function (k) ### para cada sitio
	site_list[[j]][k,][which (site_list[[j]][k,] > 0.1)]))

list_pool_sites <- lapply (as.list(seq (1, length (nm_teste))), function (j)
	lapply (as.list (seq (1,length (local))), function (k) ### para cada sitio
	list_pool_sites [[j]][[k]][,which (colnames(list_pool_sites [[j]][[k]]) %in% rownames (list_pool_traits [[j]][[k]]))]))

list_pool_sites <- lapply (as.list(seq (1, length (nm_teste))), function (j)
	lapply (as.list (seq (1,length (local))), function (k) ### para cada sitio
	list_pool_sites [[j]][[k]] [,order (colnames (list_pool_sites [[j]][[k]][which (list_pool_sites [[j]][[k]] > 0.1)]),decreasing=F)]))

### removendo sitios sem sp ou com uma especie
list_pool_sites <- lapply (as.list(seq (1, length (nm_teste))), function (j)
	list_pool_sites[[j]] [which (local >= 4)])

list_pool_traits <- lapply (as.list(seq (1, length (nm_teste))), function (j)
	list_pool_traits[[j]] [which (local >= 4)])

## conferindo se o numero de especies ? similar
unlist(lapply (list_pool_traits[[1]],function (i) dim(i)[1]))==unlist(lapply (list_pool_sites[[1]],function (i) dim(i)[2]))

### aplicar pra todas os sitios e simulacoes
#require(parallel)
#cl <- makeCluster(6) ## number of cores = generally ncores -1
#clusterEvalQ(cl, library(FD))
# export your data and function
#clusterExport(cl, c("local","nsamples","nm_teste", "traits", "list_pool_traits","list_pool_sites","list_null_communities","list_local_traits"))

#resultados_500 <- parLapply(cl, 

resultados_500_local<- lapply (as.list(seq (1, length (nm_teste))), function (j)
	lapply (as.list (seq (1,length (local [which (local >= 4)]))), function (k) ### para cada sitio
	unlist (lapply (as.list (seq(1,nsamples)), function (i) ## para cada simula??o
	dbFD (x= list_local_traits[[j]][[k]][[i]],
	a=list_null_communities[[j]][[k]][[i]],stand.x=T,calc.FRic = TRUE,calc.FGR=F,calc.CWM=F,calc.FDiv=F, 
	w.abun=F,m = "max",corr="cailliez",messages=T)$FRic))))
	

resultados_500_regional<-lapply (as.list(seq (1, length (nm_teste))), function (j)
	lapply (as.list (seq (1,length (local [which (local >= 4)]))), function (k) ### para cada sitio
	dbFD (x=list_pool_traits[[j]][[k]], a=list_pool_sites [[j]][[k]],stand.x=T,calc.FRic = TRUE,calc.FGR=F,calc.CWM=F,
	calc.FDiv=F, w.abun=F,m = "max",corr="cailliez",messages=T)$FRic))#
		
#stopCluster (cl)
#saveRDS (resultados_500, "resultados_rich_500_func.rda")
resultados_500<- readRDS ("resultados_rich_500_func.rda")

# OBSERVADOS
local_obs <- data.matrix(data_pred4 [which (local >= 4),order (colnames(data_pred4), decreasing=F)])
local_obs <- local_obs [,-which (colSums (local_obs)==0)]
list_local_obs <- apply (local_obs, 1, list)

### FRic observado na escala local
FD_local<-dbFD(traits[which (rownames(traits) %in% colnames (local_obs)),c(33,36,21,27,29)], 
	local_obs, stand.x=T,calc.FRic = TRUE,calc.FGR=F,calc.CWM=F,calc.FDiv=F,corr="cailliez",
	w.abun=F, m = "max",messages=T)$FRic

pool_obs <- lapply (site_list, function (i) i [which (local >=4),order(colnames(i), decreasing=F)])
traits_pool<-traits[which(rownames(traits) %in% colnames (pool_obs[[1]])),c(33,36,21,27,29)]
pool_obs <- lapply (pool_obs, function (i) i[,which (colnames(i) %in% rownames(traits_pool))])
pool_obs <- lapply (pool_obs, function (i) apply (i, 1, list))

pool_obs <- lapply (seq (1, length (pool_obs)), function (j)
	lapply (seq(1, length (pool_obs[[1]])), function (k)
	pool_obs [[j]][[k]][[1]][which (pool_obs [[j]][[k]][[1]] >0.1)]))

traits_pool <- lapply (seq (1, length (pool_obs)), function (j)
	lapply (seq(1, length (pool_obs[[1]])), function (k)
	traits_pool [which (rownames(traits_pool) %in% names (pool_obs[[j]][[k]])),]))

#require(parallel)
#cl <- makeCluster(6) ## number of cores = generally ncores -1
#clusterEvalQ(cl, library(FD))
# export your data and function
#clusterExport(cl, c("pool_obs","traits_pool"))
#parLapply (cl, 

### FRic observada para o pool
FD_pool<- lapply (seq (1, length (pool_obs)), function (j)
	lapply (seq(1, length (pool_obs[[1]])), function (k)
	dbFD (traits_pool[[j]][[k]],pool_obs[[j]][[k]], stand.x=T,calc.FRic = TRUE,calc.FGR=F,calc.CWM=F,calc.FDiv=F,corr="cailliez",
	w.abun=F, m = "max",messages=T)$FRic))

#stopCluster(cl)

save(FD_pool, FD_local, file="composition_FD_OBS_500_forest.RData")

################# 
## selecionar florestais 1000 traps
for_sp_composition_pool_1000<-lapply  (site_list_1000, function (i) i [colnames (i) %in% iucn_hab [iucn_hab$forest >= 1,]$sp])
for_sp_composition_local_1000 <- data_pred5 [, which (colnames (data_pred5) %in% iucn_hab [iucn_hab$forest >= 1,]$sp)]

### 1000 traps
nsamples <- 1
local_1000 <- rowSums (for_sp_composition_local_1000)

#cl <- makeCluster(6) ## number of cores = generally ncores -1
#clusterEvalQ(cl, library(phytools))
#clusterEvalQ(cl, library(phylobase))
#clusterEvalQ(cl, library(dplyr))
# export your data and function
#clusterExport(cl, c("local_1000","nsamples", "for_sp_composition_pool_1000", "sampling_prevalence_NM"))
#parLapply (cl,
nm_teste_1000<-lapply(for_sp_composition_pool_1000, function (i) sampling_prevalence_NM (i,local=local_1000,nsamples=nsamples))## sample with the size of the obs community
#stopCluster(cl)


#### list of null communities 
list_null_communities <- lapply (as.list(seq (1, length (nm_teste_1000))), function (j)
	lapply (as.list (seq (1,length (local_1000))), function (k) ### para cada sitio
	lapply (as.list (seq(1,nsamples)), function (i)
	data.frame (null_sp=nm_teste_1000[[j]]$species[[k]][,i], pres=as.numeric(rep(1, length(nm_teste_1000[[j]]$species[[k]][,i])))))))

list_null_communities <- lapply (as.list(seq (1, length (nm_teste_1000))), function (j)
	lapply (as.list (seq (1,length (local_1000))), function (k) ### para cada sitio
	lapply (as.list (seq(1,nsamples)), function (i)
	as.matrix(t(list_null_communities[[j]][[k]][[i]]))))) ## fazer uma matriz

list_null_communities <- lapply (as.list(seq (1, length (nm_teste_1000))), function (j)
	lapply (as.list (seq (1,length (local_1000))), function (k) ### para cada sitio
	lapply (as.list (seq(1,nsamples)), function (i)
	matrix (list_null_communities[[j]][[k]][[i]],byrow=F,nrow=2,dimnames=list(NULL, list_null_communities[[j]][[k]][[i]][1,])))))

list_null_communities <- lapply (as.list(seq (1, length (nm_teste_1000))), function (j)
	lapply (as.list (seq (1,length (local_1000))), function (k) ### para cada sitio
	lapply (as.list (seq(1,nsamples)), function (i)
	list_null_communities [[j]][[k]][[i]][-1,])))

list_null_communities <- lapply (as.list(seq (1, length (nm_teste_1000))), function (j)
	lapply (as.list (seq (1,length (local_1000))), function (k) ### para cada sitio
	lapply (as.list (seq(1,nsamples)), function (i)
	as.matrix(t(list_null_communities[[j]][[k]][[i]])))))

list_null_communities <- lapply (as.list(seq (1, length (nm_teste_1000))), function (j)
	lapply (as.list (seq (1,length (local_1000))), function (k) ### para cada sitio
	lapply (as.list (seq(1,nsamples)), function (i)
	list_null_communities[[j]][[k]][[i]] [,order (as.character(colnames(list_null_communities[[j]][[k]][[i]])), decreasing=F)])))

list_null_communities <- lapply (as.list(seq (1, length (nm_teste_1000))), function (j)
	lapply (as.list (seq (1,length (local_1000))), function (k) ### para cada sitio
	lapply (as.list (seq(1,nsamples)), function (i)
	as.matrix(t(list_null_communities[[j]][[k]][[i]])))))

list_null_communities <- lapply (as.list(seq (1, length (nm_teste_1000))), function (j)
	lapply (as.list (seq (1,length (local_1000))), function (k) ### para cada sitio
	lapply (as.list (seq(1,nsamples)), function (i)
	{class(list_null_communities[[j]][[k]][[i]])<- "numeric";list_null_communities[[j]][[k]][[i]]})))

#traits
list_local_traits <- lapply (as.list(seq (1, length (nm_teste_1000))), function (j)
	lapply (as.list (seq (1,length (local_1000))), function (k) ### para cada sitio
	lapply (as.list (seq(1,nsamples)), function (i)
	traits[which (rownames(traits) %in% colnames(list_null_communities[[j]][[k]][[i]])),c(33,36,21,27,29)])))

### removendo sitios sem sp ou com uma especie
list_null_communities <- lapply (as.list(seq (1, length (nm_teste_1000))), function (j)
	list_null_communities[[j]] [which (local_1000 >= 4)])

list_local_traits <- lapply (as.list(seq (1, length (nm_teste_1000))), function (j)
	list_local_traits [[j]] [which (local_1000 >= 4)])
	
list_null_communities <- lapply (as.list(seq (1, length (nm_teste_1000))), function (j)
	lapply (as.list (seq (1,length (local[which (local_1000 >= 4)]))), function (k) ### para cada sitio
	lapply (as.list (seq(1,nsamples)), function (i)
	list_null_communities [[j]][[k]][[i]][,which (colnames(list_null_communities [[j]][[k]][[i]]) %in% rownames(list_local_traits [[j]][[k]][[i]]))])))

### conferindo se tudo esta com as mesmas dimensoes
unlist(lapply (as.list(seq (1, length (nm_teste_1000))), function (j)
	lapply (as.list (seq (1,length (local[which (local_1000 >= 4)]))), function (k) ### para cada sitio
	lapply (as.list (seq(1,nsamples)), function (i)
	dim(list_local_traits[[j]][[k]][[i]])[1]  == length (list_null_communities[[j]][[k]][[i]])))))

####
############# traits for species of the pool
list_pool_traits <- lapply (as.list(seq (1, length (nm_teste_1000))), function (j)
	lapply (as.list (seq (1,length (local_1000))), function (k) ### para cada sitio
	#lapply (as.list (seq(1,nsamples)), function (i)
	traits[which (rownames (traits) %in% colnames (site_list_1000[[j]][k,][which (site_list_1000[[j]][k,] > 0.1)])),c(33,36,21,27,29)]))

### species in the sites matching with species with traits
list_pool_sites <- lapply (as.list(seq (1, length (nm_teste_1000))), function (j)
	lapply (as.list (seq (1,length (local_1000))), function (k) ### para cada sitio
	site_list_1000[[j]][k,][which (site_list_1000[[j]][k,] > 0.1)]))

list_pool_sites <- lapply (as.list(seq (1, length (nm_teste_1000))), function (j)
	lapply (as.list (seq (1,length (local_1000))), function (k) ### para cada sitio
	list_pool_sites [[j]][[k]][,which (colnames(list_pool_sites [[j]][[k]]) %in% rownames (list_pool_traits [[j]][[k]]))]))

list_pool_sites <- lapply (as.list(seq (1, length (nm_teste_1000))), function (j)
	lapply (as.list (seq (1,length (local_1000))), function (k) ### para cada sitio
	list_pool_sites [[j]][[k]] [,order (colnames (list_pool_sites [[j]][[k]][which (list_pool_sites [[j]][[k]] > 0.1)]),decreasing=F)]))

### removendo sitios sem sp ou com uma especie
list_pool_sites <- lapply (as.list(seq (1, length (nm_teste_1000))), function (j)
	list_pool_sites[[j]] [which (local_1000 >= 4)])

list_pool_traits <- lapply (as.list(seq (1, length (nm_teste_1000))), function (j)
	list_pool_traits[[j]] [which (local_1000 >= 4)])

## conferindo se o numero de especies ? similar
unlist(lapply (list_pool_traits[[1]],function (i) dim(i)[1]))==unlist(lapply (list_pool_sites[[1]],function (i) dim(i)[2]))

### aplicar pra todas os sitios e simulacoes
#require(parallel)
#cl <- makeCluster(6) ## number of cores = generally ncores -1
#clusterEvalQ(cl, library(FD))
# export your data and function
#clusterExport(cl, c("local","nsamples","nm_teste", "traits", "list_pool_traits","list_pool_sites","list_null_communities","list_local_traits"))

#resultados_500 <- parLapply(cl, 

resultados_1000_local<- lapply (as.list(seq (1, length (nm_teste))), function (j)
	lapply (as.list (seq (1,length (local [which (local >= 4)]))), function (k) ### para cada sitio
	unlist (lapply (as.list (seq(1,nsamples)), function (i) ## para cada simula??o
	dbFD (x= list_local_traits[[j]][[k]][[i]],
	a=list_null_communities[[j]][[k]][[i]],stand.x=T,calc.FRic = TRUE,calc.FGR=F,calc.CWM=F,calc.FDiv=F, 
	w.abun=F,m = "max",corr="cailliez",messages=T)$FRic))))
	

resultados_1000_regional<-lapply (as.list(seq (1, length (nm_teste))), function (j)
	lapply (as.list (seq (1,length (local [which (local >= 4)]))), function (k) ### para cada sitio
	dbFD (x=list_pool_traits[[j]][[k]], a=list_pool_sites [[j]][[k]],stand.x=T,calc.FRic = TRUE,calc.FGR=F,calc.CWM=F,
	calc.FDiv=F, w.abun=F,m = "max",corr="cailliez",messages=T)$FRic))#
		
#stopCluster (cl)
save(resultados_1000_local, resultados_1000_regional, file="composition_FD_OBS_1000_forest.RData")

# OBSERVADOS
local_obs <- data.matrix(data_pred5 [which (local_1000 >= 4),order (colnames(data_pred5), decreasing=F)])
local_obs <- local_obs [,-which (colSums (local_obs)==0)]
list_local_obs <- apply (local_obs, 1, list)

### FRic observado na escala local
FD_local<-dbFD(traits[which (rownames(traits) %in% colnames (local_obs)),c(33,36,21,27,29)], 
	local_obs, stand.x=T,calc.FRic = TRUE,calc.FGR=F,calc.CWM=F,calc.FDiv=F,corr="cailliez",
	w.abun=F, m = "max",messages=T)$FRic

pool_obs <- lapply (site_list_1000, function (i) i [which (local_1000 >=4),order(colnames(i), decreasing=F)])
traits_pool<-traits[which(rownames(traits) %in% colnames (pool_obs[[1]])),c(33,36,21,27,29)]
pool_obs <- lapply (pool_obs, function (i) i[,which (colnames(i) %in% rownames(traits_pool))])
pool_obs <- lapply (pool_obs, function (i) apply (i, 1, list))

pool_obs <- lapply (seq (1, length (pool_obs)), function (j)
	lapply (seq(1, length (pool_obs[[1]])), function (k)
	pool_obs [[j]][[k]][[1]][which (pool_obs [[j]][[k]][[1]] >0.1)]))

traits_pool <- lapply (seq (1, length (pool_obs)), function (j)
	lapply (seq(1, length (pool_obs[[1]])), function (k)
	traits_pool [which (rownames(traits_pool) %in% names (pool_obs[[j]][[k]])),]))

#require(parallel)
#cl <- makeCluster(6) ## number of cores = generally ncores -1
#clusterEvalQ(cl, library(FD))
# export your data and function
#clusterExport(cl, c("pool_obs","traits_pool"))
#parLapply (cl, 

### FRic observada para o pool
FD_pool<- lapply (seq (1, length (pool_obs)), function (j)
	lapply (seq(1, length (pool_obs[[1]])), function (k)
	dbFD (traits_pool[[j]][[k]],pool_obs[[j]][[k]], stand.x=T,calc.FRic = TRUE,calc.FGR=F,calc.CWM=F,calc.FDiv=F,corr="cailliez",
	w.abun=F, m = "max",messages=T)$FRic))

#stopCluster(cl)

save(FD_pool, FD_local, file="composition_FD_OBS_1000_forest.RData")

#### campo
## selecionar campestres
grass_sp_composition_pool<-lapply  (site_list, function (i) i [colnames (i) %in% iucn_hab [iucn_hab$grassland == 1 | iucn_hab$wetland == 1 | iucn_hab$rocky == 1 | iucn_hab$shurubland == 1 | iucn_hab$savanna == 1 | iucn_hab$desert == 1,]$sp])
grass_sp_composition_local<- data_pred4 [, which (colnames (data_pred4) %in% iucn_hab [iucn_hab$grassland == 1 | iucn_hab$wetland == 1 | iucn_hab$rocky == 1 | iucn_hab$shurubland == 1 | iucn_hab$savanna == 1 | iucn_hab$desert == 1,]$sp)]

## funcao modelo nulo
nsamples <- 1
local <- rowSums (grass_sp_composition_local)

require(parallel)
cl <- makeCluster(6) ## number of cores = generally ncores -1
clusterEvalQ(cl, library(phytools))
clusterEvalQ(cl, library(phylobase))
clusterEvalQ(cl, library(dplyr))
# export your data and function
clusterExport(cl, c("local","nsamples", "grass_sp_composition_pool", "sampling_prevalence_NM"))

nm_teste <- parLapply(cl,grass_sp_composition_pool, function (i) sampling_prevalence_NM (i,local=local,nsamples=nsamples))## sample with the size of the obs community
#nm_teste <- lapply(site_list, function (i) sampling_prevalence_NM (i,local=local,nsamples=nsamples))## sample with the size of the obs community

stopCluster(cl)

#### list of null communities 
list_null_communities <- lapply (as.list(seq (1, length (nm_teste))), function (j)
	lapply (as.list (seq (1,length (local))), function (k) ### para cada sitio
	lapply (as.list (seq(1,nsamples)), function (i)
	data.frame (null_sp=nm_teste[[j]]$species[[k]][,i], pres=as.numeric(rep(1, length(nm_teste[[j]]$species[[k]][,i])))))))

list_null_communities <- lapply (as.list(seq (1, length (nm_teste))), function (j)
	lapply (as.list (seq (1,length (local))), function (k) ### para cada sitio
	lapply (as.list (seq(1,nsamples)), function (i)
	as.matrix(t(list_null_communities[[j]][[k]][[i]]))))) ## fazer uma matriz

list_null_communities <- lapply (as.list(seq (1, length (nm_teste))), function (j)
	lapply (as.list (seq (1,length (local))), function (k) ### para cada sitio
	lapply (as.list (seq(1,nsamples)), function (i)
	matrix (list_null_communities[[j]][[k]][[i]],byrow=F,nrow=2,dimnames=list(NULL, list_null_communities[[j]][[k]][[i]][1,])))))

list_null_communities <- lapply (as.list(seq (1, length (nm_teste))), function (j)
	lapply (as.list (seq (1,length (local))), function (k) ### para cada sitio
	lapply (as.list (seq(1,nsamples)), function (i)
	list_null_communities [[j]][[k]][[i]][-1,])))

list_null_communities <- lapply (as.list(seq (1, length (nm_teste))), function (j)
	lapply (as.list (seq (1,length (local))), function (k) ### para cada sitio
	lapply (as.list (seq(1,nsamples)), function (i)
	as.matrix(t(list_null_communities[[j]][[k]][[i]])))))

list_null_communities <- lapply (as.list(seq (1, length (nm_teste))), function (j)
	lapply (as.list (seq (1,length (local))), function (k) ### para cada sitio
	lapply (as.list (seq(1,nsamples)), function (i)
	list_null_communities[[j]][[k]][[i]] [,order (as.character(colnames(list_null_communities[[j]][[k]][[i]])), decreasing=F)])))

list_null_communities <- lapply (as.list(seq (1, length (nm_teste))), function (j)
	lapply (as.list (seq (1,length (local))), function (k) ### para cada sitio
	lapply (as.list (seq(1,nsamples)), function (i)
	as.matrix(t(list_null_communities[[j]][[k]][[i]])))))

list_null_communities <- lapply (as.list(seq (1, length (nm_teste))), function (j)
	lapply (as.list (seq (1,length (local))), function (k) ### para cada sitio
	lapply (as.list (seq(1,nsamples)), function (i)
	{class(list_null_communities[[j]][[k]][[i]])<- "numeric";list_null_communities[[j]][[k]][[i]]})))

#traits
list_local_traits <- lapply (as.list(seq (1, length (nm_teste))), function (j)
	lapply (as.list (seq (1,length (local))), function (k) ### para cada sitio
	lapply (as.list (seq(1,nsamples)), function (i)
	traits[which (rownames(traits) %in% colnames(list_null_communities[[j]][[k]][[i]])),c(33,36,21,27,29)])))

### removendo sitios sem sp ou com uma especie
list_null_communities <- lapply (as.list(seq (1, length (nm_teste))), function (j)
	list_null_communities[[j]] [which (local >= 4)])

list_local_traits <- lapply (as.list(seq (1, length (nm_teste))), function (j)
	list_local_traits [[j]] [which (local >= 4)])
	
list_null_communities <- lapply (as.list(seq (1, length (nm_teste))), function (j)
	lapply (as.list (seq (1,length (local[which (local >= 4)]))), function (k) ### para cada sitio
	lapply (as.list (seq(1,nsamples)), function (i)
	list_null_communities [[j]][[k]][[i]][,which (colnames(list_null_communities [[j]][[k]][[i]]) %in% rownames(list_local_traits [[j]][[k]][[i]]))])))

### conferindo se tudo esta com as mesmas dimensoes
unlist(lapply (as.list(seq (1, length (nm_teste))), function (j)
	lapply (as.list (seq (1,length (local[which (local >= 4)]))), function (k) ### para cada sitio
	lapply (as.list (seq(1,nsamples)), function (i)
	dim(list_local_traits[[j]][[k]][[i]])[1]  == length (list_null_communities[[j]][[k]][[i]])))))

####
############# traits for species of the pool
list_pool_traits <- lapply (as.list(seq (1, length (nm_teste))), function (j)
	lapply (as.list (seq (1,length (local))), function (k) ### para cada sitio
	#lapply (as.list (seq(1,nsamples)), function (i)
	traits[which (rownames (traits) %in% colnames (site_list[[j]][k,][which (site_list[[j]][k,] > 0.1)])),c(33,36,21,27,29)]))

### species in the sites matching with species with traits
list_pool_sites <- lapply (as.list(seq (1, length (nm_teste))), function (j)
	lapply (as.list (seq (1,length (local))), function (k) ### para cada sitio
	site_list[[j]][k,][which (site_list[[j]][k,] > 0.1)]))

list_pool_sites <- lapply (as.list(seq (1, length (nm_teste))), function (j)
	lapply (as.list (seq (1,length (local))), function (k) ### para cada sitio
	list_pool_sites [[j]][[k]][,which (colnames(list_pool_sites [[j]][[k]]) %in% rownames (list_pool_traits [[j]][[k]]))]))

list_pool_sites <- lapply (as.list(seq (1, length (nm_teste))), function (j)
	lapply (as.list (seq (1,length (local))), function (k) ### para cada sitio
	list_pool_sites [[j]][[k]] [,order (colnames (list_pool_sites [[j]][[k]][which (list_pool_sites [[j]][[k]] > 0.1)]),decreasing=F)]))

### removendo sitios sem sp ou com uma especie
list_pool_sites <- lapply (as.list(seq (1, length (nm_teste))), function (j)
	list_pool_sites[[j]] [which (local >= 4)])

list_pool_traits <- lapply (as.list(seq (1, length (nm_teste))), function (j)
	list_pool_traits[[j]] [which (local >= 4)])

## conferindo se o numero de especies ? similar
unlist(lapply (list_pool_traits[[1]],function (i) dim(i)[1]))==unlist(lapply (list_pool_sites[[1]],function (i) dim(i)[2]))

### aplicar pra todas os sitios e simulacoes
#require(parallel)
#cl <- makeCluster(6) ## number of cores = generally ncores -1
#clusterEvalQ(cl, library(FD))
# export your data and function
#clusterExport(cl, c("local","nsamples","nm_teste", "traits", "list_pool_traits","list_pool_sites","list_null_communities","list_local_traits"))

#resultados_500 <- parLapply(cl, 

resultados_500_local<- lapply (as.list(seq (1, length (nm_teste))), function (j)
	lapply (as.list (seq (1,length (local [which (local >= 4)]))), function (k) ### para cada sitio
	unlist (lapply (as.list (seq(1,nsamples)), function (i) ## para cada simula??o
	dbFD (x= list_local_traits[[j]][[k]][[i]],
	a=list_null_communities[[j]][[k]][[i]],stand.x=T,calc.FRic = TRUE,calc.FGR=F,calc.CWM=F,calc.FDiv=F, 
	w.abun=F,m = "max",corr="cailliez",messages=T)$FRic))))
	
resultados_500_regional<-lapply (as.list(seq (1, length (nm_teste))), function (j)
	lapply (as.list (seq (1,length (local [which (local >= 4)]))), function (k) ### para cada sitio
	dbFD (x=list_pool_traits[[j]][[k]], a=list_pool_sites [[j]][[k]],stand.x=T,calc.FRic = TRUE,calc.FGR=F,calc.CWM=F,
	calc.FDiv=F, w.abun=F,m = "max",corr="cailliez",messages=T)$FRic))#
		
#stopCluster (cl)
#saveRDS (resultados_500, "resultados_rich_500_func.rda")
save (resultados_500_regional, resultados_500_local, file="grassland_FRic_NULL_500traps.RData")

# OBSERVADOS
local_obs <- data.matrix(data_pred4 [which (local >= 4),order (colnames(data_pred4), decreasing=F)])
local_obs <- local_obs [,-which (colSums (local_obs)==0)]
list_local_obs <- apply (local_obs, 1, list)

### FRic observado na escala local
FD_local<-dbFD(traits[which (rownames(traits) %in% colnames (local_obs)),c(33,36,21,27,29)], 
	local_obs, stand.x=T,calc.FRic = TRUE,calc.FGR=F,calc.CWM=F,calc.FDiv=F,corr="cailliez",
	w.abun=F, m = "max",messages=T)$FRic

pool_obs <- lapply (site_list, function (i) i [which (local >=4),order(colnames(i), decreasing=F)])
traits_pool<-traits[which(rownames(traits) %in% colnames (pool_obs[[1]])),c(33,36,21,27,29)]
pool_obs <- lapply (pool_obs, function (i) i[,which (colnames(i) %in% rownames(traits_pool))])
pool_obs <- lapply (pool_obs, function (i) apply (i, 1, list))

pool_obs <- lapply (seq (1, length (pool_obs)), function (j)
	lapply (seq(1, length (pool_obs[[1]])), function (k)
	pool_obs [[j]][[k]][[1]][which (pool_obs [[j]][[k]][[1]] >0.1)]))

traits_pool <- lapply (seq (1, length (pool_obs)), function (j)
	lapply (seq(1, length (pool_obs[[1]])), function (k)
	traits_pool [which (rownames(traits_pool) %in% names (pool_obs[[j]][[k]])),]))

#require(parallel)
#cl <- makeCluster(6) ## number of cores = generally ncores -1
#clusterEvalQ(cl, library(FD))
# export your data and function
#clusterExport(cl, c("pool_obs","traits_pool"))
#parLapply (cl, 

### FRic observada para o pool
FD_pool<- lapply (seq (1, length (pool_obs)), function (j)
	lapply (seq(1, length (pool_obs[[1]])), function (k)
	dbFD (traits_pool[[j]][[k]],pool_obs[[j]][[k]], stand.x=T,calc.FRic = TRUE,calc.FGR=F,calc.CWM=F,calc.FDiv=F,corr="cailliez",
	w.abun=F, m = "max",messages=T)$FRic))

#stopCluster(cl)

save(FD_pool, FD_local, file="composition_FD_OBS_500_grassland.RData")

### 1000 traps
################# 
## selecionar campestres
grass_sp_composition_pool_1000<-lapply  (site_list_1000, function (i) i [colnames (i) %in% iucn_hab [iucn_hab$grassland == 1 | iucn_hab$wetland == 1 | iucn_hab$rocky == 1 | iucn_hab$shurubland == 1 | iucn_hab$savanna == 1 | iucn_hab$desert == 1,]$sp])
grass_sp_composition_local_1000<- data_pred5 [, which (colnames (data_pred5) %in% iucn_hab [iucn_hab$grassland == 1 | iucn_hab$wetland == 1 | iucn_hab$rocky == 1 | iucn_hab$shurubland == 1 | iucn_hab$savanna == 1 | iucn_hab$desert == 1,]$sp)]

## funcao modelo nulo
nsamples <- 1
local_1000 <- rowSums (grass_sp_composition_local_1000)

#cl <- makeCluster(6) ## number of cores = generally ncores -1
#clusterEvalQ(cl, library(phytools))
#clusterEvalQ(cl, library(phylobase))
#clusterEvalQ(cl, library(dplyr))
# export your data and function
#clusterExport(cl, c("local_1000","nsamples", "grass_sp_composition_pool_1000", "sampling_prevalence_NM"))
#parLapply (cl,
nm_teste_1000<-lapply(grass_sp_composition_pool_1000, function (i) sampling_prevalence_NM (i,local=local_1000,nsamples=nsamples))## sample with the size of the obs community
#stopCluster(cl)


#### list of null communities 
list_null_communities <- lapply (as.list(seq (1, length (nm_teste_1000))), function (j)
	lapply (as.list (seq (1,length (local_1000))), function (k) ### para cada sitio
	lapply (as.list (seq(1,nsamples)), function (i)
	data.frame (null_sp=nm_teste_1000[[j]]$species[[k]][,i], pres=as.numeric(rep(1, length(nm_teste_1000[[j]]$species[[k]][,i])))))))

list_null_communities <- lapply (as.list(seq (1, length (nm_teste_1000))), function (j)
	lapply (as.list (seq (1,length (local_1000))), function (k) ### para cada sitio
	lapply (as.list (seq(1,nsamples)), function (i)
	as.matrix(t(list_null_communities[[j]][[k]][[i]]))))) ## fazer uma matriz

list_null_communities <- lapply (as.list(seq (1, length (nm_teste_1000))), function (j)
	lapply (as.list (seq (1,length (local_1000))), function (k) ### para cada sitio
	lapply (as.list (seq(1,nsamples)), function (i)
	matrix (list_null_communities[[j]][[k]][[i]],byrow=F,nrow=2,dimnames=list(NULL, list_null_communities[[j]][[k]][[i]][1,])))))

list_null_communities <- lapply (as.list(seq (1, length (nm_teste_1000))), function (j)
	lapply (as.list (seq (1,length (local_1000))), function (k) ### para cada sitio
	lapply (as.list (seq(1,nsamples)), function (i)
	list_null_communities [[j]][[k]][[i]][-1,])))

list_null_communities <- lapply (as.list(seq (1, length (nm_teste_1000))), function (j)
	lapply (as.list (seq (1,length (local_1000))), function (k) ### para cada sitio
	lapply (as.list (seq(1,nsamples)), function (i)
	as.matrix(t(list_null_communities[[j]][[k]][[i]])))))

list_null_communities <- lapply (as.list(seq (1, length (nm_teste_1000))), function (j)
	lapply (as.list (seq (1,length (local_1000))), function (k) ### para cada sitio
	lapply (as.list (seq(1,nsamples)), function (i)
	list_null_communities[[j]][[k]][[i]] [,order (as.character(colnames(list_null_communities[[j]][[k]][[i]])), decreasing=F)])))

list_null_communities <- lapply (as.list(seq (1, length (nm_teste_1000))), function (j)
	lapply (as.list (seq (1,length (local_1000))), function (k) ### para cada sitio
	lapply (as.list (seq(1,nsamples)), function (i)
	as.matrix(t(list_null_communities[[j]][[k]][[i]])))))

list_null_communities <- lapply (as.list(seq (1, length (nm_teste_1000))), function (j)
	lapply (as.list (seq (1,length (local_1000))), function (k) ### para cada sitio
	lapply (as.list (seq(1,nsamples)), function (i)
	{class(list_null_communities[[j]][[k]][[i]])<- "numeric";list_null_communities[[j]][[k]][[i]]})))

#traits
list_local_traits <- lapply (as.list(seq (1, length (nm_teste_1000))), function (j)
	lapply (as.list (seq (1,length (local_1000))), function (k) ### para cada sitio
	lapply (as.list (seq(1,nsamples)), function (i)
	traits[which (rownames(traits) %in% colnames(list_null_communities[[j]][[k]][[i]])),c(33,36,21,27,29)])))

### removendo sitios sem sp ou com uma especie
list_null_communities <- lapply (as.list(seq (1, length (nm_teste_1000))), function (j)
	list_null_communities[[j]] [which (local_1000 >= 4)])

list_local_traits <- lapply (as.list(seq (1, length (nm_teste_1000))), function (j)
	list_local_traits [[j]] [which (local_1000 >= 4)])
	
list_null_communities <- lapply (as.list(seq (1, length (nm_teste_1000))), function (j)
	lapply (as.list (seq (1,length (local[which (local_1000 >= 4)]))), function (k) ### para cada sitio
	lapply (as.list (seq(1,nsamples)), function (i)
	list_null_communities [[j]][[k]][[i]][,which (colnames(list_null_communities [[j]][[k]][[i]]) %in% rownames(list_local_traits [[j]][[k]][[i]]))])))

### conferindo se tudo esta com as mesmas dimensoes
unlist(lapply (as.list(seq (1, length (nm_teste_1000))), function (j)
	lapply (as.list (seq (1,length (local[which (local_1000 >= 4)]))), function (k) ### para cada sitio
	lapply (as.list (seq(1,nsamples)), function (i)
	dim(list_local_traits[[j]][[k]][[i]])[1]  == length (list_null_communities[[j]][[k]][[i]])))))

####
############# traits for species of the pool
list_pool_traits <- lapply (as.list(seq (1, length (nm_teste_1000))), function (j)
	lapply (as.list (seq (1,length (local_1000))), function (k) ### para cada sitio
	#lapply (as.list (seq(1,nsamples)), function (i)
	traits[which (rownames (traits) %in% colnames (site_list_1000[[j]][k,][which (site_list_1000[[j]][k,] > 0.1)])),c(33,36,21,27,29)]))

### species in the sites matching with species with traits
list_pool_sites <- lapply (as.list(seq (1, length (nm_teste_1000))), function (j)
	lapply (as.list (seq (1,length (local_1000))), function (k) ### para cada sitio
	site_list_1000[[j]][k,][which (site_list_1000[[j]][k,] > 0.1)]))

list_pool_sites <- lapply (as.list(seq (1, length (nm_teste_1000))), function (j)
	lapply (as.list (seq (1,length (local_1000))), function (k) ### para cada sitio
	list_pool_sites [[j]][[k]][,which (colnames(list_pool_sites [[j]][[k]]) %in% rownames (list_pool_traits [[j]][[k]]))]))

list_pool_sites <- lapply (as.list(seq (1, length (nm_teste_1000))), function (j)
	lapply (as.list (seq (1,length (local_1000))), function (k) ### para cada sitio
	list_pool_sites [[j]][[k]] [,order (colnames (list_pool_sites [[j]][[k]][which (list_pool_sites [[j]][[k]] > 0.1)]),decreasing=F)]))

### removendo sitios sem sp ou com uma especie
list_pool_sites <- lapply (as.list(seq (1, length (nm_teste_1000))), function (j)
	list_pool_sites[[j]] [which (local_1000 >= 4)])

list_pool_traits <- lapply (as.list(seq (1, length (nm_teste_1000))), function (j)
	list_pool_traits[[j]] [which (local_1000 >= 4)])

## conferindo se o numero de especies ? similar
unlist(lapply (list_pool_traits[[1]],function (i) dim(i)[1]))==unlist(lapply (list_pool_sites[[1]],function (i) dim(i)[2]))

### aplicar pra todas os sitios e simulacoes
#require(parallel)
#cl <- makeCluster(6) ## number of cores = generally ncores -1
#clusterEvalQ(cl, library(FD))
# export your data and function
#clusterExport(cl, c("local","nsamples","nm_teste", "traits", "list_pool_traits","list_pool_sites","list_null_communities","list_local_traits"))

#resultados_500 <- parLapply(cl, 

resultados_1000_local<- lapply (as.list(seq (1, length (nm_teste_1000))), function (j)
	lapply (as.list (seq (1,length (local_1000 [which (local_1000 >= 4)]))), function (k) ### para cada sitio
	unlist (lapply (as.list (seq(1,nsamples)), function (i) ## para cada simula??o
	dbFD (x= list_local_traits[[j]][[k]][[i]],
	a=list_null_communities[[j]][[k]][[i]],stand.x=T,calc.FRic = TRUE,calc.FGR=F,calc.CWM=F,calc.FDiv=F, 
	w.abun=F,m = "max",corr="cailliez",messages=T)$FRic))))
	
resultados_1000_regional<-lapply (as.list(seq (1, length (nm_teste_1000))), function (j)
	lapply (as.list (seq (1,length (local_1000 [which (local_1000 >= 4)]))), function (k) ### para cada sitio
	dbFD (x=list_pool_traits[[j]][[k]], a=list_pool_sites [[j]][[k]],stand.x=T,calc.FRic = TRUE,calc.FGR=F,calc.CWM=F,
	calc.FDiv=F, w.abun=F,m = "max",corr="cailliez",messages=T)$FRic))#
		
#stopCluster (cl)
save(resultados_1000_local, resultados_1000_regional, file="composition_FD_NULL_1000_grassland.RData")

# OBSERVADOS
local_obs <- data.matrix(data_pred5 [which (local_1000 >= 4),order (colnames(data_pred5), decreasing=F)])
local_obs <- local_obs [,-which (colSums (local_obs)==0)]
list_local_obs <- apply (local_obs, 1, list)

### FRic observado na escala local
FD_local1000<-dbFD(traits[which (rownames(traits) %in% colnames (local_obs)),c(33,36,21,27,29)], 
	local_obs, stand.x=T,calc.FRic = TRUE,calc.FGR=F,calc.CWM=F,calc.FDiv=F,corr="cailliez",
	w.abun=F, m = "max",messages=T)$FRic

pool_obs <- lapply (site_list_1000, function (i) i [which (local_1000 >=4),order(colnames(i), decreasing=F)])
traits_pool<-traits[which(rownames(traits) %in% colnames (pool_obs[[1]])),c(33,36,21,27,29)]
pool_obs <- lapply (pool_obs, function (i) i[,which (colnames(i) %in% rownames(traits_pool))])
pool_obs <- lapply (pool_obs, function (i) apply (i, 1, list))

pool_obs <- lapply (seq (1, length (pool_obs)), function (j)
	lapply (seq(1, length (pool_obs[[1]])), function (k)
	pool_obs [[j]][[k]][[1]][which (pool_obs [[j]][[k]][[1]] >0.1)]))

traits_pool <- lapply (seq (1, length (pool_obs)), function (j)
	lapply (seq(1, length (pool_obs[[1]])), function (k)
	traits_pool [which (rownames(traits_pool) %in% names (pool_obs[[j]][[k]])),]))

#require(parallel)
#cl <- makeCluster(6) ## number of cores = generally ncores -1
#clusterEvalQ(cl, library(FD))
# export your data and function
#clusterExport(cl, c("pool_obs","traits_pool"))
#parLapply (cl, 

### FRic observada para o pool
FD_pool1000<- lapply (seq (1, length (pool_obs)), function (j)
	lapply (seq(1, length (pool_obs[[1]])), function (k)
	dbFD (traits_pool[[j]][[k]],pool_obs[[j]][[k]], stand.x=T,calc.FRic = TRUE,calc.FGR=F,calc.CWM=F,calc.FDiv=F,corr="cailliez",
	w.abun=F, m = "max",messages=T)$FRic))

#stopCluster(cl)

save(FD_pool1000, FD_local1000, file="composition_FD_OBS_1000_grassland.RData")

################### ARTIFICIAL HABITATS
## selecionar campestres
artif_sp_composition_pool<-lapply  (site_list, function (i) i [colnames (i) %in% iucn_hab [iucn_hab$artificial >= 1,]$sp])
artif_sp_composition_local <- data_pred4 [, which (colnames (data_pred4) %in% iucn_hab [iucn_hab$artificial >= 1,]$sp)]

## funcao modelo nulo
nsamples <- 1
local <- rowSums (artif_sp_composition_local)

require(parallel)
cl <- makeCluster(6) ## number of cores = generally ncores -1
clusterEvalQ(cl, library(phytools))
clusterEvalQ(cl, library(phylobase))
clusterEvalQ(cl, library(dplyr))
# export your data and function
clusterExport(cl, c("local","nsamples", "artif_sp_composition_pool", "sampling_prevalence_NM"))

nm_teste <- parLapply(cl,artif_sp_composition_pool, function (i) sampling_prevalence_NM (i,local=local,nsamples=nsamples))## sample with the size of the obs community
#nm_teste <- lapply(site_list, function (i) sampling_prevalence_NM (i,local=local,nsamples=nsamples))## sample with the size of the obs community

stopCluster(cl)

#### list of null communities 
list_null_communities <- lapply (as.list(seq (1, length (nm_teste))), function (j)
	lapply (as.list (seq (1,length (local))), function (k) ### para cada sitio
	lapply (as.list (seq(1,nsamples)), function (i)
	data.frame (null_sp=nm_teste[[j]]$species[[k]][,i], pres=as.numeric(rep(1, length(nm_teste[[j]]$species[[k]][,i])))))))

list_null_communities <- lapply (as.list(seq (1, length (nm_teste))), function (j)
	lapply (as.list (seq (1,length (local))), function (k) ### para cada sitio
	lapply (as.list (seq(1,nsamples)), function (i)
	as.matrix(t(list_null_communities[[j]][[k]][[i]]))))) ## fazer uma matriz

list_null_communities <- lapply (as.list(seq (1, length (nm_teste))), function (j)
	lapply (as.list (seq (1,length (local))), function (k) ### para cada sitio
	lapply (as.list (seq(1,nsamples)), function (i)
	matrix (list_null_communities[[j]][[k]][[i]],byrow=F,nrow=2,dimnames=list(NULL, list_null_communities[[j]][[k]][[i]][1,])))))

list_null_communities <- lapply (as.list(seq (1, length (nm_teste))), function (j)
	lapply (as.list (seq (1,length (local))), function (k) ### para cada sitio
	lapply (as.list (seq(1,nsamples)), function (i)
	list_null_communities [[j]][[k]][[i]][-1,])))

list_null_communities <- lapply (as.list(seq (1, length (nm_teste))), function (j)
	lapply (as.list (seq (1,length (local))), function (k) ### para cada sitio
	lapply (as.list (seq(1,nsamples)), function (i)
	as.matrix(t(list_null_communities[[j]][[k]][[i]])))))

list_null_communities <- lapply (as.list(seq (1, length (nm_teste))), function (j)
	lapply (as.list (seq (1,length (local))), function (k) ### para cada sitio
	lapply (as.list (seq(1,nsamples)), function (i)
	list_null_communities[[j]][[k]][[i]] [,order (as.character(colnames(list_null_communities[[j]][[k]][[i]])), decreasing=F)])))

list_null_communities <- lapply (as.list(seq (1, length (nm_teste))), function (j)
	lapply (as.list (seq (1,length (local))), function (k) ### para cada sitio
	lapply (as.list (seq(1,nsamples)), function (i)
	as.matrix(t(list_null_communities[[j]][[k]][[i]])))))

list_null_communities <- lapply (as.list(seq (1, length (nm_teste))), function (j)
	lapply (as.list (seq (1,length (local))), function (k) ### para cada sitio
	lapply (as.list (seq(1,nsamples)), function (i)
	{class(list_null_communities[[j]][[k]][[i]])<- "numeric";list_null_communities[[j]][[k]][[i]]})))

#traits
list_local_traits <- lapply (as.list(seq (1, length (nm_teste))), function (j)
	lapply (as.list (seq (1,length (local))), function (k) ### para cada sitio
	lapply (as.list (seq(1,nsamples)), function (i)
	traits[which (rownames(traits) %in% colnames(list_null_communities[[j]][[k]][[i]])),c(33,36,21,27,29)])))

### removendo sitios sem sp ou com uma especie
list_null_communities <- lapply (as.list(seq (1, length (nm_teste))), function (j)
	list_null_communities[[j]] [which (local >= 4)])

list_local_traits <- lapply (as.list(seq (1, length (nm_teste))), function (j)
	list_local_traits [[j]] [which (local >= 4)])
	
list_null_communities <- lapply (as.list(seq (1, length (nm_teste))), function (j)
	lapply (as.list (seq (1,length (local[which (local >= 4)]))), function (k) ### para cada sitio
	lapply (as.list (seq(1,nsamples)), function (i)
	list_null_communities [[j]][[k]][[i]][,which (colnames(list_null_communities [[j]][[k]][[i]]) %in% rownames(list_local_traits [[j]][[k]][[i]]))])))

### conferindo se tudo esta com as mesmas dimensoes
unlist(lapply (as.list(seq (1, length (nm_teste))), function (j)
	lapply (as.list (seq (1,length (local[which (local >= 4)]))), function (k) ### para cada sitio
	lapply (as.list (seq(1,nsamples)), function (i)
	dim(list_local_traits[[j]][[k]][[i]])[1]  == length (list_null_communities[[j]][[k]][[i]])))))

####
############# traits for species of the pool
list_pool_traits <- lapply (as.list(seq (1, length (nm_teste))), function (j)
	lapply (as.list (seq (1,length (local))), function (k) ### para cada sitio
	#lapply (as.list (seq(1,nsamples)), function (i)
	traits[which (rownames (traits) %in% colnames (site_list[[j]][k,][which (site_list[[j]][k,] > 0.1)])),c(33,36,21,27,29)]))

### species in the sites matching with species with traits
list_pool_sites <- lapply (as.list(seq (1, length (nm_teste))), function (j)
	lapply (as.list (seq (1,length (local))), function (k) ### para cada sitio
	site_list[[j]][k,][which (site_list[[j]][k,] > 0.1)]))

list_pool_sites <- lapply (as.list(seq (1, length (nm_teste))), function (j)
	lapply (as.list (seq (1,length (local))), function (k) ### para cada sitio
	list_pool_sites [[j]][[k]][,which (colnames(list_pool_sites [[j]][[k]]) %in% rownames (list_pool_traits [[j]][[k]]))]))

list_pool_sites <- lapply (as.list(seq (1, length (nm_teste))), function (j)
	lapply (as.list (seq (1,length (local))), function (k) ### para cada sitio
	list_pool_sites [[j]][[k]] [,order (colnames (list_pool_sites [[j]][[k]][which (list_pool_sites [[j]][[k]] > 0.1)]),decreasing=F)]))

### removendo sitios sem sp ou com uma especie
list_pool_sites <- lapply (as.list(seq (1, length (nm_teste))), function (j)
	list_pool_sites[[j]] [which (local >= 4)])

list_pool_traits <- lapply (as.list(seq (1, length (nm_teste))), function (j)
	list_pool_traits[[j]] [which (local >= 4)])

## conferindo se o numero de especies ? similar
unlist(lapply (list_pool_traits[[1]],function (i) dim(i)[1]))==unlist(lapply (list_pool_sites[[1]],function (i) dim(i)[2]))

### aplicar pra todas os sitios e simulacoes
#require(parallel)
#cl <- makeCluster(6) ## number of cores = generally ncores -1
#clusterEvalQ(cl, library(FD))
# export your data and function
#clusterExport(cl, c("local","nsamples","nm_teste", "traits", "list_pool_traits","list_pool_sites","list_null_communities","list_local_traits"))

#resultados_500 <- parLapply(cl, 

resultados_500_local<- lapply (as.list(seq (1, length (nm_teste))), function (j)
	lapply (as.list (seq (1,length (local [which (local >= 4)]))), function (k) ### para cada sitio
	unlist (lapply (as.list (seq(1,nsamples)), function (i) ## para cada simula??o
	dbFD (x= list_local_traits[[j]][[k]][[i]],
	a=list_null_communities[[j]][[k]][[i]],stand.x=T,calc.FRic = TRUE,calc.FGR=F,calc.CWM=F,calc.FDiv=F, 
	w.abun=F,m = "max",corr="cailliez",messages=T)$FRic))))
	
resultados_500_regional<-lapply (as.list(seq (1, length (nm_teste))), function (j)
	lapply (as.list (seq (1,length (local [which (local >= 4)]))), function (k) ### para cada sitio
	dbFD (x=list_pool_traits[[j]][[k]], a=list_pool_sites [[j]][[k]],stand.x=T,calc.FRic = TRUE,calc.FGR=F,calc.CWM=F,
	calc.FDiv=F, w.abun=F,m = "max",corr="cailliez",messages=T)$FRic))#
		
#stopCluster (cl)
#saveRDS (resultados_500, "resultados_rich_500_func.rda")
save (resultados_500_regional, resultados_500_local, file="artif_FRic_NULL_500traps.RData")

# OBSERVADOS
local_obs <- data.matrix(data_pred4 [which (local >= 4),order (colnames(data_pred4), decreasing=F)])
local_obs <- local_obs [,-which (colSums (local_obs)==0)]
list_local_obs <- apply (local_obs, 1, list)

### FRic observado na escala local
FD_local<-dbFD(traits[which (rownames(traits) %in% colnames (local_obs)),c(33,36,21,27,29)], 
	local_obs, stand.x=T,calc.FRic = TRUE,calc.FGR=F,calc.CWM=F,calc.FDiv=F,corr="cailliez",
	w.abun=F, m = "max",messages=T)$FRic

pool_obs <- lapply (site_list, function (i) i [which (local >=4),order(colnames(i), decreasing=F)])
traits_pool<-traits[which(rownames(traits) %in% colnames (pool_obs[[1]])),c(33,36,21,27,29)]
pool_obs <- lapply (pool_obs, function (i) i[,which (colnames(i) %in% rownames(traits_pool))])
pool_obs <- lapply (pool_obs, function (i) apply (i, 1, list))

pool_obs <- lapply (seq (1, length (pool_obs)), function (j)
	lapply (seq(1, length (pool_obs[[1]])), function (k)
	pool_obs [[j]][[k]][[1]][which (pool_obs [[j]][[k]][[1]] >0.1)]))

traits_pool <- lapply (seq (1, length (pool_obs)), function (j)
	lapply (seq(1, length (pool_obs[[1]])), function (k)
	traits_pool [which (rownames(traits_pool) %in% names (pool_obs[[j]][[k]])),]))

#require(parallel)
#cl <- makeCluster(6) ## number of cores = generally ncores -1
#clusterEvalQ(cl, library(FD))
# export your data and function
#clusterExport(cl, c("pool_obs","traits_pool"))
#parLapply (cl, 

### FRic observada para o pool,
FD_pool<- lapply (seq (1, length (pool_obs)), function (j)
	lapply (seq(1, length (pool_obs[[1]])), function (k)
	dbFD (traits_pool[[j]][[k]],pool_obs[[j]][[k]], stand.x=T,calc.FRic = TRUE,calc.FGR=F,calc.CWM=F,calc.FDiv=F,corr="cailliez",
	w.abun=F, m = "max",messages=T)$FRic))

#stopCluster(cl)

save(FD_pool, FD_local, file="composition_FD_OBS_500_artif.RData")

### 1000 traps
################# 
## selecionar ARTIFICIAIS
## selecionar florestais
artif_sp_composition_pool_1000<-lapply  (site_list_1000, function (i) i [colnames (i) %in% iucn_hab [iucn_hab$artificial >= 1,]$sp])
artif_sp_composition_local_1000 <- data_pred5 [, which (colnames (data_pred5) %in% iucn_hab [iucn_hab$artificial >= 1,]$sp)]

### 1000 traps
local_1000 <- rowSums (artif_sp_composition_local_1000)

nsamples <- 1

#cl <- makeCluster(6) ## number of cores = generally ncores -1
#clusterEvalQ(cl, library(phytools))
#clusterEvalQ(cl, library(phylobase))
#clusterEvalQ(cl, library(dplyr))
# export your data and function
#clusterExport(cl, c("local_1000","nsamples", "artif_sp_composition_pool_1000", "sampling_prevalence_NM"))
#parLapply (cl,
nm_teste_1000<-lapply(artif_sp_composition_pool_1000, function (i) sampling_prevalence_NM (i,local=local_1000,nsamples=nsamples))## sample with the size of the obs community
#stopCluster(cl)


#### list of null communities 
list_null_communities <- lapply (as.list(seq (1, length (nm_teste_1000))), function (j)
	lapply (as.list (seq (1,length (local_1000))), function (k) ### para cada sitio
	lapply (as.list (seq(1,nsamples)), function (i)
	data.frame (null_sp=nm_teste_1000[[j]]$species[[k]][,i], pres=as.numeric(rep(1, length(nm_teste_1000[[j]]$species[[k]][,i])))))))

list_null_communities <- lapply (as.list(seq (1, length (nm_teste_1000))), function (j)
	lapply (as.list (seq (1,length (local_1000))), function (k) ### para cada sitio
	lapply (as.list (seq(1,nsamples)), function (i)
	as.matrix(t(list_null_communities[[j]][[k]][[i]]))))) ## fazer uma matriz

list_null_communities <- lapply (as.list(seq (1, length (nm_teste_1000))), function (j)
	lapply (as.list (seq (1,length (local_1000))), function (k) ### para cada sitio
	lapply (as.list (seq(1,nsamples)), function (i)
	matrix (list_null_communities[[j]][[k]][[i]],byrow=F,nrow=2,dimnames=list(NULL, list_null_communities[[j]][[k]][[i]][1,])))))

list_null_communities <- lapply (as.list(seq (1, length (nm_teste_1000))), function (j)
	lapply (as.list (seq (1,length (local_1000))), function (k) ### para cada sitio
	lapply (as.list (seq(1,nsamples)), function (i)
	list_null_communities [[j]][[k]][[i]][-1,])))

list_null_communities <- lapply (as.list(seq (1, length (nm_teste_1000))), function (j)
	lapply (as.list (seq (1,length (local_1000))), function (k) ### para cada sitio
	lapply (as.list (seq(1,nsamples)), function (i)
	as.matrix(t(list_null_communities[[j]][[k]][[i]])))))

list_null_communities <- lapply (as.list(seq (1, length (nm_teste_1000))), function (j)
	lapply (as.list (seq (1,length (local_1000))), function (k) ### para cada sitio
	lapply (as.list (seq(1,nsamples)), function (i)
	list_null_communities[[j]][[k]][[i]] [,order (as.character(colnames(list_null_communities[[j]][[k]][[i]])), decreasing=F)])))

list_null_communities <- lapply (as.list(seq (1, length (nm_teste_1000))), function (j)
	lapply (as.list (seq (1,length (local_1000))), function (k) ### para cada sitio
	lapply (as.list (seq(1,nsamples)), function (i)
	as.matrix(t(list_null_communities[[j]][[k]][[i]])))))

list_null_communities <- lapply (as.list(seq (1, length (nm_teste_1000))), function (j)
	lapply (as.list (seq (1,length (local_1000))), function (k) ### para cada sitio
	lapply (as.list (seq(1,nsamples)), function (i)
	{class(list_null_communities[[j]][[k]][[i]])<- "numeric";list_null_communities[[j]][[k]][[i]]})))

#traits
list_local_traits <- lapply (as.list(seq (1, length (nm_teste_1000))), function (j)
	lapply (as.list (seq (1,length (local_1000))), function (k) ### para cada sitio
	lapply (as.list (seq(1,nsamples)), function (i)
	traits[which (rownames(traits) %in% colnames(list_null_communities[[j]][[k]][[i]])),c(33,36,21,27,29)])))

### removendo sitios sem sp ou com uma especie
list_null_communities <- lapply (as.list(seq (1, length (nm_teste_1000))), function (j)
	list_null_communities[[j]] [which (local_1000 >= 4)])

list_local_traits <- lapply (as.list(seq (1, length (nm_teste_1000))), function (j)
	list_local_traits [[j]] [which (local_1000 >= 4)])
	
list_null_communities <- lapply (as.list(seq (1, length (nm_teste_1000))), function (j)
	lapply (as.list (seq (1,length (local[which (local_1000 >= 4)]))), function (k) ### para cada sitio
	lapply (as.list (seq(1,nsamples)), function (i)
	list_null_communities [[j]][[k]][[i]][,which (colnames(list_null_communities [[j]][[k]][[i]]) %in% rownames(list_local_traits [[j]][[k]][[i]]))])))

### conferindo se tudo esta com as mesmas dimensoes
unlist(lapply (as.list(seq (1, length (nm_teste_1000))), function (j)
	lapply (as.list (seq (1,length (local[which (local_1000 >= 4)]))), function (k) ### para cada sitio
	lapply (as.list (seq(1,nsamples)), function (i)
	dim(list_local_traits[[j]][[k]][[i]])[1]  == length (list_null_communities[[j]][[k]][[i]])))))

####
############# traits for species of the pool
list_pool_traits <- lapply (as.list(seq (1, length (nm_teste_1000))), function (j)
	lapply (as.list (seq (1,length (local_1000))), function (k) ### para cada sitio
	#lapply (as.list (seq(1,nsamples)), function (i)
	traits[which (rownames (traits) %in% colnames (site_list_1000[[j]][k,][which (site_list_1000[[j]][k,] > 0.1)])),c(33,36,21,27,29)]))

### species in the sites matching with species with traits
list_pool_sites <- lapply (as.list(seq (1, length (nm_teste_1000))), function (j)
	lapply (as.list (seq (1,length (local_1000))), function (k) ### para cada sitio
	site_list_1000[[j]][k,][which (site_list_1000[[j]][k,] > 0.1)]))

list_pool_sites <- lapply (as.list(seq (1, length (nm_teste_1000))), function (j)
	lapply (as.list (seq (1,length (local_1000))), function (k) ### para cada sitio
	list_pool_sites [[j]][[k]][,which (colnames(list_pool_sites [[j]][[k]]) %in% rownames (list_pool_traits [[j]][[k]]))]))

list_pool_sites <- lapply (as.list(seq (1, length (nm_teste_1000))), function (j)
	lapply (as.list (seq (1,length (local_1000))), function (k) ### para cada sitio
	list_pool_sites [[j]][[k]] [,order (colnames (list_pool_sites [[j]][[k]][which (list_pool_sites [[j]][[k]] > 0.1)]),decreasing=F)]))

### removendo sitios sem sp ou com uma especie
list_pool_sites <- lapply (as.list(seq (1, length (nm_teste_1000))), function (j)
	list_pool_sites[[j]] [which (local_1000 >= 4)])

list_pool_traits <- lapply (as.list(seq (1, length (nm_teste_1000))), function (j)
	list_pool_traits[[j]] [which (local_1000 >= 4)])

## conferindo se o numero de especies ? similar
unlist(lapply (list_pool_traits[[1]],function (i) dim(i)[1]))==unlist(lapply (list_pool_sites[[1]],function (i) dim(i)[2]))

### aplicar pra todas os sitios e simulacoes
#require(parallel)
#cl <- makeCluster(6) ## number of cores = generally ncores -1
#clusterEvalQ(cl, library(FD))
# export your data and function
#clusterExport(cl, c("local","nsamples","nm_teste", "traits", "list_pool_traits","list_pool_sites","list_null_communities","list_local_traits"))

#resultados_500 <- parLapply(cl, 

resultados_1000_local<- lapply (as.list(seq (1, length (nm_teste_1000))), function (j)
	lapply (as.list (seq (1,length (local_1000 [which (local_1000 >= 4)]))), function (k) ### para cada sitio
	unlist (lapply (as.list (seq(1,nsamples)), function (i) ## para cada simula??o
	dbFD (x= list_local_traits[[j]][[k]][[i]],
	a=list_null_communities[[j]][[k]][[i]],stand.x=T,calc.FRic = TRUE,calc.FGR=F,calc.CWM=F,calc.FDiv=F, 
	w.abun=F,m = "max",corr="cailliez",messages=T)$FRic))))
	
resultados_1000_regional<-lapply (as.list(seq (1, length (nm_teste_1000))), function (j)
	lapply (as.list (seq (1,length (local_1000 [which (local_1000 >= 4)]))), function (k) ### para cada sitio
	dbFD (x=list_pool_traits[[j]][[k]], a=list_pool_sites [[j]][[k]],stand.x=T,calc.FRic = TRUE,calc.FGR=F,calc.CWM=F,
	calc.FDiv=F, w.abun=F,m = "max",corr="cailliez",messages=T)$FRic))#
		
#stopCluster (cl)
save(resultados_1000_local, resultados_1000_regional, file="composition_FD_NULL_1000_artif.RData")

# OBSERVADOS
local_obs <- data.matrix(data_pred5 [which (local_1000 >= 4),order (colnames(data_pred5), decreasing=F)])
local_obs <- local_obs [,-which (colSums (local_obs)==0)]
list_local_obs <- apply (local_obs, 1, list)

### FRic observado na escala local
FD_local1000<-dbFD(traits[which (rownames(traits) %in% colnames (local_obs)),c(33,36,21,27,29)], 
	local_obs, stand.x=T,calc.FRic = TRUE,calc.FGR=F,calc.CWM=F,calc.FDiv=F,corr="cailliez",
	w.abun=F, m = "max",messages=T)$FRic

pool_obs <- lapply (site_list_1000, function (i) i [which (local_1000 >=4),order(colnames(i), decreasing=F)])
traits_pool<-traits[which(rownames(traits) %in% colnames (pool_obs[[1]])),c(33,36,21,27,29)]
pool_obs <- lapply (pool_obs, function (i) i[,which (colnames(i) %in% rownames(traits_pool))])
pool_obs <- lapply (pool_obs, function (i) apply (i, 1, list))

pool_obs <- lapply (seq (1, length (pool_obs)), function (j)
	lapply (seq(1, length (pool_obs[[1]])), function (k)
	pool_obs [[j]][[k]][[1]][which (pool_obs [[j]][[k]][[1]] >0.1)]))

traits_pool <- lapply (seq (1, length (pool_obs)), function (j)
	lapply (seq(1, length (pool_obs[[1]])), function (k)
	traits_pool [which (rownames(traits_pool) %in% names (pool_obs[[j]][[k]])),]))

#require(parallel)
#cl <- makeCluster(6) ## number of cores = generally ncores -1
#clusterEvalQ(cl, library(FD))
# export your data and function
#clusterExport(cl, c("pool_obs","traits_pool"))
#parLapply (cl, 

### FRic observada para o pool
FD_pool1000<- lapply (seq (1, length (pool_obs)), function (j)
	lapply (seq(1, length (pool_obs[[1]])), function (k)
	dbFD (traits_pool[[j]][[k]],pool_obs[[j]][[k]], stand.x=T,calc.FRic = TRUE,calc.FGR=F,calc.CWM=F,calc.FDiv=F,corr="cailliez",
	w.abun=F, m = "max",messages=T)$FRic))

#stopCluster(cl)

save(FD_pool1000, FD_local1000, file="composition_FD_OBS_1000_artif.RData")



