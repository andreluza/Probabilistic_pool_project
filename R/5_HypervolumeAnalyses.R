#########################
# FUNCTIONAL HYPERVOLUME
# for datasets of 100, 500 and 1,000 trap-nights
############################


# --------------------- #
## load packages
source ("R/Packages.R")
source ("R/Functions.R")

# ------------------------------ #
# load trait data

load (here ("data","traits","imputed_traits_to_hypervolume.RData"))

# ------------------------------ #
# load pool rasters
# list of rasters
list_pool <- list.files(here ("data","pool_data", "ProbPool_sp_specific"))
# apply the raster function to them
probabilistic_pool <- lapply (list_pool, function (i) # for each tif
	raster(here ("data","pool_data", "ProbPool_sp_specific",i))) # create raster

# stack pool rasters
probabilistic_pool <- stack (probabilistic_pool)
# set names
spp_names <- read.csv (here("data","pool_data","names_997spp.csv"))
names(probabilistic_pool) <- spp_names$x

# functional diversity of the pool ( FD for each cell - WHOLE WORLD)
## niterations to sampling
niterations <- 3 # I advise you to do a toy example with 3 iterations
LPSP_composition <- values (probabilistic_pool)

# define the minimum number of species per site to be analyzed
min_spp <- 5
LPSP_composition_higher_five_spp <- LPSP_composition[which(round (rowSums (LPSP_composition,na.rm=T))>= min_spp),]

# run a pool sample to mapping functional diversity
pool_sample <- replicate (niterations,
                          
                  	lapply (seq(1,dim(LPSP_composition_higher_five_spp)[1]), function (cell) ## for each cell with more than 5 spp
                  	  # run a random sampling
                  	  sample (LPSP_composition_higher_five_spp [cell,], 
                  	          # with size equal to probabilistic pool size (rounded to next decimal)
                  	          size= round(sum(LPSP_composition_higher_five_spp [cell,],
                  	                          na.rm=T)), 
                  	          # probability-weighted sampling
                            	prob= LPSP_composition_higher_five_spp [cell,])), 
                  	# simplify replicate by showing us an array
                  	simplify = "array" )

# array have these dims:
### pool_sample [1,1] ## [species specific pool n, sample/iteration n]

## get trait values for each species in the pool
traits_pool <- traits[which(rownames (traits) %in% colnames (LPSP_composition_higher_five_spp)),]

# trait composition per sample
traits_pool <- lapply (seq(1,nrow(pool_sample)), function (cell) # for each cell
			              # and iteration
                    lapply (seq(1,ncol(pool_sample)), function (iteration)
                      # get trait composition  
  			              traits [which (rownames(traits) %in% names(pool_sample [cell,iteration][[1]])),]
  			              )
			              )

# observed pool hypervolume 
## 10 ITERATIONS TOOK 2.5 DAYS TO RUN -USING 3 OF 4 CORES OF MY PC
cl <- makeCluster(7) ## number of cores = generally ncores -1
clusterEvalQ(cl, library(hypervolume))
# export your data and function
clusterExport(cl, c("traits_pool"))

hypervolume_pool_obs <- parLapply(cl,traits_pool, function (cell)# for each sample
	lapply(cell, function (iteration)# for each community ##   
      	hypervolume_gaussian (iteration,
		kde.bandwidth=estimate_bandwidth(iteration,method="silverman"))@Volume))

stopCluster(cl)

## mean pool hypervolume (averaged across n interations)

mean_hypervolume_pool <-lapply (hypervolume_pool_obs, function (cell) # for each cell
	# calculate the mean pool hypervolume for each community (averaged across iterations)
	rowMeans (matrix (unlist(lapply (cell, unlist)),ncol=niterations,byrow=T))) 

# save
save(hypervolume_pool_obs ,mean_hypervolume_pool,
	file=here("output","mean_hypervolume_pool.RData")) 

# Here I ran in several steps of 10 or 20 iterations
# saving different files along the way 

# open these files saved at different steps
list.FD <- list.files (here ("output"),
                       pattern = "mean_hypervolume_pool_*")
# load into a list
list_FD_Res <- lapply (list.FD, function (i) {
  
  load(here ("output", i))
  res <-  do.call(rbind,hypervolume_pool_obs)
  
})

# melt and cbind
list_FD_Res <- do.call(cbind,list_FD_Res)
# mean functional hypervolume per cell
FD_res<-lapply(seq (1,nrow(list_FD_Res)), function (i)
  mean(unlist(list_FD_Res[i,]))
)
# melt
FD_res <- unlist(FD_res)

# map of FD
FD_POOL <- FD_res
all_cells <- data.frame ( # cells with more than 5 spp
	m5 = round (rowSums (LPSP_composition,na.rm=T))>=min_spp)
# matching cells
all_cells$hyper <- ifelse (all_cells$m5 == FALSE,
		0,1)
# which cells have => 5 spp
all_cells$hyper[which(round (rowSums (LPSP_composition,na.rm=T))>=min_spp)] <- FD_POOL
## cf
# which(all_cells$hyper > 0) == which (round (rowSums (LPSP_composition,na.rm=T))>=min_spp)
# cells into into raster
cp_pool <- raster(probabilistic_pool) # just a shell to receive values
values (cp_pool ) <- all_cells$hyper  # atribute FD values to raster obj
writeRaster(cp_pool,filename= here("output","FD_pool_final"),  # write it out
	bylayer=TRUE, format="GTiff")

# richness map
SR_pool<- cp_pool
values(SR_pool)<-rowSums(LPSP_composition) # richness per cell
writeRaster(SR_pool,filename= here("output","SR_pool_final"),  # write it out
	bylayer=TRUE, format="GTiff")

########################################
########## null model
# pool = matrix with pool probabilities
# local = a vector with local richness (rowSums of the local matrix)

load(here("data","community_data","community_matrices_effort.RData"))# load comm data
load(here ("data","pool_data","pool_probabilistico_add_especies.RData")) # load pool data

#########################
# NULL MODEL AND HYPERVOLUME ANALYSIS
##########################

# define the sample size per community (richness to the null model; data of 100 trap-nights)
local_100 <- rowSums (data_100_traps_var)
# define the number of random samples per community
# define null model
nsamples <- 3 # the number of samples

### run the null model to sort species from the pool (the same number of specis as observed in the local communities)

cl <- makeCluster(7) ## number of cores = generally ncores -1
clusterEvalQ(cl, library(phytools))
clusterEvalQ(cl, library(phylobase))
clusterEvalQ(cl, library(dplyr))

# export your data and function
clusterExport(cl, c("local_100",
	"nsamples", 
	"site_list_100", 
	"sampling_prevalence_NM"))

nm_teste <- parLapply(cl,site_list_100, function (i) 
	sampling_prevalence_NM (i,local=local_100,nsamples=nsamples))## sample with the size of the obs community

stopCluster(cl)

#### took the traits of the sampled species
#traits
list_local_traits <- lapply (as.list(seq (1, length (nm_teste))), function (j)## for each pool
	lapply (as.list (seq (1,length (local_100))), function (k) ### site
	lapply (as.list (seq(1,nsamples)), function (i) ## and sample
	traits[which (rownames(traits) %in% nm_teste[[j]]$species[[k]] [,i]),]))) ## took the trait

### removing sites with less than min species
list_local_traits <- lapply (as.list(seq (1, length (nm_teste))), function (j)
	list_local_traits [[j]] [which (local_100 >= min_spp)])
	
## NULL DELTA FD

### calculate null hypervolume for that random communities
cl <- makeCluster(7) ## number of cores = generally ncores -1
clusterEvalQ(cl, library(hypervolume))
#export your data and function
clusterExport(cl, c("nsamples","nm_teste", "list_local_traits"))

# working just with species specific dispersal (list_local_traits[[1]])
hypervolume_local_null <- parLapply (cl, list_local_traits[[1]], function (k)
			lapply (k, function (j) 
				hypervolume_gaussian (j,
					kde.bandwidth=estimate_bandwidth(j,
					method="silverman"))@Volume))
stopCluster (cl)

#save
save(hypervolume_local_null , 
	file=here ("output","hypervolume_null_deltaFD_100traps.RData"))

# OBSERVED DELTA FD

### data for observed communities
local_obs <- data.matrix(data_100_traps_var [which (rowSums(data_100_traps_var) >= min_spp),# rm sites with less than 5 spp
							order (colnames(data_100_traps_var), decreasing=F)]) # ordering cols
local_obs <- local_obs [,-which (colSums (local_obs)==0)] # removing absent species
list_local_obs <- apply (local_obs, 1, list) # it's faster to work with lists

## remove absent species
## get the list of detected spp
list_local_obs <- lapply (as.list (seq(1,length (list_local_obs))), function (i) 
	list_local_obs [[i]] [[1]] [which (list_local_obs [[i]] [[1]] == 1)])

### took the traits of species present/occurring in local communities
list_local_traits_obs <- lapply (as.list (seq(1,length (list_local_obs))), function (i)
	traits[which (rownames(traits) %in% names (list_local_obs [[i]])),])

## calculate observed hypervolume

cl <- makeCluster(7) ## number of cores = generally ncores -1
clusterEvalQ(cl, library(hypervolume))
#export your data and function
clusterExport(cl, c("list_local_traits_obs"))

hypervolume_local_obs<- parLapply (cl,list_local_traits_obs, function (i)
	hypervolume_gaussian(i,	
				kde.bandwidth=estimate_bandwidth(i,method="silverman"))@Volume) # calculate and extract the object "Volume"

stopCluster (cl)

# bandwidths - importnt parameters to check hypervolume quality
cl <- makeCluster(7) ## number of cores = generally ncores -1
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

## save obs hypervolume and bandwidths
save (bandwidths_local,
	hypervolume_local_obs,
	file=here("output", "hypervolume_obs_deltaFD_100traps.RData"))

## calculate hypervolume for the pool
pool_obs <- lapply (site_list_100, function (i) i [which (local_100 >=min_spp),# rm communities with <5 spp
									order(colnames(i), decreasing=F)])# order cols
traits_pool<-traits[which(rownames(traits) %in% colnames (pool_obs[[1]])),]
pool_obs <- lapply (pool_obs, function (i) i[,which (colnames(i) %in% rownames(traits_pool))])
pool_obs <- lapply (pool_obs, function (i) apply (i, 1, list))

# run null model with sampling by prevalence
niterations <- 3
pool_sample <- replicate (niterations , lapply (seq(1,length (pool_obs)), function (pool) ## for each pool
	lapply (seq(1,length (pool_obs[[1]])), function (comm) ## for each community
	sample (pool_obs[[pool]][[comm]][[1]], size=round(sum(pool_obs[[pool]][[comm]][[1]])), 
	prob=(pool_obs[[pool]][[comm]][[1]])))), simplify = "array" )

# obtain the traits of random samp
traits_pool <- lapply (seq(1,nrow(pool_sample)), function (pool)
			lapply (seq(1,ncol(pool_sample)), function (sample)
			lapply (seq(1,length (pool_sample[1,1][[1]])), function (comm)
				traits [which (rownames(traits) %in% names(pool_sample [pool,sample][[1]][[comm]])),])))

# pool hypervolume 
cl <- makeCluster(7) ## number of cores = generally ncores -1
clusterEvalQ(cl, library(hypervolume))
# export your data and function
clusterExport(cl, c("traits_pool"))

hypervolume_pool_comm <- parLapply(cl, traits_pool [[1]], function (sample)# for each sample
	lapply (sample, function (comm)# for each community ##   
      	hypervolume_gaussian (comm,# estimate hypervolume
		kde.bandwidth=estimate_bandwidth(comm,method="silverman"))@Volume))

stopCluster(cl)

## mean pool hypervolume (averaged across n interations)
mean_hypervolume_pool_comm<-lapply (hypervolume_pool_obs, function (pool) # for each pool
	# calculate the mean pool hypervolume for each community (averaged across iterations)
	rowMeans (matrix (unlist(lapply (pool, unlist)),ncol=niterations,byrow=T))) 
# save
save (hypervolume_pool_comm,mean_hypervolume_pool_comm,
	 file=here ("output","hypervolumePOOL_100traps_observed.RData"))

# bandwidths used to calculate the pool hypervolume
cl <- makeCluster(7) ## number of cores = generally ncores -1
clusterEvalQ(cl, library(hypervolume))
# export your data and function
clusterExport(cl, c("traits_pool"))

bandwidths_pool_obs <- parLapply(cl, traits_pool [[1]], function (sample)# for each sample
  lapply (sample, function (comm)# for each community ##   
    hypervolume_gaussian (comm,# estimate hypervolume
                          kde.bandwidth=estimate_bandwidth(comm,method="silverman"))@Parameters))
stopCluster(cl)

# Other parameters important to report
bandwidth_kde_pool<-lapply (bandwidths_pool_obs, function (i) # each iteration
  do.call(rbind, lapply (i, function (k)  # each community
      k$kde.bandwidth)))
bandwidth_kde_pool <- do.call(rbind,bandwidth_kde_pool)
# mean bandwidths over samples
apply(bandwidth_kde_pool,2,mean) 
apply(bandwidth_kde_pool,2,sd) 

## SD - ideal value would be 3
sd_count_pool <- lapply(bandwidths_pool_obs, function (i)
      
                      lapply(i, function (k) 
                          
                          k$sd.count))
  
# ideal 3
range(unlist(sd_count_pool))

## samples (shoots) per point
sample_per_point_pool<- lapply(bandwidths_pool_obs, function (i)
  
  do.call(rbind,lapply(i, function (k) 
    
    k$samples.per.point))
)
## mean and sd per pool
mean(do.call(cbind,sample_per_point_pool))
sd(do.call(cbind,sample_per_point_pool))

#save (bandwidths_pool, bandwidths_local, file="bandwidth_100.RData")


#######################
######    500 traps  ##
#######################

local_500 <- rowSums (data_500_traps_var) # local richness (sample size)

##########
# run null model
require(parallel)
cl <- makeCluster(7) ## number of cores = generally ncores -1
clusterEvalQ(cl, library(phytools))
clusterEvalQ(cl, library(phylobase))
clusterEvalQ(cl, library(dplyr))
# export your data and function
clusterExport(cl, c("local_500","nsamples", "site_list_500", "sampling_prevalence_NM"))

nm_teste <- parLapply(cl,site_list_500, function (i) 
				sampling_prevalence_NM (i,local=local_500,
				nsamples=nsamples))## sample with the size of the obs community

stopCluster(cl)

#### took the traits of the sampled species
#traits
list_local_traits <- lapply (as.list(seq (1, length (nm_teste))), function (j)# for each pool
	lapply (as.list (seq (1,length (local_500))), function (k) ### for each community
	lapply (as.list (seq(1,nsamples)), function (i) # for each iteration
	traits[which (rownames(traits) %in% nm_teste[[j]]$species[[k]] [,i]),]))) ## took the traits

### removing the sites with less than 5 species
list_local_traits <- lapply (as.list(seq (1, length (nm_teste))), function (j)
	list_local_traits [[j]] [which (local_500 >= min_spp)])
	
### apply to all communities and iterations
#require(parallel)
cl <- makeCluster(7) ## number of cores = generally ncores -1
clusterEvalQ(cl, library(hypervolume))
#export your data and function
clusterExport(cl, c("nsamples","nm_teste", "list_local_traits"))

# working just with species specific dispersal (list_local_traits[[1]])
hypervolume_local_null_500traps <- parLapply (cl, list_local_traits[[1]], function (k)
			lapply (k, function (j) 
				hypervolume_gaussian (j,
					kde.bandwidth=estimate_bandwidth(j,
					method="silverman"))@Volume))

stopCluster (cl)

save(hypervolume_local_null_500traps , file=here("output","hypervolume_null_deltaFD_500traps.RData"))

# OBSERVED
### data of local species composition
local_obs <- data.matrix(data_500_traps_var[which (local_500 >= min_spp),
							order (colnames(data_500_traps_var), decreasing=F)])
local_obs <- local_obs [,-which (colSums (local_obs)==0)]
list_local_obs <- apply (local_obs, 1, list)

### removing absent species
list_local_obs <- lapply (as.list (seq(1,length (list_local_obs))), function (i) 
	list_local_obs [[i]] [[1]] [which (list_local_obs [[i]] [[1]] == 1)])

### took the traits of the species of each communities
list_local_traits_obs <- lapply (as.list (seq(1,length (list_local_obs))), function (i)
	traits[which (rownames(traits) %in% names (list_local_obs [[i]])),])

## calculate local community hypervolume
cl <- makeCluster(7) ## number of cores = generally ncores -1
clusterEvalQ(cl, library(hypervolume))
#export your data and function
clusterExport(cl, c("list_local_traits_obs"))

hypervolume_local_obs_500<- parLapply (cl,list_local_traits_obs, function (i)
	hypervolume_gaussian(i,	kde.bandwidth=estimate_bandwidth(i,method="silverman"))@Volume)

stopCluster (cl)

## bandwidth
cl <- makeCluster(7) ## number of cores = generally ncores -1
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

## save obs hypervolume and bandwidths
save (bandwidths_local_500,
	hypervolume_local_obs_500,
	file=here("output", "hypervolume_obs_deltaFD_500traps.RData"))

## pool composition data
pool_obs <- lapply (site_list_500, function (i) i [which (local_500 >=min_spp),order(colnames(i), decreasing=F)])
traits_pool<-traits[which(rownames(traits) %in% colnames (pool_obs[[1]])),]
pool_obs <- lapply (pool_obs, function (i) i[,which (colnames(i) %in% rownames(traits_pool))])
pool_obs <- lapply (pool_obs, function (i) apply (i, 1, list))

# Run null model (sampling by prevalence
pool_sample <- replicate (niterations , lapply (seq(1,length (pool_obs)), function (pool) ## for each pool
	lapply (seq(1,length (pool_obs[[1]])), function (comm) ## for each community
	sample (pool_obs[[pool]][[comm]][[1]], 
			size=round(sum(pool_obs[[pool]][[comm]][[1]])), 
			prob=(pool_obs[[pool]][[comm]][[1]])))), simplify = "array" )

## pool_sample [1,1] ## [species specific pool 1, sample 1]

# obtain the traits of random samp
traits_pool <- lapply (seq(1,nrow(pool_sample)), function (pool)
			lapply (seq(1,ncol(pool_sample)), function (sample)
			lapply (seq(1,length (pool_sample[1,1][[1]])), function (comm)
				traits [which (rownames(traits) %in% names(pool_sample [pool,sample][[1]][[comm]])),])))


# pool hypervolume 
cl <- makeCluster(7) ## number of cores = generally ncores -1
clusterEvalQ(cl, library(hypervolume))
# export your data and function
clusterExport(cl, c("traits_pool"))

hypervolume_pool_comm_500 <- parLapply(cl, traits_pool [[1]], function (sample)# for each sample
	lapply (sample, function (comm)# for each community ##   
      	hypervolume_gaussian (comm,# estimate hypervolume
		kde.bandwidth=estimate_bandwidth(comm,method="silverman"))@Volume))

stopCluster(cl)

## mean pool hypervolume (averaged across n interations)
mean_hypervolume_pool_comm_500<-lapply (hypervolume_pool_comm_500, function (pool) # for each pool
	# calculate the mean pool hypervolume for each community (averaged across iterations)
	rowMeans (matrix (unlist(lapply (pool, unlist)),ncol=niterations,byrow=T))) 
# save
save (hypervolume_pool_comm_500,mean_hypervolume_pool_comm_500,
	 file=here ("output","hypervolumePOOL_500traps_observed.RData"))

# bandwidths we used to calculate the pool hypervolume
cl <- makeCluster(7) ## number of cores = generally ncores -1
clusterEvalQ(cl, library(hypervolume))
# export your data and function
clusterExport(cl, c("traits_pool"))

bandwidths_pool_obs <- parLapply(cl, traits_pool [[1]], function (sample)# for each sample
  lapply (sample, function (comm)# for each community ##   
    hypervolume_gaussian (comm,# estimate hypervolume
                          kde.bandwidth=estimate_bandwidth(comm,method="silverman"))@Parameters))

stopCluster(cl)

# Other parameters important to report
bandwidth_kde_pool<-lapply (bandwidths_pool_obs, function (i) # each iteration
  do.call(rbind, lapply (i, function (k)  # each community
    k$kde.bandwidth)))
bandwidth_kde_pool <- do.call(rbind,bandwidth_kde_pool)
# mean bandwidths over samples
apply(bandwidth_kde_pool,2,mean) 
apply(bandwidth_kde_pool,2,sd) 

## SD - ideal value would be 3
sd_count_pool <- lapply(bandwidths_pool_obs, function (i)
  
  lapply(i, function (k) 
    
    k$sd.count))

# ideal 3
range(unlist(sd_count_pool))

## samples (shoots) per point
sample_per_point_pool<- lapply(bandwidths_pool_obs, function (i)
  
  do.call(rbind,lapply(i, function (k) 
    
    k$samples.per.point))
)
## mean and sd per pool
mean(do.call(cbind,sample_per_point_pool))
sd(do.call(cbind,sample_per_point_pool))


#save (bandwidths_pool_obs, bandwidths_local_500, file="bandwidth_500.RData")

######################
####  1,000 traps  ###
######################

local_1000 <- rowSums (data_1000_traps_var) ## local richness

########### run null model
require(parallel)
cl <- makeCluster(7) ## number of cores = generally ncores -1
clusterEvalQ(cl, library(phytools))
clusterEvalQ(cl, library(phylobase))
clusterEvalQ(cl, library(dplyr))
# export your data and function
clusterExport(cl, c("local_1000","nsamples", "site_list_1000", "sampling_prevalence_NM"))

nm_teste <- parLapply(cl,site_list_1000, function (i) 
	sampling_prevalence_NM (i,local=local_1000,nsamples=nsamples))## sample with the size of the obs community

stopCluster(cl)

#### took the traits of the sampled species
#traits
list_local_traits <- lapply (as.list(seq (1, length (nm_teste))), function (j)
	lapply (as.list (seq (1,length (local_1000))), function (k) ### para cada sitio
	lapply (as.list (seq(1,nsamples)), function (i)
	traits[which (rownames(traits) %in% nm_teste[[j]]$species[[k]] [,i]),])))

### removing sites with less than 5 sp
list_local_traits <- lapply (as.list(seq (1, length (nm_teste))), function (j)
	list_local_traits [[j]] [which (local_1000 >= min_spp)])
	
### apply hypervolume function to all sites and iterations
#require(parallel)
cl <- makeCluster(7) ## number of cores = generally ncores -1
clusterEvalQ(cl, library(hypervolume))
#export your data and function
clusterExport(cl, c("nsamples","nm_teste", "list_local_traits"))

# working just with species specific dispersal (list_local_traits[[1]])

hypervolume_local_null_1000traps <- parLapply (cl, list_local_traits[[1]], function (k)
			lapply (k, function (j) 
				hypervolume_gaussian (j,
					kde.bandwidth=estimate_bandwidth(j,
					method="silverman"))@Volume))

stopCluster (cl)

save(hypervolume_local_null_1000traps, 
	file=here ("output","hypervolume_local_null_1000traps.RData"))

# OBSERVED
### observed local composition data
local_obs <- data.matrix(data_1000_traps_var [which (local_1000 >= min_spp),order (colnames(data_1000_traps_var), decreasing=F)])
local_obs <- local_obs [,-which (colSums (local_obs)==0)]
list_local_obs <- apply (local_obs, 1, list)

### removing absent species
list_local_obs <- lapply (as.list (seq(1,length (list_local_obs))), function (i) 
	list_local_obs [[i]] [[1]] [which (list_local_obs [[i]] [[1]] == 1)])

lapply (lapply (list_local_traits_obs, function (i) 
	apply (i, 2, duplicated)), function (k) apply (k,2,table))

### took the traits of species present in the communities
list_local_traits_obs <- lapply (as.list (seq(1,length (list_local_obs))), function (i)
	traits[which (rownames(traits) %in% names (list_local_obs [[i]])),])

## estimate local community hypervolume
cl <- makeCluster(8) ## number of cores = generally ncores -1
clusterEvalQ(cl, library(hypervolume))
#export your data and function
clusterExport(cl, c("list_local_traits_obs"))

hypervolume_local_obs_1000<- parLapply (cl,list_local_traits_obs, function (i)
	hypervolume_gaussian(i,	kde.bandwidth=estimate_bandwidth(i,method="silverman"))@Volume)
stopCluster (cl)

## bandwidths
cl <- makeCluster(8) ## number of cores = generally ncores -1
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

## save obs hypervolume and bandwidths
save (bandwidth_local_1000,
	hypervolume_local_obs_1000,
	file=here("output", "hypervolume_obs_deltaFD_1000traps.RData"))

## pool composition data
pool_obs <- lapply (site_list_1000, function (i) i [which (local_1000 >=min_spp),order(colnames(i), decreasing=F)])
traits_pool<-traits[which(rownames(traits) %in% colnames (pool_obs[[1]])),]
pool_obs <- lapply (pool_obs, function (i) i[,which (colnames(i) %in% rownames(traits_pool))])
pool_obs <- lapply (pool_obs, function (i) apply (i, 1, list))

# run null model
pool_sample <- replicate (niterations , lapply (seq(1,length (pool_obs)), function (pool) ## for each pool
	lapply (seq(1,length (pool_obs[[1]])), function (comm) ## for each community
	sample (pool_obs[[pool]][[comm]][[1]], size=round(sum(pool_obs[[pool]][[comm]][[1]])), 
	prob=(pool_obs[[pool]][[comm]][[1]])))), simplify = "array" )

## pool_sample [1,1] ## [species specific pool 1, sample 1]

## took the traits of each species in the pool

cl <- makeCluster(8) ## number of cores = generally ncores -1
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
cl <- makeCluster(8) ## number of cores = generally ncores -1
clusterEvalQ(cl, library(hypervolume))
# export your data and function
clusterExport(cl, c("traits_pool"))


hypervolume_pool_comm_1000 <- parLapply(cl, traits_pool [[1]], function (sample)# for each sample
	lapply (sample, function (comm)# for each community ##   
      	hypervolume_gaussian (comm,# estimate hypervolume
		kde.bandwidth=estimate_bandwidth(comm,method="silverman"))@Volume))


stopCluster(cl)

## mean pool hypervolume (averaged across n interations)

mean_hypervolume_pool_comm_1000<-lapply (hypervolume_pool_comm_1000, function (pool) # for each pool
	# calculate the mean pool hypervolume for each community (averaged across iterations)
	rowMeans (matrix (unlist(lapply (pool, unlist)),ncol=niterations,byrow=T))) 

#save (hypervolume_local_obs, hypervolume_pool_obs, file="hypervolume_1000traps_observed.RData")
save (hypervolume_pool_comm_1000,
	mean_hypervolume_pool_comm_1000, 
	file=here ("output","hypervolumePOOL_1000traps_observed.RData"))

# bandwidths used to calculate the pool hypervolume
cl <- makeCluster(7) ## number of cores = generally ncores -1
clusterEvalQ(cl, library(hypervolume))
# export your data and function
clusterExport(cl, c("traits_pool"))

bandwidths_pool_1000 <- parLapply(cl, traits_pool [[1]], function (sample)# for each sample
  lapply (sample, function (comm)# for each community ##   
    hypervolume_gaussian (comm,# estimate hypervolume
                          kde.bandwidth=estimate_bandwidth(comm,method="silverman"))@Parameters))

stopCluster(cl)

# Other parameters important to report
bandwidth_kde_pool<-lapply (bandwidths_pool_1000, function (i) # each iteration
  do.call(rbind, lapply (i, function (k)  # each community
    k$kde.bandwidth)))
bandwidth_kde_pool <- do.call(rbind,bandwidth_kde_pool)
# mean bandwidths over samples
apply(bandwidth_kde_pool,2,mean) 
apply(bandwidth_kde_pool,2,sd) 

## SD - ideal value would be 3
sd_count_pool <- lapply(bandwidths_pool_1000, function (i)
  
  lapply(i, function (k) 
    
    k$sd.count))

# ideal 3
range(unlist(sd_count_pool))

## samples (shoots) per point
sample_per_point_pool<- lapply(bandwidths_pool_1000, function (i)
  
  do.call(rbind,lapply(i, function (k) 
    
    k$samples.per.point))
)
## mean and sd per pool
mean(do.call(cbind,sample_per_point_pool))
sd(do.call(cbind,sample_per_point_pool))

# end


