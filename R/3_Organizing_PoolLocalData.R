# ------------------
# Using this script, we calculated the probabilistic species pool,
# based on estimates of species dispersal distances, and site environmental suitability.
# Here we also extracted, for each local community, the pool composition (we did it here because R hosts raster files in its memory, so that we cannot easily save and load rasters; thus we rather prefer to do extraction in this script).
# At the end we calculated DeltaSR for each community and pool design.

# ------------------
## load packages
source ("R/Packages.R")
## load functions
source ("R/Functions.R")

############################################
#### OCCURRENCES AND HABITAT SUITABILITY  ## 
############################################

#---------------------------------------# 
# load occurrences        

# get names of occ files
# 'occ_all' can be download from here: 10.5281/zenodo.5201760
list_occ <- list.files(here ("data","pool_data","occ_all"))
# apply the raster function to them
OCC <- lapply (list_occ, function (i) # for each tif
	raster(here ("data","pool_data","occ_all",i))) # create raster

# stack occ rasters
OCC <- stack (OCC)
# change spp names
names(OCC) <- gsub ("_occ", "", names(OCC))
names(OCC) <- gsub ("_", ".", names(OCC))

#---------------------------------------# 
# load ensemble        
# the ensemble can be download from here: 10.5281/zenodo.5201760
ensemble <- readRDS (here ("data","pool_data","ensemble", "prob.ensemble_ch_ts_sp.rds"))
ensemble <- ensemble [[order (names(ensemble))]]

# matching occ and ensemble
OCC2 <- OCC [[which (names(OCC) %in% names(ensemble))]] 
ENV2 <- ensemble[[which (names(ensemble) %in% names(OCC))]]

# --------------------------------
# load imputed dispersal data
imputed_traits_mean <- read.csv (here ("data","traits","log_dispersal_ability_updated.csv"),h=T)

#######################
#                    ##
#    DISPERSAL POOL  ##
#                    ##
#######################

# PS: due to computation limitation, decrease resolution to 2 x 2 degree cell size
OCC2 <- OCC2 [[which (names (OCC2) %in% gsub("_",".",imputed_traits_mean$species))]]
OCC3 <-lapply (unstack(OCC2), function (i) aggregate (i, fact=2, fun=mean))

##  replace NA by zero, and higher than zero by one
OCC4 <- lapply (OCC3, function (i)  {i [is.nan (values(i))] <- 0; i})
OCC4 <- lapply (OCC4, function (i)  {i [which (values(i)>0)] <- 1; i})
OCC4 <- stack (OCC4)

# decrease the resolution of ensemble too (2x2 degree cell size)
ENV2 <- ENV2 [[which (names (ENV2) %in%  gsub("_",".",imputed_traits_mean$species))]]
ENV3 <- lapply (unstack(ENV2), function (i) aggregate (i, fact=4, fun=mean))

# replacing NA by zero suitability
ENV4 <- lapply (ENV3, function (i)  {i [is.na (values(i))] <- 0; i})

# resample to matching cells and dims
ENV4 <- lapply (ENV4,  function (i) resample (i, OCC4[[1]], method="bilinear"))
ENV4 <- stack (ENV4)

# the way we calculate dispersal over time
#exp(meters/year)/exp(gen_length)/1000*40/110
# (((exp(11.07287)/1000)*40)/110) = 4.326889

# calculate spp dispersal ability
# dispersal in lat long degrees, over 40 years
meters_year <- exp(imputed_traits_mean$dispersal.mean)/exp(imputed_traits_mean$gen_length)
# into km 
km_year <- meters_year/1000
# km in 40 years
km_40years <- km_year * 40
# and then to degrees along 40 years
deg_40years <- (km_40years/110 )

# calculate 
sp_dispersal <-(((exp(imputed_traits_mean$dispersal)/exp(imputed_traits_mean$gen_length)/1000)*40)/110)
names(sp_dispersal)<-imputed_traits_mean$species

# statistics
mean(sp_dispersal)
sd(sp_dispersal)

# calculate the dispersal pool
# species specific dispersal
dispersal_sp<- disppool(sp_dispersal, OCC4,
	                    cost.surfaces=NULL, method = ("negexp"), longlat = F)# "fattail"

# plot(sum(dispersal_sp))

# average dispersal across species
dispersal_mean <-disppool(mean(sp_dispersal), OCC4, 
	                    cost.surfaces=NULL, method = ("negexp"), longlat = F)# "fattail"
# dispersal of half degree
dispersal_05<-disppool(0.5, OCC4, 
	                    cost.surfaces=NULL, method = ("negexp"), longlat = F)# "fattail"

# dispersal of one degree
dispersal_1<-disppool(1, OCC4, 
	                      cost.surfaces=NULL, method = ("negexp"), longlat = F)# "fattail"
#all.equal(names(ENV4), names(dispersal_sp))

## calculate the local probabilistic species pool (LPSP)
## dispersal multiplied by environmental suitability
probabilistic_pool <- lapply (list(dispersal_sp, dispersal_mean,
	dispersal_05, dispersal_1), function (pool)
	ENV4 * pool)

# set names
probabilistic_pool <- lapply (probabilistic_pool, function (i)
	{names(i) <- names (ENV4); i})

# -----------------------------------
# plot dispersal pool
par (mfrow=c(2,2))
lapply (list(dispersal_sp, dispersal_mean,
	dispersal_05, dispersal_1), function (i)
	
	plot(sum(i))
)

# plot prob pool
lapply (probabilistic_pool, function (i) 
	plot (sum(i)))

par (mfrow=c(1,1))

# save probabilistic pool
save (probabilistic_pool, file = here ("data","pool_data",
                                       "probabilistic_pool.RData"))

# write out these rasters in the respective folder
# species specific dispersal pool
dir.create(here ("data", "pool_data","ProbPool_sp_specific"))
writeRaster(probabilistic_pool[[1]],
	here ("data","pool_data","ProbPool_sp_specific","ProbPool.tif"), 
	bylayer=TRUE)
# mean-dispersal ability pool
dir.create(here ("data", "pool_data","ProbPool_mean"))
writeRaster(probabilistic_pool[[2]],
	here ("data","pool_data","ProbPool_mean","ProbPool.tif"), 
	bylayer=TRUE)
# half-degree dispersal ability pool
dir.create(here ("data", "pool_data","ProbPool_05deg"))
writeRaster(probabilistic_pool[[3]],
	here ("data","pool_data","ProbPool_05deg","ProbPool.tif"), 
	bylayer=TRUE)
# one-degree dispersal ability pool
dir.create(here ("data", "pool_data","ProbPool_1deg"))
writeRaster(probabilistic_pool[[4]],
	here ("data","pool_data","ProbPool_1deg","ProbPool.tif"), 
	bylayer=TRUE)

# save names
sp_names <- names(probabilistic_pool[[1]])
write.csv (sp_names,
	file= here ("data","pool_data","names_997spp.csv"))

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

# ----------------------
# get coordinates of local sites
# 100 trap nights
coords_100 <- data.frame(LAT=data_100_traps$LAT, LONG=data_100_traps$LONG)
coordinates(coords_100)<- ~ LONG + LAT
proj4string(coords_100) <- "+proj=longlat +datum=WGS84 +no_defs"

# extract pool probabilities for each site
extract_100_traps <-lapply (probabilistic_pool, function (i) # for each pool 
  extract (i, coords_100, # extract
           method='simple',
           fun=mean, small=T,df=T,na.rm=T) [,-1]) 

# 500 trap nights
coords_500 <- data.frame(LAT=data_500_traps$LAT, LONG=data_500_traps$LONG)
coordinates(coords_500)<- ~ LONG + LAT
proj4string(coords_500) <- "+proj=longlat +datum=WGS84 +no_defs"

# extract
extract_500_traps <-lapply (probabilistic_pool, function (i) # for each pool 
  extract (i, coords_500, # extract
           method='simple',
           fun=mean, small=T,df=T,na.rm=T) [,-1]) 

# 1000 trap nights
coords_1000 <- data.frame(LAT=data_1000_traps$LAT, LONG=data_1000_traps$LONG)
coordinates(coords_1000)<- ~ LONG + LAT
proj4string(coords_1000) <- "+proj=longlat +datum=WGS84 +no_defs"
# extract
extract_1000_traps <-lapply (probabilistic_pool, function (i) # for each pool 
  extract (i, coords_1000, # extract
           method='simple',
           fun=mean, small=T,df=T,na.rm=T) [,-1]) 

# Let's bind species found in local communities that are not found in the pool

site_list_100<-lapply (extract_100_traps, function (i) 
	cbind (data_100_traps_var [which (colnames (data_100_traps_var) %notin% colnames (extract_100_traps[[1]]))], i))
site_list_500<-lapply (extract_500_traps, function (i) 
	cbind (data_500_traps_var [which (colnames (data_500_traps_var) %notin% colnames (extract_500_traps[[1]]))], i))
site_list_1000<-lapply (extract_1000_traps, function (i) 
	cbind (data_1000_traps_var [which (colnames (data_1000_traps_var) %notin% colnames (extract_1000_traps[[1]]))], i))

# save this
save (site_list_100,site_list_500,site_list_1000, 
      file=here ("data","pool_data","pool_probabilistico_add_especies.RData"))

# proportion of species in each mammal order
##  presented in the methods
require(taxize)
list_spp <-  gsub ("\\."," ",colnames (site_list_1000[[1]]))
x <-  classification(list_spp,db = 'ncbi')

# rm not found
x <-x[which(unlist(lapply (x, length)) > 1)]

table(unlist(lapply (x, function (i) 
  i[which(i$rank == "order"),"name"]
)))/length(x)

# end