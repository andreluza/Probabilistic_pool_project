

##################################

#       DISPERSAL DATA      
# Imputing missing dispersal data 

##################################

## load packages
source ("R/Packages.R")
source ("R/Functions.R")


# load data taken from Pacifici et al. (2013), and Whitmee & Orme (2013)
# 'incomplete' because there is a lot of NAs in Whitmee & Orme (2013) dispersal data
disp_incomplete <- read.csv (here ("data","traits","disp_year.csv"),h=T,sep=";")

# replacing "_" by "." dot in species name (to match with occ)
disp_incomplete$species<-gsub ("_",".",disp_incomplete$species)

# transforming into number of days (originally in log scale)
disp_incomplete$gen_length <- exp(disp_incomplete$gen_length) 

# then into year (and not days)
disp_incomplete$gen_length <- disp_incomplete$gen_length/365 

# and back to log
disp_incomplete$gen_length <- log(disp_incomplete$gen_length)

# data only for species belonging to small mammal orders
disp_incomplete_small <- subset (disp_incomplete , 
                                 subset= disp_incomplete$Order == "Rodentia" | 
                                          disp_incomplete$Order == "Diprotodontia" | 
                                          disp_incomplete$Order == "Afrosoricida" | 
                                          disp_incomplete$Order == "Eulipotyphla" |
	                                  disp_incomplete$Order == "Dasyuromorphia" | 
                                          disp_incomplete$Order == "Lagomorpha" | 
                                          disp_incomplete$Order == "Didelphimorphia" | 
                                          disp_incomplete$Order == "Peramelemorphia")

# just species weighing less than 5 kg
disp_incomplete_small <- disp_incomplete_small[which(exp(disp_incomplete_small$body_mass) <= 5000),]

#disp_incomplete_small[which(exp(disp_incomplete_small$gen_length) == max(exp(disp_incomplete_small$gen_length))),]

# which species are large (not in small mammal subset)
large <- disp_incomplete [which(disp_incomplete$species %in% disp_incomplete_small$species == F),]

# -------------------------------------
# get the names of small mammal spp
# need to load occurrences and ensemble
# get names of occ files
# you can download 'occ_all' from here: 10.5281/zenodo.5201760
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

ensemble <- readRDS (here ("data","pool_data","ensemble", "prob.ensemble_ch_ts_sp.rds"))
ensemble <- ensemble [[order (names(ensemble))]]

# matching occ and ensemble
OCC2 <- OCC [[which (names(OCC) %in% names(ensemble))]] 
ENV2 <- ensemble[[which (names(ensemble) %in% names(OCC))]]

# then get the names
small_in_occurrence <- names(OCC2) [which (names(OCC2) %in% large$species  == F)]
# small_in_occurrence  [which (small_in_occurrence %in% disp_incomplete_small$species == F)]

# data to imput
disp_to_imput <- data.frame(matrix (NA, nrow=length(small_in_occurrence  [which (small_in_occurrence %in% disp_incomplete_small$species == F)]),
                                    ncol=ncol(disp_incomplete)))
colnames(disp_to_imput ) <- colnames(disp_incomplete)
disp_to_imput$species <- small_in_occurrence  [which (small_in_occurrence %in% disp_incomplete_small$species == F)]
# bind all data
disp_to_imput <- rbind (disp_incomplete_small,disp_to_imput) [,-c(1,3,4,5)]
disp_to_imput <- cbind (disp_to_imput, genus= gsub("_.*","",disp_to_imput$species))
disp_to_imput$genus <- as.numeric (disp_to_imput$genus)

# require missForest for imputing procedure
require(missForest)
# run imputing
test4<-missForest(disp_to_imput[,2:5], maxiter = 10000, ntree = 500, variablewise = T,
	decreasing = FALSE, verbose = T,replace=T)#,

# correlation before and after imputing
cor_before<- cor.test (disp_incomplete_small$body_mass,
          disp_incomplete_small$dispersal)
cor_after <- cor.test (test4$ximp[,1],
          test4$ximp[,3])
# PS: cor_after likely won't  be equal to the correlation reported in the MS
# because of the iterative process of imputation

# create a folder to host the results

dir.create ("output")
dir.create (here ("output", "figures"))

# save imputing result and comparison
png(file=here ("output","figures","imputing_FigS1_1.png"), width=18, height=18, units="cm", res=600)
par (mfrow=c(2,2))
plot(disp_incomplete_small$body_mass, 
     disp_incomplete_small$dispersal, 
     pch=19,cex=0.5, xlab="Log adult body mass (g)", 
     ylab="Log natal dispersal distance (m/yr.)",
	cex.lab=1,cex.axis=0.9, main="Original data")

abline (lm(as.numeric (disp_incomplete_small$dispersal) ~ disp_incomplete_small$body_mass),
        lwd=2)
legend ("bottomright", legend=paste ("r=",
                                     round (cor_before$estimate,2),
                                    "***"),
        bty="n")

plot(test4$ximp[,1],
     test4$ximp [,3],
     pch=19,cex=0.5, xlab="Log adult body mass (g)", ylab="Log natal dispersal distance (m/yr.)",
	cex.lab=1, cex.axis=0.9, main="Imputed data")
abline (lm(test4$ximp[,3] ~ test4$ximp[,1]),lwd=2)
legend ("bottomright", legend=paste ("r=",
                                     round (cor_after$estimate,2),
                                     "***"),
        bty="n")
dev.off()

# getw estimates in log
imputed_traits <- test4$ximp
# bind with othe trait data
imputed_traits <- cbind(imputed_traits, species=disp_to_imput$species)

# summary statistics
require(doBy)
imputed_traits_mean <- summaryBy (body_mass+gen_length+ dispersal ~species,imputed_traits, FUN=mean)
imputed_traits_mean <- imputed_traits_mean [order (as.character(imputed_traits_mean$species), decreasing=F),]

# mean value of dispersal to report in the results
meters_year <- exp(imputed_traits_mean$dispersal.mean)/exp(imputed_traits_mean$gen_length)
mean(meters_year)
sd(meters_year)
# range
imputed_traits_mean[which(meters_year == min (meters_year)),]
imputed_traits_mean[which(meters_year == max (meters_year)),]

# subset to have data of species with occurrence/ensemble data
imputed_traits_mean <- imputed_traits_mean [which (imputed_traits_mean$species %in% names(OCC2)),]

# save
write.csv (imputed_traits_mean, here ("data","traits","log_dispersal_ability_updated.csv"))

#end
