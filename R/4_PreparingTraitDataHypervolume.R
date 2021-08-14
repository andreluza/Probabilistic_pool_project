
# -------------------------------------------- # 
# Preparing traits to calculate DELTA FD
# -------------------------------------------- # 
## load packages
source ("R/Packages.R")
source ("R/Functions.R")

# help in hypervolume: https://benjaminblonder.org/hypervolume_faq.html
# load data from Penone
penone<-read.csv(here ("data","traits","Penone_et_al_2016_mammal_trait_data_imputed.csv"),
                 h=T,sep=",")
# Load Elton traits (Wilman et al 2014)
elton<-read.table (here ("data","traits","elton1.0.csv"),h=T, sep=",")

# ordering,  and subsetting small mammals
penone1 <- penone [order (penone$IUCN.binomial, decreasing=F),]
#penone1$IUCN.binomial<-droplevels (penone1$IUCN.binomial)
penone1 <-subset (penone1, subset = penone1$Body.mass.g <= 5000)
elton1 <- elton [order (elton$binomial, decreasing=F),]
#elton1$binomial<-droplevels (elton1$binomial)
elton1 <-subset (elton1, subset = elton1$BodyMass_Value <= 5000)

# matching trait data bases
penone2<-penone1 [which (penone1$IUCN.binomial %in% elton1$binomial),]
#penone2$IUCN.binomial<-droplevels (penone2$IUCN.binomial)
elton2 <- elton1 [which (elton1$binomial %in% penone2$IUCN.binomial),]
#elton2$binomial<-droplevels (elton2$binomial)

pairs(data.frame (as.numeric(as.factor(elton2$ForStrat_Value)),elton2$Activity_Nocturnal,
	elton2$Diet_Inv, elton2$Diet_PlantO, elton2$Diet_Fruit, elton2$Diet_Seed, penone2$LitPerYear,
	penone2$LitSz), 
	panel=panel.smooth, lwd=2, pch=19)

range(
	cor(data.frame (elton2$ForStrat_Value,elton2$Activity_Nocturnal,
		elton2$Diet_Inv, elton2$Diet_PlantO, elton2$Diet_Fruit, elton2$Diet_Seed,penone2$LitPerYear,penone2$LitSz)[,-1]))

# cbind these two data bases
traits_complete <- cbind (penone2, elton2)

# subsetting based on species in the probabilistic pool
# load pool data
load (here ("data","pool_data","pool_probabilistico_add_especies.RData"))

# subset
traits <- subset (traits_complete, 
	gsub ("_",".",traits_complete$IUCN.binomial) %in% colnames(site_list_100[[1]]))
# adjusting forest strata use
traits$ForStrat_Value<-as.numeric (traits$ForStrat_Value)
rownames (traits)<-traits$binomial

# select trait we will use in hypervolume analysis
sel_cols <- c("Body.mass.g","LitSz","GestLen.d","WeanAge.d","SMA.d",       
              "PopDen.n.km2")
traits <- traits [, which(colnames(traits) %in% sel_cols)]

# correlation and duplication
# cor(traits)
# range(cor(traits))
# apply (apply (traits, 2, duplicated),2,table)

### procedimento de imputing dos dados faltantes
### pegar as especies faltantes e estimar

sp_para_imputar <- colnames(site_list_100[[1]]) [which(colnames(site_list_100[[1]]) %in% gsub ("_",".",traits_complete$IUCN.binomial) == F)]
sp_para_imputar %in%  gsub ("_",".",traits_complete$IUCN.binomial)

# Imput by hand
# set the congenera average
# or instead the information for the sister species (Dasyurus_albopunctatus, D. maculatus)

# "Dasycercus.blythi"
d_blythi <- traits_complete [grep("Dasycercus", traits_complete$IUCN.binomial),c(3,4,6,9,10,13)]
d_blythi<- colMeans (d_blythi)

#  "Dasyurus.maculatus" 
d_maculatus <- traits_complete [grep("Dasyurus", traits_complete$IUCN.binomial),c(3,4,6,9,10,13)]
d_maculatus <- colMeans (d_maculatus)

# Fukomys nao ha
# https://www.researchgate.net/publication/226224782_African_Mole-rats_Bathyergidae_A_Complex_Radiation_in_Tropical_Soils/figures?lo=1
f_boca <- traits_complete [grep("Cryptomys", traits_complete$IUCN.binomial),c(3,4,6,9,10,13)]
f_boca <- colMeans (f_boca)

# hystrix sumatrae em Hystrix.crassispinis
h_criss <- traits_complete [grep("Hystrix", traits_complete$IUCN.binomial),c(3,4,6,9,10,13)]
h_criss <-  colMeans (h_criss)

# flaviventris em marmota caligata
m_caligata <- traits_complete [grep("Marmota", traits_complete$IUCN.binomial),c(3,4,6,9,10,13)]
m_caligata <-  colMeans (m_caligata)

# maxomys rajah em m taju
# https://academic.oup.com/jmammal/article/94/6/1412/907529
m_taju <- traits_complete [grep("Maxomys", traits_complete$IUCN.binomial),c(3,4,6,9,10,13)]
m_taju <- colMeans(m_taju)

# microtus  no lugar de hyperboleus
# parece ocorrer numa regiao similar
m_hyper <- traits_complete [grep("Microtus", traits_complete$IUCN.binomial),c(3,4,6,9,10,13)]
m_hyper <-  colMeans(m_hyper)

# monodelphis kunsi no lugar de peruviana
# de acordo com proximidade filogenetic
m_peru <- traits_complete [grep("Monodelphis", traits_complete$IUCN.binomial),c(3,4,6,9,10,13)]
m_peru <- colMeans(m_peru)

# nanospalax - nao ha
# https://www.researchgate.net/publication/272830161_How_can_scientific_researches_change_conservation_priorities_A_review_of_decade-long_research_on_blind_mole-rats_Rodentia_Spalacinae_in_the_Carpathian_Basin/figures?lo=1
n_ehre <- traits_complete [grep("Spalax", traits_complete$IUCN.binomial),c(3,4,6,9,10,13)]
n_ehre <- colMeans(n_ehre)

#ochotona
o_mant <- traits_complete [grep("Ochotona", traits_complete$IUCN.binomial),c(3,4,6,9,10,13)]
o_mant <- colMeans(o_mant)

# pronolagus
p_saun <- traits_complete [grep("Pronolagus", traits_complete$IUCN.binomial),c(3,4,6,9,10,13)]
p_saun <- colMeans(p_saun)

# Spilocuscus.maculatus - not available
# https://academic.oup.com/jmammal/article/85/5/825/858786
# https://www.researchgate.net/publication/302954153_Phylogenetic_Relationship_of_Cuscuses_Marsupialia_Phalangeridae_from_Papua_and_Maluku_Based_on_Mitochondrial_Sequences_of_NADH_Dehydrogenase_Sub-unit_1_Gene
s_macu <- traits_complete [grep("Phalanger", traits_complete$IUCN.binomial),c(3,4,6,9,10,13)]
s_macu <- colMeans(s_macu )

# Tamiops
t_macc <- traits_complete [grep("Tamiops", traits_complete$IUCN.binomial),c(3,4,6,9,10,13)]
t_macc <- colMeans(t_macc)

# Thylamys
t_pulc <- traits_complete [grep("Thylamys", traits_complete$IUCN.binomial),c(3,4,6,9,10,13)]
t_pulc <- colMeans(t_pulc)

#Imput mean values following previous calculations
imput_these <- rbind(d_blythi, # Dasycercus.blythi
                       d_maculatus, # Dasyurus.maculatus
                       # repeat the same information for 4 Fukomys spp with missing data
                       f_boca, # Fukomys.bocagei
                       f_boca, # Fukomys.damarensis
                       f_boca, # Fukomys.mechowi
                       f_boca, # Fukomys.ochraceocinereus
				               h_criss, # Hystrix.crassispinis
				               m_caligata, # Marmota.caligata
				               m_taju, # Maxomys.tajuddinii
				               m_hyper, # Microtus.hyperboreus
				               m_peru, # Monodelphis.peruviana
				               # repeat the same info for 2 Nannospalax spp
				               n_ehre, # Nannospalax.ehrenbergi
				               n_ehre, # Nannospalax.xanthodon
				               o_mant, # Ochotona.mantchurica
				               p_saun, # Pronolagus.saundersiae
				               s_macu, # Spilocuscus.maculatus
				               t_macc, # Tamiops.macclellandii
				               t_pulc) # Thylamys.pulchellus

rownames(imput_these) <- c(
			"Dasycercus.blythi",
			"Dasyurus.maculatus" ,
			"Fukomys.bocagei",
			"Fukomys.damarensis",
			"Fukomys.mechowi",
			"Fukomys.ochraceocinereus",
			"Hystrix.crassispinis",
			"Marmota.caligata",
			"Maxomys.tajuddinii",
			"Microtus.hyperboreus",
			"Monodelphis.peruviana",
			"Nannospalax.ehrenbergi",
			"Nannospalax.xanthodon",
			"Ochotona.mantchurica",
			"Pronolagus.saundersiae" ,
			"Spilocuscus.maculatus", 
			"Tamiops.macclellandii",
 			"Thylamys.pulchellus")

traits_imput <- rbind(traits,imput_these)
# standardize trait values

traits <- decostand(traits_imput, "standardize")
rownames(traits) <- gsub ("_",".",rownames(traits))
colnames(site_list_100[[1]]) [which(colnames(site_list_100[[1]]) %in% gsub ("_",".",rownames(traits)) == F)]

# save these data 
save (traits, file=here ("data","traits","imputed_traits_to_hypervolume.RData"))

# end