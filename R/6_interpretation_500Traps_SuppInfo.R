# -----------------------------------------------
## script interpretacao dos resultados
## diferenca de riqueza entre o que poderia ser encontrado no pool e o que eh
## observado nas comunidades locais (sem modelo nulo)
## diferenca entre o hypervolume que poderia ser encontrado no pool e o que eh observado localmetne
## Este script contém os codigos para sumarizar e comparar os dados observados e os gerados pelo modelo nulo


## load packages
source ("R/Packages.R")

# load community data
# load("tabela_comunidades_fatores.RData")
load(here ("data","community_data", "community_matrices_effort.RData"))
# load("coordenadas.RData")

# ----------------------------------------- #
#                 DeltaSR
# ----------------------------------------- #

# observed richness
obsSR_500 <- rowSums (data_500_traps_var, na.rm=T)

# load pool richness
load(here ("data","pool_data","pool_probabilistico_add_especies.RData"))
potential_500 <- apply (site_list_500[[1]],1,sum) # 500 trap-nights effort

# calculate deltaSR
deltaSR_500 <- potential_500 - obsSR_500

# subset with more than 5 spp
min_spp <- 5
deltaSR_500 <- deltaSR_500[which(obsSR_500 >= min_spp)]

# in a few cases deltaSR was lower than zero (likely due to species missing in the pool)
# set these values to zero
# table(deltaSR_100<0); table(deltaSR_500<0); table(deltaSR_500<0)
deltaSR_500 [which (deltaSR_500 < 0)] <- 0

# ----------------------------------------- #
#                  DeltaFD
# ----------------------------------------- #

#  500 TRAPS 

# FD observed for local communities
load(here ("output","hypervolume_obs_deltaFD_500traps.RData"))

# null FD for each community
load(here ("output","hypervolume_local_null_deltaFD_500traps.RData"))

# load pool FD (and then calculate average across iterations)
load(here ("output","hypervolumePOOL_500traps_observed.RData"))
hypervolume_pool_obs_500 <- lapply(hypervolume_pool_comm_500,unlist)# melt the list
hypervolume_pool_obs_500 <-do.call(cbind,hypervolume_pool_obs_500) # melt the list
hypervolume_pool_obs_average_500 <- apply(hypervolume_pool_obs_500,1,mean,na.rm=T)# average pool hypervolume

# calculate observed deltaFD
deltaFD_500 <- hypervolume_pool_obs_average_500 - unlist(hypervolume_local_obs_500)

# calculate null deltaFD
hypervolume_comm_null_500 <- lapply(hypervolume_local_null_500traps,unlist)# melt the list
hypervolume_comm_null_500 <-do.call(rbind,hypervolume_comm_null_500)

# pool FD minus local null FD
null_deltaFD_500 <- lapply (seq (1,ncol (hypervolume_comm_null_500)), function (it) # for each iteration
  
  hypervolume_pool_obs_500 [,it] - hypervolume_comm_null_500 [,it] # calculate the difference = deltaFD
  
)
null_deltaFD_500 <- do.call(cbind,null_deltaFD_500)# melt the list

# SES deltaFD
mean_deltaFD_500 <- apply(null_deltaFD_500,1,mean) # mean across iterations, per community
sd_deltaFD_500<- apply(null_deltaFD_500,1,sd) # standard deviation across iterations, per community

## standardized effect size
SES_deltaFD_500 <- (deltaFD_500 - mean_deltaFD_500)/sd_deltaFD_500

# data frame to analysis and plot
data_to_analysis_500 <- data_500_traps[which(obsSR_500>=min_spp), # sites with 5 or more spp
                                       1:5] # covariates (coordinates and habitat)
data_to_analysis_500 <- cbind (data_to_analysis_500,
                               deltaSR = deltaSR_500,# bind deltaSR
                               deltaFD = SES_deltaFD_500# bind deltaFD
)
# create a factor to describe major type of habitat (natural or modified)
# (state of human modification)
data_to_analysis_500$type <- data_to_analysis_500$hab_cut_edge
# there is a NA in hab_cut_edge (it is a forest habitat as described in the rowName)
data_to_analysis_500$type [which(is.na(data_to_analysis_500$type))] <- "forest.0.0.0"
# adjust factor
levels (data_to_analysis_500$type) [which (levels (data_to_analysis_500$type) == "forest.0.0.0")] <- "Natural"
levels (data_to_analysis_500$type) [which (levels (data_to_analysis_500$type) == "grassland.0.0.0")] <- "Natural"
levels (data_to_analysis_500$type) [which (levels (data_to_analysis_500$type) != "Natural")] <- "Human-modified"

# changing the order habitat factor to have a good plot
data_to_analysis_500$hab_cut_edge <- factor(data_to_analysis_500$hab_cut_edge, levels = c("forest.0.0.0", "grassland.0.0.0",
                                                                        "open.0.0.0", "open.1.0.0","edge.0.0.0",
                                                                        "edge.0.1.0", "tree_plantation.0.0.0"))

# Figure 1 (main text)

pdf(here ("output","figures", "Fig2_500T.pdf"),width = 4,height=4)
# scatter histogram
plot500<- ggscatterhist(
      data_to_analysis_500, x = "deltaSR", y = "deltaFD",
      color = "type", size = 2, shape = 21,alpha = 0.8,
      palette = c("red", "blue"),
      title=NULL,
      margin.plot = "boxplot",ylim=c(-2.5,2.5),
      ggtheme = theme_classic()+ theme (legend.direction = "horizontal",
                                        legend.title=element_blank()),
      
      xlab =expression (paste(Delta, "SR",sep="")),
      ylab=expression (paste(Delta, "FD"["SES"],sep=""))
    )
# lines
plot500$sp <- plot500$sp+ 
      geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
      geom_vline(xintercept = mean(data_to_analysis_500$deltaSR), linetype = "dashed", color = "black") 
plot500    
dev.off()   
   

# average values to report in the results
# global average of deltaSR
mean(data_to_analysis_500$deltaSR) 
sd(data_to_analysis_500$deltaSR) # sd
# natural habtiat average
mean(data_to_analysis_500$deltaSR [which(data_to_analysis_500$type == "Natural")])
sd(data_to_analysis_500$deltaSR [which(data_to_analysis_500$type == "Natural")])
# human modified
mean(data_to_analysis_500$deltaSR [which(data_to_analysis_500$type == "Human-modified")])
sd(data_to_analysis_500$deltaSR [which(data_to_analysis_500$type == "Human-modified")])

# global average of deltaFD
mean(data_to_analysis_500$deltaFD) 
sd(data_to_analysis_500$deltaFD) # sd
# natural habtiat average
mean(data_to_analysis_500$deltaFD [which(data_to_analysis_500$type == "Natural")])
sd(data_to_analysis_500$deltaFD [which(data_to_analysis_500$type == "Natural")])
# human modified
mean(data_to_analysis_500$deltaFD [which(data_to_analysis_500$type == "Human-modified")])
sd(data_to_analysis_500$deltaFD [which(data_to_analysis_500$type == "Human-modified")])

#---------------------------------------
# linear regression (or anova) to test variation between / among factors
# major type (state of human modification)
summary(lm_500_deltaSR <- (lm (deltaSR~type,data=data_to_analysis_500))) # deltaSR
summary(lm_500_deltaFD <- (lm (deltaFD~type,data=data_to_analysis_500))) # deltaFD

#---------------------------------------
# differences between categories of habitat type
lm_500_deltaSR_habitats <- (lm (deltaSR~hab_cut_edge,data=data_to_analysis_500)) # deltaSR
lm_500_deltaFD_habitats <- (lm (deltaFD~hab_cut_edge,data=data_to_analysis_500)) # deltaFD

# aov object to run Tukey Honest Significance test
summary (lm_500_deltaSR_habitats_reaj <- (aov (deltaSR~hab_cut_edge-1,data=data_to_analysis_500)))
summary (lm_500_deltaFD_habitats_reaj <- (aov (deltaFD~hab_cut_edge-1,data=data_to_analysis_500)))

# run contrast test
TukeyHSD (lm_500_deltaSR_habitats_reaj, "hab_cut_edge", ordered=T)
TukeyHSD (lm_500_deltaFD_habitats_reaj, "hab_cut_edge", ordered=T)

#-----------------
# delta SR
# contrast plot
## extracting parameters 
nms <- unique (names (lm_500_deltaSR_habitats_reaj$coefficients))
# adjusting names
nms <- gsub ("hab_cut_edge","",nms)
nms <- sub("forest.0.0.0", "Forest", nms)
nms <- sub("grassland.0.0.0", "Grassland", nms)
nms <- sub("open.0.0.0", "Crop field", nms)
nms <- sub("open.1.0.0", "Clear-cut", nms)
nms <- sub("edge.0.0.0", "Forest edge", nms)
nms <- sub("edge.0.1.0", "Grassland edge", nms)
nms <- sub("tree_plantation.0.0.0", "Tree plantation", nms)
# datafarme with estimates  
df <- data.frame(term     = rep(nms, 1),
                 estimate = unname(lm_500_deltaSR_habitats_reaj$coefficients))
# confidence intervals
df <- transform(df,
                lower = confint (lm_500_deltaSR_habitats_reaj)[,1],
                upper = confint (lm_500_deltaSR_habitats_reaj)[,2])
# ordering factor levels
df$term <- factor(df$term, levels = c("Tree plantation",
                                      "Forest edge",
                                      "Grassland edge", 
                                      "Clear-cut",
                                      "Crop field", 
                                      "Grassland",
                                      "Forest"))
df$tipo <- c(rep("N",2),rep("HM",5)) # state of human modification to set colors
# plot
pd<-position_dodge(0.5)
p1 <- ggplot (df,  aes (x=term, y=estimate, fill=NULL,group = NULL,colour=tipo))+
  scale_x_discrete() + #scale_colour_grey() +
  geom_errorbar(aes(ymin=lower, ymax=upper), width=.2,position=pd) +
  geom_line(position=pd) +
  geom_point(position=pd, aes(colour=tipo),size=4) +
  theme (legend.position = "none") + coord_fixed(ratio = 1)
# flip coords  
p2 <- p1 + coord_flip() +  theme_classic() +
  labs(x = NULL,
       y = expression (paste(Delta, "SR",sep="")),
       #title = expression (paste ("Estimate of ",Delta, "SR per habitat")),
       subtitle =NULL)  + 
  scale_color_manual(values=c('red','blue'))

# annotate
p3_tukey_deltaSR <- p2 + geom_hline (yintercept = 0,linetype = "dashed") + 
  theme (legend.position="none",
         axis.title.x = element_text(size=12),
         axis.text.x =  element_text(size=11),
         axis.text.y = element_text(angle = 35, hjust =1,size=12,vjust = 0))
  
#-----------------
# delta FD
# contrast plot
## extracting parameters 
nms <- unique (names (lm_500_deltaFD_habitats_reaj$coefficients))
# adjusting names
nms <- gsub ("hab_cut_edge","",nms)
nms <- sub("forest.0.0.0", "Forest", nms)
nms <- sub("grassland.0.0.0", "Grassland", nms)
nms <- sub("open.0.0.0", "Crop field", nms)
nms <- sub("open.1.0.0", "Clear-cut", nms)
nms <- sub("edge.0.0.0", "Forest edge", nms)
nms <- sub("edge.0.1.0", "Grassland edge", nms)
nms <- sub("tree_plantation.0.0.0", "Tree plantation", nms)
# datafarme with estimates  
df <- data.frame(term     = rep(nms, 1),
                 estimate = unname(lm_500_deltaFD_habitats_reaj$coefficients))
# confidence intervals
df <- transform(df,
                lower = confint (lm_500_deltaFD_habitats_reaj)[,1],
                upper = confint (lm_500_deltaFD_habitats_reaj)[,2])
# ordering factor levels
df$term <- factor(df$term, levels = c("Tree plantation",
                                      "Forest edge",
                                      "Grassland edge", 
                                      "Clear-cut",
                                      "Crop field", 
                                      "Grassland",
                                      "Forest"))
df$tipo <- c(rep("N",2),rep("HM",5)) # state of human modification to set colors
# plot
pd<-position_dodge(0.5)
p1 <- ggplot (df,  aes (x=term, y=estimate, fill=NULL,group = NULL,colour=tipo))+
  scale_x_discrete() + #scale_colour_grey() +
  geom_errorbar(aes(ymin=lower, ymax=upper), width=.2,position=pd) +
  geom_line(position=pd) +
  geom_point(position=pd, aes(colour=tipo),size=4) +
  theme (legend.position = "none") + coord_fixed(ratio = 1)
# flip coords  
p2 <- p1 + coord_flip() +  theme_classic() +
  labs(x = NULL,
       y = expression (paste(Delta, "FD",sep="")),
       #title = expression (paste ("Estimate of ",Delta, "SR per habitat")),
       subtitle =NULL)  + 
  scale_color_manual(values=c('red','blue'))
# annotate
p3_tukey_deltaFD <- p2 + geom_hline (yintercept = 0,linetype = "dashed") + 
  theme (legend.position="none",
         axis.title.x = element_text(size=12),
         axis.text.x =  element_text(size=11),
         axis.text.y = element_text(angle = 35, hjust =1,size=12,vjust = 0))


# arrange grid

pdf (here ("output","figures", "Fig3_500T.pdf"),height=4,width=4,onefile=T)
p3_tukey_deltaSR
p3_tukey_deltaFD
dev.off()

## Figure 4

# scenarios
data_to_analysis_500 <- cbind (data_to_analysis_500,
                                scenarioSR = ifelse (data_to_analysis_500$deltaSR >=mean(data_to_analysis_500$deltaSR),
                                                     "Higher-SR","Lower-SR"), # S1 and S3
                                scenarioFD = ifelse (data_to_analysis_500$deltaFD >= 0,
                                                     "Higher-FD","Lower-FD")) # S2 and S4

# interaction to have scenarios
data_to_analysis_500 <- cbind (data_to_analysis_500, 
                                IntDeltaSR_FD=interaction(data_to_analysis_500$scenarioSR,
                                                          data_to_analysis_500$scenarioFD))

# edit levels of the factor 'scenarios'
data_to_analysis_500$IntDeltaSR_FD <- factor (data_to_analysis_500$IntDeltaSR_FD, 
                                  levels=c("Lower-SR.Lower-FD",
                                           "Lower-SR.Higher-FD",
                                           "Higher-SR.Lower-FD",
                                           "Higher-SR.Higher-FD"))


# Table 1
## Pearson's chi-square analysis (states of human modification)
with(data_to_analysis_500, chisq.test (type,IntDeltaSR_FD,simulate.p.value = F))# to have degrees of freedom
with(data_to_analysis_500, chisq.test (type,IntDeltaSR_FD,simulate.p.value = T))
with(data_to_analysis_500, chisq.test (type,IntDeltaSR_FD,simulate.p.value = T))$expected
with(data_to_analysis_500, chisq.test (type,IntDeltaSR_FD,simulate.p.value = T))$observed

# Table 2
# extract the biome of each community
biomes <- readOGR (dsn= path.expand(here ("data","terr-ecoregions-TNC")),
                                  "tnc_terr_ecoregions")

# coords into spatial df
sp_df_coord <- data_to_analysis_500 [,c("LONG", "LAT")]
coordinates (sp_df_coord) <- ~ LONG + LAT
crs (sp_df_coord) <- crs(biomes)
# overlap communities on wwf biomes
overlap_pts_biomes <- over (sp_df_coord,biomes)
# bind into the data frame
data_to_analysis_500$biomes <- as.factor (overlap_pts_biomes$WWF_MHTNAM)
data_to_analysis_500$biomes <-  droplevels(data_to_analysis_500$biomes)

# biomes into four main groups
# Tropical forests = Tropical and Subtropical Moist Broadleaf Forests, Tropical and Subtropical Dry Broadleaf Forests;
tropical_forests <- c ("Tropical and Subtropical Moist Broadleaf Forests", 
                       "Tropical and Subtropical Dry Broadleaf Forests",
                       "Tropical and Subtropical Coniferous Forests")
# Temperate Forests = Boreal Forests/Taiga, Mediterranean Forests, Woodlands and Scrub, Temperate Broadleaf and Mixed Forests, Temperate Conifer Forests; 
temperate_forests <- c ("Boreal Forests/Taiga", 
                        "Mediterranean Forests, Woodlands and Scrub", 
                        "Temperate Broadleaf and Mixed Forests", 
                        "Temperate Conifer Forests")
# Tropical grasslands = Montane Grasslands and Shrublands, Tropical and Subtropical Grasslands, Savannas and Shrublands; 
tropical_grasslands <- c("Montane Grasslands and Shrublands", 
                         "Tropical and Subtropical Grasslands, Savannas and Shrublands",
                         "Deserts and Xeric Shrublands",
                         "Flooded Grasslands and Savannas")
# Temperate Grasslands = Temperate Grasslands, Savannas and Shrublands.
temperate_grasslands <- c ("Temperate Grasslands, Savannas and Shrublands",
                           "Tundra")

# transforming
levels (data_to_analysis_500$biomes) [which(levels (data_to_analysis_500$biomes) %in% tropical_forests)] <- "trop.for"
levels (data_to_analysis_500$biomes) [which(levels (data_to_analysis_500$biomes) %in% temperate_forests)] <- "temp.for"
levels (data_to_analysis_500$biomes) [which(levels (data_to_analysis_500$biomes) %in% tropical_grasslands)] <- "trop.grass"
levels (data_to_analysis_500$biomes) [which(levels (data_to_analysis_500$biomes) %in% temperate_grasslands)] <- "temp.grass"

# interaction between biome and state of human modification
data_to_analysis_500$IntStateBiome <- interaction (data_to_analysis_500$type,
                                                    data_to_analysis_500$biomes)

# Chi-squared analysis
with (data_to_analysis_500, chisq.test(IntStateBiome, 
                                        IntDeltaSR_FD,
                                        simulate.p.value = F)) # to get the degrees of freedom
with (data_to_analysis_500, chisq.test(IntStateBiome, 
                                        IntDeltaSR_FD,
                                        simulate.p.value = T))
with (data_to_analysis_500, chisq.test(IntStateBiome, 
                                        IntDeltaSR_FD,
                                        simulate.p.value = T))$observed
# totals
rowSums(with (data_to_analysis_500, chisq.test(IntStateBiome, 
                                                IntDeltaSR_FD,
                                                simulate.p.value = T))$observed
)
colSums(with (data_to_analysis_500, chisq.test(IntStateBiome, 
                                                IntDeltaSR_FD,
                                                simulate.p.value = T))$observed
)
with (data_to_analysis_500, chisq.test(IntStateBiome, 
                                        IntDeltaSR_FD,
                                        simulate.p.value = T))$expected

# -------------------------------
# spatial distribution of scenarios
# Figure 4

# plotting
require("rnaturalearth")
require("rnaturalearthdata")

world <- ne_countries(scale = "medium", returnclass = "sf")

require(ggplot2)
# world map
wm <- ggplot() + geom_sf (data=world, size = 0.1, 
                          fill="gray95", colour="gray95") +
  coord_sf (xlim = c(-200, 200),  ylim = c(-57, 100)) +
  theme_bw() + xlab(NULL) + ylab(NULL) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_blank(),axis.ticks.x=element_blank(),
        axis.text.y = element_blank(),axis.ticks.y=element_blank(),
        title = element_text(size=8)) 
# subtly jittering points
jitter <- position_jitter(width = 3.5, height = 3.5)
# show points based on scenarios
wm_b <- wm + 
  geom_point (data=data_to_analysis_500 [which(data_to_analysis_500$type == "Human-modified"),], 
              aes (x=LONG, y=LAT,
                    colour=data_to_analysis_500$IntDeltaSR_FD[which(data_to_analysis_500$type == "Human-modified")]),
              size=1.7,shape=1, stroke = 0.9, position = jitter,
              alpha = ifelse(as.numeric(data_to_analysis_500$IntDeltaSR_FD[which(data_to_analysis_500$type == "Human-modified")])>=3,
                             1,0.7)) + #
  
  scale_color_manual(expression (paste ("Relationship between ", Delta,"SR and ","FD"["SES"], " (Scenarios)",sep="")),
                     values = c("Lower-SR.Lower-FD" = "#0571b0",
                                "Lower-SR.Higher-FD" = "#80cdc1",#"springgreen3",
                                "Higher-SR.Lower-FD" = "#ca0020",
                                "Higher-SR.Higher-FD" = "#dfc27d"),
                     
                     labels=c(expression (paste ("Lower than ",bar(paste (Delta,"SR")), ", ","lower than ",paste (Delta,"FD"["SES"]),"=0 (S1)")),
                              expression (paste ("Lower than ",bar(paste (Delta,"SR")), ", ","higher than ",paste (Delta,"FD"["SES"]),"=0 (S2)")),
                              expression (paste ("Higher than ",bar(paste (Delta,"SR")), ", ","lower than ",paste (Delta,"FD"["SES"]),"=0 (S3)")),
                              expression (paste ("Higher than ",bar(paste (Delta,"SR")), ", ","higher than ",paste (Delta,"FD"["SES"]),"=0 (S4)"))
                     ))
# set title
wm_b <- wm_b + ggtitle ("Human-modified habitats")
# theme
wm_c_legend<-wm_b+theme(legend.position = "top",
                        legend.title.align=0.5,
                        legend.direction="horizontal",
                        legend.justification="top",
                        legend.box="vertical",
                        legend.margin = margin (0,0,0,0),
                        legend.box.margin = margin(-2,-30,-30,-2),
                        legend.title = element_text(color = "black", size = 9),
                        legend.text = element_text(color = "black", size = 7.5),
                        strip.text = element_text(size=7.5),
                        strip.background = element_rect(fill="white",colour="gray80"),
                        panel.grid.minor = element_blank(),               #removes minor grid lines
                        panel.grid.major = element_blank()) +
  guides(colour = guide_legend(title.position="top",
                               ncol=2, bycol=T,override.aes = list(size=2))) 



# capture the legend
require(gridExtra)
get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

common_legend <- get_legend(wm_c_legend) ## get the legend of a map with legend
# remove legend of the first plot (we will show the legend in a composition of maps)
wm_c<-wm_b+theme(plot.title=element_text(hjust=0.5,size=8,
                                         margin=margin(10,-10,-10,-100)),
                 legend.position = "none",
                 legend.title.align=0.5,
                 legend.direction="horizontal",
                 legend.justification="top",
                 legend.box="vertical",
                 legend.margin = margin (0,0,0,0),
                 legend.box.margin = margin(-2,-30,-30,-2),
                 legend.title = element_text(color = "black", size = 9),
                 legend.text = element_text(color = "black", size = 7.5),
                 strip.text = element_text(size=7.5),
                 strip.background = element_rect(fill="white",colour="gray80"),
                 panel.grid.minor = element_blank(),               #removes minor grid lines
                 panel.grid.major = element_blank(),
                 plot.margin = unit(c(-6,-0.7,-5, -0.5), "lines")) +
  guides(colour = guide_legend(title.position="top",
                               ncol=2, bycol=T,override.aes = list(size=2))) 
# final plot
(plot1 <- wm_c) # human- modified habitats

# ----------------------------
# plot 2, natural habitats
wm_b2 <- wm + 
  geom_point (data=data_to_analysis_500 [which(data_to_analysis_500$type == "Natural"),], 
              aes (x=LONG, y=LAT,
                   colour=data_to_analysis_500$IntDeltaSR_FD[which(data_to_analysis_500$type == "Natural")]),
              size=1.7,shape=1, stroke = 0.9, position = jitter,
              alpha = ifelse(as.numeric(data_to_analysis_500$IntDeltaSR_FD[which(data_to_analysis_500$type == "Natural")])>=3,
                             1,0.7)) + #
  
  scale_color_manual(expression (paste ("Relationship between ", Delta,"SR and ","FD"["SES"], " (Scenarios)",sep="")),
                     values = c("Lower-SR.Lower-FD" = "#0571b0",
                                "Lower-SR.Higher-FD" = "#80cdc1",#"springgreen3",
                                "Higher-SR.Lower-FD" = "#ca0020",
                                "Higher-SR.Higher-FD" = "#dfc27d"),
                     
                     labels=c(expression (paste ("Lower than ",bar(paste (Delta,"SR")), ", ","lower than ",paste (Delta,"FD"["SES"]),"=0 (S1)")),
                              expression (paste ("Lower than ",bar(paste (Delta,"SR")), ", ","higher than ",paste (Delta,"FD"["SES"]),"=0 (S2)")),
                              expression (paste ("Higher than ",bar(paste (Delta,"SR")), ", ","lower than ",paste (Delta,"FD"["SES"]),"=0 (S3)")),
                              expression (paste ("Higher than ",bar(paste (Delta,"SR")), ", ","higher than ",paste (Delta,"FD"["SES"]),"=0 (S4)"))
                     ))
# set title
wm_b2 <- wm_b2 + ggtitle ("Natural habitats")

# plot without a legend (shared with the previous plot)
wm_b2<-wm_b2+theme(plot.title=element_text(hjust=0.5,size=8,
                                         margin=margin(10,-10,-10,-100)),
                 legend.position = "none",
                 legend.title.align=0.5,
                 legend.direction="horizontal",
                 legend.justification="top",
                 legend.box="vertical",
                 legend.margin = margin (0,0,0,0),
                 legend.box.margin = margin(-2,-30,-30,-2),
                 legend.title = element_text(color = "black", size = 9),
                 legend.text = element_text(color = "black", size = 7.5),
                 strip.text = element_text(size=7.5),
                 strip.background = element_rect(fill="white",colour="gray80"),
                 panel.grid.minor = element_blank(),               #removes minor grid lines
                 panel.grid.major = element_blank(),
                 plot.margin = unit(c(-6,-0.7,-5, -0.5), "lines")) +
                                # top, right, bottom, and left margins
  guides(colour = guide_legend(title.position="top",
                               ncol=2, bycol=T,override.aes = list(size=2))) 
# final plot
(plot1b <- wm_b2) # human- modified habitats

# ----------------------------------------
# now we need to open global maps of non-volant small mammal species richness

#setwd ("C:/Users/topoa/OneDrive/capII/trait_data")
require(raster)
# species richness
SR_POOL <- raster (here ("output","SR_pool_final.tif"))
# Functional diversity
# listing FD files (we ran hypervolume analysis in six different steps
# because it is very demanding in terms of computation time and memory))
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

# FD (previous map)
FD_POOL <- raster (here("output","FD_pool_final.tif"))
FD_POOL[which(values(FD_POOL) >0 )] <- FD_res
# world map to mask
world <- ne_countries(scale = "medium", returnclass = "sf")
masked_SR_POOL <- mask (SR_POOL,world)
masked_FD_POOL <- mask (FD_POOL,world)

#require("fasterize")
#world_raster <- fasterize (world, masked_SR_POOL)

library(rasterVis)
library(rgdal)
library(viridis)

# POOL RICHNESS
plot2 <- gplot(masked_SR_POOL) +
  geom_tile(aes(x=x, y=y, fill=value), alpha=1) + 
  coord_fixed (xlim = c(-180, 180),  ylim = c(-51, 80), ratio = 1) +
  scale_fill_viridis(direction=-1,begin=0.1,
                     breaks = seq (0,50,10),
                     limits=c(0,50),
                     na.value=NA) +
  ggtitle ("Species richness") + 
  theme_classic() 

plot2

plot2b<- plot2 +   theme (axis.title = element_blank(),
                          axis.text = element_blank(),
                          axis.line = element_blank(),
                          axis.ticks = element_blank(),
                          legend.title = element_blank(),
                          legend.position = "bottom",
                          legend.direction="horizontal",
                          legend.justification="top",
                          legend.box="vertical",
                          legend.margin = margin (0,0,0,0),
                          legend.box.margin = margin(1,-10,-10,-10),
                          #legend.title = element_text(color = "black", size = 9,hjust=.5),
                          legend.text = element_text(color = "black", size = 7),
                          strip.text = element_text(size=5),
                          plot.title = element_blank(),#element_text(size=9,hjust = 0.5),
                          legend.key.width=unit(0.6,"cm"),
                          legend.key.size = unit(0.15,"cm")) + 
  guides(fill = guide_colourbar(ticks.colour = "white",
                                ticks.linewidth = 2,
                                frame.colour = "white",
                                frame.linewidth = 0.5))

plot2b

# POOL FUNCCTIONAL DIVERSITY
plot3 <-  gplot(sqrt(masked_FD_POOL)) +
  geom_tile(aes(x=x, y=y, fill=value), alpha=1) + 
  coord_fixed (xlim = c(-180, 180),  ylim = c(-51, 80), ratio = 1) +
  scale_fill_viridis(direction=-1,begin=0.1,
                     breaks = seq (0,
                                   85,
                                   20),
                     na.value=NA) +
  
  ggtitle ("Functional diversity") + 
  theme_classic() 

plot3 

plot3b <- plot3 +   theme (axis.title = element_blank(),
                           axis.text = element_blank(),
                           axis.line = element_blank(),
                           axis.ticks = element_blank(),
                           legend.title = element_blank(),
                           legend.position = "bottom",
                           legend.direction="horizontal",
                           legend.justification="top",
                           legend.box="vertical",
                           legend.margin = margin (0,0,0,0),
                           legend.box.margin = margin(1,-10,-10,-10),
                           #legend.title = element_text(color = "black", size = 9,hjust=.5),
                           legend.text = element_text(color = "black", size = 7),
                           strip.text = element_text(size=5),
                           plot.title = element_blank(),#element_text(size=9,hjust = 0.5),
                           legend.key.width=unit(0.6,"cm"),
                           legend.key.size = unit(0.15,"cm")) + 
  guides(fill = guide_colourbar(ticks.colour = "white",
                                ticks.linewidth = 2,
                                frame.colour = "white",
                                frame.linewidth = 0.5))

plot3b

# arrange maps and scenarios

library("gridExtra")

pdf (here ("output","figures","fig4_500T.pdf"), width=7,height=4.5,family="serif")
grid.arrange(plot2b, plot3b,                               # global maps
             common_legend,
             plot1,plot1b,
             ncol = 8, nrow = 5, 
             layout_matrix = rbind(c(1,1,1,1,2,2,2,2),
                                   c(1,1,1,1,2,2,2,2),
                                   c(4,4,4,4,5,5,5,5),
                                   c(4,4,4,4,5,5,5,5),
                                   c(NA,NA,3,3,3,3,NA,NA) 
             ))

dev.off()


# end