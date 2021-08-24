# load packages

# maps
require (raster)
require(rgdal)
require(maps)
require(sp)
require(rgeos)

# data formatting and transformation
require(here)
require(vegan)
require(reshape2)
require (reshape)
require(clipr)# clear_clip()

# analyses
require(parallel)
require(hypervolume)
# install "probpool" Available here: https://rdrr.io/github/ChrKoenig/probpool/
#install.packages("remotes") 
#remotes::install_github("ChrKoenig/probpool")
require(probpool)

# plots
library(rasterVis)
library(viridis)
require(ggplot2)
require("rnaturalearth")
require("rnaturalearthdata")
library("gridExtra")
library(ggpubr)
require(grid)