source("utils.R") # load the harmonic function

library(raster)
library(ggplot2)
library(ncdf4) #?
library(dplyr)
library(tidyverse)
library(splines)
#library(lmerTest)

# IMPORT COPEPOD DATABASE BIOMASS GLMM
glm1 <- readRDS(file.path("Output","glm1.rds"))

## IMPORT BATHYMETRY DATA AND ORIENT LATITUDES TO BE SOUTH TO NORTH
bathy_data <- readRDS(file.path("Data","Bathy_raster_oneDeg.rds"))
bathy_matrix <- t(as.matrix(bathy_data$Bathy))
bathy_matrix <- bathy_matrix[,180:1]


## SOME CODE TO EXPLORE THE RASTER FILE
plot(bathy_data) #check out bathy map 
#prep your raster for ggplot by converting to a dataframe
bathy.df <- as.data.frame(bathy_data, xy = TRUE) 
ggplot() +  #gglot version of map 
  geom_raster(data = bathy.df, aes(x = x, y = y, fill = Bathy)) +
  scale_fill_viridis_c(na.value="#000000")

length(bathy.df$Bathy)
sum(is.na(bathy.df$Bathy)) #NAs: 22079 - because earth has land...

