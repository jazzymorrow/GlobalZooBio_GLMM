############ Preliminaries ################
library(tidyverse)
library(lme4)
library(ggplot2)
library(splines)
library(DHARMa)
library(MuMIn)


###########################################
##        Data and modifications         ##
###########################################

dat <- readRDS("GlobalBiomassData.rds") ## need to relocate this 

## Reduce (but don't remove) some extreme values: as in GlobalZooBio_01_Model.R

dat <- dat %>% 
  mutate(
    HarmTOD = (TimeLocal/24)*2*pi, # Convert to radians
    HarmDOY = (DOY2/365)*2*pi, # Convert to radians
    Latitude2 = abs(Latitude),
    Mesh = replace(Mesh, Mesh > 1000, 1000),
    Depth = replace(Depth, Depth > 1500, 1500),
    Bathy = replace(Bathy, Bathy > 7000, 7000),
    SST = replace(SST, SST > 31, 31),
    Biomass = replace(Biomass, Biomass > 10000, 10000)) %>% 
  unite(Gear_Mesh, Gear, Mesh, remove = FALSE) %>%  # Calculate Gear_Mesh
  mutate(Gear_Mesh = as.factor(Gear_Mesh)) %>% 
  droplevels()
##add NHemis column