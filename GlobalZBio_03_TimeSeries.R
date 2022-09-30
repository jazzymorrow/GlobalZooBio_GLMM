############ Preliminaries ################
library(tidyverse)
library(lme4)
library(ggplot2)
library(splines)
#library(DHARMa)
#library(MuMIn)
library(visreg)
library(mgcv)

source("utils.R") #harmonic and visreg functions

dat <- readRDS("Data/GlobalBiomassData.rds") 

dat <- dat %>% 
  mutate(
    HarmTOD = (TimeLocal/24)*2*pi, # Convert to radians
    HarmDOY = (DOY/365)*2*pi, # Convert to radians
    Mesh = replace(Mesh, Mesh > 1000, 1000),
    Depth = replace(Depth, Depth > 1500, 1500),
    Depth2 = Depth/1000, #scaled depth variable 
    Bathy = replace(Bathy, Bathy > 7000, 7000),
    SST = replace(SST, SST > 31, 31),
    NorthHemis = as.factor(ifelse(Latitude >= 0, 1,0)), #hemisphere variable
    Biomass = replace(Biomass, Biomass > 10000, 10000)) 

#data below -50 degrees Latitude
dat_50S <- dat %>%
  filter(Latitude < -50)

ggplot() + geom_point(data = dat_50S, aes(x = Year, y = Biomass)) +
  scale_y_log10()

gam1 <- gam(Biomass ~ s(Year) + fHarmonic(HarmDOY, k = 1),
            family = Gamma("log"),
            data = dat_50S)


