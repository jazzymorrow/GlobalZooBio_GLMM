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
    Latitude2 = abs(Latitude), #may not need this 
    Mesh = replace(Mesh, Mesh > 1000, 1000),
    Depth = replace(Depth, Depth > 1500, 1500),
    Bathy = replace(Bathy, Bathy > 7000, 7000),
    SST = replace(SST, SST > 31, 31),
    NorthHemis <- ifelse(dat$Latitude >= 0, 1,0), #binary hemisphere variable
    Biomass = replace(Biomass, Biomass > 10000, 10000)) 



###############################################################
##                  fitting log10 lmm                        ##
###############################################################

## first a log10 LMM
StartTime <- Sys.time()
m_linear <- lmer(log10(Biomass) ~ BiomassMethod + Mesh + 
                   exp(-Depth/1000)*fHarmonic(HarmTOD, k = 1) +
                   log10(Chl) + ns(Bathy, df = 3)+
                   ns(Bathy, df = 3) +
                   ns(SST, df = 3)*fHarmonic(HarmDOY, k = 1) + 
                   #add hemisphere indicator variable??
                   (1|Gear) + 
                   (1|Institution),
                 data = dat)
EndTime <- Sys.time()
lmm_time <- StartTime - EndTime

## produces warning messages about variable scale 
summary(m_linear)

