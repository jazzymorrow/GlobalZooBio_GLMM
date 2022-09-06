############ Preliminaries ################
library(tidyverse)
library(lme4)
library(ggplot2)
library(splines)
library(DHARMa)
library(MuMIn)
library(visreg)

source("fHarmonic.R") #harmonic functions

###########################################
##        Data and modifications         ##
###########################################

dat <- readRDS("Data/GlobalBiomassData.rds") 

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

## remove zero biomass measures 
length(dat$Biomass[dat$Biomass==0]) #5 zero biomass
dat <- dat %>%
  filter(Biomass > 0)

###############################################################
##                  fitting lmms                             ##
###############################################################
## linear mixed model
StartTime <- Sys.time()
m_linear <- lmer(Biomass ~ BiomassMethod + Mesh + 
                   exp(-Depth/1000)*fHarmonic(HarmTOD, k = 1) +
                   log10(Chl) + ns(Bathy, df = 3)+
                   ns(Bathy, df = 3) +
                   ns(SST, df = 3)*fHarmonic(HarmDOY, k = 1) + 
                   #add hemisphere indicator variable??
                   (1|Gear) + 
                   (1|Institution),
                 data = dat)
EndTime <- Sys.time()

# Residual normality check
qqnorm(residuals(m_linear))
qqline(residuals(m_linear))
# variance homogeneity check
plot(m_linear)
## plot of effects 
fPlotBiomassLM(m_linear, "BiomassLM", Y_transform = 0)


## linear mixed model with log10 response variable 
StartTime <- Sys.time()
m_loglinear <- lmer(log10(Biomass) ~ BiomassMethod + Mesh + 
                   exp(-Depth/1000)*fHarmonic(HarmTOD, k = 1) +
                   log10(Chl) + ns(Bathy, df = 3)+
                   ns(Bathy, df = 3) +
                   ns(SST, df = 3)*fHarmonic(HarmDOY, k = 1) + 
                   #add hemisphere indicator variable??
                   (1|Gear) + 
                   (1|Institution),
                 data = dat)
EndTime <- Sys.time()
EndTime - StartTime ## 17 sec

## produces warning messages about variable scale 
summary(m_loglinear)
anova(m_loglinear)

# Residual normality check
qqnorm(residuals(m_loglinear))
qqline(residuals(m_loglinear))
# variance homogeneity check
plot(m_linear)
## plot of effects 
fPlotBiomassLM(m_loglinear, "BiomassLogLM", Y_transform = 1)

###############################################################
##                fitting gamma glmms                        ##
###############################################################
StartTime <- Sys.time()
m1 <- glmer(Biomass ~ BiomassMethod + Mesh +
              exp(-Depth/1000)*fHarmonic(HarmTOD, k = 1) + 
              log10(Chl) + ns(Bathy, df = 3) +
              fHarmonic(HarmDOY, k = 1) * ns(SST, df = 3) +
             (1|Gear) + (1|Institution),
           data = dat,
           family = Gamma(link = "log"), nAGQ = 0)
EndTime <- Sys.time()
EndTime - StartTime 
## n = 30000: 15.7 minutes, 6.5 secs if nAGQ = 0 ???
## full data set runs in 1 min if nAGQ = 0


##??outliers
r.squaredGLMM(m1)
summary(m1)
plot(m1)

## DHARMa diagnostics 
#res <- simulateResiduals(m1, plot = T)

# visualise effect of variables 
fPlotBiomassGLM(m1, "Biomass_glmm1")

#analysis of deviance - iteratively drop each predictor
drop1(m1, test = "Chi")
