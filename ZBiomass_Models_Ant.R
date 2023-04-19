# Script for Z biomass models
# 19/4/2023

library(tidyverse)
library(lme4)

source("utils.R") # Harmonic and visreg functions

#### Data wrangling ####
dat <- readRDS("Data/GlobalBiomassDataUpdated.rds")

dat <- dat %>% 
  mutate(
    HarmTOD = (TimeLocal/24)*2*pi, # Convert to radians
    HarmDOY = (DOY/365)*2*pi, # Convert to radians
    Mesh = replace(Mesh, Mesh > 1000, 1000),
    Depth = replace(Depth, Depth > 1500, 1500),
    Depth2 = Depth/1000,  #scaled depth variable
    Bathy = replace(Bathy, Bathy > 7000, 7000),
    Thetao = replace(Thetao, Thetao > 31, 31),
    Chl = replace(Chl, Chl > 3*10^-4, 3*10^-4),
    NorthHemis = as.factor(ifelse(Latitude >= 0, 1,0)), #hemis indic
    Biomass = replace(Biomass, Biomass > 10000, 10000), 
    Biomass = if_else(DatasetID == "IESSNS", Biomass * 1000, Biomass), # IESSNS should be in mg/m3 but was in g/m3
    Biomass = if_else(DatasetID == "NCEI Accession 9500090", Biomass / 1000, Biomass) # NCEI should be in ml/m3 but was in ml/1000m3
    )

# Is N vs S Hemi correct? YES
ggplot(dat, 
       aes(x = NorthHemis, y = Latitude)) +
  geom_point()


#### Models ####
# Random effects  partially crossed partially nested gear within datasetID

m1 <- glmer(Biomass ~ BiomassMethod + 
                 Mesh + 
                 exp(-Depth2)*fHarmonic(HarmTOD, k = 1) + 
                 log10(Chl) + 
                 ns(Bathy, df = 3) +
                 fHarmonic(HarmDOY, k = 1)*ns(Latitude, df = 3) +
                 ns(Thetao, df = 3) + 
                 (1|Gear) +
                 (1|DatasetID),
               data = dat,
               family = Gamma(link = "log"), nAGQ = 0)

fPlotBiomassGLM(m1, "Biomass_m1")

r.squaredGLMM(m1)
summary(m1)

# ## diagnostics: look good
# res <- simulateResiduals(m1)
# graphics.off()
# plotQQunif(res, testDispersion = F, testUniformity = F, 
#            testOutliers = F)
# 
# plot(residuals(m1) ~ predict(m1,type="link"),
#      xlab=expression(hat(eta)),ylab="Deviance residuals",
#      pch=20,col="blue")
# 
# car::vif(m1, type = 'predictor') # All vifs good except for variables with interactions 

saveRDS(m1, "Output/m1_Ant.rds") # Save output for mapping 
