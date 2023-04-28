### Preliminaries ----
library(tidyverse)
library(lme4)
library(ggplot2)
library(splines)
library(DHARMa)
library(MuMIn)
library(visreg)

source("utils.R") #harmonic and visreg functions


## Data and modifications ----

dat <- readRDS("Data/GlobalBiomassDataUpdated.rds") 


## Reduce (but don't remove) some extreme values: 
# as in (Heneghan et al., 2020)

dat <- dat %>% 
  mutate(
    HarmTOD = (TimeLocal/24)*2*pi, # Convert to radians
    HarmDOY = (DOY/365)*2*pi, # Convert to radians
    Mesh = replace(Mesh, Mesh > 1000, 1000),
    Depth = replace(Depth, Depth > 1500, 1500),
    Depth2 = Depth/1000,  #scaled depth variable 
    Bathy = replace(Bathy, Bathy > 7000, 7000),
    Thetao = replace(Thetao, Thetao > 31, 31),
    NorthHemis = as.factor(ifelse(Latitude >= 0, 1,0)), #hemis indic
    Biomass = replace(Biomass, Biomass > 10000, 10000), 
    # IESSNS should be in mg/m3 but was in g/m3
    Biomass = if_else(DatasetID == "IESSNS", 
                      Biomass * 1000, Biomass), 
    # NCEI should be in ml/m3 but was in ml/1000m3
    Biomass = if_else(DatasetID == "NCEI Accession 9500090", 
                      Biomass / 1000, Biomass) 
    ) 


  
### Prelim exploration of data ----

## check for missing measures 
sum(is.na(dat$Institution)) ##14547 missing institutions
sum(is.na(dat$Project)) ## 56811
sum(is.na(dat$Tow)) #3194 missing 
sum(is.na(dat$DatasetID)) # nothing missing 


## check out institutions/projects/gear
nlevels(dat$ShpCruise)
nlevels(dat$Gear)
nlevels(dat$Project)
nlevels(dat$Institution)
nlevels(dat$DatasetID)

sum(table(dat$Gear, dat$DatasetID)>0) / 
  (nlevels(dat$Gear) * nlevels(dat$DatasetID)) *100 #2.11%
sum(table(dat$Gear, dat$Institution)>0) / 
  (nlevels(dat$Gear) * nlevels(dat$Institution)) *100 #2.22%
sum(table(dat$Gear, dat$Project)>0) / 
  (nlevels(dat$Gear) * nlevels(dat$Project)) * 100 #1.95%

## check out number of gears used in each dataset
temp <- dat %>% 
  group_by(DatasetID) %>% 
  summarise(N = length(unique(Gear)))



### LMMs ----

#### No transformation ----
m_linear <- lmer(Biomass ~ BiomassMethod + Mesh + 
                   exp(-Depth/1000)*fHarmonic(HarmTOD, k = 1) +
                   log10(Chl) + ns(Bathy, df = 3)+
                   ns(Bathy, df = 3) +
                   ns(SST, df = 3)*fHarmonic(HarmDOY, k = 1) + 
                   (1|Gear) + 
                   (1|Institution),
                 data = dat)


# Residual normality check
qqnorm(residuals(m_linear))
qqline(residuals(m_linear))

# variance homogeneity check
plot(m_linear)

## plot of effects 
fPlotBiomassLM(m_linear, "BiomassLM", Y_transform = 0)

## Random Effects ##
RE <- ranef(m_linear)
dotplot.ranef.mer(RE)$Institution ##check plots of Random effects 

saveRDS(m_linear, "Output/m_linear.rds")


#### Log10 response variable ----

m_loglinear <- lmer(log10(Biomass) ~ BiomassMethod + Mesh + 
                   exp(-Depth/1000)*fHarmonic(HarmTOD, k = 1) +
                   log10(Chl) + ns(Bathy, df = 3)+
                   ns(Bathy, df = 3) +
                   ns(Latitude, df = 3)*fHarmonic(HarmDOY, k = 1) + 
                  ns(SST, df = 3) +
                   (1|Gear) + 
                   (1|DatasetID),
                 data = dat)


## produces warning messages about variable scale 
summary(m_loglinear)
anova(m_loglinear)

# Residual normality check
qqnorm(residuals(m_loglinear))
qqline(residuals(m_loglinear))

# variance homogeneity check
plot(m_loglinear)
## plot of effects 
fPlotBiomassLM(m_loglinear, "BiomassLogLM", Y_transform = 1)

## Random Effects ##
RE <- ranef(m_loglinear)
qqnorm(RE$Institution$`(Intercept)`)
qqnorm(RE$Gear$`(Intercept)`)
dotplot.ranef.mer(RE) ##check plots of Random effects 

saveRDS(m_loglinear, "Output/m_loglinear.rds")



### FINAL: Gamma GLMM  ----

#random effects are partially crossed partially nested gear within datasetID

glmFinal <- glmer(Biomass ~ BiomassMethod + 
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

## visreg plots of glmm
fPlotBiomassGLM(glmFinal, "Biomass_glmmFinal")

## model summary
r.squaredGLMM(glmFinal)
summary(glmFinal)

## diagnostics: look good
res <- simulateResiduals(glmFinal)
plotQQunif(res, testDispersion = F, testUniformity = F, 
           testOutliers = F)

plot(residuals(glmFinal) ~ predict(glmFinal,type="link"),
     xlab=expression(hat(eta)),ylab="Deviance residuals",
     pch=20,col="blue")

## check for collinearity issues 
car::vif(glmFinal, type = 'predictor') 
#all vifs good except for variables with interactions 

## save model for mapping 
saveRDS(glmFinal, "Output/glmFinal.rds") 
