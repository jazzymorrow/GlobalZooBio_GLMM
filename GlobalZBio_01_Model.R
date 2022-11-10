############ Preliminaries ################
library(tidyverse)
library(lme4)
library(ggplot2)
library(splines)
library(DHARMa)
library(MuMIn)
library(visreg)

source("utils.R") #harmonic and visreg functions

#####################################################
##             Data and modifications              ##
#####################################################

dat <- readRDS("Data/GlobalBiomassData.rds") 

## Reduce (but don't remove) some extreme values: 
### as in GlobalZooBio_01_Model.R
dat <- dat %>% 
  mutate(
    HarmTOD = (TimeLocal/24)*2*pi, # Convert to radians
    HarmDOY = (DOY/365)*2*pi, # Convert to radians
    Mesh = replace(Mesh, Mesh > 1000, 1000),
    Depth = replace(Depth, Depth > 1500, 1500),
    Depth2 = Depth/1000, #scaled depth variable 
    Bathy = replace(Bathy, Bathy > 7000, 7000),
    SST = replace(SST, SST > 31, 31),
    NorthHemis = as.factor(ifelse(Latitude >= 0, 1,0)), #hemis indic
    Biomass = replace(Biomass, Biomass > 10000, 10000)) %>% 
  filter(Biomass > 0) ## remove 5 zero biomass measures 

  
########### Prelim exploration of data #############
## check for missing measures 
sum(is.na(dat$Institution)) ##14274 missing institutions
sum(is.na(dat$Project)) ## 55544
sum(is.na(dat$Tow)) #2132 missing 
sum(is.na(dat$DatasetID)) # nothing missing 


## check out institutions/projects/gear
nlevels(dat$ShpCruise)
nlevels(dat$Gear)
nlevels(dat$Project)
nlevels(dat$Institution)
nlevels(dat$DatasetID)

sum(table(dat$Gear, dat$DatasetID)>0) / 
  (nlevels(dat$Gear) * nlevels(dat$DatasetID)) *100 #1.74%
sum(table(dat$Gear, dat$Institution)>0) / 
  (nlevels(dat$Gear) * nlevels(dat$Institution)) *100 #1.91%
sum(table(dat$Gear, dat$Project)>0) / 
  (nlevels(dat$Gear) * nlevels(dat$Project)) * 100 #1.67%

# check out number of gears used in each dataset
tempory <- dat %>% 
  group_by(DatasetID) %>% 
  summarise(N = length(unique(Gear)))

###############################################################
##                  fitting lmms                             ##
###############################################################
## linear mixed model
#StartTime <- Sys.time()
m_linear <- lmer(Biomass ~ BiomassMethod + Mesh + 
                   exp(-Depth/1000)*fHarmonic(HarmTOD, k = 1) +
                   log10(Chl) + ns(Bathy, df = 3)+
                   ns(Bathy, df = 3) +
                   ns(SST, df = 3)*fHarmonic(HarmDOY, k = 1) + 
                   (1|Gear) + 
                   (1|Institution),
                 data = dat)
#EndTime <- Sys.time()

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

#############################################################
## linear mixed model with log10 response variable 
#StartTime <- Sys.time()
m_loglinear <- lmer(log10(Biomass) ~ BiomassMethod + Mesh + 
                   exp(-Depth/1000)*fHarmonic(HarmTOD, k = 1) +
                   log10(Chl) + ns(Bathy, df = 3)+
                   ns(Bathy, df = 3) +
                   ns(Latitude, df = 3)*fHarmonic(HarmDOY, k = 1) + 
                  ns(SST, df = 3) +
                   (1|Gear) + 
                   (1|DatasetID),
                 data = dat)
#EndTime <- Sys.time()
#EndTime - StartTime ## 17 sec

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


###############################################################
##                fitting gamma glmms                        ##
###############################################################
#StartTime <- Sys.time()
glm1 <- glmer(Biomass ~ BiomassMethod + Mesh +
              exp(-Depth/1000)*fHarmonic(HarmTOD, k = 1) + 
              log10(Chl) + ns(Bathy, df = 3) +
              fHarmonic(HarmDOY, k = 1) * ns(SST, df = 3) +
             (1|Gear) + (1|Institution),
           data = dat,
           family = Gamma(link = "log"), nAGQ = 0)
#EndTime <- Sys.time()
#EndTime - StartTime 
## n = 30000: 15.7 minutes, 6.5 secs if nAGQ = 0

#model assessment
r.squaredGLMM(glm1)
summary(glm1)

# DHARMa diagnostics: QQplot 
res <- simulateResiduals(glm1)
plotQQunif(res)

# residuals v fitted in link scale 
plot(residuals(glm1) ~ predict(glm1,type="link"),
     xlab=expression(hat(eta)),ylab="Deviance residuals",
     pch=20,col="blue")

# visualise effect of variables 
fPlotBiomassGLM(glm1, "Biomass_glmm1")

#analysis of deviance - iteratively drop each predictor
drop1(glm1, test = "Chi")

saveRDS(glm1, "Output/glm1.rds")


## random effects ## 
RE <- ranef(glm1)
qqnorm(RE$Institution$`(Intercept)`)
qqnorm(RE$Gear$`(Intercept)`)

############ Find outliers, fit new model to compare ##################
hist(dat$Biomass[dat$Biomass>4000]) #58 measures
hist(dat$Biomass[dat$Biomass>8000]) #15

i_n <- influence(glm1)$hat
halfnorm((i_n))

str(which(residuals(glm1)>16))
length(residuals(glm1)) #182426 total deviance residuals? "88650","141078"

glm2 <- glmer(Biomass ~ BiomassMethod + Mesh +
                exp(-Depth/1000)*fHarmonic(HarmTOD, k = 1) + 
                log10(Chl) + ns(Bathy, df = 3) +
                fHarmonic(HarmDOY, k = 1) * ns(SST, df = 3) +
                (1|Gear) + (1|Institution),
              data = dat[-c(88650, 141078),],
              family = Gamma(link = "log"), nAGQ = 0)

plot(residuals(glm2) ~ predict(glm2,type="link"),
     xlab=expression(hat(eta)),ylab="Deviance residuals",pch=20,col="blue")
#looks pretty good

#compare model fit of glm1 and glm2
summary(glm1)
summary(glm2) 

#estimates all look pretty similar
#outliers clearly aren't overly influential 

#################### GLMM without interactions ######################
glm3 <- glmer(Biomass ~ BiomassMethod + Mesh +
              exp(-Depth/1000) + fHarmonic(HarmTOD, k = 1) + 
              log10(Chl) + ns(Bathy, df = 3) +
              fHarmonic(HarmDOY, k = 1) + ns(SST, df = 3) +
              (1|Gear) + (1|Institution),
            data = dat,
            family = Gamma(link = "log"), nAGQ = 0)

summary(glm3)
anova(glm1,glm3)
## model with interactions has lower deviance chisq = 4596.3
## glm1 best model so far 


################### GLMM with depth:TOD, SST:DOY ###############
glm4 <- glmer(Biomass ~ BiomassMethod + Mesh +
                exp(-Depth/1000):fHarmonic(HarmTOD, k = 1) + 
                log10(Chl) + ns(Bathy, df = 3) +
                fHarmonic(HarmDOY, k = 1):ns(SST, df = 3) +
                (1|Gear) + (1|Institution),
              data = dat,
              family = Gamma(link = "log"), nAGQ = 0)
car::vif(glm4)
car::vif(glm1)
anova(glm1, glm4) #AIC
## no multicollinearity but glmm1 is still better 

################ add tow as a predictor #####################
glm5 <- glmer(Biomass ~ BiomassMethod + Mesh + Tow +
                exp(-Depth/1000)*fHarmonic(HarmTOD, k = 1) + 
                log10(Chl) + ns(Bathy, df = 3) +
                fHarmonic(HarmDOY, k = 1) * ns(SST, df = 3) +
                (1|Gear) + (1|Institution),
              data = dat,
              family = Gamma(link = "log"), nAGQ = 0)
summary(glm5) #Tow coefficient estimates not significant 
summary(glm1) 

################ add hemisphere factor #####################
glm6 <- glmer(Biomass ~ BiomassMethod + Mesh +
                exp(-Depth/1000)*fHarmonic(HarmTOD, k = 1) + 
                log10(Chl) + ns(Bathy, df = 3) +
                fHarmonic(HarmDOY, k = 1) * ns(SST, df = 3) * NorthHemis +
                (1|Gear) + (1|Institution),
              data = dat,
              family = Gamma(link = "log"), nAGQ = 0)
summary(glm6) #BIC: 259994.4, Deviance: 259473.5 
summary(glm1)
anova(glm1, glm6) #Lots more parameters but has lower deviance/BIC???


################ Test different random effects ################
glm7 <- glmer(Biomass ~ BiomassMethod + Mesh +
                exp(-Depth/1000)*fHarmonic(HarmTOD, k = 1) + 
                log10(Chl) + ns(Bathy, df = 3) +
                fHarmonic(HarmDOY, k = 1) *NorthHemis* ns(SST, df = 3) +
                (1|Gear)  + (1|DatasetID),
              data = dat,
              family = Gamma(link = "log"), nAGQ = 0)

summary(glm7) #BIC: 367553.6 Deviance: 367029.5

glm8 <- glmer(Biomass ~ BiomassMethod + Mesh +
                exp(-Depth/1000)*fHarmonic(HarmTOD, k = 1) + 
                log10(Chl) + ns(Bathy, df = 3) +
                fHarmonic(HarmDOY, k = 1) * NorthHemis*ns(SST, df = 3) +
                (1|Institution)  + (1|DatasetID),
              data = dat,
              family = Gamma(link = "log"), nAGQ = 0)
summary(glm8)

r.squaredGLMM(glm1)
r.squaredGLMM(glm8)
r.squaredGLMM(glm7)

## glm7 has lower residual variance, higher deviance, more samples used??
# glm1 R^2: 0.4249192 trigamma, 0.8945343 lognormal
# glm7 has highest R^2: 0.477 trigamma, 0.91 lognormal 

## check random effects of glm7
RE <- ranef(glm7)
qqnorm(RE$DatasetID$`(Intercept)`)
qqnorm(RE$Gear$`(Intercept)`)

fPlotBiomassGLM(glm7, "Biomass_glmm7")


####################### Add lat*lon interaction ########################
glm9 <- glmer(Biomass ~ BiomassMethod + Mesh + 
                ns(Latitude, df = 3)*ns(Longitude, df = 3) +
                exp(-Depth/1000)*fHarmonic(HarmTOD, k = 1) + 
                log10(Chl) + ns(Bathy, df = 3) +
                fHarmonic(HarmDOY, k = 1) * ns(SST, df = 3)* NorthHemis +
                (1|Gear)  + (1|DatasetID),
              data = dat,
              family = Gamma(link = "log"), nAGQ = 0)

summary(glm9)

## SST and latitude are corellated, NorthHemis not required now 


################# Lat * Lon and SST ######################

glm10 <- glmer(Biomass ~ BiomassMethod + Mesh + 
                ns(Latitude, df = 5)*ns(Longitude, df = 5) +
                exp(-Depth2)*fHarmonic(HarmTOD, k = 1) + 
                log10(Chl) + ns(Bathy, df = 3) +
                fHarmonic(HarmDOY, k = 1):ns(Latitude, df = 3) +
                fHarmonic(HarmDOY, k = 1) +
                ns(SST, df = 3) +
                (1|Gear)  + (1|DatasetID),
              data = dat,
              family = Gamma(link = "log"), nAGQ = 0)

summary(glm10)
r.squaredGLMM(glm10)
r.squaredGLMM(glm7)
anova(glm7, glm10)

## glm10 diagnostics 
res <- simulateResiduals(glm10)
plotQQunif(res)
car::vif(glm10)

## plot predictors 
fPlotBiomassGLM(glm10, "Biomass_glmm10")
lat_lon(glm10, "LatLon_glmm10")


## Now increase the df 
glm11 <- glmer(Biomass ~ BiomassMethod + Mesh + 
                 ns(Latitude, df = 7)*ns(Longitude, df = 7) +
                 exp(-Depth2)*fHarmonic(HarmTOD, k = 1) + 
                 log10(Chl) + ns(Bathy, df = 3) +
                 fHarmonic(HarmDOY, k = 1):ns(Latitude, df = 3) +
                 fHarmonic(HarmDOY, k = 1) +
                 ns(SST, df = 3) +
                 (1|Gear)  + (1|DatasetID),
               data = dat,
               family = Gamma(link = "log"), nAGQ = 0)

## plot predictors 
fPlotBiomassGLM(glm11, "Biomass_glmm11")
lat_lon(glm11, "LatLon_glmm11")

anova(glm10, glm11)


#################### remove LatxLon surface ###################
glm12 <- glmer(Biomass ~ BiomassMethod + Mesh + 
                 exp(-Depth2)*fHarmonic(HarmTOD, k = 1) + 
                 log10(Chl) + ns(Bathy, df = 3) +
                 fHarmonic(HarmDOY, k = 1)*ns(Latitude, df = 3) +
                 ns(SST, df = 3) + 
                 (1|Gear) +
                 (1|DatasetID),
               data = dat,
               family = Gamma(link = "log"), nAGQ = 0)

fPlotBiomassGLM(glm12, "Biomass_glmm12")
#lat_lon(glm12, "LatLon_glmm12")

r.squaredGLMM(glm12)
summary(glm12)

## diagnostics: look good
res <- simulateResiduals(glm12)
plotQQunif(res, testDispersion = F, testUniformity = F, 
           testOutliers = F)

plot(residuals(glm12) ~ predict(glm12,type="link"),
     xlab=expression(hat(eta)),ylab="Deviance residuals",
     pch=20,col="blue")

car::vif(glm12, type = 'predictor') #all vifs good except for variables with interactions 

saveRDS(glm12, "Output/glm12.rds") #save output for mapping 
