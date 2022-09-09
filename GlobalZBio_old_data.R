############ Preliminaries ################
library(tidyverse)
library(lme4)
library(ggplot2)
library(splines)
library(DHARMa)
library(MuMIn)
library(visreg)

source("utils.R") #harmonic functions

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
    NorthHemis = as.factor(ifelse(Latitude >= 0, 1,0)), #hemisphere variable
    Biomass = replace(Biomass, Biomass > 10000, 10000)) 

## remove zero biomass measures 
length(dat$Biomass[dat$Biomass==0]) #5 zero biomass
dat <- dat %>%
  filter(Biomass > 0)

## check for missing measures 
sum(is.na(dat$Institution)) ##14274 missing institutions
sum(is.na(dat$Project)) ## 55544
sum(is.na(dat$Tow)) #2132 missing 

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
                   #add hemisphere indicator variable??
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
                   ns(SST, df = 3)*fHarmonic(HarmDOY, k = 1) + 
                   #add hemisphere indicator variable??
                   (1|Gear) + 
                   (1|Institution),
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
plot(m_linear)
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
StartTime <- Sys.time()
glm1 <- glmer(Biomass ~ BiomassMethod + Mesh +
              exp(-Depth/1000)*fHarmonic(HarmTOD, k = 1) + 
              log10(Chl) + ns(Bathy, df = 3) +
              fHarmonic(HarmDOY, k = 1) * ns(SST, df = 3) +
             (1|Gear) + (1|Institution),
           data = dat,
           family = Gamma(link = "log"), nAGQ = 0)
EndTime <- Sys.time()
EndTime - StartTime 
## n = 30000: 15.7 minutes, 6.5 secs if nAGQ = 0

#model assessment
r.squaredGLMM(glm1)
summary(glm1)

# DHARMa diagnostics: QQplot 
res <- simulateResiduals(glm1)
plotQQunif(res)
#qqnorm(residuals(glm1))
qqline(residuals(glm1))

# visualise effect of variables 
fPlotBiomassGLM(glm1, "Biomass_glmm1")

#analysis of deviance - iteratively drop each predictor
drop1(glm1, test = "Chi")

saveRDS(glm1, "Output/glmm1.rds")

# residuals v fitted in link scale 
plot(residuals(glm1) ~ predict(glm1,type="link"),
     xlab=expression(hat(eta)),ylab="Deviance residuals",pch=20,col="blue")

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
summary(glm2) #estimates all look pretty similar

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
glm6 <- glmer(Biomass ~ BiomassMethod + Mesh + NorthHemis +
                exp(-Depth/1000)*fHarmonic(HarmTOD, k = 1) + 
                log10(Chl) + ns(Bathy, df = 3) +
                fHarmonic(HarmDOY, k = 1) * ns(SST, df = 3) +
                (1|Gear) + (1|Institution),
              data = dat,
              family = Gamma(link = "log"), nAGQ = 0)
summary(glm6)
summary(glm1)
anova(glm1, glm6) #significant but not a huge change in AIC/DEV
