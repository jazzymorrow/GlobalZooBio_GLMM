############ Preliminaries ################
library(xgboost)
library(caret)
library(tidyverse)

dat <- readRDS("Data/GlobalBiomassData.rds") 

# not sure if I should still do this reduction/transformation of predictors 
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
    Biomass = replace(Biomass, Biomass > 10000, 10000)) %>%
  filter(Biomass > 0)

#-------------------------------------------------------
             # TRAIN TEST SPLIT
#-------------------------------------------------------
#make this example reproducible
set.seed(0)

#split into training (80%) and testing set (20%)
#by default, the split uses percentiles of y and has well balanced pred
parts = createDataPartition(dat$Biomass, p = .8, list = F)

# only include a few predictors for now 
train = as.data.frame(dat[parts,
                          c("Biomass","BiomassMethod","DOY","Depth",
                            "Bathy","SST","Chl")]) 
test = as.data.frame(dat[-parts, 
                         c("Biomass","BiomassMethod","DOY","Depth",
                           "Bathy","SST","Chl")])  

#define predictor and response variables in training set
train_x = data.matrix(train[, -1]) #column 1 is biomass 
train_y = train[,1]

#define predictor and response variables in testing set
test_x = data.matrix(test[, -1])
test_y = test[, 1]

#define final training and testing sets
xgb_train = xgb.DMatrix(data = train_x, label = train_y)
xgb_test = xgb.DMatrix(data = test_x, label = test_y)


#----------------------------------------------------------
            # TRAIN MODEL - a few predictors to start 
#----------------------------------------------------------
#define watchlist
watchlist = list(train=xgb_train, test=xgb_test)

#fit XGBoost model and display training and testing data at each round
model = xgb.train(data = xgb_train, max.depth = 3,
                  objective = "reg:gamma",
                  watchlist=watchlist, nrounds = 500)

# plot training error through time 
ggplot(data = model$evaluation_log) + 
  geom_line(aes(x = 1:length(test_gamma_nloglik), 
                y = test_gamma_nloglik), col = "red") +
  geom_line(aes(x = 1:length(test_gamma_nloglik), 
                y = train_gamma_nloglik), col = "green") +
  ylim(0,10)
# seems like gradient is very small after 200 rounds 
# way better with gamma specified - 0.99 in first 50, 0.92 at 500


#-------------------------------------------------------
                # check accuracy metrics?
#-------------------------------------------------------
mean((test_y - pred_y)^2) #mse
MAE(test_y, pred_y) #mae
RMSE(test_y, pred_y) #rmse


#-------------------------------------------------------
                  # try dismo package 
#-------------------------------------------------------
dismo::gbm.step(data = train,
  gbm.x = -1, gbm.y = 1,
                family = "gaussian", n.folds = 5,
                learning.rate = 0.01)


#-------------------------------------------------------
                 # tune hyperparameters 
#-------------------------------------------------------
# prepare a grid of parameters for tuning
# this is taken from Drago et al. 2022
#https://github.com/dlaetitia/Global_zooplankton_biomass_distribution/blob/main/Scripts/1.model-fit_habitat_models.R

eta_values = c(0.05, 0.075, 0.1)
max_depth_values = c(2, 4, 6)
min_child_weight_values = c(1, 3, 5)
rs_grid <- param_grid(rs,
# the CV resamples we will run for each parameter combination
                      eta = eta_values,
                      # learning rate
                      max_depth = max_depth_values,
                      # maximum depth of trees
                      min_child_weight = min_child_weight_values)