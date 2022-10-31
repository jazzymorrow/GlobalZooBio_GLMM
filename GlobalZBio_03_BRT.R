############ Preliminaries ################
library(xgboost)
library(caret)
library(tidyverse)

dat <- readRDS("Data/GlobalBiomassData.rds") 

# Reduce some extreme values, and remove 5 zero values  
dat <- dat %>% 
  mutate(
    Mesh = replace(Mesh, Mesh > 1000, 1000),
    Depth = replace(Depth, Depth > 1500, 1500),
    Depth2 = Depth/1000, #scaled depth variable 
    Bathy = replace(Bathy, Bathy > 7000, 7000),
    SST = replace(SST, SST > 31, 31),
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
                            "Bathy","SST","Chl","TimeLocal","Mesh",
                            "DatasetID","Gear")]) 
test = as.data.frame(dat[-parts, 
                         c("Biomass","BiomassMethod","DOY","Depth",
                           "Bathy","SST","Chl","TimeLocal","Mesh",
                           "DatasetID","Gear")])  

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

#pred_y = predict()

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

# need to look into nrounds and what I should select 

cvmodel <- xgb.cv(params = list(objective = "reg:gamma",
                     eta = 0.1,
                     max_depth = 4),
       data = xgb_train,
       nrounds = 50, 
       nfold = 5)

cvmodel$evaluation_log$test_gamma_nloglik_mean


#function that returns parameter grid - could be generalised a bit better 

eta_vals = c(0.01, 0.05, 0.1)
max_depth_vals = c(2, 4, 6)
min_child_weight_vals = c(1, 3, 5)

eval <- data.frame(expand.grid(eta_vals,max_depth_vals,min_child_weight_vals))
colnames(eval) <- c("eta_vals","max_depth_vals","min_child_weight_vals")
  



#might set up an array 
#eval <- array(NA, dim = c(length(eval[,1]), 600, 2))
#dimnames(save_array)[[1]] 


for (i in 1:length(eval[, 1])) {
  #train model with given parameter combination
  print(i)
  cvmodel <- xgb.cv(
    params = list(
      objective = "reg:gamma",
      eta = eval$eta[i],
      max_depth = eval$max_depth[i],
      min_child_weight = eval$min_child[i]),
    data = xgb_train,
    nrounds = 1000,
    early_stopping_rounds = 50,
    nfold = 5)
  #save model performance
  eval$train_dev[i] <- cvmodel$evaluation_log$train_gamma_nloglik_mean[cvmodel$niter]
  eval$test_dev[i] <- cvmodel$evaluation_log$test_gamma_nloglik_mean[cvmodel$niter]
}

eval
which.min(eval$train_dev)
which.min(eval$test_dev)
best <- eval[which.min(eval$test_dev),]

# first round with trees = 100
#eta max_depth min_child train_dev  test_dev
# 0.1         6         3 0.8851848 0.9214311
 
# second round with 1000 trees 
#0.05         6         5 0.7594031 0.8465268

# third round with 1000 trees and early stopping 
#0.10         6         5 0.7135428 0.8435824

saveRDS(eval, "BRT_eval.RDS")
#---------------------------------------------------------
    # now fix learning rate and vary rounds just to check 
#---------------------------------------------------------
eval2 <- data.frame(ntrees = c(500, 1000, 3000))


for (i in 1:length(eval2$ntrees)) {
  #train model with given parameter combination
  print(i)
  cvmodel <- xgb.cv(
    params = list(
      objective = "reg:gamma",
      eta = best$eta_vals,
      max_depth = best$max_depth_vals,
      min_child_weight = best$min_child_weight_vals),
    data = xgb_train,
    nrounds = eval2$ntrees[i],
    early_stopping_rounds = 50,
    nfold = 5)
  #save model performance
  eval2$train_dev1[i] <- cvmodel$evaluation_log$train_gamma_nloglik_mean[cvmodel$niter]
  eval2$train_dev2[i] <- cvmodel$evaluation_log$train_gamma_nloglik_mean[cvmodel$best_iteration]
  eval2$test_dev2[i] <- cvmodel$evaluation_log$test_gamma_nloglik_mean[cvmodel$niter]
  eval2$test_dev2[i] <- cvmodel$evaluation_log$test_gamma_nloglik_mean[cvmodel$best_iteration]
}

#---------------------------------------------------------
     # retrain with final parameters on test set  
#---------------------------------------------------------
#fit XGBoost model and display training and testing data at each round
model = xgb.train(data = xgb_train, max.depth = 3,
                  objective = "reg:gamma",
                  watchlist=watchlist, nrounds = 500)