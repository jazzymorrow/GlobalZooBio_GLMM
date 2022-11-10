# Estimating Global Zooplankton Biomass With a Generalised Linear Mixed Model
This repository contains generalised linear mixed models to estimate the global distribution of zooplankton biomass. This work is supervised by Professor Anthony Richardson and is a research project within my Masters in Quantitative Biology. 
 
Zooplankton biomass data used in this project were sourced from the global Coastal and Oceanic Plankton Ecology, Production, and Observation Database (COPEPOD) and from the Australian Zooplankton Biomass database. The complete dataset can be found in `Data/` folder.

The statistical model used in this research was generated in an iterative process whereby the fit, error structure and biological meaning were evaluated at each step. Model 12 (glm12) was selected as the best model and used in this analysis. This selection process can be recreated with following the code in `GlobalZBio_01_Model.R`. The output from each model is stored in `Output/` and the corresponding plots are in `Figures/`. Some residual plots for assumptions are stored in `Misc/`. 

The global zooplankton biomass estimates was produced using `GlobalZBio_02_SurfaceMap.Rmd`. Contour plots of biomass predictions at a fixed longitude can be created using `GlobalZBio_03_ContourPlot.Rmd`. The one-degree gridded satellite SST, satellite Chlorophyll a and bathymetry data used in the mapping can be found in `Data/`. The output map of the prediction is stored in `Figures/` and data can be found in `Output/`.

