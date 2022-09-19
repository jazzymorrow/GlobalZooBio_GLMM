# Supplementary functions used in GlobalZBio_old_data.R

##### Harmonic function to fit TOD and DOY #####
fHarmonic <- function (theta, k = 4) {
  X <- matrix(0, length(theta), 2 * k)
  nam <- as.vector(outer(c("c", "s"), 1:k, paste, sep = ""))
  dimnames(X) <- list(names(theta), nam)
  m <- 0
  for (j in 1:k) {
    X[, (m <- m + 1)] <- cos(j * theta)
    X[, (m <- m + 1)] <- sin(j * theta)
  }
  X
}

#### Notes ####
# creates empty matrix (length(theta) x 2k)
# columns are c1 and s1 and contain:
#c1 = cos(theta), s1 = sin(theta)
# k be greater than 1 extends the wave, multiple cycles 



#################################################################
# Function to plot the linear models using visreg: 
# from CMIP6_ZoopStatisticalModel repo
fPlotBiomassLM <- function (mdl, Name, Y_transform = 0) {
  
  # extract terms 
  Terms <- as.character(mdl@call)[2]
  
  # Y_transform: 1 if log10 on the response
  # Fix y label for each plot
  if (Y_transform == 0){
  Y_lab <- expression("Biomass")
  }
  else{Y_lab <- expression("log"[10]*"(Biomass)")}
  
  #set up figure 
  x11(width = 12, height = 6)
  if (length(mdl@frame) <= 8){r <- 2} else{r <- 3}
  par(mfrow = c(r,4), mar = c(4,4,2,2))
  
  ##PLOTS##
  
  if(grepl("BiomassMethod", Terms, fixed = TRUE)) {
    visreg(mdl, "BiomassMethod", rug = FALSE, 
           scale = "response", xlab = "Method", 
           ylab = Y_lab)}
  
  
  if(grepl("Mesh", Terms, fixed = TRUE)) { 
    visreg(mdl, "Mesh", scale = "response", 
           xlab = "Mesh  (microns)", 
           ylab = Y_lab)}
  
  
  if(grepl("Chl", Terms, fixed = TRUE)) {
    visreg(mdl, "Chl", scale = "response", 
           xlab = expression("Chl-a (mg m"^-3*")"), 
           ylab = Y_lab)
  } 
  
  
  if(grepl("Bathy", Terms, fixed = TRUE)) {
    visreg(mdl, "Bathy", scale = "response", 
           xlab = "Bathy (m)", 
           ylab = Y_lab)}
  
  if(grepl("HarmTOD", Terms, fixed = TRUE)) {
    visreg(mdl, "HarmTOD", rug = FALSE, scale = "response", 
           xlab = "Time of Day", xaxt = 'n', 
           ylab = Y_lab)
    axis(side=1, at=c(0, pi/2 , pi, pi + pi/2, pi*2), 
         labels=c("00:00","06:00","12:00","18:00","00:00"))
  }
  
  if(grepl("Depth", Terms, fixed = TRUE)) {
    visreg(mdl, "Depth", scale = "response", 
           xlab = "Depth", 
           ylab = Y_lab)}
  
  
  if(grepl("HarmDOY", Terms, fixed = TRUE)) {
    visreg(mdl, "HarmDOY", scale = "response", 
           xlab = "Day of Year", xaxt = 'n', 
           ylab = Y_lab)
    axis(side=1, at=c(0, pi/2 , pi, pi + pi/2, pi*2), 
         labels=c("1","91","182","273","365"))
  }
  
  if(grepl("SST", Terms, fixed = TRUE)) {
    visreg(mdl, "SST", scale = "response", 
           xlab = "SST (ºC)", ylab = Y_lab)}
  
  
  if(grepl("fHarmonic\\(HarmDOY, k = \\d\\) \\* ns\\(SST, df = \\d\\)", 
           Terms)|
     grepl("ns\\(SST, df = \\d\\) \\* fHarmonic\\(HarmDOY, k = \\d\\)", 
           Terms)){
    visreg(mdl, "HarmDOY", by = "SST",
           type = "conditional",
           scale = "response",
           overlay = TRUE, rug = 0, 
           breaks = c(2, 15, 30), 
           xlab = "Day of Year", 
           ylab = Y_lab,
           strip.names = c("2 ºC", "15 ºC", "30 ºC"),
           xaxt = 'n')
    axis(side=1, at=c(0, pi/2 , pi, pi + pi/2, pi*2), 
         labels=c("1","91","182","273","365"))
  }
  
  
  if(grepl('exp\\(-Depth\\/1000\\) \\* fHarmonic\\(HarmTOD, k = \\d\\)',
           Terms) |
     grepl('fHarmonic\\(HarmTOD, k = \\d\\) \\* exp\\(-Depth\\/1000)', 
           Terms)) {
    visreg(mdl, "HarmTOD", by = "Depth", breaks = c(0, 100, 500), 
           xlab = "Time of Day", ylab = Y_lab,
           type = "conditional", scale = "response", 
           overlay = TRUE, rug = 0,
           strip.names = c("Depth=0","Depth=100", "Depth=500"), 
           xaxt = 'n')
    axis(side=1, at=c(0, pi/2 , pi, pi + pi/2, pi*2), 
         labels=c("00:00","06:00","12:00","18:00","00:00"))
  }
  
  
  ### DO RANDOM EFFECTS
  if(grepl('(1 | DatasetId)', Terms, fixed=TRUE )){
    RE <- ranef(mdl, whichel = "DatasetId", condVar = TRUE)
    dotchart(RE$DatasetId$`(Intercept)`, ylab = "DatasetId")
  }
  
  if(grepl('(1 | Gear)', Terms, fixed=TRUE )){
    RE <- ranef(mdl, whichel = "Gear", condVar = TRUE)
    dotchart(RE$Gear$`(Intercept)`, ylab = "Gear")
    # labels = sort(unique(dat$Gear)))
  }
  
  if(grepl('(1 | Institution)', Terms, fixed=TRUE )){
    RE <- ranef(mdl, whichel = "Institution", condVar = TRUE)
    dotchart(RE$Institution$`(Intercept)`, ylab = "Institution")
    # labels = sort(unique(dat$Institution))
    
  }
  
  if(grepl('(1 | ShpCruise)', Terms, fixed=TRUE )){
    RE <- ranef(mdl, whichel = "ShpCruise", condVar = TRUE)
    dotchart(RE$ShpCruise$`(Intercept)`, ylab = "ShpCruise")
  }
  
  dev.print(pdf, paste0("Figures/", Name, ".pdf"))
}

#####################################################
## function to plot effects of glmm

fPlotBiomassGLM <- function (mdl, Name) {
  
  #Extract model terms
  Terms <- as.character(mdl@call)[2]
  # set figure dimensions/row number 
  x11(width = 12, height = 6)
  r <- ceiling((length(mdl@frame)+1)/4) 
  #plus one because 3-way interactions will have too many plots
  par(mfrow = c(r,4), mar = c(4,4,2,2))
  
  #Biomass method plot
  if(grepl("BiomassMethod", Terms, fixed = TRUE)) {
    visreg(mdl, "BiomassMethod", rug = FALSE, 
           scale = "response", xlab = "Method", 
           ylab = expression("Biomass"))}
  
  ## Mesh size plot 
  if(grepl("Mesh", Terms, fixed = TRUE)) { 
    visreg(mdl, "Mesh", scale = "response", 
           xlab = "Mesh  (microns)", 
           ylab = expression("Biomass"))}
  
  ## Chl plot 
  if(grepl("Chl", Terms, fixed = TRUE)) {
    visreg(mdl, "Chl", scale = "response", 
           xlab = expression("Chl-a (mg m"^-3*")"), 
           ylab = expression("Biomass"))} 
  
  ## Bathymetry plot 
  if(grepl("Bathy", Terms, fixed = TRUE)) {
    visreg(mdl, "Bathy", scale = "response", 
           xlab = "Bathy (m)", 
           ylab = expression("Biomass"))}
  
  ## time of day plot
  if(grepl("HarmTOD", Terms, fixed = TRUE)) {
    visreg(mdl, "HarmTOD", rug = FALSE, 
           scale = "response", xlab = "Time of Day", 
           xaxt = 'n', ylab = expression("Biomass"))
    axis(side=1, at=c(0, pi/2 , pi, pi + pi/2, pi*2), 
         labels=c("00:00","06:00","12:00","18:00","00:00"))
  }
  
  # Depth plot
  if(grepl("Depth", Terms, fixed = TRUE)) {
    visreg(mdl, "Depth", scale = "response", 
           xlab = "Depth (m)", ylab = expression("Biomass"))}
  
  # day of year plot
  if(grepl("HarmDOY", Terms, fixed = TRUE)) {
    visreg(mdl, "HarmDOY", scale = "response", 
           xlab = "Day of Year", xaxt = 'n', 
           ylab = expression("Biomass"))
    axis(side=1, at=c(0, pi/2 , pi, pi + pi/2, pi*2), 
         labels=c("1","91","182","273","365"))
  }
  
  #SST plot
  if(grepl("SST", Terms, fixed = TRUE)) {
    visreg(mdl, "SST", scale = "response", 
           xlab = "SST (ºC)", ylab = expression("Biomass"))}
  
  # DOY x SST plot
  if(grepl("fHarmonic\\(HarmDOY, k = \\d\\) \\* ns\\(SST, df = \\d\\)", Terms) |
     grepl("ns\\(SST, df = \\d\\) \\* fHarmonic\\(HarmDOY, k = \\d\\)", Terms)){
    visreg(mdl, "HarmDOY", by = "SST",
           type = "conditional",
           scale = "response",
           overlay = TRUE, rug = 0, 
           breaks = c(2, 15, 30), 
           xlab = "Day of Year", 
           ylab = expression("Biomass"),
           strip.names = c("2 ºC", "15 ºC", "30 ºC"),
           xaxt = 'n')
    axis(side=1, at=c(0, pi/2 , pi, pi + pi/2, pi*2), 
         labels=c("1","91","182","273","365"))
  }
  
  
  if(grepl('exp\\(-Depth\\/1000\\) \\* fHarmonic\\(HarmTOD, k = \\d\\)', 
           Terms) |
     grepl('fHarmonic\\(HarmTOD, k = \\d\\) \\* exp\\(-Depth\\/1000', 
           Terms)) {
    visreg(mdl, "HarmTOD", 
           by = "Depth", 
           breaks = c(0, 100, 500), 
           xlab = "Time of Day", 
           ylab = expression("Biomass"),
           type = "conditional", 
           scale = "response", 
           overlay = TRUE, 
           rug = 0,
           strip.names = c("Depth=0","Depth=100", "Depth=500"), 
           xaxt = 'n')
    axis(side=1, at=c(0, pi/2 , pi, pi + pi/2, pi*2), 
         labels=c("00:00","06:00","12:00","18:00","00:00"))
  }
  
  # DOY x SST x NorthHemis plot
  #hemisphere factor must be MIDDLE in interaction 
  if(grepl("fHarmonic\\(HarmDOY, k = \\d\\) \\* NorthHemis \\* ns\\(SST, df = \\d\\) ", 
           Terms) |
     grepl("ns\\(SST, df = \\d\\) \\* NorthHemis \\* fHarmonic\\(HarmDOY, k = \\d\\)", 
           Terms)){
    visreg(mdl, "HarmDOY", by = "SST",
           cond = list(NorthHemis = 0),
           type = "conditional",
           scale = "response",
           overlay = TRUE, rug = 0, 
           breaks = c(2, 15, 30), 
           xlab = "Day of Year", 
           ylab = expression("Biomass: Southern Hemisphere"),
           strip.names = c("2 ºC", "15 ºC", "30 ºC"),
           xaxt = 'n')
    axis(side=1, at=c(0, pi/2 , pi, pi + pi/2, pi*2), 
         labels=c("1","91","182","273","365"))
    
    visreg(mdl, "HarmDOY", by = "SST",
           cond = list(NorthHemis = 1),
           type = "conditional",
           scale = "response",
           overlay = TRUE, rug = 0, 
           breaks = c(2, 15, 30), 
           xlab = "Day of Year", 
           ylab = expression("Biomass: Northern Hemisphere"),
           strip.names = c("2 ºC", "15 ºC", "30 ºC"),
           xaxt = 'n')
    axis(side=1, at=c(0, pi/2 , pi, pi + pi/2, pi*2), 
         labels=c("1","91","182","273","365"))
  }
  
  
  ## DO RANDOM EFFECTS
  if(grepl('(1 | DatasetID)', Terms, fixed=TRUE )){
    RE <- ranef(mdl, whichel = "DatasetID", condVar = TRUE)
    dotchart(RE$DatasetID$`(Intercept)`, ylab = "DatasetID")
  }
  
  if(grepl('(1 | Gear)', Terms, fixed=TRUE )){
    RE <- ranef(mdl, whichel = "Gear", condVar = TRUE)
    dotchart(RE$Gear$`(Intercept)`, ylab = "Gear")
    # labels = sort(unique(dat$Gear)))
  }
  
  if(grepl('(1 | Institution)', Terms, fixed=TRUE )){
    RE <- ranef(mdl, whichel = "Institution", condVar = TRUE)
    dotchart(RE$Institution$`(Intercept)`, ylab = "Institution")
    # labels = sort(unique(dat$Institution))
    
  }
  
  dev.print(pdf, paste0("Figures/", Name, ".pdf"))
}


############# Plot three way interaction for N/S hemisphere ##########
visreg(glm9, "HarmDOY", by = "SST",
       cond = list(NorthHemis = 0),
       type = "conditional",
       scale = "response",
       overlay = TRUE, rug = 0, 
       breaks = c(2, 15, 30), 
       xlab = "Day of Year", 
       ylab = expression("Biomass: Southern Hemisphere"),
       strip.names = c("2 ºC", "15 ºC", "30 ºC"),
       xaxt = 'n')

visreg(glm9, "HarmDOY", by = "SST",
       cond = list(NorthHemis = 1),
       type = "conditional",
       scale = "response",
       overlay = TRUE, rug = 0, 
       breaks = c(2, 15, 30), 
       xlab = "Day of Year", 
       ylab = expression("Biomass: Northern Hemisphere"),
       strip.names = c("2 ºC", "15 ºC", "30 ºC"),
       xaxt = 'n')
