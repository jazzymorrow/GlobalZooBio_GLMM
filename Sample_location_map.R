## Produce a map of all observations 

library(tidyverse)
library(ggplot2)
library(maps)

## Z biomass data
dat <- readRDS("Data/GlobalBiomassData.rds") 
## world map
world <- map_data("world")

## create plot of sample locations in dataset
ggplot() +
  geom_map(data = world, map = world,
    aes(long, lat, map_id = region),
    color = "white", fill = "lightgray", 
    size = 0.1) +
  geom_point(data = dat,
    aes(Longitude, Latitude, color = DatasetID),
    alpha = 0.5,size = 0.1, show.legend = FALSE) +
  theme_minimal()


dev.print(pdf, paste0("Figures/", "sample_map", ".pdf"))
