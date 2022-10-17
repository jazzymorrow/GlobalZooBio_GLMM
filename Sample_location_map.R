## Produce a map of all observations 
library(tidyverse)
library(ggplot2)
library(maps)
library(sf)

## Z biomass data
dat <- readRDS("Data/GlobalBiomassData.rds")

## world map
world <- map_data("world")

## create plot of sample locations in dataset: default projection
ggplot() +
  geom_map(data = world, map = world,
    aes(long, lat, map_id = region),
    color = "lightgray", fill = "lightgray", 
    size = 0.1) +
  geom_point(data = dat,
    aes(Longitude, Latitude, fill = "black"),
    alpha = 0.4,size = 0.1, show.legend = FALSE) +
  labs(x = "Latitude", "Longitude")
  theme_minimal()


dev.print(pdf, paste0("Figures/", "sample_map", ".pdf"))


### try with projection
# Convert into a sf
sf <- dat %>% 
  sf::st_as_sf(coords=c("Longitude", "Latitude"))

# longitude-latitude projection
lonlat <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0" 
sf::st_crs(sf) <- lonlat

landmass <- rnaturalearth::ne_countries(scale = "large") %>% 
  # get the landmass like this
  sf::st_as_sf(crs = latlon)

#Mollweide projection 
moll <- "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m no_defs"

#transform dat using new projection
sf_transformed <- sf %>% 
  sf::st_transform(crs = moll)

# Transform the landmass to the same projection...
landmass <- landmass %>% 
  sf::st_transform(crs = moll)


ggplot() +
  geom_sf(data = landmass, fill = "grey64", 
          color = NA, size = 0.01) +
  geom_sf(data = sf_transformed, aes(col = BiomassMethod), 
          size = 0.01, alpha = 0.4) + 
  guides(colour = guide_legend(override.aes = list(size=3))) +
  theme_bw()
