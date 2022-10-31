## Produce a map of all observations 
library(tidyverse)
library(ggplot2)
library(sf)
library(patchwork)

## Z biomass data - adjust as we did in the model 
dat <- readRDS("Data/GlobalBiomassData.rds")
dat <- dat %>% 
  mutate(
    Mesh = replace(Mesh, Mesh > 1000, 1000),
    Depth = replace(Depth, Depth > 1500, 1500),
    Depth2 = Depth/1000, #scaled depth variable 
    Bathy = replace(Bathy, Bathy > 7000, 7000),
    SST = replace(SST, SST > 31, 31),
    Biomass = replace(Biomass, Biomass > 10000, 10000)) %>%
  filter(Biomass > 0)


#dat %>%
#  group_by(BiomassMethod) %>%
#  summarise(count = n())

# create new column for figure colour coding 
dat$BiomassMethod2 <- as.character(dat$BiomassMethod)
dat[dat$BiomassMethod2 %in% c("AshfreeDry","Carbon","CarbonCHN"),
    "BiomassMethod2"] <- "Other"


#-------------------------------------------
### PLOT SAMPLE LOCATION WITH PROJECTION 
#-------------------------------------------
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


#axes values when projected???
ggplot() +
  geom_sf(data = landmass, fill = "grey64", 
          color = NA, size = 0.01) +
  geom_sf(data = sf_transformed, aes(col = BiomassMethod2), 
          size = 0.01, alpha = 0.2) + 
  scale_color_manual("Biomass Method", 
                       breaks=c('Displacement','Dry','Settled',
                                'Wet','Other'),
                       values = c("#E69F00", "#56B4E9", "#009E73",
                                  "#0072B2", "#F0E442")) +
  guides(colour = guide_legend(override.aes = list(size = 2, 
                                                   alpha = 1))) +
  theme_bw()

dev.print(pdf, paste0("Figures/", "SampleMap", ".pdf"))
#wesanderson::wes_palette("Zissou1", n = 5)

#---------------------------------------------------------
                # data distribution plots
#---------------------------------------------------------
#Add patchwork once plots are complete 
#should this be on the transformed data 
# DOY vs Latitude
p1 <- ggplot() + geom_point(data = dat, aes(x = DOY,y = Latitude),
                      size = 0.3, alpha = 0.2) +
  labs(x = "Day of year", y = "Latitude") +
  theme_bw() +
  scale_y_continuous(breaks = c(-50,0,50),
                     labels = c(expression(50~degree~S),
                                expression(0~degree),
                                expression(50~degree~N)))

# observation count per year 
p2 <- ggplot() + geom_bar(data = dat, aes(x = Year),
                          col = "grey8", fill = "grey40") +
  labs(y = "Count") + 
  theme_bw()

# mesh size hist 
p3 <- ggplot() + geom_histogram(data = dat, aes(x = Mesh), 
                                bins = 40, fill = "grey40", 
                                col = "grey8") +
  labs(x = expression(paste("Net mesh size (", mu,M, ")")), y = "Count") +
  xlim(0,1000) +
  theme_bw() 

# time of day hist
p4 <- ggplot() + geom_histogram(data = dat, aes(x = TimeLocal), 
                                bins = 24, fill = "grey40", 
                                col = "grey8") +
  labs(x = "Time of day", y = "Count") +
  theme_bw()

# sampling depth hist 
p5 <- ggplot() + geom_histogram(data = dat, aes(x = Depth), 
                                bins = 40, fill = "grey40", 
                                col = "grey8") +
  labs(x = "Depth (m)", y = "Count") +
  xlim(0,1600) +
  theme_bw()

(p1 | p2 )/( p3 | p4 | p5) + plot_annotation(tag_levels = 'A')

ggsave("./Figures/biomass_dat_plots.png")
