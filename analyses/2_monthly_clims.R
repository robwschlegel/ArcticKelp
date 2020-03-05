# analyses/2_monthly_clims.R
# The purpose of this script is to download and prep Bio-Oracle data
# This should only need to be run once at the outset of the project.
# If any sites change/are added, run "analyses/3_study_sites_clims.R" 
# to get the updated values, not this script

# The Bio-Oracle (http://www.bio-oracle.org/code.php)  are downloaded directly 
# from their server and so may be accessed from anywhere


# Setup -------------------------------------------------------------------

# The study sites and bounding box
source("analyses/1_study_sites.R")

# Bio-Oracle access
library(sdmpredictors)

# Other packages
library(tidync)
library(stringr)
library(data.table)
library(FNN)

# Set cores
doParallel::registerDoParallel(cores = 50)

# A rainbow colour palette was explicitly requested
rainbow_palette <- c("#fefefe", "#f963fa", "#020135", "#00efe1", "#057400", "#fcfd00", "#ed0000", "#3d0000")


# Download Bio-ORACLE data ------------------------------------------------

# Explore datasets in the package
list_datasets()

# Explore layers in a dataset
BO_layers <- list_layers(datasets = "Bio-ORACLE")
MAR_layers <- list_layers(datasets = "MARSPEC")

# Check layer statistics
layer_stats()

# Check Pearson correlation coefficient between layers
layers_correlation() 

# Download bathy layers
  # NB: Not used as these values are suspicious
# bathy <- load_layers(c("BO_bathymin", "BO_bathymean", "BO_bathymax"))
# Arctic_bathy <- as.data.frame(bathy, xy = T) %>% 
#   dplyr::rename(lon = x, lat = y) %>% 
#   filter(lon >= bbox_arctic[1], lon <= bbox_arctic[2],
#          lat >= bbox_arctic[3], lat <= bbox_arctic[4])
# save(Arctic_bathy, file = "data/Arctic_bathy.RData")

## The layes currently chosen for use in this study
  # NB: Many of the variables below also have surface values
  # NB: Max depth is the deepest, we checked
# Calcite - mean - BO_calcite
# Diffuse attenuation coefficient at 490 nm - mean - BO_damean
# Photosynthetically available radiation - mean - BO_parmean
# pH - mean - BO_ph
# Chl con. - mean at max depth - BO2_chlomean_bdmax
# Current velocity - mean at max depth - BO2_curvelmean_bdmax
# Dissolved oxygen - mean at max depth - BO2_dissoxmean_bdmax
# Iron con. - mean at max depth -	BO2_ironmean_bdmax
# Phos con. - mean at max depth - BO2_phosphatemean_bdmax
# Light at bottom - mean at max depth - BO2_lightbotmean_bdmax
# Nitr con. - mean at max depth - BO2_nitratemean_bdmax
# sea temp. - mean at max depth - BO2_tempmean_bdmax
# Carbon phytoplankton biomass - mean at max depth - BO2_carbonphytomean_bdmax
# Primary production - mean at max depth - BO2_ppmean_bdmax
# sea salinity - mean at max depth - BO2_salinitymean_bdmax
# Silicate con. - mean at max depth - BO2_silicatemean_bdmax
# Ice con. - mean - BO2_icecovermean_ss
# Ice thickness - mean + range - BO2_icethickmean_ss + BO2_icethickrange_ss

## Download the chosen layers
  # NB: Don't run this if nothing has changed as there is no need to ping their servers
                              # Bottom temperature
BO_layers_dl <- load_layers(c("BO2_templtmin_bdmax", "BO2_tempmean_bdmax", "BO2_templtmax_bdmax",
                              # Surface temperature
                              "BO2_templtmin_ss", "BO2_tempmean_ss", "BO2_templtmax_ss", 
                              # Bottom salinity
                              "BO2_salinityltmin_bdmax", "BO2_salinitymean_bdmax", "BO2_salinityltmax_bdmax", 
                              # Surface salinity
                              "BO2_salinityltmin_ss", "BO2_salinitymean_ss", "BO2_salinityltmax_ss", 
                              # Ice cover
                              "BO2_icecoverltmin_ss", "BO2_icecovermean_ss", "BO2_icecoverltmax_ss", 
                              # Ice thickness
                              "BO2_icethickltmin_ss", "BO2_icethickmean_ss", "BO2_icethickltmax_ss", 
                              # Photosynthetically active radiation
                              "BO_parmean", "BO_parmax", 
                              # Dissolve oxygen
                              "BO2_dissoxltmin_bdmax", "BO2_dissoxmean_bdmax", "BO2_dissoxltmax_bdmax", 
                              # pH
                              "BO_ph", 
                              # Calcite
                              "BO_calcite",
                              # Iron
                              "BO2_ironltmin_bdmax", "BO2_ironmean_bdmax", "BO2_ironltmax_bdmax", 
                              # Nitrate
                              "BO2_nitrateltmin_bdmax", "BO2_nitratemean_bdmax", "BO2_nitrateltmax_bdmax", 
                              # Phosphate
                              "BO2_phosphateltmin_bdmax", "BO2_phosphatemean_bdmax", "BO2_phosphateltmax_bdmax", 
                              # Silicate
                              "BO2_silicateltmin_bdmax", "BO2_silicatemean_bdmax", "BO2_silicateltmax_bdmax", 
                              # Current velocity
                              "BO2_curvelltmin_bdmax", "BO2_curvelmean_bdmax", "BO2_curvelltmax_bdmax"
                              ))
BO_layers_df <- as.data.frame(BO_layers_dl, xy = T)
Arctic_BO <- BO_layers_df %>%
  dplyr::rename(lon = x, lat = y) %>%
  filter(lon >= bbox_arctic[1], lon <= bbox_arctic[2],
         lat >= bbox_arctic[3], lat <= bbox_arctic[4])
save(Arctic_BO, file = "data/Arctic_BO.RData")

# Visualise
ggplot(Arctic_BO, aes(x = lon, y = lat)) +
  geom_tile(aes(fill = BO2_icecovermean_ss)) +
  borders(fill = "grey70", colour = "black") +
  scale_fill_viridis_c(option = "D") +
  coord_cartesian(xlim = c(bbox_arctic[1], bbox_arctic[2]),
                  ylim = c(bbox_arctic[3], bbox_arctic[4])) +
  theme(legend.position = "bottom")


# Download future Bio-ORACLE data -----------------------------------------

# Future scenario conversions
# https://ar5-syr.ipcc.ch/topic_futurechanges.php
# RCP8.5 ~= A2
# RCP6.0 ~= B2
# RCP4.5 ~= B1

list_layers_future()

get_future_layers()


# Load GMED ASCII files ---------------------------------------------------

# The .asc files loaded below were downloaded from:
# http://gmed.auckland.ac.nz/download.html

# Depth
depth <- read.asciigrid("data/gb_depth.asc")
depth_df <- as.data.frame(depth, xy = T)
Arctic_depth <- depth_df %>%
  dplyr::rename(bathy = data.gb_depth.asc,
                lon = s1, lat = s2) %>%
  filter(lon >= bbox_arctic[1], lon <= bbox_arctic[2],
         lat >= bbox_arctic[3], lat <= bbox_arctic[4])

# Slope
slope <- read.asciigrid("data/gb_slope.asc")
slope_df <- as.data.frame(slope, xy = T)
Arctic_slope <- slope_df %>%
  dplyr::rename(slope = data.gb_slope.asc,
                lon = s1, lat = s2) %>%
  filter(lon >= bbox_arctic[1], lon <= bbox_arctic[2],
         lat >= bbox_arctic[3], lat <= bbox_arctic[4])

# Distance to land
land_distance <- read.asciigrid("data/gb_land_distance.asc")
land_distance_df <- as.data.frame(land_distance, xy = T)
Arctic_land_distance <- land_distance_df %>%
  dplyr::rename(land_distance = data.gb_land_distance.asc,
                lon = s1, lat = s2) %>%
  filter(lon >= bbox_arctic[1], lon <= bbox_arctic[2],
         lat >= bbox_arctic[3], lat <= bbox_arctic[4])

# Combine into single file and save
Arctic_GMED <- left_join(Arctic_land_distance, Arctic_depth, by = c("lon", "lat")) %>% 
  left_join(Arctic_slope, by = c("lon", "lat")) %>% 
  dplyr::select(lon, lat, everything())
save(Arctic_GMED, file = "data/Arctic_GMED.RData")


# Combine BO and GMED -----------------------------------------------------

# Add an index for ease of joining
Arctic_BO$BO_index <- 1:nrow(Arctic_BO)

# Nearest neighbour search of GMED against BO
Arctic_env <- Arctic_GMED %>%
  mutate(BO_index = as.vector(knnx.index(as.matrix(Arctic_BO[,c("lon", "lat")]),
                                         as.matrix(Arctic_GMED[,c("lon", "lat")]), k = 1))) %>% 
  left_join(Arctic_BO, by = "BO_index") %>% 
  dplyr::select(-lon.x, -lat.x, -BO_index) %>% 
  dplyr::rename(lon = lon.y, lat = lat.y) %>% 
  group_by(lon, lat) %>% 
  summarise_all(mean, na.rm = T) %>% 
  ungroup() %>% 
  mutate(bathy = replace_na(bathy, NA),
         slope = replace_na(slope, NA),
         bathy = -bathy)
save(Arctic_env, file = "data/Arctic_env.RData")

# Visualise
ggplot(Arctic_env, aes(x = lon, y = lat)) +
  geom_tile(aes(fill = bathy)) +  
  borders(fill = "grey70", colour = "black") +
  # scale_fill_viridis_c(option = "E") +
  coord_cartesian(xlim = c(bbox_arctic[1], bbox_arctic[2]),
                  ylim = c(bbox_arctic[3], bbox_arctic[4])) +
  theme(legend.position = "bottom")

