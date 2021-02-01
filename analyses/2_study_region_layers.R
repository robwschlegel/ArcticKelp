# analyses/2_study_region_layers.R
# The purpose of this script is to download and prep Bio-Oracle data
# This should only need to be run once at the outset of the project.
# If any sites change/are added, run "analyses/3_study_site_layers.R" 
# to get the updated values, not this script

# The Bio-Oracle (http://www.bio-oracle.org/code.php) are downloaded directly 
# from their server and so may be accessed from anywhere


# Setup -------------------------------------------------------------------

# The study sites and bounding box
source("analyses/1_study_region_sites.R")

# Bio-Oracle access
library(sdmpredictors)

# Other packages
library(tidync)
library(stringr)
library(data.table)
library(FNN)
library(correlation)

# Set cores
doParallel::registerDoParallel(cores = 50)

# Disable scientific notation
options(scipen = 999)

# A rainbow colour palette was explicitly requested
rainbow_palette <- c("#fefefe", "#f963fa", "#020135", "#00efe1", "#057400", "#fcfd00", "#ed0000", "#3d0000")

# Convenience function for loading .asc files
load_asc <- function(file_name, col_name){
  df <- as.data.frame(raster(file_name), xy = T) %>% 
    `colnames<-`(c("lon", "lat", col_name))  #%>% 
    # mutate(lon = round(lon, 4),
           # lat = round(lat, 4))
}


# Download Bio-ORACLE data ------------------------------------------------

# Explore datasets in the package
list_datasets()

# Explore layers in a dataset
BO_layers <- list_layers(datasets = "Bio-ORACLE")
MAR_layers <- list_layers(datasets = "MARSPEC")

# Check layer statistics
layer_stats()

# Download bathy layers
  # NB: Not used as these values are suspicious
# bathy <- load_layers(c("BO_bathymin", "BO_bathymean", "BO_bathymax"))
# Arctic_bathy <- as.data.frame(bathy, xy = T) %>% 
#   dplyr::rename(lon = x, lat = y) %>% 
#   filter(lon >= bbox_arctic[1], lon <= bbox_arctic[2],
#          lat >= bbox_arctic[3], lat <= bbox_arctic[4])
# save(Arctic_bathy, file = "data/Arctic_bathy.RData")

# The layers currently chosen for use in this study
  # NB: Many of the variables below also have surface values
  # NB: Max depth is the deepest, we checked

# Download the chosen layers
  # NB: Don't run this if nothing has changed as there is no need to ping their servers
                              # Bottom temperature
BO_layers_dl <- load_layers(c("BO2_templtmin_bdmax", "BO2_tempmean_bdmax", "BO2_templtmax_bdmax", 
                              # Surface temperature
                              "BO2_templtmin_ss", "BO2_tempmean_ss", "BO2_templtmax_ss", 
                              # Bottom salinity
                              "BO2_salinityltmin_bdmax", "BO2_salinitymean_bdmax", "BO2_salinityltmax_bdmax", 
                              # Surface salinity
                              "BO2_salinityltmin_ss", "BO2_salinitymean_ss", "BO2_salinityltmax_ss", 
                              # Ice thickness
                              "BO2_icethickltmin_ss", "BO2_icethickmean_ss", "BO2_icethickltmax_ss", 
                              # Current velocity
                                # Note that the v2.1 files are needed and downloaded here: https://www.bio-oracle.org/downloads-to-email.php
                              # "BO2_curvelltmin_bdmax", "BO2_curvelmean_bdmax", "BO2_curvelltmax_bdmax",
                              # Photosynthetically active radiation
                              "BO_parmean", "BO_parmax", 
                              # Dissolve oxygen
                              "BO2_dissoxltmin_bdmax", "BO2_dissoxmean_bdmax", "BO2_dissoxltmax_bdmax", 
                              # Iron
                              "BO2_ironltmin_bdmax", "BO2_ironmean_bdmax", "BO2_ironltmax_bdmax", 
                              # Nitrate
                              "BO2_nitrateltmin_bdmax", "BO2_nitratemean_bdmax", "BO2_nitrateltmax_bdmax", 
                              # Phosphate
                              "BO2_phosphateltmin_bdmax", "BO2_phosphatemean_bdmax", "BO2_phosphateltmax_bdmax"))

# Convert to dataframe
BO_layers_df <- as.data.frame(BO_layers_dl, xy = T) %>% 
  dplyr::rename(lon = x, lat = y) %>% 
  filter(BO2_icethickmean_ss >= 0,
         lat >= min(Arctic_boundary$lat)) #%>%
  # mutate(lon = round(lon, 4), 
         # lat = round(lat, 4))
rm(BO_layers_dl); gc()

# Join the current velocity v2.1 layers
BO_layers_present <- BO_layers_df %>% 
  left_join(load_asc("data/Present.Benthic.Max.Depth.Current.Velocity.Lt.min.asc.BOv2_1.asc", "BO21_curvelltmin_bdmax")) %>% 
  left_join(load_asc("data/Present.Benthic.Max.Depth.Current.Velocity.Mean.asc.BOv2_1.asc", "BO21_curvelmean_bdmax")) %>% 
  left_join(load_asc("data/Present.Benthic.Max.Depth.Current.Velocity.Lt.max.asc.BOv2_1.asc", "BO21_curvelltmax_bdmax"))
rm(BO_layers_df); gc()

# Clip to Arctic study region
Arctic_BO <- BO_layers_present %>%
  mutate(in_grid = sp::point.in.polygon(point.x = BO_layers_present[["lon"]], point.y = BO_layers_present[["lat"]], 
                                        pol.x = Arctic_boundary[["lon"]], pol.y = Arctic_boundary[["lat"]])) %>% 
  filter(in_grid >= 1) %>% 
  dplyr::select(-in_grid)
save(Arctic_BO, file = "data/Arctic_BO.RData")

# Visualise
ggplot(Arctic_BO, aes(x = lon, y = lat)) +
  # geom_raster(aes(fill = BO2_templtmax_bdmax)) +
  geom_raster(aes(fill = BO21_curvelltmax_bdmax)) +
  borders(fill = "grey70", colour = "black") +
  scale_fill_viridis_c(option = "D") +
  coord_quickmap(xlim = c(bbox_arctic[1], bbox_arctic[2]),
                 ylim = c(bbox_arctic[3], bbox_arctic[4]),
                 expand = F) +
  theme(legend.position = "bottom")


# Load future Bio-ORACLE data ---------------------------------------------

## NB: All of the v2.1 future layers had issues so had to be downloaded manually from the BO website:
# https://www.bio-oracle.org/downloads-to-email.php

# Future scenario conversions
# https://ar5-syr.ipcc.ch/topic_futurechanges.php
# RCP8.5 ~= A2
# RCP6.0 ~= B2
# RCP4.5 ~= B1

# Look at possible layers
BO_layers_future <- list_layers_future(datasets = "Bio-ORACLE") %>% 
  filter(scenario == "RCP85")

# Load present data for easier joining
load("data/Arctic_BO.RData")

# Load and process 2050 data layers
  # Note that the layer names must be the same as the present data even though these are v2.1 data
Arctic_BO_2050 <- Arctic_BO %>% 
  dplyr::select(lon, lat, BO_parmean:BO2_phosphateltmax_bdmax) %>% 
  # Add bottom temperature
  left_join(load_asc("data/2050AOGCM.RCP85.Benthic.Max.Depth.Temperature.Lt.min.asc.BOv2_1.asc", "BO2_templtmin_bdmax")) %>%
  left_join(load_asc("data/2050AOGCM.RCP85.Benthic.Max.Depth.Temperature.Mean.asc.BOv2_1.asc", "BO2_tempmean_bdmax")) %>% 
  left_join(load_asc("data/2050AOGCM.RCP85.Benthic.Max.Depth.Temperature.Lt.max.asc.BOv2_1.asc", "BO2_templtmax_bdmax")) %>% 
  # Add SST
  left_join(load_asc("data/2050AOGCM.RCP85.Surface.Temperature.Lt.min.asc.BOv2_1.asc", "BO2_templtmin_ss")) %>%
  left_join(load_asc("data/2050AOGCM.RCP85.Surface.Temperature.Mean.asc.BOv2_1.asc", "BO2_tempmean_ss")) %>% 
  left_join(load_asc("data/2050AOGCM.RCP85.Surface.Temperature.Lt.max.asc.BOv2_1.asc", "BO2_templtmax_ss")) %>% 
  # Add bottom salinity
  left_join(load_asc("data/2050AOGCM.RCP85.Benthic.Max.Depth.Salinity.Lt.min.asc.BOv2_1.asc", "BO2_salinityltmin_bdmax")) %>%
  left_join(load_asc("data/2050AOGCM.RCP85.Benthic.Max.Depth.Salinity.Mean.asc.BOv2_1.asc", "BO2_salinitymean_bdmax")) %>%
  left_join(load_asc("data/2050AOGCM.RCP85.Benthic.Max.Depth.Salinity.Lt.max.asc.BOv2_1.asc", "BO2_salinityltmax_bdmax")) %>%
  # Add SSS
  left_join(load_asc("data/2050AOGCM.RCP85.Surface.Salinity.Lt.min.asc.BOv2_1.asc", "BO2_salinityltmin_ss")) %>%
  left_join(load_asc("data/2050AOGCM.RCP85.Surface.Salinity.Mean.asc.BOv2_1.asc", "BO2_salinitymean_ss")) %>%
  left_join(load_asc("data/2050AOGCM.RCP85.Surface.Salinity.Lt.max.asc.BOv2_1.asc", "BO2_salinityltmax_ss")) %>%
  # Add ice thickness
  left_join(load_asc("data/2050AOGCM.RCP85.Surface.Ice.thickness.Lt.min.asc.BOv2_1.asc", "BO2_icethickltmin_ss")) %>%
  left_join(load_asc("data/2050AOGCM.RCP85.Surface.Ice.thickness.Mean.asc.BOv2_1.asc", "BO2_icethickmean_ss")) %>%
  left_join(load_asc("data/2050AOGCM.RCP85.Surface.Ice.thickness.Lt.max.asc.BOv2_1.asc", "BO2_icethickltmax_ss")) %>%
  # Add current velocities
  left_join(load_asc("data/2050AOGCM.RCP85.Benthic.Max.Depth.Current.Velocity.Lt.min.asc.BOv2_1.asc", "BO21_curvelltmin_bdmax")) %>%
  left_join(load_asc("data/2050AOGCM.RCP85.Benthic.Max.Depth.Current.Velocity.Mean.asc.BOv2_1.asc", "BO21_curvelmean_bdmax")) %>%
  left_join(load_asc("data/2050AOGCM.RCP85.Benthic.Max.Depth.Current.Velocity.Lt.max.asc.BOv2_1.asc", "BO21_curvelltmax_bdmax"))
save(Arctic_BO_2050, file = "data/Arctic_BO_2050.RData")

# Load and process 2100 data layers
# Note that the layer names must be the same as the present data even though these are v2.1 data
Arctic_BO_2100 <- Arctic_BO %>% 
  dplyr::select(lon, lat, BO_parmean:BO2_phosphateltmax_bdmax) %>% 
  # Add bottom temperature
  left_join(load_asc("data/2100AOGCM.RCP85.Benthic.Max.Depth.Temperature.Lt.min.asc.BOv2_1.asc", "BO2_templtmin_bdmax")) %>%
  left_join(load_asc("data/2100AOGCM.RCP85.Benthic.Max.Depth.Temperature.Mean.asc.BOv2_1.asc", "BO2_tempmean_bdmax")) %>% 
  left_join(load_asc("data/2100AOGCM.RCP85.Benthic.Max.Depth.Temperature.Lt.max.asc.BOv2_1.asc", "BO2_templtmax_bdmax")) %>% 
  # Add SST
  left_join(load_asc("data/2100AOGCM.RCP85.Surface.Temperature.Lt.min.asc.BOv2_1.asc", "BO2_templtmin_ss")) %>%
  left_join(load_asc("data/2100AOGCM.RCP85.Surface.Temperature.Mean.asc.BOv2_1.asc", "BO2_tempmean_ss")) %>% 
  left_join(load_asc("data/2100AOGCM.RCP85.Surface.Temperature.Lt.max.asc.BOv2_1.asc", "BO2_templtmax_ss")) %>% 
  # Add bottom salinity
  left_join(load_asc("data/2100AOGCM.RCP85.Benthic.Max.Depth.Salinity.Lt.min.asc.BOv2_1.asc", "BO2_salinityltmin_bdmax")) %>%
  left_join(load_asc("data/2100AOGCM.RCP85.Benthic.Max.Depth.Salinity.Mean.asc.BOv2_1.asc", "BO2_salinitymean_bdmax")) %>%
  left_join(load_asc("data/2100AOGCM.RCP85.Benthic.Max.Depth.Salinity.Lt.max.asc.BOv2_1.asc", "BO2_salinityltmax_bdmax")) %>%
  # Add SSS
  left_join(load_asc("data/2100AOGCM.RCP85.Surface.Salinity.Lt.min.asc.BOv2_1.asc", "BO2_salinityltmin_ss")) %>%
  left_join(load_asc("data/2100AOGCM.RCP85.Surface.Salinity.Mean.asc.BOv2_1.asc", "BO2_salinitymean_ss")) %>%
  left_join(load_asc("data/2100AOGCM.RCP85.Surface.Salinity.Lt.max.asc.BOv2_1.asc", "BO2_salinityltmax_ss")) %>%
  # Add ice thickness
  left_join(load_asc("data/2100AOGCM.RCP85.Surface.Ice.thickness.Lt.min.asc.BOv2_1.asc", "BO2_icethickltmin_ss")) %>%
  left_join(load_asc("data/2100AOGCM.RCP85.Surface.Ice.thickness.Mean.asc.BOv2_1.asc", "BO2_icethickmean_ss")) %>%
  left_join(load_asc("data/2100AOGCM.RCP85.Surface.Ice.thickness.Lt.max.asc.BOv2_1.asc", "BO2_icethickltmax_ss")) %>%
  # Add current velocities
  left_join(load_asc("data/2100AOGCM.RCP85.Benthic.Max.Depth.Current.Velocity.Lt.min.asc.BOv2_1.asc", "BO21_curvelltmin_bdmax")) %>%
  left_join(load_asc("data/2100AOGCM.RCP85.Benthic.Max.Depth.Current.Velocity.Mean.asc.BOv2_1.asc", "BO21_curvelmean_bdmax")) %>%
  left_join(load_asc("data/2100AOGCM.RCP85.Benthic.Max.Depth.Current.Velocity.Lt.max.asc.BOv2_1.asc", "BO21_curvelltmax_bdmax"))
save(Arctic_BO_2100, file = "data/Arctic_BO_2100.RData")

# Visualise
ggplot(Arctic_BO_2100, aes(x = lon, y = lat)) +
  # geom_raster(aes(fill = BO_parmax)) +
  # geom_raster(aes(fill = BO2_templtmax_bdmax)) +
  geom_raster(aes(fill = BO21_curvelltmax_bdmax)) +
  borders(fill = "grey70", colour = "black") +
  scale_fill_viridis_c(option = "D") +
  coord_cartesian(xlim = c(bbox_arctic[1], bbox_arctic[2]),
                  ylim = c(bbox_arctic[3], bbox_arctic[4]),
                  expand = F) +
  theme(legend.position = "bottom")


# Load GMED ASCII files ---------------------------------------------------

# The .asc files loaded below were sent by Jesi
# They are the Aqua Maps layers but have been 
# regridded to match the BioOracle data

# Depth
depth <- read.asciigrid("data/depthclip.asc")
depth_df <- as.data.frame(depth, xy = T)
Arctic_depth <- depth_df %>%
  dplyr::rename(bathy = data.depthclip.asc,
                lon = s1, lat = s2) %>%
  filter(lon >= bbox_arctic[1], lon <= bbox_arctic[2],
         lat >= bbox_arctic[3], lat <= bbox_arctic[4])

# Distance to land
land_distance <- read.asciigrid("data/landdistclip.asc")
land_distance_df <- as.data.frame(land_distance, xy = T)
Arctic_land_distance <- land_distance_df %>%
  dplyr::rename(land_distance = data.landdistclip.asc,
                lon = s1, lat = s2) %>%
  filter(lon >= bbox_arctic[1], lon <= bbox_arctic[2],
         lat >= bbox_arctic[3], lat <= bbox_arctic[4])

# Combine into single file and save
Arctic_AM <- left_join(Arctic_land_distance, Arctic_depth, by = c("lon", "lat")) %>% 
  dplyr::select(lon, lat, everything()) %>% 
  mutate(lon = round(lon, 5),
         lat = round(lat, 5))
save(Arctic_AM, file = "data/Arctic_AM.RData")

# Visualise
ggplot(Arctic_AM, aes(x = lon, y = lat)) +
  geom_tile(aes(fill = land_distance)) +
  borders(fill = "grey70", colour = "black") +
  scale_fill_viridis_c(option = "D") +
  coord_cartesian(xlim = c(bbox_arctic[1], bbox_arctic[2]),
                  ylim = c(bbox_arctic[3], bbox_arctic[4]),
                  expand = F) +
  theme(legend.position = "bottom")


# Establish the correlation matrix ----------------------------------------

# Check Pearson correlation coefficient between layers
  # These warnings are fine, we don't 
BO_cor_matrix <- Arctic_BO %>% 
  dplyr::select(-lon, -lat) %>% 
  correlation(redundant = T) %>% 
  dplyr::select(Parameter1:r) %>% 
  pivot_wider(names_from = Parameter2, values_from = r)
save(BO_cor_matrix, file = "data/BO_cor_matrix.RData")
write_csv(BO_cor_matrix, "data/BO_cor_matrix.csv")

BO_cor_groups <- correlation_groups(layers_correlation(colnames(dplyr::select(Arctic_env, BO2_templtmin_bdmax:BO2_phosphateltmax_bdmax))))
BO_cor_groups

# ANother method of screening vartiables based on correlations
usdm::vifcor()

