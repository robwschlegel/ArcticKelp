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
BO_layers_dl <- load_layers(c("BO21_templtmin_bdmax", "BO21_tempmean_bdmax", "BO21_templtmax_bdmax", 
                              # Surface temperature
                              "BO21_templtmin_ss", "BO21_tempmean_ss", "BO21_templtmax_ss", 
                              # Bottom salinity
                              "BO21_salinityltmin_bdmax", "BO21_salinitymean_bdmax", "BO21_salinityltmax_bdmax", 
                              # Surface salinity
                              "BO21_salinityltmin_ss", "BO21_salinitymean_ss", "BO21_salinityltmax_ss", 
                              # Ice thickness
                              "BO21_icethickltmin_ss", "BO21_icethickmean_ss", "BO21_icethickltmax_ss", 
                              # Current velocity
                              "BO21_curvelltmin_bdmax", "BO21_curvelmean_bdmax", "BO21_curvelltmax_bdmax",
                              # Photosynthetically active radiation
                              "BO_parmean", "BO_parmax", 
                              # Dissolve oxygen
                              "BO21_dissoxltmin_bdmax", "BO21_dissoxmean_bdmax", "BO21_dissoxltmax_bdmax", 
                              # Iron
                              "BO21_ironltmin_bdmax", "BO21_ironmean_bdmax", "BO21_ironltmax_bdmax", 
                              # Nitrate
                              "BO21_nitrateltmin_bdmax", "BO21_nitratemean_bdmax", "BO21_nitrateltmax_bdmax", 
                              # Phosphate
                              "BO21_phosphateltmin_bdmax", "BO21_phosphatemean_bdmax", "BO21_phosphateltmax_bdmax"))

# Convert to dataframe
BO_layers_present <- as.data.frame(BO_layers_dl, xy = T) %>% 
  dplyr::rename(lon = x, lat = y) %>% 
  filter(BO21_icethickmean_ss >= 0,
         lat >= min(Arctic_boundary$lat)) #%>%
# mutate(lon = round(lon, 5), 
# lat = round(lat, 5))
rm(BO_layers_dl); gc()

# Clip to Arctic study region
Arctic_BO <- BO_layers_present %>%
  mutate(in_grid = sp::point.in.polygon(point.x = BO_layers_present[["lon"]], point.y = BO_layers_present[["lat"]], 
                                        pol.x = Arctic_boundary[["lon"]], pol.y = Arctic_boundary[["lat"]])) %>% 
  filter(in_grid >= 1) %>% 
  dplyr::select(-in_grid)
save(Arctic_BO, file = "data/Arctic_BO.RData")
rm(BO_layers_present); gc()

# Visualise
ggplot(Arctic_BO, aes(x = lon, y = lat)) +
  geom_raster(aes(fill = BO_parmean)) +
  borders(fill = "grey70", colour = "black") +
  scale_fill_viridis_c(option = "D") +
  # coord_quickmap(xlim = c(bbox_arctic[1], bbox_arctic[2]),
                 # ylim = c(bbox_arctic[3], bbox_arctic[4]),
                 # expand = F) +
  theme(legend.position = "bottom")


# Load future Bio-ORACLE data ---------------------------------------------

# Future scenario conversions
# https://ar5-syr.ipcc.ch/topic_futurechanges.php
# RCP8.5 ~= A2
# RCP6.0 ~= B2
# RCP4.5 ~= B1

# Load present data for easier joining
load("data/Arctic_BO.RData")

# Look at possible layers
BO_layers_future <- list_layers_future(datasets = "Bio-ORACLE") %>% 
  filter(scenario == "RCP85")

# Download as similar of layers as possible to present data
# Bottom temperature
BO_layers_future_dl <- get_future_layers(c("BO21_templtmin_bdmax", "BO21_tempmean_bdmax", "BO21_templtmax_bdmax", 
                                           # Surface temperature
                                           "BO21_templtmin_ss", "BO21_tempmean_ss", "BO21_templtmax_ss", 
                                           # Bottom salinity
                                           "BO21_salinityltmin_bdmax", "BO21_salinitymean_bdmax", "BO21_salinityltmax_bdmax", 
                                           # Surface salinity
                                           "BO21_salinityltmin_ss", "BO21_salinitymean_ss", "BO21_salinityltmax_ss", 
                                           # Ice thickness
                                           "BO21_icethickltmin_ss", "BO21_icethickmean_ss", "BO21_icethickltmax_ss",
                                           # Current velocity
                                           "BO21_curvelltmin_bdmax", "BO21_curvelmean_bdmax", "BO21_curvelltmax_bdmax"), 
                                         scenario = "RCP85", year = c(2050, 2100))
BO_layers_future_dl <- load_layers(BO_layers_future_dl$layer_code)

# Convert to data.frame
BO_layers_future_df <- as.data.frame(BO_layers_future_dl, xy = T) %>% 
  dplyr::rename(lon = x, lat = y) %>% 
  filter(BO21_RCP85_2100_icethickmean_ss >= 0,
         lat >= min(Arctic_boundary$lat))
rm(BO_layers_future_dl); gc()

# Create 2050 data.frame
Arctic_BO_2050 <- Arctic_BO %>% 
  dplyr::select(lon, lat, BO_parmean:BO21_phosphateltmax_bdmax) %>% 
  left_join(dplyr::select(BO_layers_future_df, lon, lat, BO21_RCP85_2050_curvelltmax_bdmax:BO21_RCP85_2050_tempmean_ss), 
            by = c("lon", "lat"))
colnames(Arctic_BO_2050) <- sub("RCP85_2050_", "", colnames(Arctic_BO_2050))
save(Arctic_BO_2050, file = "data/Arctic_BO_2050.RData")

# Create 2100 data.frame
Arctic_BO_2100 <- Arctic_BO %>% 
  dplyr::select(lon, lat, BO_parmean:BO21_phosphateltmax_bdmax) %>% 
  left_join(dplyr::select(BO_layers_future_df, lon, lat, BO21_RCP85_2100_curvelltmax_bdmax:BO21_RCP85_2100_tempmean_ss), 
            by = c("lon", "lat"))
colnames(Arctic_BO_2100) <- sub("RCP85_2100_", "", colnames(Arctic_BO_2100))
save(Arctic_BO_2100, file = "data/Arctic_BO_2100.RData")
rm(BO_layers_future_df); gc()

# Visualise
ggplot(Arctic_BO_2100, aes(x = lon, y = lat)) +
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
                lon = s1, lat = s2) #%>%
# filter(lon >= bbox_arctic[1], lon <= bbox_arctic[2],
# lat >= bbox_arctic[3], lat <= bbox_arctic[4])

# Distance to land
land_distance <- read.asciigrid("data/landdistclip.asc")
land_distance_df <- as.data.frame(land_distance, xy = T)
Arctic_land_distance <- land_distance_df %>%
  dplyr::rename(land_distance = data.landdistclip.asc,
                lon = s1, lat = s2) #%>%
# filter(lon >= bbox_arctic[1], lon <= bbox_arctic[2],
# lat >= bbox_arctic[3], lat <= bbox_arctic[4])

# Combine into single file and save
Arctic_AM <- left_join(Arctic_land_distance, Arctic_depth, by = c("lon", "lat")) %>% 
  dplyr::select(lon, lat, everything()) %>% 
  mutate(in_grid = sp::point.in.polygon(point.x = Arctic_land_distance[["lon"]], point.y = Arctic_land_distance[["lat"]], 
                                        pol.x = Arctic_boundary[["lon"]], pol.y = Arctic_boundary[["lat"]])) %>% 
  filter(in_grid >= 1) %>% 
  dplyr::select(-in_grid) #%>% 
  # mutate(lon = round(lon, 5),
         # lat = round(lat, 5))
save(Arctic_AM, file = "data/Arctic_AM.RData")

# Visualise
ggplot(Arctic_AM, aes(x = lon, y = lat)) +
  geom_tile(aes(fill = land_distance)) +
  borders(fill = "grey70", colour = "black") +
  scale_fill_viridis_c(option = "D") +
  # coord_cartesian(xlim = c(bbox_arctic[1], bbox_arctic[2]),
  # ylim = c(bbox_arctic[3], bbox_arctic[4]),
  # expand = F) +
  theme(legend.position = "bottom")


# Create coastal data.frame -----------------------------------------------

# The present data
load("data/Arctic_BO.RData")
Arctic_BO <- Arctic_BO %>%
  mutate(lon_match = round(lon, 5),
         lat_match = round(lat, 5))

# The depth/distance info
load("data/Arctic_AM.RData")
Arctic_AM <- Arctic_AM %>% 
  dplyr::rename(depth = bathy) %>% 
  mutate(lon = round(lon, 5),
         lat = round(lat, 5))

# Merge to extract only "coastal" pixels
Arctic_coast <- left_join(Arctic_BO , Arctic_AM, 
                          by = c("lon_match" = "lon", "lat_match" = "lat")) %>% 
  filter(land_distance <= 50 | depth <= 100)
save(Arctic_coast, file = "data/Arctic_coast.RData")

# Plot coastal pixels
ggplot(data = Arctic_coast, aes(x = lon, y = lat)) +
# ggplot(data = filter(Arctic_AM, depth > 100), aes(x = lon, y = lat)) +
  geom_tile(aes(fill = BO21_templtmin_bdmax)) +
  # geom_tile(aes(fill = depth)) +
  borders(fill = "grey90", colour = "black") +
  coord_quickmap(ylim = c(50, 90), expand = F)

# Coastal coordinates
coastal_coords <- dplyr::select(Arctic_coast, lon, lat) %>%
  mutate(env_index = 1:n())
save(coastal_coords, file = "metadata/coastal_coords.RData")


# Establish the correlation matrix ----------------------------------------

# Check Pearson correlation coefficient between layers
BO_cor_matrix <- Arctic_BO %>% 
  dplyr::select(-lon, -lat) %>% 
  correlation(redundant = T) %>% 
  dplyr::select(Parameter1:r) %>% 
  pivot_wider(names_from = Parameter2, values_from = r)
save(BO_cor_matrix, file = "data/BO_cor_matrix.RData")
write_csv(BO_cor_matrix, "data/BO_cor_matrix.csv")

# Another method
# BO_cor_groups <- correlation_groups(layers_correlation(colnames(dplyr::select(Arctic_env, BO21_templtmin_bdmax:BO21_phosphateltmax_bdmax))))
# BO_cor_groups

# Yet another method of screening variables based on correlations
# usdm::vifcor()

