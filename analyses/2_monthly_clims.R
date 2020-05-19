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

# Download bathy layers
  # NB: Not used as these values are suspicious
# bathy <- load_layers(c("BO_bathymin", "BO_bathymean", "BO_bathymax"))
# Arctic_bathy <- as.data.frame(bathy, xy = T) %>% 
#   dplyr::rename(lon = x, lat = y) %>% 
#   filter(lon >= bbox_arctic[1], lon <= bbox_arctic[2],
#          lat >= bbox_arctic[3], lat <= bbox_arctic[4])
# save(Arctic_bathy, file = "data/Arctic_bathy.RData")

# The layes currently chosen for use in this study
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
                              # Diffuse attenuation coefficient at 490 nm 
                              "BO_damin", "BO_damean", "BO_damax", 
                              # Primary productivity
                              # "BO2_ppltmin_ss", "BO2_ppmean_ss", "BO2_ppltmax_ss", 
                              # Calcite
                              # "BO_calcite",
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

# NB: Add the diffuse attenuation value e-mailed from Jorge

# Convert to dataframe
BO_layers_df <- as.data.frame(BO_layers_dl, xy = T)

# Clip to Arctic study region
Arctic_BO <- BO_layers_df %>%
  dplyr::rename(lon = x, lat = y) %>%
  filter(lon >= bbox_arctic[1], lon <= bbox_arctic[2],
         lat >= bbox_arctic[3], lat <= bbox_arctic[4],
         BO2_icecovermean_ss >= 0) %>% 
  mutate(lon = round(lon, 5),
         lat = round(lat, 5))
save(Arctic_BO, file = "data/Arctic_BO.RData")

# Visualise
ggplot(Arctic_BO, aes(x = lon, y = lat)) +
  geom_tile(aes(fill = BO_damax)) +
  borders(fill = "grey70", colour = "black") +
  scale_fill_viridis_c(option = "D") +
  coord_cartesian(xlim = c(bbox_arctic[1], bbox_arctic[2]),
                  ylim = c(bbox_arctic[3], bbox_arctic[4]),
                  expand = F) +
  theme(legend.position = "bottom")


# Download future Bio-ORACLE data -----------------------------------------

# Future scenario conversions
# https://ar5-syr.ipcc.ch/topic_futurechanges.php
# RCP8.5 ~= A2
# RCP6.0 ~= B2
# RCP4.5 ~= B1

# Look at possible layers
BO_layers_future <- list_layers_future(datasets = "Bio-ORACLE") %>% 
  filter(scenario == "RCP85")

# Download as similar of layers as possible to present data
                                           # Bottom temperature
BO_layers_future_dl <- get_future_layers(c("BO2_templtmin_bdmax", "BO2_tempmean_bdmax", "BO2_templtmax_bdmax", 
                                           # Surface temperature
                                           "BO2_templtmin_ss", "BO2_tempmean_ss", "BO2_templtmax_ss", 
                                           # Bottom salinity
                                           "BO2_salinityltmin_bdmax", "BO2_salinitymean_bdmax", "BO2_salinityltmax_bdmax", 
                                           # Surface salinity
                                           "BO2_salinityltmin_ss", "BO2_salinitymean_ss", "BO2_salinityltmax_ss", 
                                           # Ice thickness
                                           "BO2_icethickltmin_ss", "BO2_icethickmean_ss", "BO2_icethickltmax_ss", 
                                           # Current velocity
                                           "BO2_curvelltmin_bdmax", "BO2_curvelmean_bdmax", "BO2_curvelltmax_bdmax"), 
                                         scenario = "RCP85", year = c(2050, 2100))
BO_layers_future_dl <- load_layers(BO_layers_future_dl$layer_code)

# Convert to dataframe
BO_layers_future_df <- as.data.frame(BO_layers_future_dl, xy = T)

# Clip to Arctic study region
Arctic_BO_future <- BO_layers_future_df %>%
  dplyr::rename(lon = x, lat = y) %>%
  filter(lon >= bbox_arctic[1], lon <= bbox_arctic[2],
         lat >= bbox_arctic[3], lat <= bbox_arctic[4],
         BO2_RCP85_2100_icethickmean_ss >= 0) %>% 
  mutate(lon = round(lon, 5),
         lat = round(lat, 5))
save(Arctic_BO_future, file = "data/Arctic_BO_future.RData")

# Visualise
ggplot(Arctic_BO_future, aes(x = lon, y = lat)) +
  geom_tile(aes(fill = BO2_RCP85_2100_tempmean_ss)) +
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


# Combine BO and GMED -----------------------------------------------------

# Load data from previous steps as necessary
load("data/Arctic_BO.RData")
load("data/Arctic_BO_future.RData")
load("data/Arctic_AM.RData")

# Merge and save
Arctic_env <- left_join(Arctic_BO, Arctic_AM, by = c("lon", "lat")) %>% 
  left_join(Arctic_BO_future, by = c("lon", "lat")) %>% 
  na.omit() # There are some pixels from the BO data that dont' match with BO2
save(Arctic_env, file = "data/Arctic_env.RData")
rm(Arctic_AM, Arctic_BO, Arctic_BO_future); gc()

# Visualise
ggplot(Arctic_env, aes(x = lon, y = lat)) +
  geom_tile(aes(fill = BO2_tempmean_ss)) +  
  borders(fill = "grey70", colour = "black") +
  scale_fill_viridis_c(option = "E") +
  coord_cartesian(xlim = c(bbox_arctic[1], bbox_arctic[2]),
                  ylim = c(bbox_arctic[3], bbox_arctic[4]),
                  expand = F) +
  theme(legend.position = "bottom")


# Establish the correlation matrix ----------------------------------------

# Check Pearson correlation coefficient between layers
  # These warnings are fine, we don't 
BO_cor_matrix <- layers_correlation(colnames(dplyr::select(Arctic_env, 
                                                           BO2_templtmin_bdmax:BO2_curvelltmax_bdmax))) %>% 
  mutate(var = row.names(.)) %>% 
  dplyr::select(var, everything())
save(BO_cor_matrix, file = "data/BO_cor_matrix.RData")
