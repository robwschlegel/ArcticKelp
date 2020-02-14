# analyses/2_monthly_clims.R
# The purpose of this script is to create monthly climatologies at each pixel
# in the arctic from the NAPA model and Bio-Oracle data
# This should only need to be run once at the outset of the project.
# If any sites change/are added, run "analyses/3_study_sites_clims.R" 
# to get the updated values, not this script

# NB: This script will only run on Eric Oliver's tikoraluk server
# This is because it is utilising files that are up to 25GB in size
# This shouldn't be an issue though as it only needs to be run once

# The first round of data compiled come from the NAPA model created at BIO in Halifax
  # These data are only available through a specific server and so much of the code used below
  # to access them won't work outside of that environment.
  # For this reason much of the NAPA code is commented out
# The second round of data are from Bio-Oracle (http://www.bio-oracle.org/code.php)
  # They are downloaded directly from their server and so may be accessed from anywhere

# Because nearly all of this code is only meant to be run once,
# much of it has been commented out to prevent it accidentally being re-run
# This wouldn't create erors etc., it would just be time consuming


# Setup -------------------------------------------------------------------

# The study sites and bounding box
source("analyses/1_study_sites.R")

# Bio-Oracle access
library(sdmpredictors)

# Other packages
library(tidync)
library(stringr)
library(data.table)

# Set cores
doParallel::registerDoParallel(cores = 50)


# NAPA NetCDF information -------------------------------------------------

# NB: The code in this section is greyed out as it will only run correctly on the server where the model data live
  # It also only needs to be run once, so is also just stored here for future convenience

# Check the grid data
# info_grid <- ncdump::NetCDF("../../data/NAPA025/mesh_grid/bathy_creg025_extended_5m.nc")$variable %>%
#   dplyr::select(-addOffset, -scaleFact)

# Check the surface NetCDF information
# info_surface <- ncdump::NetCDF("../../data/NAPA025/1d_grid_T_2D/CREG025E-GLSII_1d_grid_T_2D_19931001-19931005.nc")$variable

# Check the ice NetCDF information
# info_ice <- ncdump::NetCDF("../../data/NAPA025/5d_icemod/CREG025E-GLSII_5d_icemod_19931001-19931005.nc")$variable#$longname

# Check the depth NetCDF information
# info_depth_T <- ncdump::NetCDF("../../data/NAPA025/5d_grid_T/CREG025E-GLSII_5d_grid_T_19931001-19931005.nc")$variable#$longname
# info_depth_U <- ncdump::NetCDF("../../data/NAPA025/5d_grid_U/CREG025E-GLSII_5d_grid_U_19931001-19931005.nc")$variable#$longname
# info_depth_V <- ncdump::NetCDF("../../data/NAPA025/5d_grid_V/CREG025E-GLSII_5d_grid_V_19931001-19931005.nc")$variable#$longname
# info_depth_W <- ncdump::NetCDF("../../data/NAPA025/5d_grid_W/CREG025E-GLSII_5d_grid_W_19931001-19931005.nc")$variable#$longname

# Combine and clean-up
# info_ALL <- rbind(info_grid, info_surface, info_ice, info_depth_T, info_depth_U, info_depth_V, info_depth_W) %>%
  # dplyr::select(name, ndims, units, longname) %>%
  # unique() %>%
  # arrange(name)

# Save
# write_csv(info_ALL, "metadata/model_info.csv")
model_info <- read_csv("metadata/model_info.csv")


# NAPA data ---------------------------------------------------------------

# Extract NAPA grid from NetCDF and save as an RData file
# NAPA_bathy <- tidync("../../data/NAPA025/mesh_grid/bathy_creg025_extended_5m.nc") %>%
#   hyper_tibble() %>%
#   select(-time_counter) %>%
#   dplyr::rename(bathy = Bathymetry)
# NAPA_lon <- tidync("../../data/NAPA025/mesh_grid/bathy_creg025_extended_5m.nc") %>%
#   activate(nav_lon) %>%
#   hyper_tibble()
# NAPA_lat <- tidync("../../data/NAPA025/mesh_grid/bathy_creg025_extended_5m.nc") %>%
#   activate(nav_lat) %>%
#   hyper_tibble()
# NAPA_grid <- left_join(NAPA_lon, NAPA_lat, by = c("x", "y")) %>%
#   left_join(NAPA_bathy, by = c("x", "y")) %>%
#   select(x, y, nav_lon, nav_lat, bathy)
# NAPA_arctic <- NAPA_grid %>%
#   filter(nav_lon >= bbox_arctic[1], nav_lon <= bbox_arctic[2],
#          nav_lat >= bbox_arctic[3], nav_lat <= bbox_arctic[4],
#          bathy > 0) %>%
#   mutate(nav_lon = round(nav_lon, 4),
#          nav_lat = round(nav_lat, 4),
#          bathy = round(bathy))
# save(NAPA_arctic, file = "metadata/NAPA_arctic.RData")
load("metadata/NAPA_arctic.RData")

# NB: Should not use data before 1998-01-01 as this is model spin-up

# The surface files location
NAPA_surface_files <- dir("../../data/NAPA025/1d_grid_T_2D", full.names = T)[-c(1:311)]
# head(NAPA_surface_files)

# The ice file location
NAPA_ice_files <- dir("../../data/NAPA025/5d_icemod", full.names = T)[-c(1:311)]
# head(NAPA_ice_files)

# The depth files location
NAPA_depth_T_files <- dir("../../data/NAPA025/5d_grid_T", full.names = T)[-c(1:311)]
NAPA_depth_U_files <- dir("../../data/NAPA025/5d_grid_U", full.names = T)[-c(1:311)]
NAPA_depth_V_files <- dir("../../data/NAPA025/5d_grid_V", full.names = T)[-c(1:311)]
NAPA_depth_W_files <- dir("../../data/NAPA025/5d_grid_W", full.names = T)[-c(1:311)]
# head(NAPA_depth_T_files)

# A custom rainbow palette
# A rainbow colour palette was explicitly requested
rainbow_palette <- c("#fefefe", "#f963fa", "#020135", "#00efe1", "#057400", "#fcfd00", "#ed0000", "#3d0000")


# Loading functions -------------------------------------------------------

# Function for loading the individual NAPA NetCDF surface files
# testers...
# file_name <- NAPA_surface_files[1]
# file_name <- NAPA_ice_files[1]
# layer_name <- "D2,D1,D0"
load_NAPA_surface <- function(file_name, layer_name = "D2,D1,D0"){
  res <- tidync(file_name, what = layer_name) %>%
    # Quick filter before loading into memory
    hyper_filter(x = dplyr::between(x, min(NAPA_arctic$x), max(NAPA_arctic$x)),
                 y = dplyr::between(y, min(NAPA_arctic$y), max(NAPA_arctic$y))) %>%
    hyper_tibble() %>%
    # More precise filter to the study area before further processing
    # This is necessary anyway to attach lon/lat/bathy info
    right_join(NAPA_arctic, by = c("x", "y")) %>% 
    dplyr::rename(t = time_counter) %>%
    mutate(t = as.Date(as.POSIXct(t, origin = "1900-01-01")),
           month = lubridate::month(t, label = T),
           bathy = round(bathy)) %>% 
    dplyr::select(x, y, nav_lon, nav_lat, t, month, everything())
  if("ncatice" %in% colnames(res)){
    res <- res %>% 
      dplyr::select(-ncatice) %>%
      unique() %>% 
      mutate(iceconc_cat = round(iceconc_cat, 2),
             icethic_cat = round(icethic_cat, 2)) %>% 
      unique()
  }
  return(res)
}
# system.time(test <- load_NAPA_surface(NAPA_surface_files[111], "D2,D1,D0")) # 1 second
# system.time(test <- load_NAPA_surface(NAPA_ice_files[111], "D2,D1,D3,D0")) # 2 second
# system.time(test <- load_NAPA_surface(NAPA_ice_files[111])) # 2 second

# Function for loading the individual NAPA NetCDF depth files
# testers...
# file_name <- NAPA_depth_T_files[1]
# file_name <- NAPA_depth_U_files[1]
# file_name <- NAPA_depth_V_files[1]
# file_name <- NAPA_depth_W_files[1]
load_NAPA_depth <- function(file_name){
  res <- tidync(file_name) %>%
    hyper_filter(x = dplyr::between(x, min(NAPA_arctic$x), max(NAPA_arctic$x)),
                 y = dplyr::between(y, min(NAPA_arctic$y), max(NAPA_arctic$y))) %>%
    hyper_tibble() %>%
    dplyr::select(-starts_with("e3")) %>%
    dplyr::rename(depth = starts_with("depth")) %>% 
    filter(depth < 32.5) %>%
    mutate(depth = plyr::round_any(depth, 5)) %>%
    right_join(NAPA_arctic, by = c("x", "y")) %>% 
    dplyr::rename(t = time_counter) %>%
    mutate(t = as.Date(as.POSIXct(t, origin = "1900-01-01"))) %>% 
    group_by(nav_lon, nav_lat, t, depth) %>% 
    summarise_all(mean, na.rm = T) %>% 
    ungroup() %>% 
    na.omit() %>% 
    mutate(bathy = round(bathy),
           month = lubridate::month(t, label = T)) %>% 
    select(x, y, nav_lon, nav_lat, t, month, depth, bathy, everything())
  return(res)
}
# system.time(test <- load_NAPA_depth(NAPA_depth_U_files)) # 4 seconds


# Create base data.frames -------------------------------------------------

# NB: Can only run one of these at a time
# Must purge memory and restart RStudio before processing a new one

# Surface data
# system.time(Arctic_surface <- plyr::ldply(NAPA_surface_files,
#                                           .fun = load_NAPA_surface,
#                                           .parallel = TRUE)) # 268 seconds
# # It would be better to save these files as .Rds format for the compression
# # but this causes RStudio to hang terribly
# system.time(save(Arctic_surface, file = "data/Arctic_surface.RData")) # 37 seconds, 18.5 GB

# Ice data
# system.time(Arctic_ice_1 <- plyr::ldply(NAPA_ice_files,
#                                         .fun = load_NAPA_surface,
#                                         .parallel = TRUE,
#                                         layer_name = "D2,D1,D3,D0")) # 157 seconds
# system.time(Arctic_ice_2 <- plyr::ldply(NAPA_ice_files,
#                                         .fun = load_NAPA_surface,
#                                         .parallel = TRUE)) # 127 seconds
# system.time(Arctic_ice <- left_join(Arctic_ice_1, Arctic_ice_2, by = c("x", "y", "nav_lon", "nav_lat", "t", "month", "bathy"))) # 90 seconds
# system.time(save(Arctic_ice, file = "data/Arctic_ice.RData")) # 45 seconds, 24.6 GB

# Depth T data
# system.time(Arctic_depth_T <- plyr::ldply(NAPA_depth_T_files,
#                                           .fun = load_NAPA_depth,
#                                           .parallel = TRUE)) # 548 seconds
# system.time(save(Arctic_depth_T, file = "data/Arctic_depth_T.RData")) # 27 seconds, 13.2 GB

# Depth U data
# system.time(Arctic_depth_U <- plyr::ldply(NAPA_depth_U_files,
#                                           .fun = load_NAPA_depth,
#                                           .parallel = TRUE)) # 372 seconds
# system.time(save(Arctic_depth_U, file = "data/Arctic_depth_U.RData")) # 23 seconds, 10.7 GB

# Depth V data
# system.time(Arctic_depth_V <- plyr::ldply(NAPA_depth_V_files,
#                                           .fun = load_NAPA_depth,
#                                           .parallel = TRUE)) # 358 seconds
# system.time(save(Arctic_depth_V, file = "data/Arctic_depth_V.RData")) # 23 seconds, 10.7 GB

# Depth W data
# system.time(Arctic_depth_W <- plyr::ldply(NAPA_depth_W_files,
#                                           .fun = load_NAPA_depth,
#                                           .parallel = TRUE)) # 389 seconds
# system.time(save(Arctic_depth_W, file = "data/Arctic_depth_W.RData")) # 24 seconds, 12.0 GB


# Monthly climatologies ---------------------------------------------------

monthly_clims <- function(df, depth = F){
  # system.time(
  df_clim <- data.table(dplyr::select(df, -t))
  #) # 56 seconds
  if(!depth){
    # system.time(
    setkey(df_clim, nav_lon, nav_lat, month)
    # ) # xxx seconds
    # system.time(
    df_clim <- df_clim[, lapply(.SD, mean), by = list(nav_lon, nav_lat, month)]
    # ) # xxx seconds
  } else if(depth){
    # system.time(
    setkey(df_clim, nav_lon, nav_lat, month, depth)
    # ) # 81 seconds
    # system.time(
    df_clim <- df_clim[, lapply(.SD, mean), by = list(nav_lon, nav_lat, month, depth)]
    # ) # 31 seconds
  } else{
    stop("Something has gone wrong...")
  }
  return(df_clim)
}

# Calculate monthly clims for surface data
# if(!exists("Arctic_surface")) load("data/Arctic_surface.RData")
# system.time(Arctic_surface_clim <- monthly_clims(Arctic_surface, depth = F)) # 47 seconds
# save(Arctic_surface_clim, file = "data/Arctic_surface_clim.RData") # 34.6 MB

# Calculate monthly clims for ice data
# if(!exists("Arctic_ice")) load("data/Arctic_ice.RData")
# system.time(Arctic_ice_clim <- monthly_clims(Arctic_ice, depth = F)) # 55 seconds
# save(Arctic_ice_clim, file = "data/Arctic_ice_clim.RData") # 66.6 MB

# Calculate monthly clims for depth_T data
# if(!exists("Arctic_depth_T")) load("data/Arctic_depth_T.RData")
# system.time(Arctic_depth_T_clim <- monthly_clims(Arctic_depth_T, depth = T)) # 38 seconds
# save(Arctic_depth_T_clim, file = "data/Arctic_depth_T_clim.RData") # 112.1 MB

# Calculate monthly clims for depth_U data
# if(!exists("Arctic_depth_U")) load("data/Arctic_depth_U.RData")
# system.time(Arctic_depth_U_clim <- monthly_clims(Arctic_depth_U, depth = T)) # 35 seconds
# save(Arctic_depth_U_clim, file = "data/Arctic_depth_U_clim.RData") # 88.5 MB

# Calculate monthly clims for depth_V data
# if(!exists("Arctic_depth_V")) load("data/Arctic_depth_V.RData")
# system.time(Arctic_depth_V_clim <- monthly_clims(Arctic_depth_V, depth = T)) # 34 seconds
# save(Arctic_depth_V_clim, file = "data/Arctic_depth_V_clim.RData") # 88.5 MB

# Calculate monthly clims for depth_W data
# if(!exists("Arctic_depth_W")) load("data/Arctic_depth_W.RData")
# system.time(Arctic_depth_W_clim <- monthly_clims(Arctic_depth_W, depth = T)) # 37 seconds
# save(Arctic_depth_W_clim, file = "data/Arctic_depth_W_clim.RData") # 100.3 MB

# NB: Note that the climatology files above are too large to push to GitHub


# Overall means per pixel -------------------------------------------------

overall_means <- function(df, depth = F){
  # system.time(
  df_mean <- data.table(dplyr::select(df, -t, -month))
  #) # 56 seconds
  if(!depth){
    # system.time(
    setkey(df_mean, nav_lon, nav_lat)
    # ) # xxx seconds
    # system.time(
    df_mean <- df_mean[, lapply(.SD, mean), by = list(nav_lon, nav_lat)]
    # ) # xxx seconds
  } else if(depth){
    # system.time(
    setkey(df_mean, nav_lon, nav_lat, depth)
    # ) # 81 seconds
    # system.time(
    df_mean <- df_mean[, lapply(.SD, mean), by = list(nav_lon, nav_lat, depth)]
    # ) # 31 seconds
  } else{
    stop("Something has gone wrong...")
  }
  return(df_mean)
}

# Calculate monthly clims for surface data
# if(!exists("Arctic_surface")) load("data/Arctic_surface.RData")
# system.time(Arctic_surface_mean <- overall_means(Arctic_surface, depth = F)) # 50 seconds
# save(Arctic_surface_mean, file = "data/Arctic_surface_mean.RData") # 2.8 MB

# Calculate monthly clims for ice data
# if(!exists("Arctic_ice")) load("data/Arctic_ice.RData")
# system.time(Arctic_ice_mean <- overall_means(Arctic_ice, depth = F)) # 53 seconds
# save(Arctic_ice_mean, file = "data/Arctic_ice_mean.RData") # 5.5 MB

# Calculate monthly clims for depth_T data
# if(!exists("Arctic_depth_T")) load("data/Arctic_depth_T.RData")
# system.time(Arctic_depth_T_mean <- overall_means(Arctic_depth_T, depth = T)) # 37 seconds
# save(Arctic_depth_T_mean, file = "data/Arctic_depth_T_mean.RData") # 8.8 MB

# Calculate monthly clims for depth_U data
# if(!exists("Arctic_depth_U")) load("data/Arctic_depth_U.RData")
# system.time(Arctic_depth_U_mean <- overall_means(Arctic_depth_U, depth = T)) # 29 seconds
# save(Arctic_depth_U_mean, file = "data/Arctic_depth_U_mean.RData") # 6.9 MB

# Calculate monthly clims for depth_V data
# if(!exists("Arctic_depth_V")) load("data/Arctic_depth_V.RData")
# system.time(Arctic_depth_V_mean <- overall_means(Arctic_depth_V, depth = T)) # 29 seconds
# save(Arctic_depth_V_mean, file = "data/Arctic_depth_V_mean.RData") # 6.9 MB

# Calculate monthly clims for depth_W data
# if(!exists("Arctic_depth_W")) load("data/Arctic_depth_W.RData")
# system.time(Arctic_depth_W_mean <- overall_means(Arctic_depth_W, depth = T)) # 33 seconds
# save(Arctic_depth_W_mean, file = "data/Arctic_depth_W_mean.RData") # 7.9 MB


# Download Bio-ORACLE data ------------------------------------------------

# Explore datasets in the package
list_datasets()

# Explore layers in a dataset
BO_layers <- list_layers(datasets="Bio-ORACLE")

# Check layer statistics
layer_stats()

# Check Pearson correlation coefficient between layers
layers_correlation() 

# Download bathy layers
bathy <- load_layers(c("BO_bathymin", "BO_bathymean", "BO_bathymax"))
bathy_df <- as.data.frame(bathy, xy = T)

## The layes currently chosen for use in this study
  # NB: Many of the varuables below also have surface values
  # NB: Min depth is the deepest as the values are in -m
# Calcite - mean - BO_calcite
# Diffuse attenuation coefficient at 490 nm - mean - BO_damean
# Photosynthetically available radiation - mean - BO_parmean
# pH - mean - BO_ph
# Chl con. - mean at min depth - BO2_chlomean_bdmin
# Current velocity - passing for now - BO2_curvelmean_bdmin
# Dissolved oxygen - mean at min depth - BO2_dissoxmean_bdmin
# Iron con. - mean at min depth -	BO2_ironmean_bdmin
# Phos con. - mean at min depth - BO2_phosphatemean_bdmin
# Light at bottom - mean at min depth - BO2_lightbotmean_bdmin
# Nitr con. - mean at min depth - BO2_nitratemean_bdmin
# sea temp. - mean at min depth - BO2_tempmean_bdmin
# Carbon phytoplankton biomass - mean at min depth - BO2_carbonphytomean_bdmin
# Primary production - mean at min depth - BO2_ppmean_bdmin
# sea salinity - mean at min depth - BO2_salinitymean_bdmin
# Silicate con. - mean at min depth - BO2_silicatemean_bdmin
# Ice con. - mean - BO2_icecovermean_ss
# Ice thickness - mean + range - BO2_icethickmean_ss + BO2_icethickrange_ss

## Download the chosen layers
# BO_layers <- load_layers(c("BO_calcite", "BO_damean", "BO_parmean", "BO_ph",
#                            "BO2_chlomean_bdmin", "BO2_curvelmean_bdmin", "BO2_dissoxmean_bdmin",
#                            "BO2_ironmean_bdmin", "BO2_phosphatemean_bdmin", "BO2_lightbotmean_bdmin",
#                            "BO2_nitratemean_bdmin", "BO2_tempmean_bdmin", "BO2_carbonphytomean_bdmin",
#                            "BO2_ppmean_bdmin", "BO2_salinitymean_bdmin", "BO2_silicatemean_bdmin",
#                            "BO2_icecovermean_ss", "BO2_icethickmean_ss", "BO2_icethickrange_ss"))
# BO_layers_df <- as.data.frame(BO_layers, xy = T)
# Arctic_BO <- BO_layers_df %>% 
#   dplyr::rename(lon = x, lat = y) %>% 
#   filter(lon >= bbox_arctic[1], lon <= bbox_arctic[2],
#          lat >= bbox_arctic[3], lat <= bbox_arctic[4])
# save(Arctic_BO, file = "data/Arctic_BO.RData")

# Visualise
# ggplot(Arctic_BO, aes(x = lon, y = lat)) +
#   geom_tile(aes(fill = BO2_icecovermean_ss))


# Download future Bio-ORACLE data -----------------------------------------

sdmpredictors::get_future_layers()
