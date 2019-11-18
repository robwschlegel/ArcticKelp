# analyses/2_monthly_clims.R
# The purpose of this script is to create monthly climatologies at each pixel
# in the arctic from the NAPA model
# This should only need to be run once at the outset of the project.
# If any sites change/are added, run "analyses/study_sites_clims.R" 
# to get the updated values, not this script

# NB: This script will only run on Eric Oliver's tikoraluk server
# This is because it is utilising files that are up to 25GB in size
# This shouldn't be an issue though as it only needs to be run once


# Libraries ---------------------------------------------------------------

# The study sites and bounding box
source("analyses/1_study_sites.R")

# Other packages
library(tidync)
library(stringr)
library(data.table)

# Set cores
doMC::registerDoMC(cores = 50)


# NetCDF information ------------------------------------------------------

# NB: The code in this section is greyed out as it will only run correctly on the server where the model data live
  # It also only needs to be run once, so is also just stored here for future convenience

# Check the grid data
# info_grid <- ncdump::NetCDF("../../data/NAPA025/mesh_grid/bathy_creg025_extended_5m.nc")$variable %>% 
#   dplyr::select(-addOffset, -scaleFact)

# Check the surface NetCDF information
# info_surface <- ncdump::NetCDF("../../data/NAPA025/1d_grid_T_2D/CREG025E-GLSII_1d_grid_T_2D_19931001-19931005.nc")$variable

# Check the ice NetCDF information
# info_ice <- ncdump::NetCDF("../../data/NAPA025/5d_icemod/CREG025E-GLSII_5d_icemod_19931001-19931005.nc")$variable#$longname
# info_ice_2 <- info$variable

# Check the depth NetCDF information
# info_depth_T <- ncdump::NetCDF("../../data/NAPA025/5d_grid_T/CREG025E-GLSII_5d_grid_T_19931001-19931005.nc")$variable#$longname
# info_depth_U <- ncdump::NetCDF("../../data/NAPA025/5d_grid_U/CREG025E-GLSII_5d_grid_U_19931001-19931005.nc")$variable#$longname
# info_depth_V <- ncdump::NetCDF("../../data/NAPA025/5d_grid_V/CREG025E-GLSII_5d_grid_V_19931001-19931005.nc")$variable#$longname
# info_depth_W <- ncdump::NetCDF("../../data/NAPA025/5d_grid_W/CREG025E-GLSII_5d_grid_W_19931001-19931005.nc")$variable#$longname

# Combine and clean-up
# info_ALL <- rbind(info_grid, info_surface, info_ice, info_depth_T, info_depth_U, info_depth_V, info_depth_W) %>%
#   dplyr::select(name, ndims, units, longname) %>%
#   unique() %>% 
#   arrange(name)

# Save
# write_csv(info_ALL, "metadata/model_info.csv")


# Data --------------------------------------------------------------------

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

# The depth files location
NAPA_depth_T_files <- dir("../../data/NAPA025/5d_grid_T/", full.names = T)[-c(1:311)]
# head(NAPA_depth_files)
NAPA_depth_U_files <- dir("../../data/NAPA025/5d_grid_U/", full.names = T)[-c(1:311)]
NAPA_depth_V_files <- dir("../../data/NAPA025/5d_grid_V/", full.names = T)[-c(1:311)]
NAPA_depth_W_files <- dir("../../data/NAPA025/5d_grid_W/", full.names = T)[-c(1:311)]

# The ice file location
NAPA_ice_files <- dir("../../data/NAPA025/5d_icemod", full.names = T)[-c(1:311)]
# head(NAPA_ice_files)

# A custom rainbow palette
# A rainbow colour palette was explicitly requested
rainbow_palette <- c("#fefefe", "#f963fa", "#020135", "#00efe1", "#057400", "#fcfd00", "#ed0000", "#3d0000")


# Loading functions -------------------------------------------------------

# Function for loading the individual NAPA NetCDF surface files
# testers...
# file_name <- NAPA_surface_files[1]
load_NAPA_surface <- function(file_name){
  res <- tidync(file_name) %>%
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
    select(x, y, nav_lon, nav_lat, t, month, everything())
  if("ncatice" %in% colnames(res)){
    res <- res %>% 
      dplyr::select(-ncatice) %>% 
      unique() %>% 
      mutate(iceconc_cat = round(iceconc_cat, 2),
             icethic_cat = round(icethic_cat, 2))
  }
  return(res)
}
# system.time(test <- load_NAPA_surface(NAPA_surface_files[111])) # 1 second
# system.time(test <- load_NAPA_surface(NAPA_ice_files[111])) # 1 second

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
# system.time(test <- load_NAPA_depth(NAPA_depth_T_files)) # 4 seconds


# Create base data.frames -------------------------------------------------

# NB: Can only run one of these at a time
# Must purge memory and restart RStudio before processing a new one

# Surface data
# system.time(Arctic_surface <- plyr::ldply(NAPA_surface_files,
#                                           .fun = load_NAPA_surface,
#                                           .parallel = TRUE)) # 276 seconds
# It would be better to save these files as .Rds format for the compression
# but this causes RStudio to hang terribly
# system.time(saveRDS(Arctic_surface, file = "data/Arctic_surface.Rds")) # xxx seconds, xxx GB
# system.time(save(Arctic_surface, file = "data/Arctic_surface.RData")) # 46 seconds, 24.9 GB

# Depth T data
# system.time(Arctic_depth_T <- plyr::ldply(NAPA_depth_T_files,
#                                           .fun = load_NAPA_depth,
#                                           .parallel = TRUE)) # 451 seconds
# system.time(save(Arctic_depth_T, file = "data/Arctic_depth_T.RData")) # 58 seconds, 17.9 GB

# Ice data
# system.time(Arctic_ice <- plyr::ldply(NAPA_ice_files,
#                                       .fun = load_NAPA_surface,
#                                       .parallel = TRUE)) # 134 seconds
# system.time(save(Arctic_ice, file = "data/Arctic_ice.RData")) # 12 seconds, 6.2 GB


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
# system.time(Arctic_surface_clim <- monthly_clims(Arctic_surface, depth = F)) # 186 seconds
# save(Arctic_surface_clim, file = "data/Arctic_surface_clim.RData") # 46.6 MB

# Calculate monthly clims for depth_T data
# if(!exists("Arctic_depth_T")) load("data/Arctic_depth_T.RData")
# system.time(Arctic_depth_T_clim <- monthly_clims(Arctic_depth_T, depth = T)) # 74 seconds
# save(Arctic_depth_T_clim, file = "data/Arctic_depth_T_clim.RData") # 151 MB

# Calculate monthly clims for ice data
# if(!exists("Arctic_ice")) load("data/Arctic_ice.RData")
# system.time(Arctic_ice_clim <- monthly_clims(Arctic_ice, depth = F)) # 31 seconds
# save(Arctic_ice_clim, file = "data/Arctic_ice_clim.RData") # 17 MB

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
# system.time(Arctic_surface_mean <- overall_means(Arctic_surface, depth = F)) # 67 seconds
# save(Arctic_surface_mean, file = "data/Arctic_surface_mean.RData") # 3.8 MB

# Calculate monthly clims for depth_T data
# if(!exists("Arctic_depth_T")) load("data/Arctic_depth_T.RData")
# system.time(Arctic_depth_T_mean <- overall_means(Arctic_depth_T, depth = T)) # 47 seconds
# save(Arctic_depth_T_mean, file = "data/Arctic_depth_T_mean.RData") # 11.9 MB

# Calculate monthly clims for ice data
# if(!exists("Arctic_ice")) load("data/Arctic_ice.RData")
# system.time(Arctic_ice_mean <- overall_means(Arctic_ice, depth = F)) # 15 seconds
# save(Arctic_ice_mean, file = "data/Arctic_ice_mean.RData") # 1.3 MB

