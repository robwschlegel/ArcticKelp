# monthly_clims.R
# The purpose of this script is to create monthly climatologies at each pixel
# in the arctic from the NAPA model


# Libraries ---------------------------------------------------------------

library(tidyverse)
library(lubridate)
library(tidync, lib.loc = "../R-packages/")
library(stringr)
library(data.table)

# Set cores
doMC::registerDoMC(cores = 50)


# NetCDF information ------------------------------------------------------

# Check the grid data
# ncdump::NetCDF("../../data/NAPA025/mesh_grid/bathy_creg025_extended_5m.nc")

# Check the surface NetCDF information
# ncdump::NetCDF("../../data/NAPA025/1d_grid_T_2D/CREG025E-GLSII_1d_grid_T_2D_19931001-19931005.nc")

# Check the depth NetCDF information
# ncdump::NetCDF("../../data/NAPA025/5d_grid_T/CREG025E-GLSII_5d_grid_T_19931001-19931005.nc")$variable$longname
# ncdump::NetCDF("../../data/NAPA025/5d_grid_U/CREG025E-GLSII_5d_grid_U_19931001-19931005.nc")#$variable$longname
# ncdump::NetCDF("../../data/NAPA025/5d_grid_V/CREG025E-GLSII_5d_grid_V_19931001-19931005.nc")#$variable$longname
# ncdump::NetCDF("../../data/NAPA025/5d_grid_W/CREG025E-GLSII_5d_grid_W_19931001-19931005.nc")#$variable$longname

# Check the ice NetCDF information
# ncdump::NetCDF("../../data/NAPA025/5d_icemod/CREG025E-GLSII_5d_icemod_19931001-19931005.nc")#$variable$longname


# Data --------------------------------------------------------------------

# The arctic study area extent
bbox_arctic <- c(-100, -40, 50, 80)

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
#          bathy > 0)
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
           month = lubridate::month(t, label = T))
  return(res)
}
# system.time(test <- load_NAPA_surface(NAPA_surface_files[111])) # 1 second

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
# system.time(test <- load_NAPA_depth(NAPA_depth_T_files)) # 5 seconds

# Function for loading the individual NAPA NetCDF ice files
# NB: This hasn't been tested yet
# testers...
# file_name <- NAPA_ice_files[1]
load_NAPA_ice <- function(file_name){
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
           month = lubridate::month(t, label = T))
  return(res)
}

# Create base data.frames -------------------------------------------------

# NB: Can only run one of these at a time
# Must purge memory and restart RStudio before processing a new one

# Surface data
# system.time(Arctic_surface <- plyr::ldply(NAPA_surface_files,
#                               .fun = load_NAPA_surface, 
#                               .parallel = TRUE)) # 422 seconds
# It would be better to save these files as .Rds format for the compression
# but this causes RSstudio to hang terribly
# system.time(saveRDS(Arctic_surface, file = "data/Arctic_surface.Rds")) # xxx seconds, xxx GB
# system.time(save(Arctic_surface, file = "data/Arctic_surface.RData")) # xxx seconds, xxx GB

# Depth T data
# system.time(Arctic_depth_T <- plyr::ldply(NAPA_depth_T_files,
#                                           .fun = load_NAPA_depth,
#                                           .parallel = TRUE)) # 451 seconds
# system.time(save(Arctic_depth_T, file = "data/Arctic_depth_T.RData")) # 58 seconds, 34.9 GB


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

# Calculate monthly clims for depth_T data
system.time(Arctic_depth_T_clim <- monthly_clims(Arctic_depth_T, depth = T)) # xxx seconds
save(Arctic_depth_T_clim, file = "data/Arctic_depth_T_clim.RData")

# system.time(Arctic_surface_clim <- Arctic_surface %>% 
#               mutate(month = month(t, label = T)) %>% 
#               group_by(nav_lon, nav_lat, month) %>% 
#               summarise_all(.funs = mean, na.rm = T)) # xxx seconds

# Switch to data.table for faster means
# Arctic_surface_clim <- Arctic_surface %>% 
#   mutate(month = month(t, label = T)) %>% 
#   dplyr::select(-t) %>% 
#   data.table()
# setkey(Arctic_surface_clim, nav_lon, nav_lat, month)
# system.time(Arctic_surface_clim <- Arctic_surface[, lapply(.SD, mean), 
#                                                   by = list(nav_lon, nav_lat, month)]) # xxx seconds
# system.time(saveRDS(Arctic_surface, file = "data/Arctic_surface.Rds")) # xxx seconds, xxx MB


# Visualise ---------------------------------------------------------------

# Daily data
ggplot(filter(Arctic_depth_T, t == "1998-01-05"), 
       aes(x = nav_lon, y = nav_lat, colour = eken)) +
  geom_point(size = 0.001) +
  scale_colour_viridis_c()

 # Monthly clims
ggplot(filter(Arctic_depth_T_clim, month == "Jan"), 
       aes(x = nav_lon, y = nav_lat, colour = eken)) +
  geom_point(size = 0.001) +
  scale_colour_viridis_c()

# ggplot(NAPA_all_sub, aes(x = -nav_lon, y = -nav_lat, colour = temp)) +
#   geom_point(size = 0.001) +
#   scale_colour_viridis_c() +
#   coord_polar() +
#   scale_x_continuous(expand = c(0, 0)) +
#   scale_y_continuous(expand = c(0, 0)) +
#   labs(x = "", y = "")
