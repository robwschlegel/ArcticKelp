# monthly_clims.R
# The purpose of this script is to create monthly climatologies at each pixel
# in the arctic from the NAPA model


# Libraries ---------------------------------------------------------------

library(tidyverse)
library(lubridate)
library(tidync, lib.loc = "../R-packages/")
library(stringr)
# library(ncdf4)

# Set cores
doMC::registerDoMC(cores = 50)


# Data --------------------------------------------------------------------

# The surface files location
NAPA_surface_files <- dir("../../data/NAPA025/1d_grid_T_2D", full.names = T)

# The depth files location
NAPA_depth_files <- dir("../../data/NAPA025/1d_grid_T_2D", full.names = T)

# The ice file location
NAPA_ice_files <- dir("../../data/NAPA025/5d_icemod", full.names = T)

# The NAPA to OISST lon/lat mask
load("metadata/lon_lat_NAPA_OISST.RData")

# The OISST lon values for subsetting
load("metadata/lon_OISST.RData")


# Loading functions -------------------------------------------------------


# Function for loading the individual NAPA NetCDF files and subsetting SST accordingly 
load_NAPA_sst_sub <- function(file_name, coords){
  
  nc <- nc_open(as.character(file_name))
  
  date_start <- ymd(str_sub(basename(as.character(file_name)), start = 29, end = 36))
  date_end <- ymd(str_sub(basename(as.character(file_name)), start = 38, end = 45))
  date_seq <- seq(date_start, date_end, by = "day")
  
  sst <- as.data.frame(ncvar_get(nc, varid = "sst")) %>% 
    mutate(lon = as.numeric(nc$dim$x$vals)) %>%
    gather(-lon, key = lat, value = temp) %>%
    mutate(t = rep(date_seq, each = 388080),
           lat = rep(rep(as.numeric(nc$dim$y$vals), each = 528), times = 5)) %>%
    select(lon, lat, t, temp) %>%
    inner_join(coords, by = c("lon", "lat")) %>% 
    select(nav_lon, nav_lat, t, temp)
  
  nc_close(nc)
  return(sst)
}

# Function for loading the individual NAPA NetCDF files and subsetting SSH accordingly 
load_NAPA_ssh_sub <- function(file_name, coords){
  
  nc <- nc_open(as.character(file_name))
  
  date_start <- ymd(str_sub(basename(as.character(file_name)), start = 29, end = 36))
  date_end <- ymd(str_sub(basename(as.character(file_name)), start = 38, end = 45))
  date_seq <- seq(date_start, date_end, by = "day")
  
  ssh <- as.data.frame(ncvar_get(nc, varid = "ssh")) %>% 
    mutate(lon = as.numeric(nc$dim$x$vals)) %>%
    gather(-lon, key = lat, value = height) %>%
    mutate(t = rep(date_seq, each = 388080),
           lat = rep(rep(as.numeric(nc$dim$y$vals), each = 528), times = 5)) %>%
    select(lon, lat, t, height) %>%
    inner_join(coords, by = c("lon", "lat")) %>% 
    select(nav_lon, nav_lat, t, height)
  
  nc_close(nc)
  return(ssh)
}

# Function for loading the individual NAPA NetCDF files and subsetting SSH accordingly 
load_NAPA_ice_sub <- function(file_name, coords){
  
  nc <- nc_open(as.character(file_name))
  
  date_start <- ymd(str_sub(basename(as.character(file_name)), start = 26, end = 33))
  date_end <- ymd(str_sub(basename(as.character(file_name)), start = 35, end = 42))
  
  ice <- as.data.frame(ncvar_get(nc, varid = "iicefrac")) %>% 
    mutate(lon = as.numeric(nc$dim$x$vals)) %>%
    gather(-lon, key = lat, value = ice) %>%
    mutate(date_start = rep(date_start, 388080),
           date_end = rep(date_end, 388080),
           lat = rep(as.numeric(nc$dim$y$vals), each = 528)) %>%
    select(lon, lat, date_start, date_end, ice) %>%
    inner_join(coords, by = c("lon", "lat")) %>% 
    select(nav_lon, nav_lat, date_start, date_end, ice)
  
  nc_close(nc)
  return(ice)
}


# Saving functions --------------------------------------------------------

# Function for combining and saving the subsetted NAPA SST data
save_NAPA_sst_sub <- function(df){
  
  coords <- lon_lat_NAPA_OISST %>% 
    filter(lon_O == df$lon)
  
  # system.time(
  NAPA_sst_sub <- plyr::ldply(NAPA_files,
                              .fun = load_NAPA_sst_sub, 
                              .parallel = TRUE, 
                              coords = coords)
  # )
  
  lon_sub_label <- str_pad(which(lon_OISST == df$lon), width = 4, pad = "0", side = "left")
  
  if(nrow(NAPA_sst_sub) > 0){
    save(NAPA_sst_sub, file = paste0("../data/NAPA_sst_sub_",lon_sub_label,".RData"))
  }
}

# Function for combining and saving the subsetted NAPA SSH data
save_NAPA_ssh_sub <- function(df){
  
  coords <- lon_lat_NAPA_OISST %>% 
    filter(lon_O == df$lon)
  
  # system.time(
  NAPA_ssh_sub <- plyr::ldply(NAPA_files,
                              .fun = load_NAPA_ssh_sub, 
                              .parallel = TRUE, 
                              coords = coords)
  # )
  
  lon_sub_label <- str_pad(which(lon_OISST == df$lon), width = 4, pad = "0", side = "left")
  
  if(nrow(NAPA_ssh_sub) > 0){
    save(NAPA_ssh_sub, file = paste0("../data/NAPA_ssh_sub_",lon_sub_label,".RData"))
  }
}


# Function for combining and saving the subsetted NAPA ice data
save_NAPA_ice_sub <- function(lon_row){
  lon_sub_label <- str_pad(lon_row, width = 4, pad = "0", side = "left")
  
  coords <- lon_lat_NAPA_OISST %>% 
    filter(lon_O == lon_OISST[lon_row])
  
  # system.time(
  NAPA_ice_sub <- map_dfr(NAPA_ice_files, load_NAPA_ice_sub, coords = coords)
  # ) # 268 seconds
  
  save(NAPA_ice_sub, file = paste0("../data/NAPA_ice_sub_",lon_sub_label,".RData"))
  rm(NAPA_ice_sub)
}

# Extract SST -------------------------------------------------------------

# lon_OISST_multi <- data.frame(lon = lon_OISST,
#                               x = 1:length(lon_OISST))
# 
# # Re-run on Thursday, October 11th, 2018
# system.time(
#   plyr::ddply(lon_OISST_multi[1,], .variables = "x",
#               .fun = save_NAPA_sst_sub)
# ) # 58 seconds at 50 cores
# system.time(
#   plyr::ddply(lon_OISST_multi[2:100,], .variables = "x",
#               .fun = save_NAPA_sst_sub)
# ) # 4850 seconds at 50 cores
# system.time(
#   plyr::ddply(lon_OISST_multi[101:400,], .variables = "x",
#               .fun = save_NAPA_sst_sub)
# ) # 15461 seconds at 50 cores
# system.time(
#   plyr::ddply(lon_OISST_multi[401:700,], .variables = "x",
#               .fun = save_NAPA_sst_sub)
# ) # 17091 seconds at 50 cores
# system.time(
#   plyr::ddply(lon_OISST_multi[701:1000,], .variables = "x",
#               .fun = save_NAPA_sst_sub)
# ) # 16013 seconds at 50 cores
# system.time(
#   plyr::ddply(lon_OISST_multi[1001:1440,], .variables = "x",
#               .fun = save_NAPA_sst_sub)
# ) # 24572 seconds at 50 cores


# Extract SSH -------------------------------------------------------------

# lon_OISST_multi <- data.frame(lon = lon_OISST,
#                               x = 1:length(lon_OISST))
# 
# # Run on Monday, October 29th 2018
# system.time(
#   plyr::ddply(lon_OISST_multi[1:100,], .variables = "x",
#               .fun = save_NAPA_ssh_sub)
# ) # 6382 seconds at 50 cores
# system.time(
#   plyr::ddply(lon_OISST_multi[101:400,], .variables = "x",
#               .fun = save_NAPA_ssh_sub)
# ) # 17828 seconds at 50 cores
# system.time(
#   plyr::ddply(lon_OISST_multi[401:700,], .variables = "x",
#               .fun = save_NAPA_ssh_sub)
# ) # 19462 seconds at 50 cores
# system.time(
#   plyr::ddply(lon_OISST_multi[701:1000,], .variables = "x",
#               .fun = save_NAPA_ssh_sub)
# ) # 17802 seconds at 50 cores
# system.time(
#   plyr::ddply(lon_OISST_multi[1001:1440,], .variables = "x",
#               .fun = save_NAPA_ssh_sub)
# ) # 25423 seconds at 50 cores


# Extract ice -------------------------------------------------------------

# system.time(
#   save_NAPA_ice_sub(1)
# ) ### 269 seconds

# # Run on Wednesday, November 4th 2018
# system.time(
#   plyr::ldply(1:1440, save_NAPA_ice_sub, .parallel = T)
# ) # xxx seconds at 50 cores
# system.time(
#   plyr::ldply(lon_OISST_multi[101:400,], .variables = "x",
#               .fun = save_NAPA_ssh_sub)
# ) # 17828 seconds at 50 cores
# system.time(
#   plyr::ldply(lon_OISST_multi[401:700,], .variables = "x",
#               .fun = save_NAPA_ssh_sub)
# ) # 19462 seconds at 50 cores
# system.time(
#   plyr::ldply(lon_OISST_multi[701:1000,], .variables = "x",
#               .fun = save_NAPA_ssh_sub)
# ) # 17802 seconds at 50 cores
# system.time(
#   plyr::ldply(lon_OISST_multi[1001:1440,], .variables = "x",
#               .fun = save_NAPA_ssh_sub)
# ) # 25423 seconds at 50 cores


# Visualise ---------------------------------------------------------------

# NAPA_saves <- dir("../data", pattern = "NAPA_sst", full.names = T)
# 
# load_NAPA <- function(file_name){
#   load(file = file_name)
#   data_sub <- NAPA_sst_sub %>% 
#     filter(t == as.Date("1993-10-01"))
#   return(data_sub)
# }
# 
# system.time(
#   NAPA_all_sub_1 <- plyr::ldply(NAPA_saves[1:360], .parallel = T,
#                           .fun = load_NAPA)
# ) # 38 seconds at 50 cores
# system.time(
#   NAPA_all_sub_2 <- plyr::ldply(NAPA_saves[361:720], .parallel = T,
#                                 .fun = load_NAPA)
# ) # 42 seconds at 50 cores
# system.time(
#   NAPA_all_sub_3 <- plyr::ldply(NAPA_saves[721:1090], .parallel = T,
#                                 .fun = load_NAPA)
# ) # 42 seconds at 50 cores
# system.time(
#   NAPA_all_sub_4 <- plyr::ldply(NAPA_saves[1091:1440], .parallel = T,
#                                 .fun = load_NAPA)
# ) # 39 seconds at 50 cores
# 
# NAPA_all_sub <- rbind(NAPA_all_sub_1, NAPA_all_sub_2, NAPA_all_sub_3, NAPA_all_sub_4)
# 
# ggplot(NAPA_all_sub, aes(x = nav_lon, y = nav_lat, colour = temp)) +
#   geom_point(size = 0.001) +
#   scale_colour_viridis_c()
# 
# ggplot(NAPA_all_sub, aes(x = -nav_lon, y = -nav_lat, colour = temp)) +
#   geom_point(size = 0.001) +
#   scale_colour_viridis_c() +
#   coord_polar() +
#   scale_x_continuous(expand = c(0, 0)) +
#   scale_y_continuous(expand = c(0, 0)) +
#   labs(x = "", y = "")
