# analyses/3_study_site_layers.R
# The purpose of this script is to pull out the data from the layers at each study site pixel.
# Whenever new sites are added re-run this script to get the data for the new pixels.
# Otherwise the analyses properly start from "4_kelp_cover.R".


# Setup -------------------------------------------------------------------

# Load study sites and base packages
source("analyses/1_study_region_sites.R")

# Additional packages
library(FNN)


# Study site BO data ------------------------------------------------------

# NB: The following code chunks require Arctic_coast.RData
# This is not hosted on GitHub as it is 30.3 MB
# E-mail robert.schlegel@dal.ca for the file
# Or create it from '2_monthly_clims.R'
load("data/Arctic_coast.RData")
Arctic_coast$env_index <- 1:nrow(Arctic_coast)

# Find the nearest BO points to each site and add bathy data
study_site_env <- study_sites %>%
  mutate(env_index = as.vector(knnx.index(as.matrix(Arctic_coast[,c("lon", "lat")]),
                                          as.matrix(study_sites[,c("lon", "lat")]), k = 1))) %>%
  left_join(Arctic_coast, by = "env_index") %>%
  dplyr::select(-Date, -Notes, -env_index) %>%
  dplyr::rename(lon = lon.x, lat = lat.x, lon_env = lon.y, lat_env = lat.y) %>% 
  dplyr::select(site:lat_env, depth, land_distance, everything())
save(study_site_env, file = "data/study_site_env.RData")
Arctic_coast$env_index <- NULL


# Future layers -----------------------------------------------------------

# NB: These layers are not hosted on GitHub, contact Robert for access
load("data/Arctic_BO_2050.RData")
load("data/Arctic_BO_2100.RData")

# Create 2050 study site data
study_site_env_2050 <- study_site_env %>% 
  dplyr::select(site:land_distance) %>% 
  left_join(Arctic_BO_2050, by = c("lon_env" = "lon", "lat_env" = "lat"))
save(study_site_env_2050, file = "data/study_site_env_2050.RData")

# Create 2100 study site data
study_site_env_2100 <- study_site_env %>% 
  dplyr::select(site:land_distance) %>% 
  left_join(Arctic_BO_2100, by = c("lon_env" = "lon", "lat_env" = "lat"))
save(study_site_env_2100, file = "data/study_site_env_2100.RData")

