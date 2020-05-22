# analyses/3_study_site_clims
# The purpose of this script is to pull out the monthly clims and mean states at each study site
# Whenever new sites are added re-run this script to get the clims for the new pixels.
# Otherwise the analyses properly start from "4_kelp_cover.R".


# Setup -------------------------------------------------------------------

# Load study sites and base packages
source("analyses/1_study_sites.R")

# Additional packages
library(FNN)


# Study site BO data ------------------------------------------------------

# NB: The following code chunks require Arctic_env.RData
# This is not hosted on GitHub as it is 17.3 MB
# E-mail robert.schlegel@dal.ca for the file
# Or create it from '2_monthly_clims.R'
load("data/Arctic_env.RData")
Arctic_env$env_index <- 1:nrow(Arctic_env)

# Find the nearest BO points to each site and add bathy data
study_site_env <- study_sites %>%
  mutate(env_index = as.vector(knnx.index(as.matrix(Arctic_env[,c("lon", "lat")]),
                                          as.matrix(study_sites[,c("lon", "lat")]), k = 1))) %>%
  left_join(Arctic_env, by = "env_index") %>%
  dplyr::select(-Date, -Notes) %>%
  dplyr::rename(lon = lon.x, lat = lat.x, lon_env = lon.y, lat_env = lat.y)
save(study_site_env, file = "data/study_site_env.RData")
Arctic_env$env_index <- NULL
