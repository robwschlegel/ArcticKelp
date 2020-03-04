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
bathy <- load_layers(c("BO_bathymin", "BO_bathymean", "BO_bathymax"))
Arctic_bathy <- as.data.frame(bathy, xy = T) %>% 
  dplyr::rename(lon = x, lat = y) %>% 
  filter(lon >= bbox_arctic[1], lon <= bbox_arctic[2],
         lat >= bbox_arctic[3], lat <= bbox_arctic[4])
save(Arctic_bathy, file = "data/Arctic_bathy.RData")

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
BO_layers_dl <- load_layers(c("BO_calcite", "BO_damean", "BO_parmean", "BO_ph",
                              "BO2_chlomin_bdmax", "BO2_chlomean_bdmax", "BO2_chlomax_bdmax", 
                              "BO2_curvelmin_bdmax", "BO2_curvelmean_bdmax", "BO2_curvelmax_bdmax", 
                              "BO2_dissoxmin_bdmax", "BO2_dissoxmean_bdmax", "BO2_dissoxmax_bdmax", 
                              "BO2_ironmin_bdmax", "BO2_ironmean_bdmax", "BO2_ironmax_bdmax", 
                              "BO2_phosphatemin_bdmax", "BO2_phosphatemean_bdmax", "BO2_phosphatemax_bdmax", 
                              "BO2_lightbotmin_bdmax", "BO2_lightbotmean_bdmax", "BO2_lightbotmax_bdmax", 
                              "BO2_nitratemin_bdmax", "BO2_nitratemean_bdmax", "BO2_nitratemax_bdmax", 
                              "BO2_tempmin_bdmax", "BO2_tempmean_bdmax", "BO2_tempmax_bdmax", 
                              "BO2_templtmin_bdmean", "BO2_templtmax_bdmean",
                              "BO2_carbonphytomin_bdmax", "BO2_carbonphytomean_bdmax", "BO2_carbonphytomax_bdmax",
                              "BO2_ppmin_bdmax", "BO2_ppmean_bdmax", "BO2_ppmax_bdmax", 
                              "BO2_salinitymin_bdmax", "BO2_salinitymean_bdmax", "BO2_salinitymax_bdmax", 
                              "BO2_silicatemin_bdmax", "BO2_silicatemean_bdmax", "BO2_silicatemax_bdmax",
                              "BO2_icecoverltmin_ss",  "BO2_icecoverltmax_ss",
                              "BO2_icecovermin_ss", "BO2_icecovermean_ss", "BO2_icecovermax_ss",
                              "BO2_icethickmin_ss", "BO2_icethickmean_ss", "BO2_icethickmax_ss", 
                              "BO2_icethickrange_ss"))
BO_layers_df <- as.data.frame(BO_layers_dl, xy = T)
Arctic_BO <- BO_layers_df %>%
  dplyr::rename(lon = x, lat = y) %>%
  filter(lon >= bbox_arctic[1], lon <= bbox_arctic[2],
         lat >= bbox_arctic[3], lat <= bbox_arctic[4])
save(Arctic_BO, file = "data/Arctic_BO.RData")

# Visualise
ggplot(Arctic_BO, aes(x = lon, y = lat)) +
  geom_tile(aes(fill = BO2_icecovermean_ss)) +
  coord_quickmap(expand = F) +
  theme(legend.position = "bottom")


# Download future Bio-ORACLE data -----------------------------------------

# Future scenario conversions
# https://ar5-syr.ipcc.ch/topic_futurechanges.php
# RCP8.5 ~= A2
# RCP6.0 ~= B2
# RCP4.5 ~= B1

list_layers_future()

get_future_layers()
