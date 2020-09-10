# tests/test_BO.R
# The purpose of this script is to test the BioOracle data

# The study sites and bounding box
source("analyses/1_study_region_sites.R")

# Bio-Oracle access
library(sdmpredictors)


# Tests of bottom current layers ------------------------------------------

# Test bottom currents
current_bdmax_layers <- load_layers(c("BO2_curvelltmin_bdmax", "BO2_curvelmean_bdmax", "BO2_curvelltmax_bdmax"))
current_bdmax_test <- as.data.frame(current_bdmax_layers, xy = T) %>% 
  dplyr::rename(lon = x, lat = y) %>% 
  mutate(lon = round(lon, 4), 
         lat = round(lat, 4)) %>% 
  na.omit() %>% 
  mutate(max_min = ifelse(BO2_curvelltmax_bdmax > BO2_curvelltmin_bdmax, TRUE, FALSE))

# Visualise pixels where the max and min values are not as expected
current_bdmax_global <- ggplot(data = current_bdmax_test, aes(x = lon, y = lat)) +
  geom_raster(aes(fill = max_min)) +
  coord_quickmap(expand = F) +
  labs(fill = "Max greater than min", x = NULL, y = NULL,
       title = "Bottom currents downloaded via R on May 27th") +
  theme(legend.position = "bottom")
ggsave(plot = current_bdmax_global, filename = "tests/current_bdmax_global_R.png", height = 5, width = 8)

# Test the first set of bottom layers Jorge sent
curvel_bdmax_min_old <- as.data.frame(raster("data/SeaWaterVelocity Benthic Mean Pred LtMin old.tif"), xy = T) %>% 
  `colnames<-`(c("lon", "lat", "curvel_bdmax_min"))  %>% 
  mutate(lon = round(lon, 4),
         lat = round(lat, 4))
curvel_bdmax_max_old <- as.data.frame(raster("data/SeaWaterVelocity Benthic Mean Pred LtMax old.tif"), xy = T) %>% 
  `colnames<-`(c("lon", "lat", "curvel_bdmax_max")) %>% 
  mutate(lon = round(lon, 4),
         lat = round(lat, 4))
curvel_bdmax_old <- left_join(curvel_bdmax_min_old, curvel_bdmax_max_old, by = c("lon", "lat")) %>% 
  na.omit() %>% 
  mutate(max_min = ifelse(curvel_bdmax_max >= curvel_bdmax_min, TRUE, FALSE)) %>% 
  ggplot(aes(x = lon, y = lat)) +
  geom_raster(aes(fill = max_min)) +
  coord_quickmap(expand = F) +
  labs(fill = "Max greater than min", x = NULL, y = NULL,
       title = "Bottom currents from Dropbox files on May 22nd") +
  theme(legend.position = "bottom")
ggsave(plot = curvel_bdmax_old , filename = "tests/current_bdmax_global_old.png", height = 5, width = 8)

# Test the second set of bottom layers Jorge sent
curvel_bdmax_min_new <- as.data.frame(raster("data/SeaWaterVelocity Benthic Mean Pred LtMin new.tif"), xy = T) %>% 
  `colnames<-`(c("lon", "lat", "curvel_bdmax_min"))  %>% 
  mutate(lon = round(lon, 4),
         lat = round(lat, 4))
curvel_bdmax_max_new <- as.data.frame(raster("data/SeaWaterVelocity Benthic Mean Pred LtMax new.tif"), xy = T) %>% 
  `colnames<-`(c("lon", "lat", "curvel_bdmax_max")) %>% 
  mutate(lon = round(lon, 4),
         lat = round(lat, 4))
curvel_bdmax_new <- left_join(curvel_bdmax_min_new, curvel_bdmax_max_new, by = c("lon", "lat")) %>% 
  na.omit() %>% 
  mutate(max_min = ifelse(curvel_bdmax_max >= curvel_bdmax_min, TRUE, FALSE)) %>% 
  ggplot(aes(x = lon, y = lat)) +
  geom_raster(aes(fill = max_min)) +
  coord_quickmap(expand = F) +
  labs(fill = "Max greater than min", x = NULL, y = NULL,
       title = "Bottom currents from Dropbox files on May 26th") +
  theme(legend.position = "bottom")
ggsave(plot = curvel_bdmax_new , filename = "tests/current_bdmax_global_new.png", height = 5, width = 8)

# Test the v2.1 layers
  # Note that these layers are not on GitHub as they are too large
  # They may be downloaded here: https://www.bio-oracle.org/downloads-to-email.php
curvel_bdmax_min_21 <- as.data.frame(raster("data/Present.Benthic.Max.Depth.Current.Velocity.Lt.min.asc.BOv2_1.asc"), xy = T) %>% 
  `colnames<-`(c("lon", "lat", "curvel_bdmax_min"))  %>% 
  mutate(lon = round(lon, 4),
         lat = round(lat, 4))
curvel_bdmax_max_21 <- as.data.frame(raster("data/Present.Benthic.Max.Depth.Current.Velocity.Lt.max.asc.BOv2_1.asc"), xy = T) %>% 
  `colnames<-`(c("lon", "lat", "curvel_bdmax_max")) %>% 
  mutate(lon = round(lon, 4),
         lat = round(lat, 4))
curvel_bdmax_21 <- left_join(curvel_bdmax_min_21, curvel_bdmax_max_21, by = c("lon", "lat")) %>% 
  na.omit() %>% 
  mutate(max_min = ifelse(curvel_bdmax_max >= curvel_bdmax_min, TRUE, FALSE)) %>% 
  ggplot(aes(x = lon, y = lat)) +
  geom_raster(aes(fill = max_min)) +
  coord_quickmap(expand = F) +
  labs(fill = "Max greater than min", x = NULL, y = NULL,
       title = "Bottom currents from v2.1") +
  theme(legend.position = "bottom")
ggsave(plot = curvel_bdmax_21 , filename = "tests/current_bdmax_global_21.png", height = 5, width = 8)


# Test surface current layer ----------------------------------------------

# Test surface currents
curvel_ss_layers <- load_layers(c("BO2_curvelltmin_ss", "BO2_curvelmean_ss", "BO2_curvelltmax_ss"))
curvel_ss_test <- as.data.frame(curvel_ss_layers, xy = T) %>% 
  dplyr::rename(lon = x, lat = y) %>% 
  mutate(lon = round(lon, 4), 
         lat = round(lat, 4)) %>% 
  na.omit() %>% 
  mutate(max_min = ifelse(BO2_curvelltmax_ss > BO2_curvelltmin_ss, TRUE, FALSE))
curvel_ss_global <- ggplot(data = curvel_ss_test, aes(x = lon, y = lat)) +
  geom_raster(aes(fill = max_min)) +
  coord_quickmap(expand = F) +
  labs(fill = "Max greater than min", x = NULL, y = NULL,
       title = "Surface currents downloaded via R on May 27th") +
  theme(legend.position = "bottom")
ggsave(plot = curvel_ss_global, filename = "tests/curvel_ss_global.png", height = 5, width = 8)

# Test surface currents for v2.1
curvel_ss_min_21 <- as.data.frame(raster("data/Present.Surface.Current.Velocity.Lt.min.asc.BOv2_1.asc"), xy = T) %>% 
  `colnames<-`(c("lon", "lat", "curvel_ss_min"))  %>% 
  mutate(lon = round(lon, 4),
         lat = round(lat, 4))
curvel_ss_max_21 <- as.data.frame(raster("data/Present.Surface.Current.Velocity.Lt.max.asc.BOv2_1.asc"), xy = T) %>% 
  `colnames<-`(c("lon", "lat", "curvel_ss_max")) %>% 
  mutate(lon = round(lon, 4),
         lat = round(lat, 4))
curvel_ss_21 <- left_join(curvel_ss_min_21, curvel_ss_max_21, by = c("lon", "lat")) %>% 
  na.omit() %>% 
  mutate(max_min = ifelse(curvel_ss_max >= curvel_ss_min, TRUE, FALSE)) %>% 
  ggplot(aes(x = lon, y = lat)) +
  geom_raster(aes(fill = max_min)) +
  coord_quickmap(expand = F) +
  labs(fill = "Max greater than min", x = NULL, y = NULL,
       title = "Surface currents from v2.1") +
  theme(legend.position = "bottom")
ggsave(plot = curvel_ss_21, filename = "tests/curvel_ss_global_21.png", height = 5, width = 8)


# Test SST layer ----------------------------------------------------------

# Test SST
SST_layers <- load_layers(c("BO2_templtmin_ss", "BO2_tempmean_ss", "BO2_templtmax_ss"))
SST_test <- as.data.frame(SST_layers, xy = T) %>% 
  dplyr::rename(lon = x, lat = y) %>% 
  mutate(lon = round(lon, 4), 
         lat = round(lat, 4)) %>% 
  na.omit() %>% 
  mutate(max_min = ifelse(BO2_templtmax_ss > BO2_templtmin_ss, TRUE, FALSE))

SST_global <-  ggplot(data = SST_test, aes(x = lon, y = lat)) +
  geom_raster(aes(fill = max_min)) +
  coord_quickmap(expand = F) +
  labs(fill = "Max greater than min", x = NULL, y = NULL) +
  theme(legend.position = "bottom")
ggsave(plot = SST_global, filename = "tests/SST_global.png", height = 5, width = 8)

# Compare data layer downloaded via R against a layer download manually
SST_mean_manual <- as.data.frame(read.asciigrid("data/Present.Surface.Temperature.Mean.asc"), xy = T) %>% 
  `colnames<-`(c("SST_mean", "lon", "lat")) %>% 
  mutate(lon = round(lon, 4), 
         lat = round(lat, 4)) %>% 
  na.omit() %>% 
  as.data.frame() %>% 
  left_join(SST_test, by = c("lon", "lat")) %>% 
  mutate(BO2_tempmean_ss = round(BO2_tempmean_ss, 6),
         SST_mean = round(SST_mean, 6),
         mean_comp = round(BO2_tempmean_ss-SST_mean, 6)) %>% 
  dplyr::select(lon, lat, mean_comp)

# Plot the difference
SST_mean_diff <- ggplot(SST_mean_manual, aes(x = lon, y = lat)) +
  geom_raster(aes(fill = mean_comp)) +
  coord_quickmap(expand = F) +
  labs(fill = "R download minus \nmanual download (°C)", x = NULL, y = NULL,
       title = "Difference in SST between layer downloaded in R vs. manually from BO website") +
  theme(legend.position = "bottom")
ggsave(plot = SST_mean_diff, filename = "tests/SST_global_mean_diff.png", height = 5, width = 8)

# Test future SST for v2.1
SST_min_21 <- as.data.frame(raster("data/2050AOGCM.RCP85.Surface.Temperature.Lt.min.asc.BOv2_1.asc"), xy = T) %>% 
  `colnames<-`(c("lon", "lat", "SST_min"))  %>% 
  mutate(lon = round(lon, 4),
         lat = round(lat, 4))
SST_max_21 <- as.data.frame(raster("data/2050AOGCM.RCP85.Surface.Temperature.Lt.max.asc.BOv2_1.asc"), xy = T) %>% 
  `colnames<-`(c("lon", "lat", "SST_max")) %>% 
  mutate(lon = round(lon, 4),
         lat = round(lat, 4))
SST_21 <- left_join(SST_min_21, SST_max_21, by = c("lon", "lat")) %>% 
  na.omit() %>% 
  mutate(max_min = ifelse(SST_max >= SST_min, TRUE, FALSE)) %>% 
  ggplot(aes(x = lon, y = lat)) +
  geom_raster(aes(fill = max_min)) +
  coord_quickmap(expand = F) +
  labs(fill = "Max greater than min", x = NULL, y = NULL,
       title = "SST 2050 RCP 8.5 from v2.1") +
  theme(legend.position = "bottom")
ggsave(plot = SST_21, filename = "tests/SST_2050_RCP85_21.png", height = 5, width = 8)


# Quick tests of Arctic area only -----------------------------------------

# Load the Arctic cropped data and check all remaining min max layers
load("data/Arctic_env.RData")

# Function for comparing max and min of a chosen variable
# chosen_var <- "BO2_templt..._bdmax"
# chosen_var <- "BO2_RCP85_2050_templt..._bdmax"
# chosen_var <- "BO2_RCP85_2050_salinitylt..._bdmax"
max_min_comp <- function(chosen_var){
  
  # Set chosen columns
  col_sub <- c("lon", "lat", colnames(Arctic_env)[grep(chosen_var, colnames(Arctic_env))])
  Arctic_sub <- Arctic_env[, col_sub]
  
  # Find max and min columns
  max_col <- grep("max_", colnames(Arctic_sub))
  min_col <- grep("min_", colnames(Arctic_sub))
  
  # Calculate if max is greater than min
  Arctic_sub %>% 
    na.omit() %>% 
    mutate(max_min = ifelse(Arctic_sub[max_col] >= Arctic_sub[min_col], TRUE, FALSE)) %>% 
    ggplot(aes(x = lon, y = lat)) +
    geom_raster(aes(fill = max_min)) +
    coord_quickmap(expand = F) +
    labs(fill = "Max greater than min", x = NULL, y = NULL,
         title = chosen_var) +
    theme(legend.position = "bottom")
}

## The tests
colnames(Arctic_env)
# Present
max_min_comp("BO2_templt..._bdmax")
max_min_comp("BO2_templt..._ss")
max_min_comp("BO2_salinitylt..._bdmax")
max_min_comp("BO2_salinitylt..._ss")
max_min_comp("BO2_icethicklt..._ss")
max_min_comp("BO2_dissoxlt..._bdmax")
max_min_comp("BO2_ironlt..._bdmax")
max_min_comp("BO2_nitratelt..._bdmax")
max_min_comp("BO2_phosphatelt..._bdmax")
max_min_comp("BO2_curvellt..._bdmax") # Problems everywhere
ggsave("graph/tests/curvel_bdmax.png")
# 2050
max_min_comp("BO2_RCP85_2050_curvellt..._bdmax") # Nearly identical to present day curvel
ggsave("graph/tests/curvel_bdmax_2050.png")
max_min_comp("BO2_RCP85_2050_salinitylt..._bdmax") # The largest min values are larger than the largest max values
ggsave("graph/tests/salinity_bdmax_2050.png")
max_min_comp("BO2_RCP85_2050_salinitylt..._ss")
max_min_comp("BO2_RCP85_2050_templt..._bdmax") # Some issues in Baffin Bay and further north, there is no apparent pattern
ggsave("graph/tests/temp_bdmax_2050.png")
max_min_comp("BO2_RCP85_2050_templt..._ss")
max_min_comp("BO2_RCP85_2050_icethicklt..._ss")
# 2100
max_min_comp("BO2_RCP85_2100_curvellt..._bdmax") # Nearly identical to present day curvel
ggsave("graph/tests/curvel_bdmax_2100.png")
max_min_comp("BO2_RCP85_2100_salinitylt..._bdmax") # Similar problems to 2050 data
ggsave("graph/tests/salinity_bdmax_2100.png")
max_min_comp("BO2_RCP85_2100_salinitylt..._ss")
max_min_comp("BO2_RCP85_2100_templt..._bdmax") # Some issues as 2050 data
ggsave("graph/tests/temp_bdmax_2100.png")
max_min_comp("BO2_RCP85_2100_templt..._ss")
max_min_comp("BO2_RCP85_2100_icethicklt..._ss")

# PAR # Some issues in the far north
Arctic_env %>% 
  na.omit() %>% 
  mutate(max_min = ifelse(BO_parmax >= BO_parmean, TRUE, FALSE)) %>% 
  ggplot(aes(x = lon, y = lat)) +
  geom_raster(aes(fill = max_min)) +
  coord_quickmap(expand = F) +
  labs(fill = "Max greater than mean", x = NULL, y = NULL,
       title = "PAR") +
  theme(legend.position = "bottom")


# Global look at problem layers from Arctic test --------------------------

# Download problematic future layers
BO_layers_future_dl <- get_future_layers(c("BO2_templtmin_bdmax", "BO2_templtmax_bdmax", 
                                           "BO2_salinityltmin_bdmax", "BO2_salinityltmax_bdmax", 
                                           "BO2_curvelltmin_bdmax", "BO2_curvelltmax_bdmax"), 
                                         scenario = "RCP85", year = c(2050, 2100))
BO_layers_future_dl <- load_layers(BO_layers_future_dl$layer_code)

# Convert to dataframe
BO_layers_future_df <- as.data.frame(BO_layers_future_dl, xy = T) %>% 
  dplyr::rename(lon = x, lat = y) %>% 
  mutate(lon = round(lon, 4), 
         lat = round(lat, 4)) %>% 
  na.omit() %>% 
  mutate(max_min_curvel_2050 = ifelse(BO2_RCP85_2050_curvelltmax_bdmax >= BO2_RCP85_2050_curvelltmin_bdmax, TRUE, FALSE),
         max_min_curvel_2100 = ifelse(BO2_RCP85_2100_curvelltmax_bdmax >= BO2_RCP85_2100_curvelltmin_bdmax, TRUE, FALSE),
         max_min_salinity_2050 = ifelse(BO2_RCP85_2050_salinityltmax_bdmax >= BO2_RCP85_2050_salinityltmin_bdmax, TRUE, FALSE),
         max_min_salinity_2100 = ifelse(BO2_RCP85_2100_salinityltmax_bdmax >= BO2_RCP85_2100_salinityltmin_bdmax, TRUE, FALSE),
         max_min_temp_2050 = ifelse(BO2_RCP85_2050_templtmax_bdmax >= BO2_RCP85_2050_templtmin_bdmax, TRUE, FALSE),
         max_min_temp_2100 = ifelse(BO2_RCP85_2100_templtmax_bdmax >= BO2_RCP85_2100_templtmin_bdmax, TRUE, FALSE))

# Convenience function for plots
test_global_fig <- function(chosen_col, chosen_title){
  ggplot(data = BO_layers_future_df, aes(x = lon, y = lat)) +
    geom_raster(aes_string(fill = chosen_col)) +
    coord_quickmap(expand = F) +
    labs(fill = "Max greater than min", x = NULL, y = NULL,
         title = chosen_title) +
    theme(legend.position = "bottom")
}

# Curvel 2050
curvel_2050_global <- test_global_fig("max_min_curvel_2050", "Bottom currents: 2050; RCP8.5")
ggsave(plot = curvel_2050_global, filename = "tests/curvel_2050_global.png", height = 5, width = 8)

# Curvel 2100
curvel_2100_global <- test_global_fig("max_min_curvel_2100", "Bottom currents: 2100; RCP8.5")
ggsave(plot = curvel_2100_global, filename = "tests/curvel_2100_global.png", height = 5, width = 8)

# salinity 2050
salinity_2050_global <- test_global_fig("max_min_salinity_2050", "Bottom salinity: 2050; RCP8.5")
ggsave(plot = salinity_2050_global, filename = "tests/salinity_2050_global.png", height = 5, width = 8)

# salinity 2100
salinity_2100_global <- test_global_fig("max_min_salinity_2100", "Bottom salinity: 2100; RCP8.5")
ggsave(plot = salinity_2100_global, filename = "tests/salinity_2100_global.png", height = 5, width = 8)

# temp 2050
temp_2050_global <- test_global_fig("max_min_temp_2050", "Bottom temperature: 2050; RCP8.5")
ggsave(plot = temp_2050_global, filename = "tests/temp_2050_global.png", height = 5, width = 8)

# temp 2100
temp_2100_global <- test_global_fig("max_min_temp_2100", "Bottom temperature: 2100; RCP8.5")
ggsave(plot = temp_2100_global, filename = "tests/temp_2100_global.png", height = 5, width = 8)

