# tests/test_BO.R
# The purpose of this script is to test the BioOracle data

# The study sites and bounding box
source("analyses/1_study_sites.R")

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
current_bdmax_global <- ggplot(data = current_bd_test, aes(x = lon, y = lat)) +
  geom_raster(aes(fill = max_min)) +
  coord_quickmap(expand = F) +
  labs(fill = "Max greater than min", x = NULL, y = NULL) +
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
  labs(fill = "Max greater than min", x = NULL, y = NULL) +
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
  labs(fill = "Max greater than min", x = NULL, y = NULL) +
  theme(legend.position = "bottom")
ggsave(plot = curvel_bdmax_new , filename = "tests/current_bdmax_global_new.png", height = 5, width = 8)


# Test for issues in the BO layers ----------------------------------------

# Test surface currents
current_ss_layers <- load_layers(c("BO2_curvelltmin_ss", "BO2_curvelmean_ss", "BO2_curvelltmax_ss"))
current_ss_test <- as.data.frame(current_ss_layers, xy = T) %>% 
  dplyr::rename(lon = x, lat = y) %>% 
  mutate(lon = round(lon, 4), 
         lat = round(lat, 4)) %>% 
  na.omit() %>% 
  mutate(max_min = ifelse(BO2_curvelltmax_ss > BO2_curvelltmin_ss, TRUE, FALSE),
         mean_min = ifelse(BO2_curvelmean_ss > BO2_curvelltmin_ss, TRUE, FALSE),
         max_min_int = as.integer(max_min))
current_ss_pixel_test <- ggplot(data = current_ss_test, aes(x = lon, y = lat)) +
  geom_raster(aes(fill = max_min)) +
  coord_quickmap(expand = F) +
  labs(fill = "Max greater than min", x = NULL, y = NULL) +
  theme(legend.position = "bottom")
ggsave(plot = current_ss_pixel_test, filename = "tests/current_ss_pixel_test.png", height = 5, width = 8)

# Test SST
SST_layers <- load_layers(c("BO2_templtmin_ss", "BO2_tempmean_ss", "BO2_templtmax_ss"))
SST_test <- as.data.frame(SST_layers, xy = T) %>% 
  dplyr::rename(lon = x, lat = y) %>% 
  mutate(lon = round(lon, 4), 
         lat = round(lat, 4)) %>% 
  na.omit() %>% 
  mutate(max_min = ifelse(BO2_templtmax_ss > BO2_templtmin_ss, TRUE, FALSE),
         mean_min = ifelse( BO2_tempmean_ss > BO2_templtmin_ss, TRUE, FALSE),
         max_min_int = as.integer(max_min))

# Visualise pixels where the max and min values are not as expected
SST_pixel_test <- ggplot(data = SST_test, aes(x = lon, y = lat)) +
  geom_raster(aes(fill = max_min)) +
  coord_quickmap(expand = F) +
  labs(fill = "Max greater than min", x = NULL, y = NULL) +
  theme(legend.position = "bottom")
ggsave(plot = SST_pixel_test, filename = "graph/SST_pixel_test.png", height = 5, width = 8)

# Compare data layer downloaded via R against a layer download manually
SST_mean_manual <- as.data.frame(read.asciigrid("data/Present.Surface.Temperature.Mean.asc"), xy = T) %>% 
  `colnames<-`(c("SST_mean_manual", "lon", "lat")) %>% 
  mutate(lon = round(lon, 4), 
         lat = round(lat, 4)) %>% 
  na.omit() 
SST_test <- left_join(SST_test, SST_mean_manual, by = c("lon", "lat")) %>% 
  mutate(mean_comp = BO2_tempmean_ss-SST_mean_manual)

# Plot the difference
SST_mean_diff <- ggplot(SST_test, aes(x = lon, y = lat)) +
  geom_raster(aes(fill = mean_comp)) +
  coord_quickmap(expand = F) +
  labs(fill = "R download minus \nmanual download (C)", x = NULL, y = NULL) +
  theme(legend.position = "bottom", legend.key.width = unit(3, "cm"))
ggsave(plot = SST_mean_diff, filename = "graph/SST_pixel_mean_diff.png", height = 5, width = 8)

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
max_min_comp("BO2_RCP85_2100_icethicklt..._ss") # Much of Hudson Bay and Labrador Sea are wrong

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