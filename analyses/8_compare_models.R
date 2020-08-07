# 8_compare_models.R
# The purpose of this script is to load results from Jesi's
# MAXENT model and compare them against the RF model


# Setup -------------------------------------------------------------------

# Load all data nad previous libraries
source("analyses/7_random_forests.R")
library(ggridges) # For ridgeplots
# library(sp) # For reading ASCII files
library(raster)

# Load Arctic data for land distance and bathy only
load("data/Arctic_AM.RData")
colnames(Arctic_AM)[4] <- "depth"

# The MAXENT lon/lat grid
# NB: Not currently necessary, and too large to save and push to GitHub; ~117 MB
# MAX_grid <- as.data.frame(sp::read.asciigrid(file_name), xy = T) %>% 
#   `colnames<-`(c("suitability", "lon", "lat")) %>%
#   dplyr::select(lon, lat)
# save(MAX_grid, file = "metadata/Max_grid.RData")
# load("metadata/Max_grid.RData")


# Predict coverage --------------------------------------------------------

# ALl kelp
pred_kelpcover <- Arctic_cover_predict(best_rf_kelpcover$choice_reg, base)

# Laminariales
pred_laminariales <- Arctic_cover_predict(best_rf_laminariales$choice_reg, base)

# Agarum
pred_agarum <- Arctic_cover_predict(best_rf_agarum$choice_reg, base) #%>% 
  # dplyr::rename(pred_base = pred_val)
# pred_agarum_2050 <- Arctic_cover_predict(best_rf_agarum$choice_reg, future_2050) %>% 
#   dplyr::rename(pred_2050 = pred_val)
# pred_agarum_2100 <- Arctic_cover_predict(best_rf_agarum$choice_reg, future_2100) %>% 
#   dplyr::rename(pred_2100 = pred_val)
# pred_agarum_ALL <- left_join(pred_agarum, pred_agarum_2050) %>% 
#   left_join(pred_agarum_2100)

# Alaria
pred_alaria <- Arctic_cover_predict(best_rf_alaria$choice_reg, base)


# Load MAXENT data --------------------------------------------------------

# Function for loading MAXENT .tif files as data.frames in the study area
load_MAX_sub <- function(file_name, binary = F){
  if(binary == F) MAX_df <- as.data.frame(sp::read.asciigrid(file_name), xy = T)
  # if(binary == T) MAX_raw <- tiff::readTIFF(file_name, native = T)
  if(binary == T){
    MAX_df <- as.data.frame(raster::raster(file_name), xy = T) %>% 
      cbind(., MAX_grid)
    
  } 
  # MAX_df <- as.data.frame(MAX_raw, xy = T)
  MAX_Arctic <- MAX_df %>%
    `colnames<-`(c("suitability", "lon", "lat")) %>%
    filter(lon >= bbox_arctic[1], lon <= bbox_arctic[2],
           lat >= bbox_arctic[3], lat <= bbox_arctic[4]) %>% 
    mutate(lon = round(lon, 5),
           lat = round(lat, 5))
  return(MAX_Arctic)
}

## Load continuous values
# Laminariales
MAX_Ldig <- load_MAX_sub("data/Ldig_avg.asc")
MAX_Lsol <- load_MAX_sub("data/Lsol_avg.asc")
MAX_Slat <- load_MAX_sub("data/Slat_avg.asc")

# Alaria
MAX_Alaria <- load_MAX_sub("data/Aesc_avg.asc")

# Agarum
MAX_Agarum <- load_MAX_sub("data/Acla_avg.asc")

## Load binary values
# Laminariales
MAX_Ldig <- load_MAX_sub("data/Ldig_binary.tif", binary = T)
MAX_Lsol <- load_MAX_sub("data/Lsol_binary.tif", binary = T)
MAX_Slat <- load_MAX_sub("data/Slat_binary.tif", binary = T)

# Alaria
MAX_Alaria <- load_MAX_sub("data/Aesc_binary.tif", binary = T)

# Agarum
MAX_Agarum <- load_MAX_sub("data/Acla_binary.tif", binary = T)

# Merge models ------------------------------------------------------------

# Used for the mapping figures
ALL_Agarum <- right_join(pred_agarum_ALL, MAX_Agarum_Arctic, by = c("lon", "lat")) %>% 
  left_join(MAX_Agarum_future_Arctic, by = c("lon", "lat")) %>% 
  dplyr::rename(suitability = suitability.x, suitability_2050 = suitability.y) %>% 
  left_join(Arctic_AM, by = c("lon", "lat", "depth", "land_distance")) %>% 
  mutate(pred_base_perc = pred_base/100,
         pred_2050_perc = pred_2050/100,
         pred_2100_perc = pred_2100/100,
         diff_model_base = pred_base_perc - suitability,
         diff_model_2050 = pred_2050_perc - suitability_2050,
         pred_diff_2050 = pred_2050 - pred_base,
         pred_diff_2100 = pred_2100 - pred_base,
         both_high =  case_when(
           pred_base_perc >= 0.25 & suitability >= 0.5 ~ 1),
         pred_cat = case_when(
           pred_base_perc >= 0.5 ~ "High",
           pred_base_perc < 0.5 & pred_base >= 0.1 ~ "Low",
           pred_base_perc <= 0.1 ~ "Scarce"),
         suit_cat = cut(suitability, breaks = c(0, 30, 50, 100), ordered_result = T))
           # pred_base <= 0.1 & pred_base > 0 ~ "Scarce"))#,
           # pred_base == 0 ~ "None"))

# Used for most of the other figures
ALL_Agarum_long <- ALL_Agarum %>% 
  dplyr::select(-both_high, -pred_cat, -suit_cat) %>% 
  pivot_longer(cols = pred_base:pred_diff_2100)

# ALL_Agarum_cut <- left_join(pred_agarum_cut, Arctic_MAX_Agarum, by = c("lon", "lat")) %>%
  # mutate(suit_cat = cut(suitability, breaks = c(0, 30, 50, 100), ordered_result = T))

# Create a third value that boths outputs can be brought closer to

# Create a yes/no threshold for both

# May be good to create four categories for MAXENT data, too


# Map figure --------------------------------------------------------------

# Data points for map
data_point <- adf %>% 
  dplyr::select(Campaign, site, depth, -c(Bedrock..:sand), kelp.cover, Laminariales, Agarum, Alaria) %>% 
  left_join(study_site_env, by = c("Campaign", "site")) %>%
  mutate(kelp.cover = ifelse(kelp.cover > 100, 100, kelp.cover)) %>%
  dplyr::select(lon, lat, Agarum) %>% 
  group_by(lon, lat) %>%
  summarise_all(mean) %>%
  ungroup()

# Projected current Agarum percent cover
ggplot(filter(ALL_Agarum, land_distance <= 50 | depth <= 50), aes(x = lon, y = lat)) +
  geom_tile(aes(fill = pred_base)) +
  borders(fill = "grey70", colour = "black") +
  geom_point(data = data_point, colour = "red", shape = 21,
             aes(x = lon, y = lat, size = Agarum)) +
  scale_fill_viridis_c(option = "D") +
  coord_cartesian(xlim = c(bbox_arctic[1], bbox_arctic[2]),
                  ylim = c(bbox_arctic[3]-5, bbox_arctic[4]),
                  expand = F) +
  theme(legend.position = "bottom") +
  labs(x = NULL, y = NULL, fill = "Agarum cover")
# ggsave("graph/map_agarum_base.png")

# Difference in 2050 projected Agarum cover
ggplot(filter(ALL_Agarum, land_distance <= 50 | depth <= 50), aes(x = lon, y = lat)) +
  geom_tile(aes(fill = pred_diff_2050)) +
  borders(fill = "grey70", colour = "black") +
  scale_fill_gradient2(low = "blue", high = "red") +
  coord_cartesian(xlim = c(bbox_arctic[1], bbox_arctic[2]),
                  ylim = c(bbox_arctic[3]-5, bbox_arctic[4]),
                  expand = F) +
  theme(legend.position = "bottom") +
  labs(x = NULL, y = NULL, fill = "Agarum cover")
# ggsave("graph/map_agarum_2050.png")

# Difference in 2100 projected Agarum cover
ggplot(filter(ALL_Agarum, land_distance <= 50 | depth <= 50), aes(x = lon, y = lat)) +
  geom_tile(aes(fill = pred_diff_2100)) +
  borders(fill = "grey70", colour = "black") +
  scale_fill_gradient2(low = "blue", high = "red") +
  coord_cartesian(xlim = c(bbox_arctic[1], bbox_arctic[2]),
                  ylim = c(bbox_arctic[3]-5, bbox_arctic[4]),
                  expand = F) +
  theme(legend.position = "bottom") +
  labs(x = NULL, y = NULL, fill = "Agarum cover")
# ggsave("graph/map_agarum_2100.png")

# Difference between model projection
ggplot(filter(ALL_Agarum, land_distance <= 50 | depth <= 50), aes(x = lon, y = lat)) +
  geom_tile(aes(fill = diff_model_base)) +
  borders(fill = "grey70", colour = "black") +
  geom_point(data = data_point, colour = "red", shape = 21,
             aes(x = lon, y = lat, size = Agarum)) +
  scale_fill_gradient2(low = "blue", high = "red") +
  coord_cartesian(xlim = c(bbox_arctic[1], bbox_arctic[2]),
                  ylim = c(bbox_arctic[3]-5, bbox_arctic[4]),
                  expand = F) +
  theme(legend.position = "bottom") +
  labs(x = NULL, y = NULL, fill = "Agarum cover")
# ggsave("graph/map_agarum_model_base.png")

# Difference between model projection for 2050
ggplot(filter(ALL_Agarum, land_distance <= 50 | depth <= 50), aes(x = lon, y = lat)) +
  geom_tile(aes(fill = diff_model_2050)) +
  borders(fill = "grey70", colour = "black") +
  geom_point(data = data_point, colour = "red", shape = 21,
             aes(x = lon, y = lat, size = Agarum)) +
  scale_fill_gradient2(low = "blue", high = "red") +
  coord_cartesian(xlim = c(bbox_arctic[1], bbox_arctic[2]),
                  ylim = c(bbox_arctic[3]-5, bbox_arctic[4]),
                  expand = F) +
  theme(legend.position = "bottom") +
  labs(x = NULL, y = NULL, fill = "Agarum cover")
# ggsave("graph/map_agarum_model_2050.png")


# Difference between model projection for 2100
## Still need 2100 projections from MAXENT


# A nice map figure that shows regions of agreement between models
# Make a figure like this for the present, future, and both

# It may be good to accompany these maps with density plots showing how the pixels relate to each other


# Barplot figure ----------------------------------------------------------

mean_agarum <- MAX_Agarum_ALL_Arctic %>% 
  group_by(model_run) %>% 
  summarise_all(mean)

ggplot(data = mean_agarum, aes(x = model_run, y = suitability)) +
  geom_bar(stat = "identity", aes(fill = model_run))


# Boxplot figure ----------------------------------------------------------

ggplot(data = MAX_Agarum_ALL_Arctic, aes(x = model_run, y = suitability)) +
  geom_boxplot(aes(fill = model_run))

ALL_Agarum_long %>% 
  dplyr::filter(name %in% c("pred_base_perc", "pred_2050_perc", 
                            "suitability", "suitability_2050")) %>% 
  ggplot(aes(x = name, y = value)) +
  geom_boxplot(aes(fill = name))


# Density plot ------------------------------------------------------------

ggplot(data = MAX_Agarum_ALL_Arctic, aes(x = suitability)) +
  geom_density(aes(fill = model_run), alpha = 0.5)
# ggsave("graph/draft_density.png")

ALL_Agarum_long %>% 
  dplyr::filter(name %in% c("pred_base_perc", "pred_2050_perc", 
                            "suitability", "suitability_2050")) %>% 
  ggplot(aes(x = value)) +
  geom_density(aes(fill = name), alpha = 0.5)


# Ridgeplot ---------------------------------------------------------------

ggplot(data = MAX_Agarum_ALL_Arctic, aes(x = suitability, y = model_run)) +
  stat_density_ridges(alpha = 0.7)
# ggsave("graph/draft_ridge.png")

ALL_Agarum_long %>% 
  dplyr::filter(name %in% c("pred_base_perc", "pred_2050_perc", 
                            "suitability", "suitability_2050")) %>% 
  ggplot(aes(x = value, y = name)) +
  stat_density_ridges(alpha = 0.7)

