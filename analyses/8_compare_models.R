# 8_compare_models.R
# The purpose of this script is to load results from Jesi's
# MAXENT model and compare them against the RF model


# Setup -------------------------------------------------------------------

# Load all data nad previous libraries
source("analyses/7_random_forests.R")
library(ggridges) # For ridgeplots
# library(sp) # For reading ASCII files
library(raster)
library(doParallel); registerDoParallel(cores = 10)

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

# Load kelp cover projections
load("data/best_rf_kelpcover.RData")
load("data/best_rf_laminariales.RData")
load("data/best_rf_agarum.RData")
load("data/best_rf_alaria.RData")

# Project the best fits
pred_kelpcover <- Arctic_cover_predict(best_rf_kelpcover$choice_reg, base)
pred_laminariales <- Arctic_cover_predict(best_rf_laminariales$choice_reg, base)
pred_agarum <- Arctic_cover_predict(best_rf_agarum$choice_reg, base)
pred_alaria <- Arctic_cover_predict(best_rf_alaria$choice_reg, base)


# Load MAXENT data --------------------------------------------------------

# Function for loading MAXENT .tif files as data.frames in the study area
load_MAX_sub <- function(file_name){
  MAX_raw <- sp::read.asciigrid(file_name)
  MAX_df <- as.data.frame(MAX_raw, xy = T)
  # MAX_df <- as.data.frame(sp::read.asciigrid(file_name), xy = T)
  MAX_Arctic <- MAX_df %>%
    `colnames<-`(c("suitability", "lon", "lat")) %>%
    filter(lon >= bbox_arctic[1], lon <= bbox_arctic[2],
           lat >= bbox_arctic[3], lat <= bbox_arctic[4]) %>% 
    mutate(lon = round(lon, 5),
           lat = round(lat, 5)) %>% 
    dplyr::select(lon, lat, suitability)
  return(MAX_Arctic)
}

## Load continuous values
# Laminariales
MAX_Laminariales_cont <- plyr::ldply(c("data/Ldig_avg.asc", "data/Lsol_avg.asc", "data/Slat_avg.asc"),
                                load_MAX_sub, .parallel = T) %>% 
  plyr::ddply(.data = ., .variables = c("lon", "lat"), .parallel = T,
              .fun = plyr::summarise, suitability = mean(suitability))

# Alaria
MAX_Alaria_cont <- load_MAX_sub("data/Aesc_avg.asc")

# Agarum
MAX_Agarum_cont <- load_MAX_sub("data/Acla_avg.asc")

## Create binary values
# Laminariales; Lsol - 0.441; Ldig - 0.222; Slat 0.199
MAX_Laminariales_binary <- plyr::ldply(c("data/Ldig_avg.asc", "data/Lsol_avg.asc", "data/Slat_avg.asc"),
                                load_MAX_sub, .parallel = T) %>% 
  plyr::ddply(.data = ., .variables = c("lon", "lat"), .parallel = T,
              .fun = plyr::summarise, suitability = mean(suitability))

# Alaria; Aesc - 0.233
MAX_Alaria_binary <- MAX_Alaria_cont %>% 
  mutate(presence = case_when())

# Agarum; Acla - 0.156
MAX_Agarum_binary <- MAX_Agarum_cont


# Merge models ------------------------------------------------------------

# Used for the mapping figures
merge_pred_MAX <- function(pred_layer, MAX_layer){
  ALL_Agarum <- left_join(pred_layer, MAX_layer, by = c("lon", "lat")) %>% 
    mutate(pred_val_perc = pred_val/100,
           pred_relative = pred_val/max(pred_val, na.rm = T),
           diff_model_base = pred_val_perc - suitability,
           diff_model_base_relative = pred_relative - suitability,
           diff_relative_cat = case_when(
             diff_model_base_relative >= -0.2 & diff_model_base_relative <= 0.2 ~ "Similar",
             diff_model_base_relative > 0.2 ~ "RF",
             diff_model_base_relative < -0.2 ~ "MAXENT"),
           both_high =  case_when(
             pred_relative >= 0.5 & suitability >= 0.5 ~ 1,
             TRUE ~ 0),
           pred_cat = case_when(
             pred_val_perc >= 0.5 ~ "High",
             pred_val_perc < 0.5 & pred_val >= 0.1 ~ "Low",
             pred_val_perc <= 0.1 ~ "Scarce"),
           suit_cat = cut(suitability, breaks = c(0, 0.30, 0.50, 1.00), ordered_result = T))
}

# Merge layers
ALL_laminariales <- merge_pred_MAX(pred_laminariales, MAX_Laminariales_cont)
ALL_agarum <- merge_pred_MAX(pred_agarum, MAX_Agarum_cont)
ALL_alaria <- merge_pred_MAX(pred_alaria, MAX_Alaria_cont)

# Convenience wrapper to create long dataframes
long_kelp_df <- function(df){
  ALL_res_long <- df %>% 
    dplyr::select(-both_high, -pred_cat, -suit_cat) %>% 
    pivot_longer(cols = pred_val:diff_model_base_relative)
}

# ALL_Agarum_cut <- left_join(pred_agarum_cut, Arctic_MAX_Agarum, by = c("lon", "lat")) %>%
  # mutate(suit_cat = cut(suitability, breaks = c(0, 30, 50, 100), ordered_result = T))

# Create a third value that both outputs can be brought closer to

# Create a yes/no threshold for both

# May be good to create four categories for MAXENT data, too


# Map figure --------------------------------------------------------------

# Data points for map
data_point <- adf %>% 
  dplyr::select(Campaign, site, depth, -c(Bedrock..:sand), kelp.cover, Laminariales, Agarum, Alaria) %>% 
  left_join(study_site_env, by = c("Campaign", "site")) %>%
  mutate(kelp.cover = ifelse(kelp.cover > 100, 100, kelp.cover)) %>%
  dplyr::select(lon, lat, Agarum, Alaria, Laminariales) %>%
  group_by(lon, lat) %>%
  summarise_all(mean) %>%
  ungroup()

# Projected current percent cover figure
project_cover_fig <- function(df, kelp_choice){
  cover_fig <- ggplot(filter(df, land_distance <= 50 | depth <= 50), aes(x = lon, y = lat)) +
    geom_tile(aes(fill = pred_val)) +
    borders(fill = "grey70", colour = "black") +
    geom_point(data = data_point, colour = "red", shape = 21,
               aes_string(x = "lon", y = "lat", size = {{kelp_choice}})) +
    scale_fill_viridis_c(option = "D") +
    coord_quickmap(xlim = bbox_arctic[1:2], ylim = bbox_arctic[3:4], expand = F) +
    # theme(legend.position = "bottom") +
    labs(x = NULL, y = NULL, fill = paste0(kelp_choice," cover"))
  ggsave(paste0("graph/map_",tolower(kelp_choice),"_base.png"), cover_fig)
}

# Create the figures
project_cover_fig(ALL_laminariales, "Laminariales")
project_cover_fig(ALL_agarum, "Agarum")
project_cover_fig(ALL_alaria, "Alaria")

# Difference between model projection
project_diff_fig <- function(df, kelp_choice){
  diff_fig <- ggplot(filter(df, land_distance <= 50 | depth <= 50), aes(x = lon, y = lat)) +
    # geom_tile(aes(fill = diff_model_base_relative)) +
    geom_tile(aes(fill = diff_relative_cat)) +
    borders(fill = "white", colour = "grey70", size = 0.1) +
    geom_point(data = data_point, colour = "black", shape = 21,
               aes_string(x = "lon", y = "lat", size = {{kelp_choice}})) +
    # scale_fill_gradient2(low = "blue", mid = "grey30", high = "red", limits = c(-1, 1)) +
    # scale_fill_brewer(palette = "Dark2") +
    scale_fill_manual(values = RColorBrewer::brewer.pal(n = 6, name = 'YlOrRd')[c(2,4,6)]) +
    coord_quickmap(xlim = bbox_arctic[1:2], ylim = bbox_arctic[3:4], expand = F) +
    # theme(legend.position = "bottom") +
    labs(x = NULL, y = NULL, size = paste0(kelp_choice," %"), 
         fill = paste0(kelp_choice," cover\n RF - MAXENT"))
  ggsave(paste0("graph/map_",tolower(kelp_choice),"_model_diff.png"), diff_fig)
}

# Create the figures
project_diff_fig(ALL_laminariales, "Laminariales")
project_diff_fig(ALL_agarum, "Agarum")
project_diff_fig(ALL_alaria, "Alaria")

# Combine this plot with the ridgeplot below

# A nice map figure that shows regions of agreement between models

# It may be good to accompany these maps with density plots showing how the pixels relate to each other


# Barplot figure ----------------------------------------------------------

# mean_agarum <- MAX_Agarum_ALL_Arctic %>% 
#   group_by(model_run) %>% 
#   summarise_all(mean)
# 
# ggplot(data = mean_agarum, aes(x = model_run, y = suitability)) +
#   geom_bar(stat = "identity", aes(fill = model_run))


# Boxplot figure ----------------------------------------------------------

ALL_agarum %>% 
  dplyr::select(-both_high, -pred_cat, -suit_cat, -pred_val) %>% 
  pivot_longer(cols = suitability:diff_model_base_relative) %>% 
  ggplot(aes(x = name, y = value)) +
  geom_boxplot(aes(fill = name))


# Density plot ------------------------------------------------------------

ALL_agarum %>% 
  dplyr::select(-both_high, -pred_cat, -suit_cat, -pred_val) %>% 
  pivot_longer(cols = suitability:diff_model_base_relative) %>% 
  ggplot(aes(x = value)) +
  geom_density(aes(fill = name), alpha = 0.5)


# Ridgeplot ---------------------------------------------------------------

ALL_agarum %>% 
  dplyr::select(-both_high, -pred_cat, -suit_cat, -pred_val, -pred_val_perc, -diff_model_base) %>% 
  pivot_longer(cols = suitability:diff_model_base_relative) %>% 
  ggplot(aes(x = value, y = name)) +
  stat_density_ridges(alpha = 0.7)
# Combine this with the diff map
# Put everything in panels


# Create categories for the kelps that are all relative to their own percent covers
# Use binary models in combination with RF categories to see how well those match
  # Create a RF figure showing the categories, low-high, and panel that next to the binary MAXENT map
# Another comparison to look at is in areas where one model may have much more optimistic projections than the other
# MAXENT can be verified again by looking at where it predicts suitability and where the field data show high percentage covers




