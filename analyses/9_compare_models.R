# 8_compare_models.R
# The purpose of this script is to load results from Jesi's
# MAXENT model and compare them against the RF model


# Setup -------------------------------------------------------------------

# Load all data and previous libraries
source("analyses/7_random_forests.R")
library(ggridges) # For ridgeplots
library(raster)
library(doParallel); registerDoParallel(cores = 50)
# library(sp) # For reading ASCII files

# Load Arctic data for land distance and bathy only
load("data/Arctic_AM.RData")
colnames(Arctic_AM)[4] <- "depth"

# Data points for maps
RF_points <- adf %>% 
  dplyr::select(Campaign, site, depth, -c(Bedrock..:sand), kelp.cover, Laminariales, Agarum, Alaria) %>% 
  left_join(study_site_env, by = c("Campaign", "site")) %>%
  mutate(kelp.cover = ifelse(kelp.cover > 100, 100, kelp.cover)) %>%
  dplyr::select(lon, lat, Agarum, Alaria, Laminariales) %>%
  group_by(lon, lat) %>%
  summarise_all(mean) %>%
  ungroup()

# The coords for the MAXENT data
MAXENT_points <- rbind(
  read_csv("metadata/Acla_MaxEnt.csv")[1:4],
  read_csv("metadata/Aesc_MaxEnt.csv")[1:4],
  read_csv("metadata/Slat_MaxEnt.csv")[1:4],
  read_csv("metadata/Ldig_MaxEnt.csv")[1:4],
  read_csv("metadata/Lsol_MaxEnt.csv")[1:4])
MAXENT_points <- MAXENT_points %>% 
  mutate(family = case_when(Sp == "Aesc" ~ "Alaria",
                            Sp == "Acla" ~ "Agarum",
                            Sp %in% c("Slat", "Ldig", "Lsol") ~ "Laminariales")) %>% 
  filter(lon >= bbox_arctic[1], lon <= bbox_arctic[2],
         lat >= bbox_arctic[3], lat <= bbox_arctic[4])

# The MAXENT lon/lat grid
# NB: Not currently necessary, and too large to save and push to GitHub; ~117 MB
# MAX_grid <- as.data.frame(sp::read.asciigrid(file_name), xy = T) %>% 
#   `colnames<-`(c("suitability", "lon", "lat")) %>%
#   dplyr::select(lon, lat)
# save(MAX_grid, file = "metadata/Max_grid.RData")
# load("metadata/Max_grid.RData")

# MAXENT can be verified again by looking at where it predicts suitability and where the field data show high percentage covers

# Consider depths of 60 to 80 metres as a filter after the models

# Run a presence absence random forest to compare to MAXENT. And also the percent cover RF afterwards when the comparisons of the models has been made.

# Build a story around model distribution, and as point B we can talk about percent cover projections.

# Consider knocking out variables and seeing if some regions agreement change based on that one missing variable.


# Predict coverage --------------------------------------------------------

# Load kelp cover projections
load("data/best_rf_kelpcover.RData")
load("data/best_rf_laminariales.RData")
load("data/best_rf_agarum.RData")
load("data/best_rf_alaria.RData")


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

## Load continuous values and create binaries
# Laminariales; Ldig - 0.222; , Lsol - 0.441; Slat 0.199
MAX_Laminariales <- plyr::ldply(c("data/Ldig_avg.asc", "data/Lsol_avg.asc", "data/Slat_avg.asc"),
                                load_MAX_sub, .parallel = T) %>% 
  mutate(species = rep(c("Ldig", "Lsol", "Slat"), each = 101257),
         presence = case_when(species == "Ldig" & suitability >= 0.222 ~ 1,
                              species == "Lsol" & suitability >= 0.441 ~ 1,
                              species == "Slat" & suitability >= 0.199 ~ 1,
                              TRUE ~ 0)) %>% 
  dplyr::select(-species) %>% 
  plyr::ddply(.data = ., .variables = c("lon", "lat"), .parallel = T,
              .fun = plyr::summarise, suitability = mean(suitability), presence = ceiling(mean(presence)))

# Alaria Aesc - 0.233
MAX_Alaria <- load_MAX_sub("data/Aesc_avg.asc" ) %>% 
  mutate(presence = case_when(suitability >= 0.233 ~ 1, TRUE ~ 0))

# Agarum; Acla - 0.156
MAX_Agarum <- load_MAX_sub("data/Acla_avg.asc") %>% 
  mutate(presence = case_when(suitability >= 0.156 ~ 1, TRUE ~ 0))


# Merge models ------------------------------------------------------------

# Used for the mapping figures
merge_pred_MAX <- function(pred_layer, MAX_layer){
  max_val <- max(pred_layer$reg_val)
  ALL_res <- left_join(pred_layer, MAX_layer, by = c("lon", "lat")) %>% 
    mutate(reg_perc = reg_val/100,
           reg_relative = reg_val/max_val,
           diff_perc = reg_perc - suitability,
           diff_relative = reg_relative - suitability,
           diff_relative_cat = case_when(
             diff_relative >= -0.2 & diff_relative <= 0.2 ~ "Similar",
             diff_relative > 0.2 ~ "RF",
             diff_relative < -0.2 ~ "MAXENT"),
           diff_cat = case_when(
             presence == 1 & cat_val != "none" ~ "both",
             presence == 0 & cat_val == "none" ~ "none",
             presence == 1 & cat_val == "none" ~ "MAXENT",
             presence == 0 & cat_val != "none" ~ "RF"))
}

# Merge layers
ALL_laminariales <- merge_pred_MAX(best_rf_laminariales$project_multi, MAX_Laminariales)
ALL_agarum <- merge_pred_MAX(best_rf_agarum$project_multi, MAX_Agarum)
ALL_alaria <- merge_pred_MAX(best_rf_alaria$project_multi, MAX_Alaria)


# Percent cover figure ----------------------------------------------------

project_reg_fig <- function(df, kelp_choice){
  cover_fig <- df %>% 
    dplyr::select(lon:land_distance, reg_relative, suitability) %>% 
    pivot_longer(cols = c(reg_relative, suitability)) %>% 
    mutate(name = case_when(name == "reg_relative" ~ "RF: relative",
                            name == "suitability" ~ "MAXENT: suitability")) %>% 
    filter(land_distance <= 100 | depth <= 100) %>% 
    ggplot(aes(x = lon, y = lat)) +
    geom_tile(aes(fill = value)) +
    borders(fill = "grey70", colour = "black") +
    geom_point(data = filter(MAXENT_points, family == kelp_choice), colour = "red", shape = 15, alpha = 0.5) +
    geom_point(data = RF_points, colour = "red", shape = 21,
               aes_string(x = "lon", y = "lat", size = {{kelp_choice}})) +
    scale_fill_viridis_c(option = "D") +
    coord_quickmap(xlim = bbox_arctic[1:2], ylim = bbox_arctic[3:4], expand = F) +
    theme(legend.position = "bottom") +#,
          # legend.box = "vertical") +
    labs(x = NULL, y = NULL, fill = paste0("cover (model)"),
         title = kelp_choice, size = paste0("cover (%)")) +
    facet_wrap(~name)
  # cover_fig
  ggsave(paste0("graph/map_",tolower(kelp_choice),"_reg.png"), cover_fig)
}

# Create the figures
project_reg_fig(ALL_laminariales, "Laminariales")
project_reg_fig(ALL_agarum, "Agarum")
project_reg_fig(ALL_alaria, "Alaria")


# Category cover figure ---------------------------------------------------

# Create categories for the kelps that are all relative to their own percent covers
# Use binary models in combination with RF categories to see how well those match
# Create a RF figure showing the categories, low-high, and panel that next to the binary MAXENT map
project_cat_fig <- function(df, kelp_choice){
  cat_val_fig <- df %>% 
    filter(land_distance <= 100 | depth <= 100) %>% 
    mutate(name = "RF: category") %>% 
    ggplot(aes(x = lon, y = lat)) +
    geom_tile(aes(fill = cat_val)) +
    borders(fill = "grey70", colour = "black") +
    geom_point(data = RF_points, colour = "purple", shape = 21,
               aes_string(x = "lon", y = "lat", size = {{kelp_choice}})) +
    scale_fill_brewer(palette = "Dark2") +
    coord_quickmap(xlim = bbox_arctic[1:2], ylim = bbox_arctic[3:4], expand = F) +
    theme(legend.position = "bottom",
          legend.box = "vertical") +
    labs(x = NULL, y = NULL, fill = paste0("category"),
         title = kelp_choice, size = paste0("cover (%)")) +
      facet_wrap(~name)
  # cat_val_fig
  
  presence_fig <- df %>% 
    filter(land_distance <= 100 | depth <= 100) %>% 
    mutate(name = "MAXENT: presence",
           presence = case_when(presence == 1 ~ "yes",
                                presence == 0 ~ "no")) %>% 
    ggplot(aes(x = lon, y = lat)) +
    geom_tile(aes(fill = presence)) +
    borders(fill = "grey70", colour = "black") +
    geom_point(data = filter(MAXENT_points, family == kelp_choice), colour = "purple", shape = 15, alpha = 0.5) +
    scale_fill_brewer(palette = "Set1") +
    coord_quickmap(xlim = bbox_arctic[1:2], ylim = bbox_arctic[3:4], expand = F) +
    theme(legend.position = "bottom",
          legend.box = "vertical") +
    labs(x = NULL, y = NULL, fill = paste0("presence"), size = paste0("cover (%)")) +
    facet_wrap(~name)
  # presence_fig
  cat_fig <- ggpubr::ggarrange(cat_val_fig, presence_fig, ncol = 2, nrow = 1, 
                               align = "hv", common.legend = F, legend = "bottom")
  ggsave(paste0("graph/map_",tolower(kelp_choice),"_cat.png"), cat_fig, width = 6.5)
}

# Create the figures
project_cat_fig(ALL_laminariales, "Laminariales")
project_cat_fig(ALL_agarum, "Agarum")
project_cat_fig(ALL_alaria, "Alaria")


# Difference regression ---------------------------------------------------

# Change the colour palette to show where they agree, but by how much.
# Also showing that they agree because they are both low or high, etc.
# This may require a figure showing how they are similar at certain levels,
# and a figure showing dissimilarity at similar levels.
project_reg_diff_fig <- function(df, kelp_choice){
  reg_diff_map <- ggplot(filter(df, land_distance <= 100 | depth <= 100), aes(x = lon, y = lat)) +
    geom_tile(aes(fill = diff_relative_cat)) +
    # geom_tile(aes(fill = diff_relative_cat)) +
    borders(fill = "white", colour = "grey70", size = 0.1) +
    geom_point(data = RF_points, colour = "black", shape = 21,
               aes_string(x = "lon", y = "lat", size = {{kelp_choice}})) +
    geom_point(data = filter(MAXENT_points, family == kelp_choice), colour = "black", shape = 15, alpha = 0.5) +
    # scale_fill_gradient2(low = "blue", mid = "grey30", high = "red", limits = c(-1, 1)) +
    # scale_fill_brewer(palette = "Dark2") +
    scale_fill_manual(values = RColorBrewer::brewer.pal(n = 6, name = 'YlOrRd')[c(2,4,6)]) +
    coord_quickmap(xlim = bbox_arctic[1:2], ylim = bbox_arctic[3:4], expand = F) +
    theme(legend.position = "bottom",
          legend.box = "vertical") +
    labs(x = NULL, y = NULL, size = paste0("cover (%)"), 
         fill = paste0("greater"), title = kelp_choice)
  # reg_diff_map
  # reg_diff_fig <- ggpubr::ggarrange(reg_diff_map, reg_diff_ridge, ncol = 2, nrow = 1, 
  #                                   align = "v", widths = c(1.5, 1))
  ggsave(paste0("graph/map_",tolower(kelp_choice),"_reg_diff.png"), reg_diff_map)
}

# Create the figures
project_reg_diff_fig(ALL_laminariales, "Laminariales")
project_reg_diff_fig(ALL_agarum, "Agarum")
project_reg_diff_fig(ALL_alaria, "Alaria")


# Difference category -----------------------------------------------------

# Could show the RF categories of low, medium, high but as increasing shades of the same colour
# And then don't show where they don't agree or both predict none

project_cat_diff_fig <- function(df, kelp_choice){
  cat_diff_map <- df %>% 
    filter(land_distance <= 100 | depth <= 100) %>% 
    mutate(diff_cat = factor(diff_cat, levels = c("both", "MAXENT", "RF", "none"))) %>% 
    ggplot(aes(x = lon, y = lat)) +
    geom_tile(aes(fill = diff_cat)) +
    # geom_tile(aes(fill = diff_relative_cat)) +
    borders(fill = "white", colour = "grey70", size = 0.1) +
    geom_point(data = RF_points, colour = "black", shape = 21,
               aes_string(x = "lon", y = "lat", size = {{kelp_choice}})) +
    geom_point(data = filter(MAXENT_points, family == kelp_choice), colour = "black", shape = 15, alpha = 0.5) +
    scale_fill_brewer(palette = "Accent") +
    coord_quickmap(xlim = bbox_arctic[1:2], ylim = bbox_arctic[3:4], expand = F) +
    theme(legend.position = "bottom",
          legend.box = "vertical") +
    labs(x = NULL, y = NULL, size = paste0("cover (%)"), 
         fill = paste0("presence"), title = kelp_choice)
  # cat_diff_map
  # cat_diff_fig <- ggpubr::ggarrange(cat_diff_map, cat_diff_ridge, ncol = 2, nrow = 1, 
  #                                   align = "v", widths = c(1.5, 1))
  ggsave(paste0("graph/map_",tolower(kelp_choice),"_cat_diff.png"), cat_diff_map)
}

# Create the figures
project_cat_diff_fig(ALL_laminariales, "Laminariales")
project_cat_diff_fig(ALL_agarum, "Agarum")
project_cat_diff_fig(ALL_alaria, "Alaria")

