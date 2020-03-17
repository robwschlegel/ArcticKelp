# 8_compare_models.R
# The purpose of this script is to load results from Jesi's
# MAXENT model and compare them against the RF model


# Setup -------------------------------------------------------------------

# Load all data nad previous libraries
source("analyses/4_kelp_cover.R")
library(ggridges) # For ridgeplots
library(sp) # For reading ASCII files

# Remove scientific notation from data.frame displays in RStudio
options(scipen = 9999)

# The base map to use for everything else
Arctic_map <- ggplot() +
  borders(fill = "grey70", colour = "black") +
  coord_cartesian(expand = F,
                  xlim = c(bbox_arctic[1], bbox_arctic[2]),
                  ylim = c(bbox_arctic[3], bbox_arctic[4])) +
  labs(x = NULL, y = NULL)

# Load top variable choices
load("data/top_var_kelpcover.RData")
load("data/top_var_laminariales.RData")
load("data/top_var_agarum.RData")
load("data/top_var_alaria.RData")

# The best random forest models
load("data/best_rf_kelpcover.RData")
load("data/best_rf_laminariales.RData")
load("data/best_rf_agarum.RData")
load("data/best_rf_alaria.RData")

# Load the Arctic data
load("data/Arctic_env.RData")
load("data/Arctic_AM.RData")


# Predict coverage --------------------------------------------------------

# Convenience function for final step before prediction
Arctic_cover_predict <- function(model_choice){
  pred_df <- data.frame(lon = Arctic_env$lon, lat = Arctic_env$lat,
                        depth = Arctic_env$bathy,
                        pred_val = predict(model_choice, Arctic_env))
}

# Predict the covers
pred_kelpcover <- Arctic_cover_predict(best_rf_kelpcover$choice_model)
pred_laminariales <- Arctic_cover_predict(best_rf_laminariales$choice_model)
pred_agarum <- Arctic_cover_predict(best_rf_agarum$choice_model)
pred_alaria <- Arctic_cover_predict(best_rf_alaria$choice_model)


# Load MAXENT data --------------------------------------------------------

# Agarum base
MAX_Agarum <- read.asciigrid("data/Acla.asc")
MAX_Agarum_df <- as.data.frame(MAX_Agarum, xy = T)
Arctic_MAX_Agarum <- MAX_Agarum_df %>%
  dplyr::rename(suitability = data.Acla.asc,
                lon = s1, lat = s2) %>%
  filter(lon >= bbox_arctic[1], lon <= bbox_arctic[2],
         lat >= bbox_arctic[3]-5, lat <= bbox_arctic[4]) %>% 
  mutate(lon = round(lon, 5),
         lat = round(lat, 5))

# Agarum sub
MAX_Agarum_sub <- read.asciigrid("data/Acla.sub.asc")
MAX_Agarum_sub_df <- as.data.frame(MAX_Agarum_sub, xy = T)
Arctic_MAX_Agarum_sub <- MAX_Agarum_sub_df %>%
  dplyr::rename(suitability = data.Acla.sub.asc,
                lon = s1, lat = s2) %>%
  filter(lon >= bbox_arctic[1], lon <= bbox_arctic[2],
         lat >= bbox_arctic[3]-5, lat <= bbox_arctic[4]) %>% 
  mutate(lon = round(lon, 5),
         lat = round(lat, 5))

# Agarum future
MAX_Agarum_future <- read.asciigrid("data/Acla.sub.fut.asc")
MAX_Agarum_future_df <- as.data.frame(MAX_Agarum_future, xy = T)
Arctic_MAX_Agarum_future <- MAX_Agarum_future_df %>%
  dplyr::rename(suitability = data.Acla.sub.fut.asc,
                lon = s1, lat = s2) %>%
  filter(lon >= bbox_arctic[1], lon <= bbox_arctic[2],
         lat >= bbox_arctic[3]-5, lat <= bbox_arctic[4]) %>% 
  mutate(lon = round(lon, 5),
         lat = round(lat, 5))

Arctic_MAX_Agarum_ALL <- rbind(Arctic_MAX_Agarum, Arctic_MAX_Agarum_sub, Arctic_MAX_Agarum_future) %>% 
  mutate(model_run = factor(c(rep("base", nrow(Arctic_MAX_Agarum)), 
                              rep("sub", nrow(Arctic_MAX_Agarum_sub)),
                              rep("future_50", nrow(Arctic_MAX_Agarum_future))),
                            levels = c("base", "sub", "future_50")))

# Merge models ------------------------------------------------------------

ALL_Agarum <- right_join(pred_agarum, Arctic_MAX_Agarum, by = c("lon", "lat")) %>% 
  left_join(Arctic_AM, by = c("lon", "lat")) %>% 
  mutate(pred_val = pred_val/100,
         diff_model = pred_val - suitability,
         both_high =  case_when(
           pred_val >= 0.25 & suitability >= 0.5 ~ 1),
         pred_cat = case_when(
           pred_val >= 0.5 ~ "High",
           pred_val < 0.5 & pred_val >= 0.1 ~ "Low",
           pred_val < 0.1 & pred_val > 0 ~ "Scarce",
           pred_val == 0 ~ "None"))

# ALL_Agarum_cut <- left_join(pred_agarum_cut, Arctic_MAX_Agarum, by = c("lon", "lat")) %>% 
  # mutate(suit_cat = cut(suitability, breaks = c(0, 30, 50, 100), ordered_result = T))

# Create a third value that boths outputs can be brought closer to

# Create a yes/no threshold for both

# May be good to create four categories for MAXENT data, too


# Map figure --------------------------------------------------------------

# Visualise
ggplot(filter(ALL_Agarum, land_distance <= 50 | depth <= 100), aes(x = lon, y = lat)) +
  # geom_tile(aes(fill = as.factor(pred_cat))) +
  # geom_tile(aes(fill = as.factor(both_high))) +
  # geom_tile(aes(fill = pred_val)) +
  geom_tile(aes(fill = pred_val)) +
  borders(fill = "grey70", colour = "black") +
  geom_point(data = study_site_env, colour = "red") +
  # geom_point(aes(size = pred_val), shape = 21, fill = NA) +
  scale_fill_viridis_c(option = "D") +
  # scale_fill_gradient2() +
  coord_cartesian(xlim = c(bbox_arctic[1], bbox_arctic[2]),
                  ylim = c(bbox_arctic[3]-5, bbox_arctic[4]),
                  expand = F) +
  theme(legend.position = "bottom") +
  labs(x = NULL, y = NULL)
ggsave("graph/draft_map.png")

# A nice map figure that shows regions of agreement between models
# Make a figure like this for the present, future, and both

# It may be good to accompany these maps with density plots showing how the pixels relate to each other


# Barplot figure ----------------------------------------------------------

mean_agarum <- Arctic_MAX_Agarum_ALL %>% 
  group_by(model_run) %>% 
  summarise_all(mean)

ggplot(data = mean_agarum, aes(x = model_run, y = suitability)) +
  geom_bar(stat = "identity", aes( fill = model_run))


# Ridgeplot ---------------------------------------------------------------


ggplot(data = Arctic_MAX_Agarum_ALL, aes(x = suitability)) +
  geom_density(aes(fill = model_run), alpha = 0.5)
ggsave("graph/draft_density.png")

ggplot(data = Arctic_MAX_Agarum_ALL, aes(x = suitability, y = model_run)) +
  stat_density_ridges(alpha = 0.7)
ggsave("graph/draft_ridge.png")
