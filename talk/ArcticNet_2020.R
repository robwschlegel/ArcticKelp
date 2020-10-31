# talk/ArcticNet_2020.R
# The purpose of this script is to house all of the code used to 
# create the figures etc. used in the talk with the same name

# 1: Setup ----------------------------------------------------------------

source("analyses/4_kelp_cover.R")

# Other libraries
library(randomForest)
library(sdmpredictors)

# The base map to use for everything else
Arctic_map <- ggplot() +
  borders(fill = "grey70", colour = "black") +
  scale_y_continuous(breaks = c(60, 70), labels = c("60°N", "70°N")) +
  scale_x_continuous(breaks = c(-80, -60), labels = c("80°W", "60°W")) +
  coord_quickmap(xlim = c(bbox_arctic[1], bbox_arctic[2]),
                 ylim = c(bbox_arctic[3], bbox_arctic[4]), expand = F) +
  labs(x = NULL, y = NULL) +
  theme(panel.border = element_rect(colour = "black", fill = NA))

# Presentation standard figure width/height
fig_width <- 5
fig_height <- 5

# Final BO variables
load("metadata/BO_vars.RData")

# Function for re-loading .RData files as necessary
loadRData <- function(fileName){
  load(fileName)
  get(ls()[ls() != "fileName"])
}

# alaria_trading_points <- data.frame(lon = c(-67.491047),
#                                    lat = c(66.552445))

# Site coordinates
adf_summary_coords <- left_join(adf_summary, study_sites, by = c("Campaign", "site"))

# Mean values per site
adf_summary_mean_coords <- filter(adf_summary_coords, family == "kelp.cover") %>% 
  group_by(lon, lat) %>%
  summarise(mean_cover = mean(mean_cover),
            range_cover = max(mean_cover)-min(mean_cover))

# 2: Arctic region figures ------------------------------------------------

# Overall regions
kelp_question <- Arctic_map +
  geom_label(aes(x = -70, y = 65, label = "?"), size = 60, alpha = 0.6)
kelp_question
ggsave(kelp_question, filename = "talk/figure/kelp_question.png", height = fig_height)

# ArcticKelp campaign map
study_site_campaign <-  Arctic_map +
  geom_point(data = study_sites,
             aes(x = lon, y = lat, colour = Campaign), size = 4) +
  guides(colour = guide_legend(nrow = 2, byrow = TRUE)) +
  labs(colour = "Campaign") +
  theme(strip.background = element_rect(colour = "white", fill = "white"),
        legend.position = c(0.55, 0.96),
        legend.direction = "horizontal",
        legend.spacing.y = unit(0, "mm"))
study_site_campaign
ggsave(study_site_campaign, filename = "talk/figure/study_site_campaign.png", height = fig_height)

# ArcticKelp percent cover map
study_site_mean_cover <-  Arctic_map +
  geom_point(data = adf_summary_mean_coords,
             aes(x = lon, y = lat, colour = mean_cover), size = 4) +
  scale_colour_viridis_c() +
  labs(colour = "Mean cover (%)", shape = "Depth (m)") +
  # guides(color = guide_colourbar(order = 1)) +
  theme(strip.background = element_rect(colour = "white", fill = "white"),
        legend.position = c(0.55, 0.96),
        legend.direction = "horizontal",
        legend.spacing.y = unit(0, "mm"))
study_site_mean_cover
ggsave(study_site_mean_cover, filename = "talk/figure/study_site_mean_cover.png", height = fig_height)


# 3: Model accuracy -------------------------------------------------------

## Ensemble


## RF
# First load the best random forest models produced above
load("data/best_rf_kelpcover.RData")
load("data/best_rf_laminariales.RData")
load("data/best_rf_agarum.RData")
load("data/best_rf_alaria.RData")

# Quick visuals
# Histogram of accuracy of predictions
ggplot(filter(best_rf_kelpcover$model_accuracy, portion == "validate"), aes(x = accuracy)) +
  geom_histogram()
# Scatterplot of accuracy with a linear model overalid
ggplot(filter(best_rf_kelpcover$model_accuracy, portion == "validate"), aes(x = original, y = pred)) +
  geom_point() +
  geom_smooth(method = "lm")
# Boxplots of accuracy
ggplot(filter(best_rf_kelpcover$model_accuracy, portion == "validate"), aes(x = as.factor(original), y = pred)) +
  geom_boxplot()

# Function for creating figure showing confidence intervals of prediction accuracy
conf_plot <- function(df, plot_title){
  
  # 90 CI around predictions per step
  conf_acc <- df %>% 
    filter(portion == "validate") %>% 
    group_by(original) %>% 
    mutate(accuracy = round(accuracy)) %>% 
    summarise(q05 = quantile(accuracy, 0.05),
              q25 = quantile(accuracy, 0.25),
              q50 = median(accuracy),
              mean = mean(accuracy),
              q75 = quantile(accuracy, 0.75),
              q95 = quantile(accuracy, 0.95)) %>% 
    ungroup()
  
  conf_mean <- df %>% 
    filter(portion == "validate") %>% 
    group_by(model_id) %>% 
    mutate(accuracy = round(accuracy)) %>% 
    summarise(mean_acc = mean(abs(accuracy)),
              sd_acc = sd(abs(accuracy)),
              r_acc = cor(x = original, y = pred)) %>% 
    ungroup()
  
  conf_mean_label <- conf_mean %>% 
    summarise(mean_acc = round(mean(abs(mean_acc))),
              sd_acc = round(mean(abs(sd_acc))),
              r_acc = round(mean(r_acc), 2))
  
  conf_best_label <- conf_mean %>% 
    filter(mean_acc == min(mean_acc)) %>% 
    mutate(mean_acc = round(abs(mean_acc)),
           sd_acc = round(sd_acc),
           r_acc = round(r_acc, 2))
  
  conf_best <- df %>% 
    filter(model_id == conf_best_label$model_id,
           portion == "validate") %>% 
    group_by(original) %>% 
    summarise(mean_acc = mean(accuracy),
              sd_acc = sd(accuracy)) %>% 
    ungroup()
  
  # Visualise
  ggplot(conf_acc, aes(x = original, y = mean)) +
    geom_hline(yintercept = 0, size = 2, colour = "red") +
    geom_crossbar(aes(y = 0, ymin = q05, ymax = q95),
                  fatten = 0, fill = "grey70", colour = NA, width = 1) +
    geom_crossbar(aes(ymin = q25, ymax = q75),
                  fatten = 0, fill = "grey50", width = 1) +
    geom_crossbar(aes(ymin = q50, ymax = q50),
                  fatten = 0, fill = NA, colour = "black", width = 1) +
    geom_segment(data = conf_best, aes(xend = original, y = mean_acc, yend = 0), 
                 colour = "purple", size = 1.2, alpha = 0.8) +
    geom_point(data = conf_best, aes(y = mean_acc), colour = "purple", size = 3, alpha = 0.8) +
    geom_label(data = conf_mean_label, 
               aes(x = 75, y = 75, label = paste0("Mean accuracy: ",mean_acc,"±",sd_acc,"; r = ",r_acc))) +
    geom_label(data = conf_best_label, colour = "purple",
               aes(x = 75, y = 60, label = paste0("Best accuracy: ",mean_acc,"±",sd_acc,"; r = ",r_acc))) +
    scale_y_continuous(limits = c(-100, 100)) +
    labs(y = "Range in accuracy of predictions", x = "Original value (% cover)", title = plot_title)
}

# Create the plots
conf_plot(best_rf_kelpcover$model_accuracy, "Total cover confidence")
ggsave("talk/figure/conf_kelpcover.png", height = fig_height, width = fig_width)
conf_plot(best_rf_laminariales$model_accuracy, "Laminariales cover confidence")
ggsave("talk/figure/conf_laminariales.png", height = fig_height, width = fig_width)
conf_plot(best_rf_agarum$model_accuracy, "Agarum cover confidence")
ggsave("talk/figure/conf_agarum.png", height = fig_height, width = fig_width)
conf_plot(best_rf_alaria$model_accuracy, "Alaria cover confidence")
ggsave("talk/figure/conf_alaria.png", height = fig_height, width = fig_width)


# 4: Variable importance --------------------------------------------------

## Ensemble


## RF

# Load top variable choices
load("data/top_var_kelpcover.RData")
load("data/top_var_laminariales.RData")
load("data/top_var_agarum.RData")
load("data/top_var_alaria.RData")

# Need to load for max cover value
load("data/best_rf_kelpcover.RData")
load("data/best_rf_laminariales.RData")
load("data/best_rf_agarum.RData")
load("data/best_rf_alaria.RData")

test <- best_rf_laminariales$accuracy_reg %>% 
  filter(model_id == "1") %>% 
  summarise(mean = round(mean(original)))

# Load long names of variables
model_info <- read_csv("metadata/model_info.csv") %>% 
  dplyr::select(-ndims) %>% 
  dplyr::rename(var = name)

# Load BIO names
BO_layers <- list_layers(datasets = "Bio-ORACLE") %>% 
  dplyr::select(layer_code, units, name) %>% 
  dplyr::rename(longname = name, var = layer_code)

# Combine
all_layers <- rbind(model_info, BO_layers)

# Convenience function to get the cleaned up variables
top_3_bottom <- function(df){
  # mean_cov <- df2$model_accuracy %>% 
  #   filter(model_id == "1") %>% 
  #   summarise(mean = round(mean(original))) %>% 
  #   as.numeric()
  res <- df %>% 
    left_join(all_layers, by = "var") %>% 
    mutate(longname = case_when(var == "depth" ~ "Depth",
                                var == "lon" ~ "Longitude",
                                var == "lat" ~ "Latitude",
                                var == "idive" ~ "Ice divergence",
                                TRUE ~ longname),
           mean_IncMSE = round(sum_IncMSE/1000)) %>% 
    # arrange(-mean_IncMSE) #%>% 
    dplyr::select(longname, mean_IncMSE) %>%
    dplyr::rename(`Data layer` = longname, `% Inc. MSE` = mean_IncMSE)
  top_bottom <- rbind(head(res, 3), tail(res, 3))
}

# Order and make pretty for presentation
top_full_kelpcover <- top_3_bottom(top_var_kelpcover)
save(top_full_kelpcover, file = "talk/data/top_full_kelpcover.RData")
top_full_laminariales <- top_3_bottom(top_var_laminariales)
save(top_full_laminariales, file = "talk/data/top_full_laminariales.RData")
top_full_agarum <- top_3_bottom(top_var_agarum)
save(top_full_agarum, file = "talk/data/top_full_agarum.RData")
top_full_alaria <- top_3_bottom(top_var_alaria)
save(top_full_alaria, file = "talk/data/top_full_alaria.RData")


# 5: Project values -------------------------------------------------------

## Ensemble

## RF

# Load the best random forest models if not already loaded
load("data/best_rf_kelpcover.RData")
load("data/best_rf_laminariales.RData")
load("data/best_rf_agarum.RData")
load("data/best_rf_alaria.RData")

# Load the Arctic data
load("data/Arctic_BO.RData")

# Prep the data for lazy joining
Arctic_mean_prep <- Arctic_mean %>% 
  # dplyr::select(-qla_oce, -qsb_oce) %>%  # These two columns have no values
  na.omit() %>% 
  dplyr::rename(lon = nav_lon, lat = nav_lat) %>% 
  mutate(lon = plyr::round_any(lon, 0.25),
         lat = plyr::round_any(lat, 0.25)) %>% 
  group_by(lon, lat) %>% 
  summarise_if(is.numeric, mean, na.rm = T) %>% 
  ungroup()

Arctic_BO_prep <- na.omit(Arctic_BO) %>% 
  mutate(lon = plyr::round_any(lon, 0.25),
         lat = plyr::round_any(lat, 0.25)) %>% 
  group_by(lon, lat) %>% 
  summarise_if(is.numeric, mean, na.rm = T) %>% 
  ungroup()

# Join the data for being fed to the model
Arctic_data <- left_join(Arctic_BO_prep, Arctic_mean_prep, by = c("lon", "lat")) %>%
  # na.omit() %>% 
  dplyr::select(-depth, -x, -y) %>% 
  dplyr::rename(depth = bathy) %>% 
  ungroup() %>% 
  tidyr::fill(lon:vtau_ice, .direction = "downup")

# Convenience function for final step before prediction
Arctic_cover_predict <- function(model_choice){
  pred_df <- data.frame(lon = Arctic_data$lon, lat = Arctic_data$lat,
                        depth = Arctic_data$depth,
                        pred_val = predict(model_choice, Arctic_data))
}

# Visualise a family of cover
cover_squiz <- function(df, legend_title, x_nudge){
  Arctic_map +
    geom_tile(data = filter(df, depth <= 50), 
                aes(x = lon, y = lat, fill = pred_val)) +
    scale_fill_viridis_c(legend_title) +
    theme(strip.background = element_rect(colour = "white", fill = "white"),
          # legend.position = "topright",
          legend.position = c(x_nudge, 0.96),
          # legend.justification = c(0.9, 0.9),
          # legend.position = c(0, 0),
          legend.direction = "horizontal",
          legend.spacing.y = unit(0, "mm"))
}

# Predictions
pred_kelpcover <- Arctic_cover_predict(best_rf_kelpcover$choice_model)
cover_squiz(pred_kelpcover, "Total cover (%)", 0.785)
ggsave("talk/figure/prediction_kelpcover.png", width = 6, height = 6)
pred_laminariales <- Arctic_cover_predict(best_rf_laminariales$choice_model)
cover_squiz(pred_laminariales, "Laminariales cover (%)", 0.745)
ggsave("talk/figure/prediction_laminariales.png", width = 6, height = 6)
pred_agarum <- Arctic_cover_predict(best_rf_agarum$choice_model)
cover_squiz(pred_agarum, "Agarum cover (%)", 0.77)
ggsave("talk/figure/prediction_agarum.png", width = 6, height = 6)
pred_alaria <- Arctic_cover_predict(best_rf_alaria$choice_model)
cover_squiz(pred_alaria, "Alaria cover (%)", 0.78)
ggsave("talk/figure/prediction_alaria.png", width = 6, height = 6)


# 6: Linear regression ----------------------------------------------------


