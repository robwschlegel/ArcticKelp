# talk/ArcticNet_2020.R
# The purpose of this script is to house all of the code used to 
# create the figures etc. used in the talk with the same name

# 1: Setup ----------------------------------------------------------------

source("analyses/4_kelp_cover.R")

# Other libraries
library(randomForest)
library(sdmpredictors)
library(biomod2)
library(doParallel)
library(FNN)

# Random forest model results
load("data/best_rf_kelpcover.RData")
load("data/best_rf_laminariales.RData")
load("data/best_rf_agarum.RData")
load("data/best_rf_alaria.RData")

# Load top variable dataframes
load("data/top_var_kelpcover.RData")
load("data/top_var_laminariales.RData")
load("data/top_var_agarum.RData")
load("data/top_var_alaria.RData")

# The base global map with some corrections
load("../MHWapp/metadata/map_base.Rdata")

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

# Function for re-loading .RData files as necessary
loadRData <- function(fileName){
  load(fileName)
  get(ls()[ls() != "fileName"])
}

# Site coordinates
adf_summary_coords <- left_join(adf_summary, study_sites, by = c("Campaign", "site"))

# Mean values per site
adf_summary_mean_coords <- filter(adf_summary_coords, family == "kelp.cover") %>% 
  group_by(lon, lat) %>%
  summarise(mean_cover = mean(mean_cover),
            range_cover = max(mean_cover)-min(mean_cover))

# The species occurrence data
sps_files <- dir("metadata", full.names = T, pattern = "rarefied")
sps_names <- str_remove(dir("metadata", full.names = F, pattern = "rarefied"), pattern = "_Arct_rarefied_points.csv")

# Holder for species full names
sps_full_names <- data.frame(sps = sps_names,
                             sps_long = c("Agarum clathratum",
                                          "Alaria esculenta",
                                          "Laminaria digitata",
                                          "Laminaria solidungula",
                                          "Saccharina latissima"))

# Load BO layer names used for the ensemble models
# NB: Somehow the ensembles were given a different salinity variable
# This will need to be corrected later
load("metadata/BO_vars.RData")
BO_vars[4] <- "BO2_salinityltmax_ss" 
var_full_names <- data.frame(var = BO_vars,
                             var_long = c("Temperature (bottom, min.)",
                                          "Temperature (bottom, max.)",
                                          "Temperature (surface, max.)",
                                          "Salinity (surface, max.)",
                                          "Ice thickness (surface, min.)",
                                          "Iron (bottom, max.)",
                                          "Phosphate (bottom, max.)",
                                          "Current velocity (bottom, min.)"))

# Load data used for maps etc.
load("data/Arctic_AM.RData")
colnames(Arctic_AM)[4] <- "depth"
Arctic_AM <- Arctic_AM %>%
  mutate(lon = round(lon, 4), lat = round(lat, 4))

# Load the Arctic data
load("data/Arctic_BO.RData")

# Coordinates only
global_coords <- dplyr::select(Arctic_BO, lon, lat) %>% 
  mutate(env_index = 1:nrow(Arctic_BO))

# Crop variables to project standard
study_site_env <- study_site_env %>% 
  dplyr::select(site:land_distance, all_of(BO_vars))

# All of the BO variables matched to sites
kelp_all <- adf %>% 
  dplyr::select(Campaign, site, depth, -c(Bedrock..:sand), kelp.cover, Laminariales, Agarum, Alaria) %>% 
  left_join(study_site_env, by = c("Campaign", "site")) %>%
  mutate(kelp.cover = ifelse(kelp.cover > 100, 100, kelp.cover)) %>% # Correct values over 100
  dplyr::select(-lon_env, -lat_env, -lon, -lat) %>%
  dplyr::select(-depth) %>% # This is not used in Jesi's model
  dplyr::select(-bathy, -land_distance) %>% # Decided against these variables
  na.omit() # No missing data


# 2: Arctic region figures ------------------------------------------------

# Load the Arctic study region shape file
Arctic_poly <- readOGR(dsn = "metadata/", layer = "amaplim_lam_poly")
plot(Arctic_poly)

# Convert to a more useful coordinate system
Arctic_flat <- spTransform(Arctic_poly, "+init=epsg:4326 +proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0")
plot(Arctic_flat)

# Convert that to a data.frame
Arctic_boundary <- as.data.frame(Arctic_flat@polygons[[1]]@Polygons[[1]]@coords) %>%
  `colnames<-`(c("lon", "lat")) %>%
  arrange(lon)
Arctic_boundary <- rbind(Arctic_boundary,
                         data.frame(lon = rev(Arctic_boundary$lon),
                                    lat = rep(90, nrow(Arctic_boundary))))

# Overall regions
Arctic_region_map <- ggplot(data = Arctic_boundary, aes(x = lon, y = lat)) +
  geom_rect(aes(xmin = bbox_arctic[1], xmax = bbox_arctic[2],
                ymin = bbox_arctic[3], ymax = bbox_arctic[4]),
            fill = "forestgreen", colour = "black", alpha = 0.2) +
  geom_polygon(data = map_base, aes(group = group)) +
  geom_polygon(fill = "cadetblue1", colour = "black", alpha = 0.2) +
  geom_point(data = CANA_kelp, aes(x = Longitude, y = Latitude), colour = "red") +
  coord_quickmap(ylim = c(40, 90), expand = F) +
  scale_y_continuous(breaks = c(60, 80), labels = c("60°N", "80°N")) +
  scale_x_continuous(breaks = c(-100, 0, 100), labels = c("100°W", "0°E", "100°E")) +
  labs(x = NULL, y = NULL, title = "Arctic boundary (blue)\nStudy region (green)\nCANA data (red)")#, subtitle = "Study region (green)\nCANA data (red)")
ggsave("talk/figure/Arctic_region_map.png", Arctic_region_map, width = 10)

# Eastern Canadian Arctic question mark map
kelp_question <- Arctic_map +
  geom_label(aes(x = -70, y = 65, label = "?"), size = 60, alpha = 0.6)
kelp_question
ggsave("talk/figure/kelp_question.png", kelp_question, height = fig_height)

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
ggsave("talk/figure/study_site_campaign.png", study_site_campaign, height = fig_height)

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
ggsave("talk/figure/study_site_mean_cover.png", study_site_mean_cover, height = fig_height)


# 3: Model accuracy -------------------------------------------------------

## Ensemble
# FUnction for processing ensemble model accuracy results
conf_plot_ensemble <- function(sps_choice){
  
  # Load chosen biomod_model and print evaluation scores
  biomod_model <- loadRData(paste0(sps_choice,"/",sps_choice,".",sps_choice,".models.out"))
  
  # Get species name
  sps_label <- sps_full_names$sps_long[sps_full_names$sps == sps_choice]
  
  # Model evaluation by algorithm
  plot_out <- models_scores_graph(biomod_model, by = "models", metrics = c('ROC','TSS'),
                      xlim = c(0.5,1), ylim = c(0.5,1)) + ggtitle(sps_label) +
    geom_hline(aes(yintercept = 0.7), colour = "red") +
    coord_fixed(xlim = c(0.5, 1), ylim = c(0.5, 1)) +
    # guides(colour = guide_legend(nrow = 2, byrow = TRUE)) +
    theme(panel.border = element_rect(fill = NA, colour = "black"),
          legend.position = "right")
  # plot_out
  ggsave(paste0("talk/figure/conf_ensemble_",sps_choice,".png"), width = 4, height = 3)
}

# Run all of them
registerDoParallel(cores = 5)
plyr::l_ply(sps_names, conf_plot_ensemble, .parallel = T)


## RF
# Function for creating figure showing confidence intervals of prediction accuracy
conf_plot_RF <- function(df, plot_title){
  
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
    # geom_segment(data = conf_best, aes(xend = original, y = mean_acc, yend = 0), 
                 # colour = "purple", size = 1.2, alpha = 0.8) +
    # geom_point(data = conf_best, aes(y = mean_acc), colour = "purple", size = 3, alpha = 0.8) +
    geom_label(data = conf_mean_label, 
               aes(x = 75, y = 75, label = paste0("Mean accuracy: ",mean_acc,"±",sd_acc,"; r = ",r_acc))) +
    # geom_label(data = conf_best_label, colour = "purple",
               # aes(x = 75, y = 60, label = paste0("Best accuracy: ",mean_acc,"±",sd_acc,"; r = ",r_acc))) +
    scale_y_continuous(limits = c(-100, 100)) +
    labs(y = "Range in accuracy of predictions", x = "Original value (% cover)", title = plot_title) +
    theme(panel.border = element_rect(fill = NA, colour = "black"))
}

# Create the plots
conf_plot_RF(best_rf_kelpcover$accuracy_reg, "Total cover confidence")
ggsave("talk/figure/conf_RF_kelpcover.png", height = fig_height-2, width = fig_width)
conf_plot_RF(best_rf_laminariales$accuracy_reg, "Laminariales cover confidence")
ggsave("talk/figure/conf_RF_laminariales.png", height = fig_height-2, width = fig_width)
conf_plot_RF(best_rf_agarum$accuracy_reg, "Agarum cover confidence")
ggsave("talk/figure/conf_RF_agarum.png", height = fig_height-2, width = fig_width)
conf_plot_RF(best_rf_alaria$accuracy_reg, "Alaria cover confidence")
ggsave("talk/figure/conf_RF_alaria.png", height = fig_height-2, width = fig_width)


# 4: Variable importance --------------------------------------------------

## Ensemble
# Function for processing ensemble variable importance
var_imp_ensemble <- function(sps_choice){
  
  # Load chosen biomod_model and print evaluation scores
  biomod_model <- loadRData(paste0(sps_choice,"/",sps_choice,".",sps_choice,".models.out"))
  
  # Get species name
  sps_label <- sps_full_names$sps_long[sps_full_names$sps == sps_choice]
  
  # Get variable importance
  models_var_import <- get_variables_importance(biomod_model)
  var_import <- data.frame(imp = round(apply(apply(models_var_import, c(1,2), mean, na.rm = T), 1, mean), 2)) %>%  # Overall mean per variable
    mutate(var = row.names(.)) %>% 
    left_join(var_full_names, by = "var") %>% 
    dplyr::select(var_long, imp) %>% 
    dplyr::rename(`Data layer` = var_long, Importance = imp) %>% 
    arrange(-Importance)
  save(var_import, file = paste0("talk/data/var_imp_ensemble_",sps_choice,".RData"))
}

# Run all
registerDoParallel(cores = 5)
plyr::l_ply(sps_names, var_imp_ensemble)


## RF
# Convenience function to get the cleaned up variables
var_imp_RF <- function(df){
  res <- df %>% 
    left_join(var_full_names, by = "var") %>% 
    dplyr::select(var_long, X.IncMSE) %>%
    mutate(X.IncMSE = round(X.IncMSE)) %>% 
    arrange(-X.IncMSE) %>% 
    dplyr::rename(`Data layer` = var_long, `% Increase MSE` = X.IncMSE) 
}

# Order and make pretty for presentation
var_imp_RF_kelpcover <- var_imp_RF(top_var_kelpcover)
save(var_imp_RF_kelpcover, file = "talk/data/var_imp_RF_kelpcover.RData")
var_imp_RF_laminariales <- var_imp_RF(top_var_laminariales)
save(var_imp_RF_laminariales, file = "talk/data/var_imp_RF_laminariales.RData")
var_imp_RF_agarum <- var_imp_RF(top_var_agarum)
save(var_imp_RF_agarum, file = "talk/data/var_imp_RF_agarum.RData")
var_imp_RF_alaria <- var_imp_RF(top_var_alaria)
save(var_imp_RF_alaria, file = "talk/data/var_imp_RF_alaria.RData")


# 5: Project values -------------------------------------------------------

## Ensemble
# Function that outputs BIOMOD projection comparison figures
plot_proj_ensemble <- function(sps_choice){
  
  # Load the species points
  sps_points <- read_csv(sps_files[str_which(sps_files,sps_choice)]) %>% 
    mutate(env_index = as.vector(knnx.index(as.matrix(global_coords[,c("lon", "lat")]),
                                            as.matrix(.[,2:3]), k = 1))) %>%
    left_join(global_coords, by = "env_index") %>% 
    dplyr::select(Sp, lon.y, lat.y) %>%
    dplyr::rename(lon = lon.y, lat = lat.y)
  
  # Get species name
  sps_label <- sps_full_names$sps_long[sps_full_names$sps == sps_choice]
  
  # Load the ensemble projections
  biomod_project_present <- loadRData(paste0(sps_choice,"/proj_present/proj_present_",sps_choice,"_ensemble_TSSbin.RData"))
  biomod_project_2050 <- loadRData(paste0(sps_choice,"/proj_2050/proj_2050_",sps_choice,"_ensemble_TSSbin.RData"))
  biomod_project_2100 <- loadRData(paste0(sps_choice,"/proj_2100/proj_2100_",sps_choice,"_ensemble_TSSbin.RData"))
  
  # Convert to data.frames
  rast_df <- function(rast){
    df_out <- as.data.frame(rast[[1]], xy = T) %>% 
      `colnames<-`(c("lon", "lat", "presence")) %>% 
      mutate(lon = round(lon, 4), lat = round(lat, 4)) %>% 
      left_join(Arctic_AM, by = c("lon", "lat")) %>% 
      na.omit() 
  }
  df_project_present <- rast_df(biomod_project_present[[1]])
  df_project_2050 <- rast_df(biomod_project_2050[[1]])
  df_project_2100 <- rast_df(biomod_project_2100[[1]])
  
  # Visualise present data
  plot_present <- df_project_present %>% 
    filter(land_distance <= 100 | depth <= 100) %>% 
    ggplot(aes(x = lon, y = lat)) +
    geom_tile(aes(fill = presence)) +
    borders(fill = "grey90", colour = "black") +
    geom_point(data = sps_points, colour = "yellow", size = 0.5) +
    scale_y_continuous(breaks = c(60, 70), labels = c("60°N", "70°N")) +
    scale_x_continuous(breaks = c(-80, -60), labels = c("80°W", "60°W")) +
    coord_quickmap(xlim = c(bbox_arctic[1], bbox_arctic[2]),
                   ylim = c(bbox_arctic[3], bbox_arctic[4]), expand = F) +
    scale_fill_manual(values = c("grey20", "forestgreen")) +
    labs(x = NULL, y = NULL, title = paste0(sps_label)) +
    theme_bw() +
    theme(legend.position = "bottom")
  
  # Function for visualising changes over time
  plot_diff <- function(df_future, year_label){
    plot_out <- left_join(df_project_present, df_future, 
                          by = c("lon", "lat", "land_distance", "depth")) %>% 
      mutate(change = factor(presence.x - presence.y, 
                             levels = c("-1", "0", "1"),
                             labels = c("increase", "same", "decrease"))) %>% 
      na.omit() %>% 
      filter(land_distance <= 100 | depth <= 100) %>% 
      ggplot(aes(x = lon, y = lat)) +
      geom_tile(aes(fill = change)) +
      borders(fill = "grey90", colour = "black") +
      scale_y_continuous(breaks = c(60, 70), labels = c("60°N", "70°N")) +
      scale_x_continuous(breaks = c(-80, -60), labels = c("80°W", "60°W")) +
      coord_quickmap(xlim = c(bbox_arctic[1], bbox_arctic[2]),
                     ylim = c(bbox_arctic[3], bbox_arctic[4]), expand = F) +
      scale_fill_brewer(palette = "Set1", direction = -1) +
      labs(x = NULL, y = NULL, title = paste0(sps_choice,": Present - ",year_label)) +
      theme_bw() +
      theme(legend.position = "bottom",
            panel.border = element_rect(fill = NA, colour = "black"),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank())
  }
  
  # Visualise present - 2050
  plot_2050 <- plot_diff(df_project_2050, "2050")
  
  # Visualise present - 2100
  plot_2100 <- plot_diff(df_project_2100, "2100")
  
  # Combine and save
  plot_ALL <- cowplot::plot_grid(
    cowplot::plot_grid(
      plot_present + theme(legend.position = "none"),
      plot_2050 + theme(legend.position = "none"),
      plot_2100 + theme(legend.position = "none"),
      ncol = 3,
      align = "v"),
    cowplot::plot_grid(
      cowplot::get_legend(plot_present),
      ggplot() + theme_void(),
      cowplot::get_legend(plot_2100),
      ncol = 3, rel_widths = c(1, 0, 1.5)),
    nrow = 2, rel_heights = c(10,1)
  )
  ggsave(paste0("talk/figure/future_diff_ensemble_",sps_choice,".png"), plot_ALL, width = 8, height = 4.7)
}

# Create all visuals
registerDoParallel(cores = 5)
plyr::l_ply(sps_names, plot_proj_ensemble, .parallel = T)

## RF
project_compare <- function(best_rf, kelp_choice){
  
  # Calculate differences
  project_diff <- best_rf$project_multi %>% 
    mutate(pred_diff_2050 = pred_2050_mean - pred_present_mean,
           pred_diff_2100 = pred_2100_mean - pred_present_mean)
  
  # Scale for difference plots
  diff_range <- range(c(project_diff$pred_diff_2050, project_diff$pred_diff_2100), na.rm = T)
  
  # Present plot
  p_present <-  Arctic_map +
    # geom_tile(data = df, # No filter
    geom_tile(data = filter(project_diff, depth <= 100 | land_distance <= 100),
              aes(x = lon, y = lat, fill = pred_present_mean)) +
    scale_fill_viridis_c(paste0(kelp_choice," cover (%)")) +
    theme(legend.position = "bottom")
  
  # 2050 plot
  p_2050 <-  Arctic_map +
    # geom_tile(data = df, # No filter
    geom_tile(data = filter(project_diff, depth <= 100 | land_distance <= 100),
              aes(x = lon, y = lat, fill = pred_diff_2050)) +
    scale_fill_gradient2("2050 - present", low = "blue", high = "red", limits = diff_range) +
    theme(legend.position = "bottom",
          axis.ticks.y = element_blank(),
          axis.text.y = element_blank())
  
  # 2050 plot
  p_2100 <-  Arctic_map +
    # geom_tile(data = df, # No filter
    geom_tile(data = filter(project_diff, depth <= 100 | land_distance <= 100),
              aes(x = lon, y = lat, fill = pred_diff_2100)) +
    scale_fill_gradient2("2100 - present", low = "blue", high = "red", limits = diff_range) +
    theme(legend.position = "bottom",
          axis.ticks.y = element_blank(),
          axis.text.y = element_blank())
  
  # Save and Print
  p_all <- ggpubr::ggarrange(p_present, p_2050, p_2100, ncol = 3, nrow = 1, align = "hv")
  ggsave(paste0("talk/figure/future_diff_RF_",kelp_choice,".png"), p_all, width = 9, height = 5)
  # p_all
}

# Visualise the comparisons
project_compare(best_rf_kelpcover, "Total_cover")
project_compare(best_rf_laminariales, "Laminariales")
project_compare(best_rf_agarum, "Agarum")
project_compare(best_rf_alaria, "Alaria")


# 6: Linear regression ----------------------------------------------------

kelp_choice <- "Laminariales"

plot_scatter_kelp <- function(kelp_choice){
  kelp_all %>% 
    pivot_longer(cols = BO2_templtmin_bdmax:BO21_curvelltmin_bdmax, names_to = "var") %>% 
    left_join(var_full_names, by = "var") %>% 
    ggplot(aes_string(x = kelp_choice, y = "value")) +
    geom_point() +
    geom_smooth(method = "lm") +
    facet_wrap(~var_long, scales = "free_y") +
    labs(x = "Percent cover", title = kelp_choice)
}

plot_scatter_kelp("kelp.cover")
ggsave("talk/figure/var_scatter_Total_cover.png", height = 5, width = 10)
plot_scatter_kelp("Laminariales")
ggsave("talk/figure/var_scatter_Laminariales.png", height = 5, width = 10)
plot_scatter_kelp("Agarum")
ggsave("talk/figure/var_scatter_Agarum.png", height = 5, width = 10)
plot_scatter_kelp("Alaria")
ggsave("talk/figure/var_scatter_Alaria.png", height = 5, width = 10)

