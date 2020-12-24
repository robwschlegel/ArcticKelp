# analyses/10_figures
# This script houses the code used to create the final figures for the manuscript


# Setup -------------------------------------------------------------------

library(tidyverse)
library(raster)

# Function for re-loading .RData files as necessary
loadRData <- function(fileName){
  load(fileName)
  get(ls()[ls() != "fileName"])
}

# Load data used for maps etc.
load("data/Arctic_AM.RData")
colnames(Arctic_AM)[4] <- "depth"
Arctic_AM <- Arctic_AM %>%
  mutate(lon = round(lon, 4), lat = round(lat, 4))

# The present data
load("data/Arctic_BO.RData")

# Coordinates only
global_coords <- dplyr::select(Arctic_BO, lon, lat) %>% 
  mutate(env_index = 1:nrow(Arctic_BO))

# The species occurrence data
sps_files <- dir("metadata", full.names = T, pattern = "rarefied")
sps_names <- str_remove(dir("metadata", full.names = F, pattern = "rarefied"), pattern = "_Arct_rarefied_points.csv")


# Table 1 -----------------------------------------------------------------
# Importance of variables by kelps and models


# Figure 1 ----------------------------------------------------------------
# The map of the study area with the sample points
# Also add colour points for abundance (%)


# Figure 2 ----------------------------------------------------------------
# Relationship with important variables and kelp PA and abundance


# Figure 3 ----------------------------------------------------------------
# The ensemble model results

# Function to convert rasters to data.frames
rast_df <- function(rast, projection_name = NULL){
  df_out <- as.data.frame(rast[[1]], xy = T) %>% 
    `colnames<-`(c("lon", "lat", "presence")) %>% 
    mutate(lon = round(lon, 4), lat = round(lat, 4)) %>% 
    left_join(Arctic_AM, by = c("lon", "lat")) %>% 
    na.omit()
  if(!is.null(projection_name)) df_out$projection <- projection_name
  return(df_out)
}

# Function for prepping the ensemble model data
ensemble_prep <- function(sps_choice){
  
  # Load the ensemble projections
  biomod_project_present <- loadRData(paste0(sps_choice,"/proj_present/proj_present_",sps_choice,"_ensemble_TSSbin.RData"))
  biomod_project_2050 <- loadRData(paste0(sps_choice,"/proj_2050/proj_2050_",sps_choice,"_ensemble_TSSbin.RData"))
  biomod_project_2100 <- loadRData(paste0(sps_choice,"/proj_2100/proj_2100_",sps_choice,"_ensemble_TSSbin.RData"))
  
  # Convert to data.frames
  df_project_present <- rast_df(biomod_project_present[[1]], "proj_pres")
  df_project_2050 <- rast_df(biomod_project_2050[[1]], "proj_2050")
  df_project_2100 <- rast_df(biomod_project_2100[[1]], "proj_2100")
  
  # Combine and exit
  df_project_all <- rbind(df_project_present, df_project_2050, df_project_2100)
  return(df_project_all)
}

# Function for visualising changes over time
ensemble_diff_plot <- function(df_future, year_label){
  plot_out <- left_join(df_project_present, df_future, 
                        by = c("lon", "lat", "land_distance", "depth")) %>% 
    mutate(change = presence.x - presence.y,
           change = ifelse(presence.x == F & presence.y == F, NA, change),
           change =  factor(change, 
                            levels = c("-1", "0", "1"),
                            labels = c("gain", "no change", "loss"))) %>% 
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
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank())
}

# Function for plotting ensemble model results
# sps_choice <- "Ldig"
ensemble_plot <- function(sps_choice){
  
  # Prep data for plotting
  if(sps_choice %in% c("Ldig", "Lsol", "Slat")){
    sps_choice <- c("Ldig", "Lsol", "Slat")
    df_project <- plyr::ddply(data.frame(sps_name = c("Ldig", "Lsol", "Slat")), c("sps_name"), ensemble_prep) %>% 
      pivot_wider(values_from = presence, names_from = sps_name) %>% 
      mutate(presence = ifelse(Ldig + Lsol + Slat >= 2, TRUE, FALSE)) %>% 
      dplyr::select(lon, lat, presence, land_distance, depth, projection)
  } else{
    df_project <- ensemble_prep(sps_choice)
  }
  
  # Load the species points
  grepl(paste(sps_choice, collapse = "|"), sps_files)
  sps_points <- map_df(sps_files[str_which(sps_files, sps_choice)]) %>% 
    mutate(env_index = as.vector(knnx.index(as.matrix(global_coords[,c("lon", "lat")]),
                                            as.matrix(.[,2:3]), k = 1))) %>%
    left_join(global_coords, by = "env_index") %>% 
    dplyr::select(Sp, lon.y, lat.y) %>%
    dplyr::rename(lon = lon.y, lat = lat.y)
  
}

# Function that outputs BIOMOD projection comparison figures
plot_biomod <- function(sps_choice){
  
  # Load the species points
  sps_points <- read_csv(sps_files[str_which(sps_files, sps_choice)]) %>% 
    mutate(env_index = as.vector(knnx.index(as.matrix(global_coords[,c("lon", "lat")]),
                                            as.matrix(.[,2:3]), k = 1))) %>%
    left_join(global_coords, by = "env_index") %>% 
    dplyr::select(Sp, lon.y, lat.y) %>%
    dplyr::rename(lon = lon.y, lat = lat.y)
  
  # Visualise present data
  plot_present <- df_project_present %>% 
    filter(land_distance <= 100 | depth <= 100) %>% 
    filter(presence == TRUE) %>% 
    ggplot(aes(x = lon, y = lat)) +
    geom_tile(aes(fill = presence)) +
    borders(fill = "grey90", colour = "black") +
    geom_point(data = sps_points, colour = "yellow", size = 0.5) +
    scale_y_continuous(breaks = c(60, 70), labels = c("60°N", "70°N")) +
    scale_x_continuous(breaks = c(-80, -60), labels = c("80°W", "60°W")) +
    coord_quickmap(xlim = c(bbox_arctic[1], bbox_arctic[2]),
                   ylim = c(bbox_arctic[3], bbox_arctic[4]), expand = F) +
    scale_fill_manual(values = c("forestgreen")) +
    labs(x = NULL, y = NULL, title = paste0(sps_choice,": Present")) +
    theme_bw() +
    theme(legend.position = "bottom")
  
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
  ggsave(paste0("graph/biomod_diff_",sps_choice,".png"), plot_ALL, width = 8, height = 4.7)
}

# Create all visuals
registerDoParallel(cores = 5)
plyr::l_ply(sps_names, plot_biomod, .parallel = T)


# Figure 4 ----------------------------------------------------------------
# The random forest model results
# Combine into one dataframe to have the same legend via facets


# Figure 5 ----------------------------------------------------------------

# Comparison of the present modelling approaches


