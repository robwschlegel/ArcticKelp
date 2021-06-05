# analyses/10_figures
# This script houses the code used to create the final figures for the manuscript

# TODO: Figures showing response curves between RF % cover and variables used
# Would be good to have for ensemble models results, too
# Calculate amount of change in cover over time


# Setup -------------------------------------------------------------------

# Load study sites and base packages
source("analyses/1_study_region_sites.R")
source("analyses/4_kelp_cover.R")

# Other libraries
library(ggOceanMaps) # https://mikkovihtakari.github.io/ggOceanMaps/index.html
library(ggtext)
library(raster)
library(FNN)
library(sdmpredictors)
library(sf)
library(sp)

# Function for re-loading .RData files as necessary
loadRData <- function(fileName){
  load(fileName)
  get(ls()[ls() != "fileName"])
}

# The species occurrence data
sps_files <- dir("metadata", full.names = T, pattern = "rarefied")
sps_names <- str_remove(dir("metadata", full.names = F, pattern = "rarefied"), pattern = "_Arct_rarefied_points.csv")

# Environmental data
load("data/Arctic_BO.RData")
Arctic_BO <- Arctic_BO %>% 
  mutate(lon = round(lon, 4), lat = round(lat, 4))
load("data/Arctic_AM.RData")
Arctic_AM <- Arctic_AM %>% 
  mutate(lon = round(lon, 4), lat = round(lat, 4)) %>% 
  dplyr::rename(depth = bathy)

# Coordinates only
global_coords <- dplyr::select(Arctic_BO, lon, lat) %>% 
  mutate(env_index = 1:nrow(Arctic_BO))

# The base map to use for everything else
Arctic_map <- ggplot() +
  borders(fill = "grey70", colour = NA) +
  scale_y_continuous(breaks = c(60, 70), labels = c("60°N", "70°N")) +
  scale_x_continuous(breaks = c(-80, -60), labels = c("80°W", "60°W")) +
  coord_quickmap(xlim = c(bbox_arctic[1], bbox_arctic[2]),
                 ylim = c(bbox_arctic[3], bbox_arctic[4]), expand = F) +
  labs(x = NULL, y = NULL) +
  theme(panel.border = element_rect(colour = "black", fill = NA))

# The base global map with some corrections
load("metadata/map_base.Rdata")

# MEOW
MEOW <- read_sf("metadata/MEOW/meow_ecos.shp") %>% 
  filter(REALM == "Arctic")


# Figure 1 ----------------------------------------------------------------

# The map of the study area with the sample points

# Load all species points
sps_data <- map_df(sps_files, read_csv) %>% 
  # mutate(Sp = case_when(Sp == "Slat" ~ "Aaa",
  #                       TRUE ~ Sp)) %>% 
  mutate(Sp = factor(Sp, levels = c("Acla", "Aesc", "Ldig", "Lsol", "Slat")))

# Site coordinates
adf_summary_coords <- left_join(adf_summary, study_sites, by = c("Campaign", "site"))

# Mean values per site
adf_summary_mean_coords <- filter(adf_summary_coords, family == "kelp.cover") %>% 
  group_by(lon, lat, Campaign) %>%
  summarise(mean_cover = mean(mean_cover),
            range_cover = max(mean_cover)-min(mean_cover), .groups = "drop")

# Create spatial polygon data frame from Arctic bounding box
bbox_top <- data.frame(lon = seq(bbox_arctic[1], bbox_arctic[2], length.out = 100), 
                       lat = bbox_arctic[3], id = "bbox")
bbox_bottom <- data.frame(lon = seq(bbox_arctic[2], bbox_arctic[1], length.out = 100), 
                          lat = bbox_arctic[4], id = "bbox")
bbox_df <- rbind(bbox_top, bbox_bottom)

# only want lon-lats in the list, not the names
bbox_list <- lapply(split(bbox_df, bbox_df$id), function(x) { x["id"] <- NULL; x })

# Convert to polygon and add id variable 
bbox_poly <- Polygons(sapply(bbox_list, Polygon), ID = 1)

# Create SpatialPolygons object
bbox_spatial <- SpatialPolygons(list(bbox_poly), 
                                proj4string = CRS("+init=epsg:4326 +proj=longlat")) 
# plot(bbox_spatial)

# Load the Arctic study region shape file
Arctic_poly <- readOGR(dsn = "metadata/", layer = "amaplim_lam_poly")
# plot(Arctic_poly)

# Convert to a square coordinate system
Arctic_flat <- spTransform(Arctic_poly, "+init=epsg:4326 +proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0")
# plot(Arctic_flat)

# Convert that to a data.frame
Arctic_boundary <- as.data.frame(Arctic_flat@polygons[[1]]@Polygons[[1]]@coords) %>%
  `colnames<-`(c("lon", "lat")) %>%
  arrange(lon)
Arctic_boundary <- rbind(Arctic_boundary,
                         data.frame(lon = rev(Arctic_boundary$lon),
                                    lat = rep(90, nrow(Arctic_boundary))))

# Overall regions
fig_1a <- basemap(limits = c(-180, 180, 40, 90), bathymetry = T) +
  annotation_spatial(bbox_spatial, fill = "forestgreen", alpha = 0.2) +
  annotation_spatial(Arctic_poly, fill = "cadetblue1", alpha = 0.2) +
  # annotation_spatial(MEOW, aes(colour = ECOREGION), fill = NA) +
  # NB: This overplotting of Slat needs to be addressed
  # This may be due to the points needing to be converted to a spatial points file first
  geom_spatial_point(data = sps_data, crs = 4326, shape = 21,
                     aes(x = lon, y = lat), fill = "hotpink", colour = "black")+
  geom_spatial_point(data = study_sites, crs = 4326, shape = 22,
                     aes(x = lon, y = lat), fill = "purple", colour = "black")
fig_1a
ggsave("figures/fig_1a.png", fig_1a, height = 6, width = 8)

# Add the names of the main regions: Hudson Bay, Hudson Strait, Lancaster Sound
label_df <- data.frame(lon = c(-85, -75.17, -82.92, -63.834153, -83.4011174, -90.0706224, -82.7015637),
                       lat = c(60, 62.544, 74.425, 72.172485, 67.5202107, 70.5938471, 53.0693041),
                       loc = c("Hudson Bay", "Hudson Strait", "Lancaster Sound", "Baffin Bay", "Foxe Basin", "Gulf of \nBoothia", "James Bay"))

# ArcticKelp campaign map
fig_1b <- ggplot() +
  # geom_tile(aes(fill = presence)) +
  borders(fill = "grey30", colour = "black") +
  geom_point(data = study_sites, shape = 21, colour = "black", 
             fill = "magenta", size = 2,
             aes(x = lon, y = lat)) +
  geom_label(data = label_df, aes(x = lon, y = lat, label = loc), alpha = 0.9) +
  scale_y_continuous(breaks = c(60, 70), labels = c("60°N", "70°N")) +
  scale_x_continuous(breaks = c(-80, -60), labels = c("80°W", "60°W")) +
  coord_quickmap(xlim = c(bbox_arctic[1], bbox_arctic[2]),
                 ylim = c(bbox_arctic[3], bbox_arctic[4]), expand = F) +
  labs(y = NULL, x = NULL) +
  # scale_fill_manual("Suitable", values = c("grey40")) +
  # labs(x = NULL, y = NULL, title = sps_title) +
  # theme_bw() +
  theme(legend.position = "bottom",
        plot.title = element_text(face = "italic"),
        panel.background = element_rect(fill = "grey95"),
        panel.border = element_rect(colour = "black", fill = NA))
fig_1b

# Combine and save
fig_1 <- ggpubr::ggarrange(fig_1a, fig_1b, ncol = 1, labels = c("A)", "B)"))
ggsave("figures/fig_1.png", fig_1, height = 12, width = 8)


# Table 1 -----------------------------------------------------------------
# Importance of variables by kelps and models


# Figure 2 ----------------------------------------------------------------
# Importance of variables per taxonomic group
# Jesi made this figure with image editing software


# Figure 3 ----------------------------------------------------------------
# The ensemble model results

# TODO: Consider using the km^2 calculations directly from the rasters

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
  
  # Calculate sqaure are in kilometres directly
  sum(biomod_project_present[] == 1, na.rm = T) * res(biomod_project_present)[1]^2
  
  # Convert to data.frames
  df_project_present <- rast_df(biomod_project_present[[1]], "proj_pres")
  df_project_2050 <- rast_df(biomod_project_2050[[1]], "proj_2050")
  df_project_2100 <- rast_df(biomod_project_2100[[1]], "proj_2100")
  
  # Combine and exit
  df_project_all <- rbind(df_project_present, df_project_2050, df_project_2100) %>% 
    mutate(presence = as.integer(presence))
  return(df_project_all)
}

# Function for visualising changes over time
ensemble_diff_plot <- function(df, year_label, sq_area_labels){
  # Prepare label
  sq_area_label_proj <- sq_area_labels[colnames(sq_area_labels) == paste0("area_",year_label)][[1]]
  sq_area_label_sub <- sq_area_label_proj-sq_area_labels$area_pres[[1]]
  sq_area_label_text <- paste0(scales::comma(sq_area_label_sub)," km<sup>2</sup>")
  if(sq_area_label_sub > 0){
    sq_area_label_text <- paste0("+",sq_area_label_text)
    lab_col <- RColorBrewer::brewer.pal(9, "Blues")[7]
  } else{
    lab_col <- RColorBrewer::brewer.pal(9, "Reds")[7]
  }
  
  
  # Create plot
  diff_plot <- df %>%
    filter(land_distance <= 50 | depth <= 30) %>% 
    ggplot(aes(x = lon, y = lat)) +
    geom_tile(aes_string(fill = paste0("change_",year_label))) +
    borders(fill = "grey50", colour = "grey90", size = 0.2) +
    scale_y_continuous(breaks = c(60, 70), labels = c("60°N", "70°N")) +
    scale_x_continuous(breaks = c(-80, -60), labels = c("80°W", "60°W")) +
    coord_quickmap(xlim = c(bbox_arctic[1], bbox_arctic[2]),
                   ylim = c(bbox_arctic[3], bbox_arctic[4]), expand = F) +
    geom_richtext(data = sq_area_labels, size = 4, hjust = 1, colour = lab_col,
                  aes(x = -50, y = 51, label = sq_area_label_text)) +
    # scale_fill_brewer(palette = "Set1", direction = -1) +
    scale_fill_manual(values = c(RColorBrewer::brewer.pal(9, "Blues")[7], "grey80",
                                 RColorBrewer::brewer.pal(9, "Reds")[7])) +
    labs(x = NULL, y = NULL, fill = "Change", title = paste0(year_label," - present")) +
    theme(legend.position = "bottom",
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          panel.background = element_rect(fill = "grey100"),
          panel.border = element_rect(colour = "black", fill = NA))
  # diff_plot
}

# Function for plotting ensemble model results
# sps_choice <- "Acla"
ensemble_plot <- function(sps_choice, add_legend = F){
  
  # Create full species name
  if(sps_choice == "Lsol"){
    sps_title <- "Laminaria solidungula"
  } else if(sps_choice == "Slat"){
    sps_title <- "Saccharina latissima"
  } else if(sps_choice == "Acla"){
    sps_title <- "Agarum clathratum"
  } else if(sps_choice == "Aesc"){
    sps_title <- "Alaria esculenta"
  } else{
    stop("*sad robot noises*")
  }
  
  # Prep data for plotting
  df_project <- ensemble_prep(sps_choice) %>% 
    pivot_wider(names_from = projection, values_from = presence) %>% 
    replace(is.na(.), 0) %>% 
    mutate(change_2050 = proj_2050 - proj_pres,
           change_2050 =  factor(change_2050, 
                                 levels = c("1", "0", "-1"),
                                 labels = c("gain", "no change", "loss")),
           change_2100 = proj_2100 - proj_pres,
           change_2100 = factor(change_2100, 
                                levels = c("1", "0", "-1"),
                                labels = c("gain", "no change", "loss")))
  
  # Calculate sq area coverage per era
  sq_area_labels <- df_project %>% 
    filter(land_distance <= 50 | depth <= 30) %>% 
    summarise(area_pres = round(sum(sq_area*proj_pres, na.rm = T)),
              area_2050 = round(sum(sq_area*proj_2050, na.rm = T)),
              area_2100 = round(sum(sq_area*proj_2100, na.rm = T)))
  pres_text <- paste0(scales::comma(sq_area_labels$area_pres), " km<sup>2</sup>")
  
  # Load the species points
  sps_points <- map_df(sps_files[grepl(paste(sps_choice, collapse = "|"), sps_files)], read_csv) %>% 
    mutate(env_index = as.vector(knnx.index(as.matrix(global_coords[,c("lon", "lat")]),
                                            as.matrix(.[,2:3]), k = 1))) %>%
    left_join(global_coords, by = "env_index") %>% 
    dplyr::select(Sp, lon.y, lat.y) %>%
    dplyr::rename(lon = lon.y, lat = lat.y)
  
  # Visualise present data
  plot_present <- df_project %>% 
    filter(land_distance <= 50 | depth <= 30,
           proj_pres == 1) %>%
    mutate(proj_pres = "") %>%
    ggplot(aes(x = lon, y = lat)) +
    geom_tile(aes(fill = proj_pres)) +
    borders(fill = "grey50", colour = "grey90", size = 0.2) +
    geom_point(data = sps_points, shape = 21, colour = "black", fill = "hotpink", size = 0.5) +
    geom_richtext(data = sq_area_labels, size = 4, hjust = 1,
                  aes(x = -50, y = 51, label = pres_text)) +
    # annotate("label", x = -58, y = 78, label = paste(pres_text, "^2", sep = ""), parse = T) +
    scale_y_continuous(breaks = c(60, 70), labels = c("60°N", "70°N")) +
    scale_x_continuous(breaks = c(-80, -60), labels = c("80°W", "60°W")) +
    coord_quickmap(xlim = c(bbox_arctic[1], bbox_arctic[2]),
                   ylim = c(bbox_arctic[3], bbox_arctic[4]), expand = F) +
    scale_fill_manual("Suitable", values = c("grey80")) +
    labs(x = NULL, y = NULL, title = sps_title) +
    # theme_bw() +
    theme(legend.position = "bottom",
          plot.title = element_text(face = "italic"),
          panel.background = element_rect(fill = "grey100"),
          panel.border = element_rect(colour = "black", fill = NA))
  # plot_present
  
  # Visualise present - 2050 and 2100
  # This is done separately to make adding the legends at the end cleaner
  plot_2050 <- df_project %>% 
    filter(proj_pres != 0 | proj_2050 != 0) %>% 
    ensemble_diff_plot("2050", sq_area_labels)
  plot_2100 <- df_project %>% 
    filter(proj_pres != 0 | proj_2100 != 0) %>% 
    ensemble_diff_plot("2100", sq_area_labels)
  
  # Combine and exit
  plot_ALL <- cowplot::plot_grid(
    plot_present + theme(legend.position = "none"),
    plot_2050 + theme(legend.position = "none"),
    plot_2100 + theme(legend.position = "none"),
    ncol = 3, align = "v")
  if(add_legend){
    plot_ALL <- cowplot::plot_grid(
      cowplot::get_legend(plot_present),
      ggplot() + theme_void(),
      cowplot::get_legend(plot_2100),
      ncol = 3, rel_widths = c(1, 0, 1.5))
  }
  return(plot_ALL)
}

# Create all visuals
ensemble_Acla <- ensemble_plot("Acla")
ensemble_Aesc <- ensemble_plot("Aesc")
ensemble_Lsol <- ensemble_plot("Lsol")
ensemble_Slat <- ensemble_plot("Slat")
ensemble_legend <- ensemble_plot("Acla", add_legend = T)

# Combine into one mecha-figure
fig_3 <- ggpubr::ggarrange(ensemble_Acla, ensemble_Aesc, ensemble_Lsol, ensemble_Slat, ensemble_legend,
                           ncol = 1, labels = c("A)", "B)", "C)", "D)", ""), heights = c(1, 1, 1, 1, 0.15))
ggsave("figures/fig_3.png", fig_3, width = 7, height = 15)

# Agarum only for demo/talk
fig_3_agarum <- ggpubr::ggarrange(ensemble_Acla, ensemble_legend, ncol = 1, heights = c(1, 0.15))
ggsave("graph/fig_3_agarum.png", fig_3_agarum, width = 7, height = 4)
ggsave("talk/figure/fig_3_agarum.png", fig_3_agarum, width = 7, height = 4)

# Alaria only for talk
fig_3_alaria <- ggpubr::ggarrange(ensemble_Aesc, ensemble_legend, ncol = 1, heights = c(1, 0.15))
ggsave("talk/figure/fig_3_alaria.png", fig_3_alaria, width = 7, height = 4)

# Laminaria only for talk
fig_3_laminaria <- ggpubr::ggarrange(ensemble_Lsol, ensemble_legend, ncol = 1, heights = c(1, 0.15))
ggsave("talk/figure/fig_3_laminaria.png", fig_3_laminaria, width = 7, height = 4)

# Saccharina only for talk
fig_3_saccharina <- ggpubr::ggarrange(ensemble_Slat, ensemble_legend, ncol = 1, heights = c(1, 0.15))
ggsave("talk/figure/fig_3_saccharina.png", fig_3_saccharina, width = 7, height = 4)

# Find current surface area of all species overlayed
proj_present_Acla <- rast_df(loadRData("Acla/proj_present/proj_present_Acla_ensemble_TSSbin.RData")[[1]], "proj_pres")
proj_present_Aesc <- rast_df(loadRData("Aesc/proj_present/proj_present_Aesc_ensemble_TSSbin.RData")[[1]], "proj_pres")
proj_present_Lsol <- rast_df(loadRData("Lsol/proj_present/proj_present_Lsol_ensemble_TSSbin.RData")[[1]], "proj_pres")
proj_present_Slat <- rast_df(loadRData("Slat/proj_present/proj_present_Slat_ensemble_TSSbin.RData")[[1]], "proj_pres")
proj_present_overlay <- rbind(proj_present_Acla, proj_present_Aesc, proj_present_Lsol, proj_present_Slat) %>% 
  filter(presence == T) %>% 
  dplyr::select(lon, lat, sq_area) %>% 
  distinct()
sum(proj_present_overlay$sq_area)


# Figure 4 ----------------------------------------------------------------
# Combination of the modeling approaches

# Join ensemble and random forest results
# model_choice <- "agarum"
model_compare_plot <- function(model_choice, add_legend = F){
  
  # Create full species name and load ensemble data
  if(model_choice == "laminariales"){
    sps_title <- "Laminariales spp."
    df_project <- left_join(ensemble_prep("Lsol"), ensemble_prep("Slat"),
                            by = c("lon", "lat", "land_distance", "depth", "projection", "sq_area")) %>% 
      mutate(presence = case_when(presence.x == 1 | presence.y == 1 ~ 1, TRUE ~ 0)) %>% 
      dplyr::select(-presence.x, -presence.y)
  } else if(model_choice == "agarum"){
    sps_title <- "Agarum clathratum"
    df_project <- ensemble_prep("Acla")
  } else if(model_choice == "alaria"){
    sps_title <- "Alaria esculenta"
    df_project <- ensemble_prep("Aesc")
  } else{
    stop("*sad robot noises*")
  }
  
  # Pivot wide and fill in NA with 0
  df_project <- df_project %>% 
    pivot_wider(names_from = projection, values_from = presence) %>% 
    replace(is.na(.), 0)
  
  # Load random forest data
  best_rf <- loadRData(paste0("data/best_rf_",model_choice,".RData"))$project_multi
  
  # Join the models
  model_join <- left_join(best_rf, df_project, by = c("lon", "lat", "depth", "land_distance")) %>% 
    mutate(pred_present_round = plyr::round_any(pred_present_mean, 20),
           pred_present_round = ifelse(pred_present_round > 60, 60, pred_present_round),
           pred_diff_2050 = plyr::round_any(pred_2050_mean - pred_present_mean, 10),
           pred_diff_2100 = plyr::round_any(pred_2100_mean - pred_present_mean, 10)) %>% 
    na.omit()
  
  # Calculate sq area coverage per era
  sq_area_labels <- model_join %>% 
    filter(land_distance <= 50 | depth <= 100) %>% 
    summarise(area_pres_mean = mean(pred_present_round, na.rm = T),
              area_2050_mean = mean(pred_diff_2050, na.rm = T),
              area_2100_mean = mean(pred_diff_2100, na.rm = T),
              area_pres_perc = mean(pred_present_round*proj_pres, na.rm = T),
              area_2050_perc = mean(pred_diff_2050*proj_2050, na.rm = T),
              area_2100_perc = mean(pred_diff_2100*proj_2100, na.rm = T)) %>% 
    mutate(area_pres_label = paste0(round(area_pres_perc),"%"),
           area_2050_label = ifelse(area_2050_perc > 0 , 
                                    paste0("+",round(area_2050_perc),"%"),
                                    paste0(round(area_2050_perc),"%")),
           area_2100_label = ifelse(area_2100_perc > 0 , 
                                    paste0("+",round(area_2100_perc),"%"),
                                    paste0(round(area_2100_perc),"%")))
  
  # Plot
  p_present <- model_join %>% 
    filter(proj_pres == 1,
           # pred_present_round >= 10,
           depth <= 100 | land_distance <= 50) %>% 
    ggplot() +
    geom_tile(aes(x = lon, y = lat, fill = pred_present_round)) +
    borders(fill = "grey50", colour = "grey90", size = 0.2) +
    geom_richtext(data = sq_area_labels, size = 6, hjust = 1, show.legend = F,
                  aes(x = -50, y = 51, label = area_pres_label, fill = area_pres_perc)) +
    scale_fill_gradientn("Cover (%)", colours = RColorBrewer::brewer.pal(9, "BuGn")[c(3,5,7,9)],
                         limits = c(0, 60), breaks = c(0, 20, 40, 60), guide = "legend") +
    # scale_fill_distiller("Cover (%)", palette = "BuGn", direction = 1, #low = "springgreen1", high = "springgreen4", 
    # limits = c(0, 70), breaks = c(0, 20, 40, 60), guide = "legend") +
    scale_y_continuous(breaks = c(60, 70), labels = c("60°N", "70°N")) +
    scale_x_continuous(breaks = c(-80, -60), labels = c("80°W", "60°W")) +
    coord_quickmap(xlim = c(bbox_arctic[1], bbox_arctic[2]),
                   ylim = c(bbox_arctic[3], bbox_arctic[4]), expand = F) +
    labs(x = NULL, y = NULL, title = paste0(sps_title)) +
    theme(legend.position = "bottom",
          axis.line = element_line(colour = NA),
          panel.background = element_rect(fill = "grey100"),
          panel.border = element_rect(colour = "black", fill = NA))
  
  # Correct title as necessary
  if(model_choice %in% c("agarum", "alaria")){
    p_present <- p_present + theme(plot.title = element_text(face = "italic"))
  }
  # p_present
  
  # 2050 plot
  p_2050 <- model_join %>% 
    filter(proj_2050 == 1,
           # pred_present_round >= 10,
           depth <= 100 | land_distance <= 50) %>% 
    ggplot() +
    geom_tile(aes(x = lon, y = lat, fill = pred_diff_2050)) +
    borders(fill = "grey50", colour = "grey90", size = 0.2) +
    geom_richtext(data = sq_area_labels, size = 6, hjust = 1, show.legend = F,
                  aes(x = -50, y = 51, label = area_2050_label, fill = area_2050_perc)) +
    # scale_fill_gradient2("Change (%)", low = "red", mid = "grey80", high = "blue",
    scale_fill_gradientn("Change (%)", 
                         colours = c(RColorBrewer::brewer.pal(9, "Reds")[c(9,8,7,6)], "grey80",
                                     RColorBrewer::brewer.pal(9, "Blues")[c(6,7,8,9)]),
                         limits = c(-40, 40), 
                         breaks = c(-40, -30, -20, -10, 0, 10, 20, 30, 40),
                         guide = "legend", na.value = NA) +
    scale_y_continuous(breaks = c(60, 70), labels = c("60°N", "70°N")) +
    scale_x_continuous(breaks = c(-80, -60), labels = c("80°W", "60°W")) +
    coord_quickmap(xlim = c(bbox_arctic[1], bbox_arctic[2]),
                   ylim = c(bbox_arctic[3], bbox_arctic[4]), expand = F) +
    labs(x = NULL, y = NULL, title = "2050 - present") +
    theme(legend.position = "bottom",
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          panel.background = element_rect(fill = "grey100"),
          panel.border = element_rect(colour = "black", fill = NA))
  # p_2050
  
  # 2050 plot
  p_2100 <- model_join %>% 
    filter(proj_2100 == 1,
           # pred_present_round >= 10,
           depth <= 100 | land_distance <= 50) %>% 
    ggplot() +
    geom_tile(aes(x = lon, y = lat, fill = pred_diff_2100)) +
    borders(fill = "grey50", colour = "grey90", size = 0.2) +
    geom_richtext(data = sq_area_labels, size = 6, hjust = 1, show.legend = F,
                  aes(x = -50, y = 51, label = area_2100_label, fill = area_2100_perc)) +
    # scale_fill_gradient2("Change (%)", low = "red", mid = "grey80", high = "blue",
    scale_fill_gradientn("Change (%)", 
                         colours = c(RColorBrewer::brewer.pal(9, "Reds")[c(9,8,7,6)], "grey80",
                                     RColorBrewer::brewer.pal(9, "Blues")[c(6,7,8,9)]),
                         limits = c(-40, 40), 
                         breaks = c(-40, -30, -20, -10, 0, 10, 20, 30, 40), 
                         guide = "legend", na.value = NA) +
    scale_y_continuous(breaks = c(60, 70), labels = c("60°N", "70°N")) +
    scale_x_continuous(breaks = c(-80, -60), labels = c("80°W", "60°W")) +
    coord_quickmap(xlim = c(bbox_arctic[1], bbox_arctic[2]),
                   ylim = c(bbox_arctic[3], bbox_arctic[4]), expand = F) +
    labs(x = NULL, y = NULL, title = "2100 - present") +
    theme(legend.position = "bottom",
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          # panel.grid = element_line(colour = "black"),
          panel.background = element_rect(fill = "grey100"),
          panel.border = element_rect(colour = "black", fill = NA))
  # p_2100
  
  # Combine and exit
  p_ALL <- cowplot::plot_grid(
    p_present + theme(legend.position = "none"),
    p_2050 + theme(legend.position = "none"),
    p_2100 + theme(legend.position = "none"),
    ncol = 3,
    align = "v")
  if(add_legend){
    p_ALL <- cowplot::plot_grid(
      cowplot::get_legend(p_present),
      ggplot() + theme_void(),
      cowplot::get_legend(p_2100),
      ncol = 3, rel_widths = c(1, 0, 1.5))
  }
  return(p_ALL)
}

# Create and save
model_compare_lam <- model_compare_plot("laminariales")
model_compare_agarum <- model_compare_plot("agarum")
model_compare_alaria <- model_compare_plot("alaria")
model_compare_legend <- model_compare_plot("laminariales", add_legend = T)

# Combine into one mecha-figure
fig_4 <- ggpubr::ggarrange(model_compare_agarum, model_compare_alaria, model_compare_lam, model_compare_legend,
                           ncol = 1, labels = c("A)", "B)", "C)", ""), heights = c(1, 1, 1, 0.15))
ggsave("figures/fig_4.png", fig_4, width = 7, height = 11)

# Agarum only for demo/talk
fig_4_agarum <- ggpubr::ggarrange(model_compare_agarum, model_compare_legend, ncol = 1, heights = c(1, 0.15))
ggsave("graph/fig_4_agarum.png", fig_4_agarum, width = 7, height = 4)
ggsave("talk/figure/fig_4_agarum.png", fig_4_agarum, width = 7, height = 4)

# Alaria only for talk
fig_4_alaria <- ggpubr::ggarrange(model_compare_alaria, model_compare_legend, ncol = 1, heights = c(1, 0.15))
ggsave("talk/figure/fig_4_alaria.png", fig_4_alaria, width = 7, height = 4)

# Laminariales only for talk
fig_4_lam <- ggpubr::ggarrange(model_compare_lam, model_compare_legend, ncol = 1, heights = c(1, 0.15))
ggsave("talk/figure/fig_4_lam.png", fig_4_lam, width = 7, height = 4)


# Figure 5 ----------------------------------------------------------------

# The model values from the ensembles; ROC vs AUC
model_stats_plot <- function(sps_choice){
  
  # Load chosen biomod_model and print evaluation scores
  biomod_model <- loadRData(paste0(sps_choice,"/",sps_choice,".",sps_choice,".models.out"))
  
  # Create full species name
  if(sps_choice == "Lsol"){
    sps_title <- "Laminaria solidungula"
  } else if(sps_choice == "Slat"){
    sps_title <- "Saccharina latissima"
  } else if(sps_choice == "Acla"){
    sps_title <- "Agarum clathratum"
  } else if(sps_choice == "Aesc"){
    sps_title <- "Alaria esculenta"
  } else{
    stop("*sad robot noises*")
  }  
  
  # Extract raw results for geom_point()
  scores <- biomod2:::get_evaluations(biomod_model, as.data.frame = T) %>% 
    separate(Model.name, into = c("model", "run", "PA"), sep = "_") %>% 
    dplyr::select(model:Testing.data) %>% 
    pivot_wider(names_from = Eval.metric, values_from = Testing.data) %>% 
    mutate(model = ifelse(model == "MAXENT.Phillips", "MAXENT", model))
  
  # Index of possible combos to catch 0 model psases
  scores_index <- data.frame(model = c("ANN", "GAM", "GLM", "MAXENT", "RF"),
                             total_count = 25)
  
  # Count of models that passed
  scores_passed <- scores %>% 
    group_by(model) %>% 
    mutate(total_count = n()) %>% 
    filter(TSS >= 0.7) %>% 
    group_by(model, total_count) %>% 
    summarise(pass_count = n(), .groups = "drop") %>% 
    right_join(scores_index, by = c("model", "total_count")) %>% 
    replace(is.na(.), 0) %>% 
    mutate(pass_label = paste0(pass_count,"/",total_count)) %>% 
    dplyr::arrange(model)
  
  # Model evaluation by algorithm
  # model_res <- biomod2::models_scores_graph(biomod_model, by = "models", metrics = c('ROC','TSS')) + 
    # geom_point(data = scores, aes(x = ROC, y = TSS, colour = model)) +
  model_res <- ggplot(scores) +
    geom_label(data = scores_passed, label.size = 1.5, size = 6, show.legend = F,
               aes(x = model, y = 0.95, label = pass_label, colour = model)) +
    geom_boxplot(aes(x = model, y = TSS, fill = model)) +
    geom_hline(aes(yintercept = 0.7), colour = "hotpink", size = 2) +
    labs(x = NULL, y = "TSS", title = sps_title, fill = "Model") +
    scale_y_continuous(limits = c(0, 1.0), expand = c(0, 0)) +
    scale_fill_brewer(palette = "Set1", aesthetics = c("colour", "fill")) +
    # guides(fill = guide_legend(override.aes = list(shape = 15))) +
    # coord_cartesian(xlim = c(0.6,1), ylim = c(0.3,1), expand = F) +
    theme(plot.title = element_text(face = "italic"),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA))
  # model_res

  return(model_res)
}

# Create plots
model_stats_Acla <- model_stats_plot("Acla")
model_stats_Aesc <- model_stats_plot("Aesc")
model_stats_Lsol <- model_stats_plot("Lsol")
model_stats_Slat <- model_stats_plot("Slat")

# Combine and save
fig_5 <- ggpubr::ggarrange(model_stats_Acla, model_stats_Aesc, model_stats_Lsol, model_stats_Slat, 
                           ncol = 2, nrow = 2, common.legend = TRUE, legend = "bottom", 
                           labels = c("A)", "B)", "C)", "D)"))
ggsave("figures/fig_5.png", fig_5, width = 10, height = 10.5)

# Agarum only for demo/talk
ggsave("talk/figure/fig_5_agarum.png", model_stats_Acla, width = 7, height = 5)

# Alaria only for talk
ggsave("talk/figure/fig_5_alaria.png", model_stats_Aesc, width = 7, height = 5)

# Laminaria only for talk
ggsave("talk/figure/fig_5_laminaria.png", model_stats_Lsol, width = 7, height = 5)

# Saccharina only for talk
ggsave("talk/figure/fig_5_saccharina.png", model_stats_Slat, width = 7, height = 5)


# Figure 6 ----------------------------------------------------------------

# Random forest model results
load("data/best_rf_kelpcover.RData")
load("data/best_rf_laminariales.RData")
load("data/best_rf_agarum.RData")
load("data/best_rf_alaria.RData")

# Confidence in RF output
# Function for creating figure showing confidence intervals of prediction accuracy
conf_plot_RF <- function(df, sps_choice){
  
  # Create full species name and load ensemble data
  if(sps_choice == "laminariales"){
    sps_title <- "Laminariales spp."
  } else if(sps_choice == "agarum"){
    sps_title <- "Agarum clathratum"
  } else if(sps_choice == "alaria"){
    sps_title <- "Alaria esculenta"
  } else{
    stop("*sad robot noises*")
  }
  
  # 90 CI around predictions per step
  conf_acc <- df %>% 
    filter(portion == "validate") %>% 
    group_by(original) %>% 
    mutate(accuracy = round(accuracy, -1)) %>% 
    summarise(count = n(),
              q05 = quantile(accuracy, 0.05),
              q25 = quantile(accuracy, 0.25),
              q50 = median(accuracy),
              mean = mean(accuracy),
              q75 = quantile(accuracy, 0.75),
              q95 = quantile(accuracy, 0.95), .groups = "drop")
  
  conf_mean <- df %>% 
    filter(portion == "validate") %>% 
    group_by(model_id) %>% 
    mutate(accuracy = round(accuracy)) %>% 
    summarise(mean_acc = mean(abs(accuracy)),
              sd_acc = sd(abs(accuracy)),
              r_acc = cor(x = original, y = pred), .groups = "drop")
  
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
              sd_acc = sd(accuracy), .groups = "drop")
  
  # Visualise
  conf_plot <- ggplot(conf_acc, aes(x = original, y = mean)) +
    geom_crossbar(aes(y = 0, ymin = q05, ymax = q95),
                  fatten = 0, fill = "grey70", colour = NA, width = 10) +
    geom_crossbar(aes(ymin = q25, ymax = q75),
                  fatten = 0, fill = "grey50", width = 10) +
    geom_crossbar(aes(ymin = q50, ymax = q50), size = 1.5,
                  fatten = 0, fill = NA, colour = "black", width = 10) +
    geom_hline(yintercept = 0, size = 2, alpha = 0.7, colour = "red") +
    geom_hline(yintercept = c(-50, 50), size = 0.5, alpha = 0.7, colour = "purple", linetype = "dashed") +
    # geom_hline(yintercept = -50, size = 1, alpha = 1, colour = "purple", linetype = "dashed") +
    geom_label(aes(label = scales::comma(count, accuracy = 1), y = q95), size = 3) +
    # geom_segment(data = conf_best, aes(xend = original, y = mean_acc, yend = 0), 
    # colour = "purple", size = 1.2, alpha = 0.8) +
    # geom_point(data = conf_best, aes(y = mean_acc), colour = "purple", size = 3, alpha = 0.8) +
    geom_label(data = conf_mean_label, 
               aes(x = 75, y = 75, label = paste0("Mean accuracy: ",mean_acc,"±",sd_acc,"; r = ",r_acc))) +
    # geom_label(data = conf_best_label, colour = "purple",
    # aes(x = 75, y = 60, label = paste0("Best accuracy: ",mean_acc,"±",sd_acc,"; r = ",r_acc))) +
    scale_y_continuous(limits = c(-100, 100), breaks = c(-50, 0 ,50), labels = c(-50, 0 ,50), 
                       minor_breaks = seq(-90, 90, 10), expand = c(0, 0)) +
    scale_x_continuous(breaks = c(0, 20, 40, 60, 80, 100), expand = c(0, 0)) +
    labs(y = "Difference from observed", x = "Observed value (% cover)", title = sps_title) +
    theme(panel.border = element_rect(fill = NA, colour = "black"))
  
  # Correct title as necessary
  if(sps_choice %in% c("agarum", "alaria")){
    conf_plot <- conf_plot + theme(plot.title = element_text(face = "italic"))
  }
  conf_plot
  return(conf_plot)
}

# Create the plots
cover_lam <- conf_plot_RF(best_rf_laminariales$accuracy_reg, "laminariales")
cover_agarum <- conf_plot_RF(best_rf_agarum$accuracy_reg, "agarum")
cover_alaria <- conf_plot_RF(best_rf_alaria$accuracy_reg, "alaria")

# Steek'em
fig_6 <- ggpubr::ggarrange(cover_agarum, cover_alaria, cover_lam, ncol = 1, 
                           labels = c("A)", "B)", "C)"))
ggsave("figures/fig_6.png", fig_6, width = 7, height = 12)

# Agarum only for demo/talk
ggsave("talk/figure/fig_6_agarum.png", cover_agarum, width = 7, height = 4)

# Alaria only for talk
ggsave("talk/figure/fig_6_alaria.png", cover_alaria, width = 7, height = 4)

# Laminariales only for talk
ggsave("talk/figure/fig_6_lam.png", cover_lam, width = 7, height = 4)


# Figure S1 ---------------------------------------------------------------
# The random forest model results
# NB: This requires functions from the Figure 3 section

# testers...
# kelp_choice <- "laminariales"
# kelp_choice <- "agarum"
rf_plot <- function(kelp_choice, add_legend = F){
  
  # Create full species name
  if(kelp_choice == "laminariales"){
    sps_title <- "Laminariales spp."
  } else if(kelp_choice == "agarum"){
    sps_title <- "Agarum clathratum"
  } else if(kelp_choice == "alaria"){
    sps_title <- "Alaria esculenta"
  } else{
    stop("*sad robot noises*")
  }
  
  # Load data
  best_rf <- loadRData(paste0("data/best_rf_",kelp_choice,".RData"))
  
  # Calculate differences
  project_diff <- best_rf$project_multi %>% 
    mutate(pred_present_round = plyr::round_any(pred_present_mean, 20),
           pred_diff_2050 = plyr::round_any(pred_2050_mean - pred_present_mean, 10),
           pred_diff_2100 = plyr::round_any(pred_2100_mean - pred_present_mean, 10))
  
  # Prep site data
  # All of the BO variables matched to sites
  kelp_mean <- adf %>% 
    dplyr::select(Campaign, site, kelp.cover, Laminariales, Agarum, Alaria) %>% 
    left_join(study_sites, by = c("Campaign", "site")) %>%
    mutate(kelp.cover = ifelse(kelp.cover > 100, 100, kelp.cover)) %>% # Correct values over 100
    # Round values to nearest 20 percent and remove uppercase characters
    group_by(Campaign, site, lon, lat) %>% 
    summarise(kelp.cover = plyr::round_any(mean(kelp.cover, na.rm = T), 20),
              laminariales = plyr::round_any(mean(Laminariales, na.rm = T), 20),
              agarum = plyr::round_any(mean(Agarum, na.rm = T), 20),
              alaria = plyr::round_any(mean(Alaria, na.rm = T), 20), .groups = "drop")
  # Add large dots to show real cover %
  
  # Calculate sq area coverage per era
  sq_area_labels <- project_diff %>% 
    filter(land_distance <= 50 | depth <= 100) %>% 
    summarise(area_pres_mean = mean(pred_present_round, na.rm = T),
              area_2050_mean = mean(pred_diff_2050, na.rm = T),
              area_2100_mean = mean(pred_diff_2100, na.rm = T)) %>% 
    mutate(area_pres_label = paste0(round(area_pres_mean),"%"),
           area_2050_label = ifelse(area_2050_mean > 0 , 
                                    paste0("+",round(area_2050_mean),"%"),
                                    paste0(round(area_2050_mean),"%")),
           area_2100_label = ifelse(area_2100_mean > 0 , 
                                    paste0("+",round(area_2100_mean),"%"),
                                    paste0(round(area_2100_mean),"%")))
  
  # Present plot
  p_present <-  ggplot() +
    # geom_tile(data = df, # No filter
    geom_tile(data = filter(project_diff,
                            pred_present_round >= 10,
                            depth <= 100 | land_distance <= 50),
              aes(x = lon, y = lat, fill = pred_present_round)) +
    borders(fill = "grey50", colour = "grey90", size = 0.2) +
    geom_richtext(data = sq_area_labels, size = 6, hjust = 1, show.legend = F,
                  aes(x = -50, y = 51, label = area_pres_label, fill = area_pres_mean)) +
    # scale_fill_viridis_c(paste0("cover (%)"), limits = c(10, 70)) +
    # scale_fill_distiller(palette = "Greens", direction = 1, limits = c(10, 70)) +
    scale_fill_gradient("Cover (%)", low = "grey90", high = "grey30", aesthetics = c("colour", "fill"),
                        limits = c(0, 60), breaks = c(20, 40, 60), guide = "legend") +
    geom_point(data = kelp_mean, aes_string(x = "lon", y = "lat"), colour = "red", size = 2) +
    geom_point(data = kelp_mean, aes_string(x = "lon", y = "lat", colour = kelp_choice), size = 1.9) +
    scale_y_continuous(breaks = c(60, 70), labels = c("60°N", "70°N")) +
    scale_x_continuous(breaks = c(-80, -60), labels = c("80°W", "60°W")) +
    coord_quickmap(xlim = c(bbox_arctic[1], bbox_arctic[2]),
                   ylim = c(bbox_arctic[3], bbox_arctic[4]), expand = F) +
    labs(x = NULL, y = NULL, title = paste0(sps_title)) +
    # theme_bw() +
    theme(legend.position = "bottom",
          panel.background = element_rect(fill = "grey95"),
          panel.border = element_rect(colour = "black", fill = NA))
  if(kelp_choice %in% c("agarum", "alaria")){
    p_present <- p_present + theme(plot.title = element_text(face = "italic"))
  }
  p_present
  
  # 2050 plot
  p_2050 <- ggplot() +
    # geom_tile(data = df, # No filter
    geom_tile(data = filter(project_diff, 
                            pred_present_round >= 10,
                            depth <= 100 | land_distance <= 50),
              aes(x = lon, y = lat, fill = pred_diff_2050)) +
    borders(fill = "grey50", colour = "grey90", size = 0.2) +
    geom_richtext(data = sq_area_labels, size = 6, hjust = 1, show.legend = F,
                  aes(x = -50, y = 51, label = area_2050_label, fill = area_2050_mean)) +
    scale_fill_gradientn("Change (%)", 
                         colours = c(RColorBrewer::brewer.pal(9, "Reds")[c(9,8,7,6)], "grey80",
                                     RColorBrewer::brewer.pal(9, "Blues")[c(6,7,8,9)]),
                         limits = c(-40, 40), 
                         breaks = c(-40, -30, -20, -10, 0, 10, 20, 30, 40),
                         guide = "legend", na.value = NA) +
    # scale_fill_discrete("2050 - present (%)") +
    scale_y_continuous(breaks = c(60, 70), labels = c("60°N", "70°N")) +
    scale_x_continuous(breaks = c(-80, -60), labels = c("80°W", "60°W")) +
    coord_quickmap(xlim = c(bbox_arctic[1], bbox_arctic[2]),
                   ylim = c(bbox_arctic[3], bbox_arctic[4]), expand = F) +
    labs(x = NULL, y = NULL, title = "2050 - present") +
    # theme_bw() +
    theme(legend.position = "bottom",
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          # panel.grid = element_line(colour = "black"),
          panel.background = element_rect(fill = "grey95"),
          panel.border = element_rect(colour = "black", fill = NA))
  # p_2050
  
  # 2050 plot
  p_2100 <-  ggplot() +
    # geom_tile(data = df, # No filter
    geom_tile(data = filter(project_diff,
                            pred_present_round >= 10,
                            depth <= 100 | land_distance <= 50),
              aes(x = lon, y = lat, fill = pred_diff_2100)) +
    borders(fill = "grey50", colour = "grey90", size = 0.2) +
    geom_richtext(data = sq_area_labels, size = 6, hjust = 1, show.legend = F,
                  aes(x = -50, y = 51, label = area_2100_label, fill = area_2100_mean)) +
    scale_fill_gradientn("Change (%)", 
                         colours = c(RColorBrewer::brewer.pal(9, "Reds")[c(9,8,7,6)], "grey80",
                                     RColorBrewer::brewer.pal(9, "Blues")[c(6,7,8,9)]),
                         limits = c(-40, 40), 
                         breaks = c(-40, -30, -20, -10, 0, 10, 20, 30, 40),
                         guide = "legend", na.value = NA) +
    # scale_fill_discrete("2100 - present (%)") +
    scale_y_continuous(breaks = c(60, 70), labels = c("60°N", "70°N")) +
    scale_x_continuous(breaks = c(-80, -60), labels = c("80°W", "60°W")) +
    coord_quickmap(xlim = c(bbox_arctic[1], bbox_arctic[2]),
                   ylim = c(bbox_arctic[3], bbox_arctic[4]), expand = F) +
    labs(x = NULL, y = NULL, title = "2100 - present") +
    # theme_bw() +
    theme(legend.position = "bottom",
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          panel.background = element_rect(fill = "grey95"),
          panel.border = element_rect(colour = "black", fill = NA))
  # p_2100
  
  # Combine and exit
  p_ALL <- cowplot::plot_grid(
    p_present + theme(legend.position = "none"),
    p_2050 + theme(legend.position = "none"),
    p_2100 + theme(legend.position = "none"),
    ncol = 3,
    align = "v")
  if(add_legend){
    p_ALL <- cowplot::plot_grid(
      cowplot::get_legend(p_present),
      ggplot() + theme_void(),
      cowplot::get_legend(p_2100),
      ncol = 3, rel_widths = c(1, 0, 1.5))
  }
  return(p_ALL)
}

rf_lam <- rf_plot("laminariales")
rf_agarum <- rf_plot("agarum")
rf_alaria <- rf_plot("alaria")
rf_legend <- rf_plot("laminariales", add_legend = T)

# Combine into one mecha-figure
fig_S1 <- ggpubr::ggarrange(rf_agarum, rf_alaria, rf_lam, rf_legend,
                            ncol = 1, labels = c("A)", "B)", "C)", ""), heights = c(1, 1, 1, 0.15))
ggsave("figures/fig_S1.png", fig_S1, width = 7, height = 11)


# Figure S2 ---------------------------------------------------------------

# Variance plot for ensemble model results


# Figure S3 ---------------------------------------------------------------

# The variance per pixel in the random forest model results
# testers...
# kelp_choice <- "laminariales"
# kelp_choice <- "agarum"
rf_var_plot <- function(kelp_choice, add_legend = F){
  
  # Create full species name
  if(kelp_choice == "laminariales"){
    sps_title <- "Laminariales spp."
  } else if(kelp_choice == "agarum"){
    sps_title <- "Agarum clathratum"
  } else if(kelp_choice == "alaria"){
    sps_title <- "Alaria esculenta"
  } else{
    stop("*sad robot noises*")
  }
  
  # Load data
  best_rf <- loadRData(paste0("data/best_rf_",kelp_choice,".RData"))
  
  # Calculate differences
  project_diff <- best_rf$project_multi
  
  # Plot
  p_present <-  ggplot() +
    geom_tile(data = filter(project_diff,
                            # proj_pres == 1,
                            # pred_present_round >= 10,
                            depth <= 100 | land_distance <= 50),
              aes(x = lon, y = lat, fill = pred_present_sd)) +
    scale_fill_distiller("Cover (SD %)", palette = "Oranges", direction = 1, limits = c(0, 8)) +
    # scale_fill_gradientn("Cover (%)", colours = RColorBrewer::brewer.pal(9, "BuGn")[c(5,6,7,8)],
    # limits = c(0, 70), breaks = c(0, 20, 40, 60), guide = "legend") +
    # scale_fill_distiller("Cover (%)", palette = "BuGn", direction = 1, #low = "springgreen1", high = "springgreen4", 
    # limits = c(0, 70), breaks = c(0, 20, 40, 60), guide = "legend") +
    borders(fill = "grey50", colour = "grey90", size = 0.2) +
    scale_y_continuous(breaks = c(60, 70), labels = c("60°N", "70°N")) +
    scale_x_continuous(breaks = c(-80, -60), labels = c("80°W", "60°W")) +
    coord_quickmap(xlim = c(bbox_arctic[1], bbox_arctic[2]),
                   ylim = c(bbox_arctic[3], bbox_arctic[4]), expand = F) +
    labs(x = NULL, y = NULL, title = paste0(sps_title)) +
    theme(legend.position = "bottom",
          axis.line = element_line(colour = NA),
          panel.background = element_rect(fill = "grey100"),
          panel.border = element_rect(colour = "black", fill = NA))
  
  # Correct title as necessary
  if(kelp_choice %in% c("agarum", "alaria")){
    p_present <- p_present + theme(plot.title = element_text(face = "italic"))
  }
  # p_present
  
  # 2050 plot
  p_2050 <- ggplot() +
    geom_tile(data = filter(project_diff,
                            # proj_2050 == 1,
                            # pred_present_round >= 10,
                            depth <= 100 | land_distance <= 50),
              aes(x = lon, y = lat, fill = pred_2050_sd)) +
    scale_fill_distiller("Cover (SD %)", palette = "Oranges", direction = 1, limits = c(0, 8)) +
    borders(fill = "grey50", colour = "grey90", size = 0.2) +
    scale_y_continuous(breaks = c(60, 70), labels = c("60°N", "70°N")) +
    scale_x_continuous(breaks = c(-80, -60), labels = c("80°W", "60°W")) +
    coord_quickmap(xlim = c(bbox_arctic[1], bbox_arctic[2]),
                   ylim = c(bbox_arctic[3], bbox_arctic[4]), expand = F) +
    labs(x = NULL, y = NULL, title = "2050 - present") +
    theme(legend.position = "bottom",
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          panel.background = element_rect(fill = "grey100"),
          panel.border = element_rect(colour = "black", fill = NA))
  # p_2050
  
  # 2050 plot
  p_2100 <-  ggplot() +
    geom_tile(data = filter(project_diff,
                            # proj_2100 == 1,
                            # pred_present_round >= 10,
                            depth <= 100 | land_distance <= 50),
              aes(x = lon, y = lat, fill = pred_2100_sd)) +
    scale_fill_distiller("Cover (SD %)", palette = "Oranges", direction = 1, limits = c(0, 8)) +
    borders(fill = "grey50", colour = "grey90", size = 0.2) +
    scale_y_continuous(breaks = c(60, 70), labels = c("60°N", "70°N")) +
    scale_x_continuous(breaks = c(-80, -60), labels = c("80°W", "60°W")) +
    coord_quickmap(xlim = c(bbox_arctic[1], bbox_arctic[2]),
                   ylim = c(bbox_arctic[3], bbox_arctic[4]), expand = F) +
    labs(x = NULL, y = NULL, title = "2100 - present") +
    theme(legend.position = "bottom",
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          panel.background = element_rect(fill = "grey100"),
          panel.border = element_rect(colour = "black", fill = NA))
  # p_2100
  
  # Combine and exit
  p_ALL <- cowplot::plot_grid(
    p_present + theme(legend.position = "none"),
    p_2050 + theme(legend.position = "none"),
    p_2100 + theme(legend.position = "none"),
    ncol = 3,
    align = "v")
  if(add_legend){
    p_ALL <- cowplot::plot_grid(
      ggplot() + theme_void(),
      cowplot::get_legend(p_present),
      ggplot() + theme_void(),
      # cowplot::get_legend(p_2100),
      ncol = 3, rel_widths = c(1, 0, 1.5))
  }
  return(p_ALL)
}

# Create figures
rf_var_lam <- rf_var_plot("laminariales")
rf_var_agarum <- rf_var_plot("agarum")
rf_var_alaria <- rf_var_plot("alaria")
rf_var_legend <- rf_var_plot("laminariales", add_legend = T)

# Combine into one mecha-figure
rf_var_ALL <- ggpubr::ggarrange(rf_var_agarum, rf_var_alaria, rf_var_lam, rf_var_legend,
                                ncol = 1, labels = c("A)", "B)", "C)", ""), heights = c(1, 1, 1, 0.15))
ggsave("figures/fig_S3.png", rf_var_ALL, width = 7, height = 11)

