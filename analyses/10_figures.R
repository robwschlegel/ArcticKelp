# analyses/10_figures
# This script houses the code used to create the final figures for the manuscript

# TODO: Slightly expand the edge of the bbox to allow the reader to be certain of the boundary of Hudson Bay


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
library(biomod2)

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
  annotation_spatial(Arctic_poly, fill = "gold", alpha = 0.2) +
  annotation_spatial(bbox_spatial, fill = "darkgreen", alpha = 0.2) +
  # annotation_spatial(MEOW, aes(colour = ECOREGION), fill = NA) +
  # NB: This overplotting of Slat needs to be addressed
  # This may be due to the points needing to be converted to a spatial points file first
  geom_spatial_point(data = sps_data, crs = 4326, shape = 21,
                     aes(x = lon, y = lat), fill = "hotpink", colour = "black") +
  geom_spatial_point(data = study_sites, crs = 4326, shape = 22,
                     aes(x = lon, y = lat), fill = "purple", colour = "black") +
  # guides(fill = guide_legend(nrow = 3, byrow = TRUE)) +
  theme(legend.position = "bottom")
# fig_1a
ggsave("figures/fig_1a.png", fig_1a, height = 8, width = 6)

# Add the names of the main regions: Hudson Bay, Hudson Strait, Lancaster Sound
label_df <- data.frame(lon = c(-85, -75.17, -82.92, -63.834153, -79.0, -90.0706224, -82.7015637, -58.16),
                       lat = c(60, 63.0, 74.425, 72.172485, 67.5202107, 70.5938471, 53.0693041, 66.0),
                       loc = c("Hudson Bay", "Hudson Strait", "Lancaster Sound", "Baffin Bay", 
                               "Foxe Basin", "Gulf of \nBoothia", "James Bay", "Davis Strait"))

# ArcticKelp campaign map
fig_1b <- ggplot() +
  # geom_tile(aes(fill = presence)) +
  borders(fill = "grey30", colour = "black") +
  geom_label(data = label_df, aes(x = lon, y = lat, label = loc), alpha = 0.9) +
  geom_point(data = study_sites, shape = 22, colour = "black", 
             fill = "purple", size = 2,
             aes(x = lon, y = lat)) +
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
# fig_1b

# Combine and save
fig_1 <- ggpubr::ggarrange(fig_1a, fig_1b, ncol = 2, labels = c("A)", "B)"), widths = c(1.5, 1))
ggsave("figures/fig_1.png", fig_1, height = 6, width = 10)


# Table 1 -----------------------------------------------------------------
# Importance of variables by kelps and models


# Figure 2 ----------------------------------------------------------------
# Importance of variables per taxonomic group
# Jesi made this figure with image editing software


# Figure 3 ----------------------------------------------------------------
# The ensemble model results

# TODO: Create a table and versions of figures of current projections that show the difference in projections/surface area km2
# Depending on which depth and land distance filter is used. This should be discussed in the paper.
# Consider creating a sub panel of a specific area where the choice of filter has a big effect
# For the area calculations a weighted average could be used by depth to determine the contribution to the total
# Consider filtering by depth south of 75°N, and by land distance above that
# Focus on the 75°N line as the most problematic shift to smaller pixels

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
  
  # Calculate square are in kilometres directly
  # This doesn't work well close to the poles
  # sum(biomod_project_present[] == TRUE, na.rm = T) * raster::res(biomod_project_present)[1]^2
  # area(biomod_project_present,  na.rm=TRUE, weights=FALSE)
  # tapply(area(biomod_project_present), biomod_project_present[], sum)
  
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
    # filter(depth <= 30) %>%
    # filter(depth <= 50) %>%
    # filter(land_distance <= 10) %>%
    filter(land_distance <= 15 | depth <= 50) %>%
    ggplot(aes(x = lon, y = lat)) +
    borders(fill = "grey50", colour = "grey90", size = 0.2) +
    geom_tile(aes_string(fill = paste0("change_",year_label))) +
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
    filter(depth <= 30) %>%
    # filter(land_distance <= 50 | depth <= 100) %>%
    summarise(area_pres = round(sum(sq_area*proj_pres, na.rm = T), -3),
              area_2050 = round(sum(sq_area*proj_2050, na.rm = T), -3),
              area_2100 = round(sum(sq_area*proj_2100, na.rm = T), -3))
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
    # filter(depth <= 30,
    # filter(depth <= 50,
    # filter(land_distance <= 10,
    filter(land_distance <= 15 | depth <= 50,
           proj_pres == 1) %>%
    mutate(proj_pres = "") %>%
    ggplot(aes(x = lon, y = lat)) +
    borders(fill = "grey50", colour = "grey90", size = 0.2) +
    geom_tile(aes(fill = proj_pres)) +
    geom_point(data = sps_points, shape = 21, colour = "black", fill = "hotpink", size = 0.5) +
    geom_richtext(data = sq_area_labels, size = 4, hjust = 1,
                  aes(x = -50, y = 51, label = pres_text)) +
    # annotate("label", x = -58, y = 78, label = paste(pres_text, "^2", sep = ""), parse = T) +
    scale_y_continuous(breaks = c(60, 70), labels = c("60°N", "70°N")) +
    scale_x_continuous(breaks = c(-80, -60), labels = c("80°W", "60°W")) +
    coord_quickmap(xlim = c(bbox_arctic[1], bbox_arctic[2]),
                   ylim = c(bbox_arctic[3], bbox_arctic[4]), expand = F) +
    scale_fill_manual("Suitable", values = c("forestgreen")) +
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
  filter(presence == T, depth <= 30) %>% 
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
    sps_title <- "Laminareacea"
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
    sps_title <- "Laminareacea"
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
    sps_title <- "Laminareacea."
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
    sps_title <- "Laminareacea"
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


# Figure S4 ---------------------------------------------------------------
# Response curves for ensemble models

# Function for extracting species response curve data
# sps_choice <- sps_names[1]
sps_response_data <- function(biomod_model, sp_dat){
  # Model scores
  model_scores <- t(data.frame(get_evaluations(biomod_model))) %>% 
    data.frame() %>% 
    mutate(idx = rownames(.),
           idx = str_replace(idx, "MAXENT.Phillips", "MAXENT_Phillips"),
           idx = str_replace(idx, "Testing.data", "Testing_data")) %>% 
    separate(idx, c("type", "model", "run", "PA"), sep = "[.]", remove = T) %>% 
    mutate(model = case_when(model == "MAXENT_Phillips" ~ "MAXENT.Phillips", TRUE ~ model)) %>% 
    filter(type == "Testing_data") %>%
    `row.names<-`(NULL) %>% 
    dplyr::select(type, model, run, PA, everything())
  
  # Clean up and exit
  sp_res <- sp_dat %>% 
    separate(pred.name, c("species", "PA", "run", "model"), sep = "_", remove = T) %>% 
    left_join(model_scores, by = c("model", "run", "PA")) %>% 
    mutate(expl.name = case_when(expl.name == "BO21_templtmax_ss" ~ "__B)__    Temperature (°C)",
                                 expl.name == "BO21_salinitymean_ss" ~ "__C)__    Salinity (PSS)",
                                 expl.name == "BO21_icethickmean_ss" ~ "__A)__    Ice thickness (m)",
                                 expl.name == "BO21_curvelmean_bdmax" ~ "__F)__    Current velocity (m-1)",
                                 expl.name == "BO21_ironmean_bdmax" ~ "__D)__    Iron (umol.m-3)",
                                 expl.name == "BO21_phosphatemean_bdmax" ~ "__E)__    Phosphate (umol.m-3)"),
           run_PA_model = paste0(run,"_",PA,"_",model))
  return(sp_res)
}

## Load model data per species
# NB: Due to BIOMOD2 internal structure, much of this needs to be run in the global environment...
# Acla
biomod_model_Acla <- loadRData(paste0(sps_names[1],"/",sps_names[1],".",sps_names[1],".models.out"))
sp_name_Acla <- BIOMOD_LoadModels(biomod_model_Acla, models = c('MAXENT.Phillips', 'GLM', 'ANN', 'RF', 'GAM'))
sp_dat_Acla <- response.plot2(models  = sp_name_Acla,
                              Data = get_formal_data(biomod_model_Acla, 'expl.var'),
                              show.variables = get_formal_data(biomod_model_Acla, 'expl.var.names'),
                              data_species = get_formal_data(biomod_model_Acla, 'resp.var'),
                              do.bivariate = FALSE, fixed.var.metric = 'median', legend = FALSE, plot = FALSE)
curve_data_Acla <- sps_response_data(biomod_model_Acla, sp_dat_Acla)
rm(list = grep("Acla_",names(.GlobalEnv),value = TRUE)); gc()
# Aesc
biomod_model_Aesc <- loadRData(paste0(sps_names[2],"/",sps_names[2],".",sps_names[2],".models.out"))
sp_name_Aesc <- BIOMOD_LoadModels(biomod_model_Aesc, models = c('MAXENT.Phillips', 'GLM', 'ANN', 'RF', 'GAM'))
sp_dat_Aesc <- response.plot2(models  = sp_name_Aesc,
                              Data = get_formal_data(biomod_model_Aesc, 'expl.var'),
                              show.variables = get_formal_data(biomod_model_Aesc, 'expl.var.names'),
                              data_species = get_formal_data(biomod_model_Aesc, 'resp.var'),
                              do.bivariate = FALSE, fixed.var.metric = 'median', legend = FALSE, plot = FALSE)
curve_data_Aesc <- sps_response_data(biomod_model_Aesc, sp_dat_Aesc)
rm(list = grep("Aesc_",names(.GlobalEnv),value = TRUE)); gc()
# Lsol
biomod_model_Lsol <- loadRData(paste0(sps_names[4],"/",sps_names[4],".",sps_names[4],".models.out"))
sp_name_Lsol <- BIOMOD_LoadModels(biomod_model_Lsol, models = c('MAXENT.Phillips', 'GLM', 'ANN', 'RF', 'GAM'))
sp_dat_Lsol <- response.plot2(models  = sp_name_Lsol,
                              Data = get_formal_data(biomod_model_Lsol, 'expl.var'),
                              show.variables = get_formal_data(biomod_model_Lsol, 'expl.var.names'),
                              data_species = get_formal_data(biomod_model_Lsol, 'resp.var'),
                              do.bivariate = FALSE, fixed.var.metric = 'median', legend = FALSE, plot = FALSE)
curve_data_Lsol <- sps_response_data(biomod_model_Lsol, sp_dat_Lsol)
rm(list = grep("Lsol_",names(.GlobalEnv),value = TRUE)); gc()
# Slat
biomod_model_Slat <- loadRData(paste0(sps_names[5],"/",sps_names[5],".",sps_names[5],".models.out"))
sp_name_Slat <- BIOMOD_LoadModels(biomod_model_Slat, models = c('MAXENT.Phillips', 'GLM', 'ANN', 'RF', 'GAM'))
# NB: Sometimes this fails on the first go. Just run it again... I don't know why...
sp_dat_Slat <- response.plot2(models  = sp_name_Slat,
                              Data = get_formal_data(biomod_model_Slat, 'expl.var'),
                              show.variables = get_formal_data(biomod_model_Slat, 'expl.var.names'),
                              data_species = get_formal_data(biomod_model_Slat, 'resp.var'),
                              do.bivariate = FALSE, fixed.var.metric = 'median', legend = FALSE, plot = FALSE)
curve_data_Slat <- sps_response_data(biomod_model_Slat, sp_dat_Slat)
rm(list = grep("Slat_",names(.GlobalEnv),value = TRUE)); gc()
# All together now
curve_data_ALL <- rbind(curve_data_Acla, curve_data_Aesc, curve_data_Lsol, curve_data_Slat)

# Plot RC for all species
response_curve_species_mean <-  curve_data_ALL %>%
  filter(TSS >= 0.7) %>% 
  group_by(id, expl.name, species) %>% 
  summarise(expl.val = mean(expl.val, na.rm = T),
            pred.val = mean(pred.val, na.rm = T), .groups = "drop") %>% 
  ggplot(aes(x = expl.val, y = pred.val, colour = species)) +
  geom_line(size = 1) +
  facet_wrap(~expl.name, scales = 'free_x') +#, strip.position = "bottom") + 
  labs(x = NULL, y = 'probability of occurence', colour = 'species') + 
  scale_color_brewer(type = 'qual', palette = 3) +
  # theme_minimal() +
  theme(legend.position = 'bottom',
        panel.border = element_rect(colour = "black", fill = NA),
        strip.text = ggtext::element_markdown(hjust = 0),
        # strip.placement = "outside", 
        strip.background = element_blank())
ggsave("figures/fig_S4.png", response_curve_species_mean, height = 8, width = 12)


# Figure S5 ---------------------------------------------------------------
# Response curves for RF

# Load RF model results
best_rf_laminariales <- loadRData("data/best_rf_laminariales.RData")$project_multi %>% mutate(species = "Lam")
best_rf_agarum <- loadRData("data/best_rf_agarum.RData")$project_multi %>% mutate(species = "Acla")
best_rf_alaria <- loadRData("data/best_rf_alaria.RData")$project_multi %>% mutate(species = "Aesc")

# Merge RF % cover predictions to BO data
curve_rf_ALL <- rbind(best_rf_laminariales, best_rf_agarum, best_rf_alaria) %>% 
  filter(depth <= 30,
         !is.na(pred_present_mean)) %>% 
  left_join(Arctic_BO, by = c("lon", "lat")) %>% 
  dplyr::select(species, pred_present_mean, BO21_templtmax_ss, BO21_salinitymean_ss, 
                BO21_icethickmean_ss, BO21_curvelmean_bdmax, BO21_ironmean_bdmax, BO21_phosphatemean_bdmax) %>% 
  group_by(species, pred_present_mean) %>% 
  # summarise_all(mean, na.rm = T, .groups = "drop") %>% 
  pivot_longer(BO21_templtmax_ss:BO21_phosphatemean_bdmax) %>% 
  mutate(name = case_when(name == "BO21_templtmax_ss" ~ "__B)__    Temperature (°C)",
                          name == "BO21_salinitymean_ss" ~ "__C)__    Salinity (PSS)",
                          name == "BO21_icethickmean_ss" ~ "__A)__    Ice thickness (m)",
                          name == "BO21_curvelmean_bdmax" ~ "__F)__    Current velocity (m-1)",
                          name == "BO21_ironmean_bdmax" ~ "__D)__    Iron (umol.m-3)",
                          name == "BO21_phosphatemean_bdmax" ~ "__E)__    Phosphate (umol.m-3)"))

# Create line plot from that with y = pred and x = var
response_curve_RF <-  curve_rf_ALL %>%
  ggplot(aes(x = value, y = pred_present_mean, colour = species)) +
  # geom_line(size = 1) +
  geom_point(alpha = 0.01) +
  geom_smooth() +
  facet_wrap(~name, scales = 'free_x') +#, strip.position = "bottom") + 
  labs(x = NULL, y = 'predicted cover (%)', colour = 'species') + 
  scale_color_brewer(type = 'qual', palette = 3) +
  # theme_minimal() +
  theme(legend.position = 'bottom',
        panel.border = element_rect(colour = "black", fill = NA),
        strip.text = ggtext::element_markdown(hjust = 0),
        # strip.placement = "outside", 
        strip.background = element_blank())
ggsave("figures/fig_S5.png", response_curve_RF, height = 8, width = 12)


# Figure S6 ---------------------------------------------------------------
# Simple maps of physical variables over time

# Load future data
Arctic_BO_present <- Arctic_BO %>% mutate(proj = "present") %>% 
  dplyr::select(lon, lat, proj, BO21_templtmax_ss, BO21_salinitymean_ss, BO21_icethickmean_ss) %>% 
  pivot_longer(BO21_templtmax_ss:BO21_icethickmean_ss)
Arctic_BO_2050 <- loadRData("data/Arctic_BO_2050.RData") %>% mutate(proj = "2050") %>% 
  dplyr::select(lon, lat, proj, BO21_templtmax_ss, BO21_salinitymean_ss, BO21_icethickmean_ss) %>% 
  pivot_longer(BO21_templtmax_ss:BO21_icethickmean_ss)
Arctic_BO_2100 <- loadRData("data/Arctic_BO_2100.RData") %>% mutate(proj = "2100") %>% 
  dplyr::select(lon, lat, proj, BO21_templtmax_ss, BO21_salinitymean_ss, BO21_icethickmean_ss) %>% 
  pivot_longer(BO21_templtmax_ss:BO21_icethickmean_ss)
Arctic_BO_futures <- rbind(Arctic_BO_present, Arctic_BO_2050, Arctic_BO_2100) %>% 
  mutate(proj = factor(proj, levels = c("present", "2050", "2100"))) %>% 
  mutate(name = case_when(name == "BO21_templtmax_ss" ~ "Temperature (°C)",
                          name == "BO21_salinitymean_ss" ~ "Salinity (PSS)",
                          name == "BO21_icethickmean_ss" ~ "Ice thickness (m)"))

# Ice plot
futures_plot_ice <- Arctic_BO_futures %>%
  filter(name == "Ice thickness (m)") %>%
  ggplot(aes(x = lon, y = lat)) +
  geom_raster(aes(fill = value)) +
  borders(fill = "grey50", colour = "grey90", size = 0.2) +
  scale_y_continuous(breaks = c(60, 70), labels = c("60°N", "70°N")) +
  scale_x_continuous(breaks = c(-80, -60), labels = c("80°W", "60°W")) +
  coord_quickmap(xlim = c(bbox_arctic[1], bbox_arctic[2]),
                 ylim = c(bbox_arctic[3], bbox_arctic[4]), expand = F) +
  scale_fill_distiller(palette = "Blues") +
  facet_grid(name ~ proj, switch = "y") +
  labs(x = NULL, y = NULL, fill = "m") +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA))
futures_plot_ice

# Salinity plot
futures_plot_sal <- Arctic_BO_futures %>%
  filter(name == "Salinity (PSS)") %>%
  ggplot(aes(x = lon, y = lat)) +
  geom_raster(aes(fill = value)) +
  borders(fill = "grey50", colour = "grey90", size = 0.2) +
  scale_y_continuous(breaks = c(60, 70), labels = c("60°N", "70°N")) +
  scale_x_continuous(breaks = c(-80, -60), labels = c("80°W", "60°W")) +
  coord_quickmap(xlim = c(bbox_arctic[1], bbox_arctic[2]),
                 ylim = c(bbox_arctic[3], bbox_arctic[4]), expand = F) +
  scale_fill_distiller(palette = "Spectral") +
  facet_grid(name ~ proj, switch = "y") +
  labs(x = NULL, y = NULL, fill = "PSS") +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        strip.text.x = element_blank(), strip.background.x = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA))
# futures_plot_sal

# Temperature plot
futures_plot_temp <- Arctic_BO_futures %>%
  filter(name == "Temperature (°C)") %>%
  ggplot(aes(x = lon, y = lat)) +
  geom_raster(aes(fill = value)) +
  borders(fill = "grey50", colour = "grey90", size = 0.2) +
  scale_y_continuous(breaks = c(60, 70), labels = c("60°N", "70°N")) +
  scale_x_continuous(breaks = c(-80, -60), labels = c("80°W", "60°W")) +
  coord_quickmap(xlim = c(bbox_arctic[1], bbox_arctic[2]),
                 ylim = c(bbox_arctic[3], bbox_arctic[4]), expand = F) +
  scale_fill_viridis_c() +
  facet_grid(name ~ proj, switch = "y") +
  labs(x = NULL, y = NULL, fill = "°C") +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        strip.text.x = element_blank(), strip.background.x = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA))
# futures_plot_temp

# Stitch them together
futures_plot_all <- ggpubr::ggarrange(futures_plot_ice, futures_plot_sal, futures_plot_temp, ncol = 1, align = "h")
ggsave("figures/fig_S6.png", futures_plot_all, width = 7, height = 10)

