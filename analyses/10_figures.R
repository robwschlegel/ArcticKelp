# analyses/10_figures
# This script houses the code used to create the final figures for the manuscript


# Setup -------------------------------------------------------------------

# Load study sites and base packages
source("analyses/1_study_region_sites.R")
source("analyses/4_kelp_cover.R")

# Other libraries
library(ggOceanMaps) # https://mikkovihtakari.github.io/ggOceanMaps/index.html
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

# Function for preparing bathymetry data for plotting in a polar projection
bbox_to_bathy <- function(coords, lon_pad = 0, lat_pad = 0,
                          bathy_file = NA, projection = NA,
                          depths = c(0, 50, 100, 200, 500, 1000, 2000, 10000)){
  
  # Get the coordinates
  if(is.data.frame(coords)){
    lon1 <- min(coords$lon1); lon2 <- max(coords$lon2)
    lat1 <- min(coords$lat1); lat2 <- max(coords$lat2)
  } else if(is.vector(coords)){
    lon1 <- coords[1]; lon2 <- coords[2]
    lat1 <- coords[3]; lat2 <- coords[4]
  } else {
    stop("Uh oh")
  }
  
  # Use the default hi-res Arctic bathy unless the user specifies something else
  # if(is.na(bathy_file)) bathy_file <- paste0(pCloud_path,"FACE-IT_data/shape_files/IBCAO_v4_200m.nc") # Super hi-res, but doesn't work...
  if(is.na(bathy_file)) bathy_file <- "~/pCloudDrive/FACE-IT_data/maps/GRIDONE_2D.nc"
  # if(is.na(bathy_file)) bathy_file <- paste0(pCloud_path,"FACE-IT_data/shape_files/GEBCO_2020.nc")
  
  # Set limits for bathy projection
  xlon <- c(lon1-lon_pad, lon2+lon_pad)
  xlat <- c(lat1-lat_pad, lat2+lat_pad)
  lims <- c(xlon, xlat)
  
  # Set projection
  if(is.na(projection)){
    # projection <- "+init=epsg:6070"
    projection <- "+init=epsg:3995" # Arctic Polar Stereographic
    # projection <- "+init=epsg:4326" # Cartesian global
    # projection <- "+init=epsg:32636"
  } 
  
  # Convert NetCDF to raster
  rb <- raster_bathymetry(bathy = bathy_file,
                          depths = depths, 
                          proj.out = projection, 
                          boundary = lims)
  
  # Convert to raster vector for plotting
  bs_bathy <- vector_bathymetry(rb)
  
  # Convert land file for use with new bathy file
  world <- rgdal::readOGR("~/pCloudDrive/FACE-IT_data/maps/ne_10m_land.shp")
  islands <- rgdal::readOGR("~/pCloudDrive/FACE-IT_data/maps/ne_10m_minor_islands.shp")
  world <- rbind(world, islands)
  # proj4string(world) <- CRS(projection)
  bs_land <- clip_shapefile(world, lims, proj.limits = projection)
  bs_land <- sp::spTransform(bs_land, CRSobj = sp::CRS(projection))
  if(!rgeos::gIsValid(bs_land)){ # Has to return TRUE, if not use rgeos::gBuffer
    bs_land <- rgeos::gBuffer(bs_land, byid = TRUE, width = 0)
  }
  
  # Create glacier shape files
  glaciers <- rgdal::readOGR("~/pCloudDrive/FACE-IT_data/maps/ne_10m_glaciated_areas.shp")
  if(!rgeos::gIsValid(glaciers)){ # Needs buffering
    glaciers <- rgeos::gBuffer(glaciers, byid = TRUE, width = 0)
  }
  bs_glacier <- clip_shapefile(glaciers, lims)
  if(dim(bs_glacier)[1] > 0){
    bs_glacier <- sp::spTransform(bs_glacier, CRSobj = sp::CRS(projection))
  } else { 
    bs_glacier <- NA
  }
  
  # Return results
  res <- list(bathy = bs_bathy, land = bs_land, glacier = bs_glacier)
  return(res)
}

# Prep Arctic bathy data
# NB: This requires too much RAM to run on a laptop
arctic_bathy <- bbox_to_bathy(c(-180, 180, 40, 90))

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
label_df <- data.frame(lon = c(-85, -75.17, -82.92),
                       lat = c(60, 62.544, 74.425),
                       loc = c("Hudson Bay", "Hudson Straight", "Lancaster Sound"))

# ArcticKelp campaign map
# fig_1b <- basemap(limits = c(bbox_arctic[1], bbox_arctic[2],
#                              bbox_arctic[3], bbox_arctic[4]),
fig_1b <- basemap(limits = bbox_arctic,
                  glaciers = TRUE, bathymetry = TRUE, rotate = TRUE) +
  geom_spatial_label(data = label_df, crs = 4326,
                     aes(x = lon, y = lat, label = loc)) +
  geom_spatial_point(data = study_sites, crs = 4326, shape = 21,
                     aes(x = lon, y = lat), colour = "black", fill = "magenta", size = 4) +
  labs(x = NULL, y = NULL)
fig_1b
fig_1b <- ggplot() +
  # geom_tile(aes(fill = presence)) +
  borders(fill = "grey30", colour = "black") +
  geom_point(data = study_sites, shape = 21, colour = "black", 
             fill = "magenta", size = 2,
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
fig_1b

# Combine and save
fig_1 <- ggpubr::ggarrange(fig_1a, fig_1b, ncol = 1, labels = c("A)", "B)"))
ggsave("figures/fig_1.png", fig_1, height = 12, width = 8)


# Table 1 -----------------------------------------------------------------
# Importance of variables by kelps and models


# Figure 2 ----------------------------------------------------------------
# Relationship with important variables and kelp PA and abundance

# Figure 3 ---------------------------------------------------------------
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
  df_project_all <- rbind(df_project_present, df_project_2050, df_project_2100) %>% 
    mutate(presence = as.integer(presence))
  return(df_project_all)
}

# Function for visualising changes over time
ensemble_diff_plot <- function(df_project, year_label){
  diff_plot <- df_project %>% 
    pivot_wider(names_from = projection, values_from = presence) %>% 
    # na.omit() %>% 
    mutate(proj_2050 = ifelse(is.na(proj_2050), FALSE, proj_2050),
           change_2050 = proj_2050 - proj_pres,
           change_2050 = ifelse(proj_pres == F & proj_2050 == F, NA, change_2050),
           change_2050 =  factor(change_2050,
                                 # levels = c("-1", "0", "1"),
                                 # labels = c("gain", "no change", "loss")),
                                 levels = c("1", "0", "-1"),
                                 labels = c("gain", "no change", "loss")),
           proj_2100 = ifelse(is.na(proj_2100), FALSE, proj_2100),
           change_2100 = proj_2100 - proj_pres,
           change_2100 = ifelse(proj_pres == F & proj_2100 == F, NA, change_2100),
           change_2100 = factor(change_2100,
                                levels = c("1", "0", "-1"),
                                labels = c("gain", "no change", "loss"))) %>%
    na.omit() %>% 
    filter(land_distance <= 50 | depth <= 100) %>% 
    ggplot(aes(x = lon, y = lat)) +
    geom_tile(aes_string(fill = paste0("change_",year_label))) +
    borders(fill = "grey50", colour = "grey90", size = 0.2) +
    scale_y_continuous(breaks = c(60, 70), labels = c("60°N", "70°N")) +
    scale_x_continuous(breaks = c(-80, -60), labels = c("80°W", "60°W")) +
    coord_quickmap(xlim = c(bbox_arctic[1], bbox_arctic[2]),
                   ylim = c(bbox_arctic[3], bbox_arctic[4]), expand = F) +
    # scale_fill_brewer(palette = "Set1", direction = -1) +
    scale_fill_manual(values = c(RColorBrewer::brewer.pal(9, "Blues")[7], "grey80",
                                 RColorBrewer::brewer.pal(9, "Reds")[7])) +
    labs(x = NULL, y = NULL, fill = "Change", title = paste0(year_label," - present")) +
    # theme_bw() +
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
  # if(sps_choice %in% c("Ldig", "Lsol", "Slat")){
  #   sps_choice <- c("Ldig", "Lsol", "Slat") # This is used by 'sps_points' below to get all three species data
  #   df_project <- plyr::ddply(data.frame(sps_name = c("Ldig", "Lsol", "Slat")), c("sps_name"), ensemble_prep) %>% 
  #     pivot_wider(values_from = presence, names_from = sps_name) %>% 
  #     mutate(presence = ifelse(Ldig + Lsol + Slat >= 2, TRUE, FALSE)) %>% 
  #     dplyr::select(lon, lat, presence, land_distance, depth, projection)
  # } else{
  df_project <- ensemble_prep(sps_choice)
  # }
  
  # test <- df_project %>% 
  #   dplyr::select(lon, lat) %>% 
  #   distinct()
  
  # Load the species points
  sps_points <- map_df(sps_files[grepl(paste(sps_choice, collapse = "|"), sps_files)], read_csv) %>% 
    mutate(env_index = as.vector(knnx.index(as.matrix(global_coords[,c("lon", "lat")]),
                                            as.matrix(.[,2:3]), k = 1))) %>%
    left_join(global_coords, by = "env_index") %>% 
    dplyr::select(Sp, lon.y, lat.y) %>%
    dplyr::rename(lon = lon.y, lat = lat.y)
  
  # Visualise present data
  plot_present <- df_project %>% 
    filter(land_distance <= 50 | depth <= 100,
           presence == 1,
           projection == "proj_pres") %>%
    # dplyr::select(lon, lat) %>% 
    # distinct()
    mutate(presence = "") %>%
    ggplot(aes(x = lon, y = lat)) +
    geom_tile(aes(fill = presence)) +
    borders(fill = "grey50", colour = "grey90", size = 0.2) +
    geom_point(data = sps_points, shape = 21, colour = "black", fill = "hotpink", size = 0.5) +
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
  plot_2050 <- ensemble_diff_plot(df_project, "2050")
  plot_2100 <- ensemble_diff_plot(df_project, "2100")
  
  # Combine and exit
  plot_ALL <- cowplot::plot_grid(
    plot_present + theme(legend.position = "none"),
    plot_2050 + theme(legend.position = "none"),
    plot_2100 + theme(legend.position = "none"),
    ncol = 3,
    align = "v")
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
ensemble_ALL <- ggpubr::ggarrange(ensemble_Acla, ensemble_Aesc, ensemble_Lsol, ensemble_Slat, ensemble_legend,
                                  ncol = 1, labels = c("A)", "B)", "C)", "D)", ""), heights = c(1, 1, 1, 1, 0.15))
ggsave("figures/fig_3.png", ensemble_ALL, width = 7, height = 14)


# Figure 4 ----------------------------------------------------------------

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
  
  # Model evaluation by algorithm
  model_res <- biomod2::models_scores_graph(biomod_model, by = "models", metrics = c('ROC','TSS')) + 
    geom_hline(aes(yintercept = 0.7), colour = "red", size = 2) +
    ggtitle(sps_title) +
    coord_cartesian(xlim = c(0.6,1), ylim = c(0.3,1), expand = F) +
    theme(plot.title = element_text(face = "italic"))
  # model_res
  
  return(model_res)
}

# Create plots
model_stats_Acla <- model_stats_plot("Acla")
model_stats_Aesc <- model_stats_plot("Aesc")
model_stats_Lsol <- model_stats_plot("Lsol")
model_stats_Slat <- model_stats_plot("Slat")

# Combine and save
fig_4 <- ggpubr::ggarrange(model_stats_Acla, model_stats_Aesc, model_stats_Lsol, model_stats_Slat, 
                           ncol = 2, nrow = 2, common.legend = TRUE, legend = "bottom", 
                           labels = c("A)", "B)", "C)", "D)"))
ggsave("figures/fig_4.png", fig_4, width = 5, height = 7)


# Figure 5 ----------------------------------------------------------------

# Combination of the modelling approaches

# Join ensemble and random forest results
model_compare_plot <- function(model_choice, add_legend = F){
  
  # Create full species name and load ensemble data
  if(model_choice == "laminariales"){
    sps_title <- "Laminariales spp."
    df_project <- left_join(ensemble_prep("Lsol"), ensemble_prep("Slat"),
                            by = c("lon", "lat", "land_distance", "depth", "projection")) %>% 
      mutate(presence = case_when(presence.x == 1 | presence.y == 1 ~ 1, TRUE ~ 0)) %>% 
      dplyr::select(-presence.x, -presence.y) %>% 
      pivot_wider(names_from = projection, values_from = presence)
  } else if(model_choice == "agarum"){
    sps_title <- "Agarum clathratum"
    df_project <- ensemble_prep("Acla") %>% 
      pivot_wider(names_from = projection, values_from = presence)
  } else if(model_choice == "alaria"){
    sps_title <- "Alaria esculenta"
    df_project <- ensemble_prep("Aesc") %>% 
      pivot_wider(names_from = projection, values_from = presence)
  } else{
    stop("*sad robot noises*")
  }
  
  # Load random forest data
  best_rf <- loadRData(paste0("data/best_rf_",model_choice,".RData"))$project_multi
  
  # Join the models
  model_join <- left_join(best_rf, df_project, by = c("lon", "lat", "depth", "land_distance")) %>% 
    mutate(pred_present_round = plyr::round_any(pred_present_mean, 20),
           pred_diff_2050 = plyr::round_any(pred_2050_mean - pred_present_mean, 10),
           pred_diff_2100 = plyr::round_any(pred_2100_mean - pred_present_mean, 10))
  
  # Plot
  p_present <-  ggplot() +
    geom_tile(data = filter(model_join,
                            proj_pres == 1,
                            # pred_present_round >= 10,
                            depth <= 100 | land_distance <= 50),
              aes(x = lon, y = lat, fill = pred_present_round)) +
    scale_fill_gradientn("Cover (%)", colours = RColorBrewer::brewer.pal(9, "BuGn")[c(5,6,7,8)],
                         limits = c(0, 70), breaks = c(0, 20, 40, 60), guide = "legend") +
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
  if(model_choice %in% c("agarum", "alaria")){
    p_present <- p_present + theme(plot.title = element_text(face = "italic"))
  }
  p_present
  
  # 2050 plot
  p_2050 <- ggplot() +
    geom_tile(data = filter(model_join,
                            proj_2050 == 1,
                            # pred_present_round >= 10,
                            depth <= 100 | land_distance <= 50),
              aes(x = lon, y = lat, fill = pred_diff_2050)) +
    # scale_fill_gradient2("Change (%)", low = "red", mid = "grey80", high = "blue",
    scale_fill_gradientn("Change (%)", 
                         colours = c(RColorBrewer::brewer.pal(9, "Reds")[c(9,8,7,6)], "grey80",
                                     RColorBrewer::brewer.pal(9, "Blues")[c(6,7,8,9)]),
                         limits = c(-40, 40), 
                         breaks = c(-40, -30, -20, -10, 0, 10, 20, 30, 40),
                         guide = "legend", na.value = NA) +
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
    geom_tile(data = filter(model_join,
                            proj_2100 == 1,
                            # pred_present_round >= 10,
                            depth <= 100 | land_distance <= 50),
              aes(x = lon, y = lat, fill = pred_diff_2100)) +
    # scale_fill_gradient2("Change (%)", low = "red", mid = "grey80", high = "blue",
    scale_fill_gradientn("Change (%)", 
                         colours = c(RColorBrewer::brewer.pal(9, "Reds")[c(9,8,7,6)], "grey80",
                                     RColorBrewer::brewer.pal(9, "Blues")[c(6,7,8,9)]),
                         limits = c(-40, 40), 
                         breaks = c(-40, -30, -20, -10, 0, 10, 20, 30, 40), 
                         guide = "legend", na.value = NA) +
    borders(fill = "grey50", colour = "grey90", size = 0.2) +
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
fig_5 <- ggpubr::ggarrange(model_compare_agarum, model_compare_alaria, model_compare_lam, model_compare_legend,
                           ncol = 1, labels = c("A)", "B)", "C)", ""), heights = c(1, 1, 1, 0.15))
ggsave("figures/fig_5.png", fig_5, width = 7, height = 11)


# Figure 6 ----------------------------------------------------------------

# Random forest model results
load("data/best_rf_kelpcover.RData")
load("data/best_rf_laminariales.RData")
load("data/best_rf_agarum.RData")
load("data/best_rf_alaria.RData")

# Confidence in RF output
# Function for creating figure showing confidence intervals of prediction accuracy
conf_plot_RF <- function(df, plot_title){
  
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
  ggplot(conf_acc, aes(x = original, y = mean)) +
    geom_hline(yintercept = 0, size = 2, colour = "red") +
    geom_crossbar(aes(y = 0, ymin = q05, ymax = q95),
                  fatten = 0, fill = "grey70", colour = NA, width = 10) +
    geom_crossbar(aes(ymin = q25, ymax = q75),
                  fatten = 0, fill = "grey50", width = 10) +
    geom_crossbar(aes(ymin = q50, ymax = q50),
                  fatten = 0, fill = NA, colour = "black", width = 10) +
    # geom_segment(data = conf_best, aes(xend = original, y = mean_acc, yend = 0), 
    # colour = "purple", size = 1.2, alpha = 0.8) +
    # geom_point(data = conf_best, aes(y = mean_acc), colour = "purple", size = 3, alpha = 0.8) +
    geom_label(data = conf_mean_label, 
               aes(x = 75, y = 75, label = paste0("Mean accuracy: ",mean_acc,"±",sd_acc,"; r = ",r_acc))) +
    # geom_label(data = conf_best_label, colour = "purple",
    # aes(x = 75, y = 60, label = paste0("Best accuracy: ",mean_acc,"±",sd_acc,"; r = ",r_acc))) +
    scale_y_continuous(limits = c(-100, 100)) +
    scale_x_continuous(breaks = c(0, 20, 40, 60, 80, 100)) +
    labs(y = "Range in accuracy of predictions", x = "Original value (% cover)", title = plot_title) +
    theme(panel.border = element_rect(fill = NA, colour = "black"))
}

# Create the plots
cover_lam <- conf_plot_RF(best_rf_laminariales$accuracy_reg, "Laminariales cover confidence")
cover_agarum <- conf_plot_RF(best_rf_agarum$accuracy_reg, "Agarum cover confidence")
cover_alaria <- conf_plot_RF(best_rf_alaria$accuracy_reg, "Alaria cover confidence")

fig_6 <- ggpubr::ggarrange(cover_agarum, cover_alaria, cover_lam, ncol = 1, 
                           labels = c("A)", "B)", "C)"))
ggsave("figures/fig_6.png", fig_6, width = 7, height = 12)


# Figure S1 ---------------------------------------------------------------
# The random forest model results
# Combine into one dataframe to have the same legend via facets

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
           pred_diff_2100 = plyr::round_any(pred_2100_mean - pred_present_mean, 10)) #%>% 
  # mutate(pred_diff_2050 = ifelse(pred_diff_2050 == 0, NA, pred_diff_2050),
  # pred_diff_2100 = ifelse(pred_diff_2100 == 0, NA, pred_diff_2100))#,
  # pred_diff_2050 = as.factor(pred_diff_2050),
  # pred_diff_2100 = as.factor(pred_diff_2100))
  # pivot_longer(cols = pred_present_mean:pred_diff_2100) %>% 
  # filter(name == "pred_diff_2050" & value != 0) %>% 
  # pivot_wider()
  # mutate(pred_diff_2050 = base::cut(pred_2050_mean - pred_present_mean, breaks = c(-0.1, 0, 0.1, 0.2)),
  #        pred_diff_2100 = base::cut(pred_2100_mean - pred_present_mean, breaks = c(-0.1, 0, 0.1, 0.2)))
  
  # Scale for difference plots
  # diff_range <- range(c(project_diff$pred_diff_2050, project_diff$pred_diff_2100), na.rm = T)
  
  # Present plot
  p_present <-  ggplot() +
    # geom_tile(data = df, # No filter
    geom_tile(data = filter(project_diff,
                            pred_present_round >= 10,
                            depth <= 100 | land_distance <= 50),
              aes(x = lon, y = lat, fill = pred_present_round)) +
    # scale_fill_viridis_c(paste0("cover (%)"), limits = c(10, 70)) +
    # scale_fill_distiller(palette = "Greens", direction = 1, limits = c(10, 70)) +
    scale_fill_gradient("Cover (%)", low = "grey90", high = "grey30", 
                        limits = c(0, 70), breaks = c(20, 40, 60), guide = "legend") +
    borders(fill = "grey50", colour = "grey90", size = 0.2) +
    scale_y_continuous(breaks = c(60, 70), labels = c("60°N", "70°N")) +
    scale_x_continuous(breaks = c(-80, -60), labels = c("80°W", "60°W")) +
    coord_quickmap(xlim = c(bbox_arctic[1], bbox_arctic[2]),
                   ylim = c(bbox_arctic[3], bbox_arctic[4]), expand = F) +
    labs(x = NULL, y = NULL, title = paste0(sps_title)) +
    # theme_bw() +
    theme(legend.position = "bottom",
          panel.background = element_rect(fill = "grey100"),
          panel.border = element_rect(colour = "black", fill = NA))
  # p_present
  
  if(kelp_choice %in% c("agarum", "alaria")){
    p_present <- p_present + theme(plot.title = element_text(face = "italic"))
  }
  
  # 2050 plot
  p_2050 <- ggplot() +
    # geom_tile(data = df, # No filter
    geom_tile(data = filter(project_diff, 
                            pred_present_round >= 10,
                            depth <= 100 | land_distance <= 50),
              aes(x = lon, y = lat, fill = pred_diff_2050)) +
    scale_fill_gradientn("Change (%)", 
                         colours = c(RColorBrewer::brewer.pal(9, "Reds")[c(9,8,7,6)], "grey80",
                                     RColorBrewer::brewer.pal(9, "Blues")[c(6,7,8,9)]),
                         limits = c(-40, 40), 
                         breaks = c(-40, -30, -20, -10, 0, 10, 20, 30, 40),
                         guide = "legend", na.value = NA) +
    # scale_fill_discrete("2050 - present (%)") +
    borders(fill = "grey50", colour = "grey90", size = 0.2) +
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
          panel.background = element_rect(fill = "grey100"),
          panel.border = element_rect(colour = "black", fill = NA))
  # p_2050
  
  # 2050 plot
  p_2100 <-  ggplot() +
    # geom_tile(data = df, # No filter
    geom_tile(data = filter(project_diff,
                            pred_present_round >= 10,
                            depth <= 100 | land_distance <= 50),
              aes(x = lon, y = lat, fill = pred_diff_2100)) +
    scale_fill_gradientn("Change (%)", 
                         colours = c(RColorBrewer::brewer.pal(9, "Reds")[c(9,8,7,6)], "grey80",
                                     RColorBrewer::brewer.pal(9, "Blues")[c(6,7,8,9)]),
                         limits = c(-40, 40), 
                         breaks = c(-40, -30, -20, -10, 0, 10, 20, 30, 40),
                         guide = "legend", na.value = NA) +
    # scale_fill_discrete("2100 - present (%)") +
    borders(fill = "grey50", colour = "grey90", size = 0.2) +
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
rf_ALL <- ggpubr::ggarrange(rf_agarum, rf_alaria, rf_lam, rf_legend,
                            ncol = 1, labels = c("A)", "B)", "C)", ""), heights = c(1, 1, 1, 0.15))
ggsave("figures/fig_S1.png", rf_ALL, width = 7, height = 11)


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
  if(model_choice %in% c("agarum", "alaria")){
    p_present <- p_present + theme(plot.title = element_text(face = "italic"))
  }
  # p_present
  
  # 2050 plot
  p_2050 <- ggplot() +
    geom_tile(data = filter(model_join,
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
    geom_tile(data = filter(model_join,
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
ggsave("figures/fig_S2.png", rf_var_ALL, width = 7, height = 11)

