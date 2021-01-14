# analyses/10_figures
# This script houses the code used to create the final figures for the manuscript


# Setup -------------------------------------------------------------------

# Load study sites and base packages
source("analyses/1_study_region_sites.R")
source("analyses/4_kelp_cover.R")

# Other libraries
library(ggOceanMaps)
library(raster)
library(FNN)
library(sdmpredictors)

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

# The base map to use for everything else
Arctic_map <- ggplot() +
  borders(fill = "grey70", colour = "black") +
  scale_y_continuous(breaks = c(60, 70), labels = c("60°N", "70°N")) +
  scale_x_continuous(breaks = c(-80, -60), labels = c("80°W", "60°W")) +
  coord_quickmap(xlim = c(bbox_arctic[1], bbox_arctic[2]),
                 ylim = c(bbox_arctic[3], bbox_arctic[4]), expand = F) +
  labs(x = NULL, y = NULL) +
  theme(panel.border = element_rect(colour = "black", fill = NA))

# The base global map with some corrections
load("metadata/map_base.Rdata")


# Figure 1 ----------------------------------------------------------------
# The map of the study area with the sample points
# Also add colour points for abundance (%)
# From Jesi: What about adding some depth contours on the figure to show 
# that region is characterized by being relatively shallow in a great extension? 
# Also add the names of the main regions: Hudson Bay, Hudson Strait, Lancaster Sound, 
# Davis Strait, etc... And I would delete the name of the stations and just leave the points
# And I can also include all of the data points from the different sources and show those as colours in a legend.
# Also add the MEOW as coloured line polygons. Or labels if too much colour is already being used.
# It may actually be better to show the full Arctic with the study area as a highlighted box 
# or standalone second panel with a reference arrow or such.
# Check the code used in the ArcticNet 2020 talk for the map figure as a starting point

# Site coordinates
adf_summary_coords <- left_join(adf_summary, study_sites, by = c("Campaign", "site"))

# Mean values per site
adf_summary_mean_coords <- filter(adf_summary_coords, family == "kelp.cover") %>% 
  group_by(lon, lat, Campaign) %>%
  summarise(mean_cover = mean(mean_cover),
            range_cover = max(mean_cover)-min(mean_cover), .groups = "drop")

# Create spatial polygon data frame from Arctic bounding box
test1 <- Polygon(cbind(bbox_arctic[c(1,2,1,2)], bbox_arctic[c(3,4,4,3)]), hole = as.logical(NA))
test2 <- Polygons(test1, ID)
SpatialPolygons(test1, pO = 1:length(test1), proj4string = CRS(as.character("+init=epsg:4326")))

bbox_df <- data.frame(lon = c(bbox_arctic[c(1,1,2,2)]), 
                      lat = bbox_arctic[c(3,4,4,3)], id = "bbox")

# make a list
# bbox_list <- split(bbox_df, bbox_df$id)

# only want lon-lats in the list, not the names
bbox_list <- lapply(split(bbox_df, bbox_df$id), function(x) { x["id"] <- NULL; x })

#  make data.frame into spatial polygon, cf. http://jwhollister.com/iale_open_science/2015/07/05/03-Spatial-Data-In-R/

# create SpatialPolygons Object, convert coords to polygon
ps <- sapply(bbox_list, Polygon)

# add id variable 
p1 <- Polygons(ps, ID = 1) 

# create SpatialPolygons object
my_spatial_polys <- SpatialPolygons(list(p1), proj4string = CRS("+proj=longlat +datum=WGS84") ) 

# let's see them
plot(my_spatial_polys)

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
fig_1a <- basemap(limits = c(-180, 180, 40, 90)) +
  annotation_spatial(my_spatial_polys, fill = "forestgreen", alpha = 0.2) +
  annotation_spatial(Arctic_poly, fill = "cadetblue1", alpha = 0.2) +
  geom_spatial_point(data = filter(CANA_kelp, Latitude > 40), 
                     aes(x = Longitude, y = Latitude), colour = "red")
fig_1a

# ArcticKelp campaign map
fig_1b <- basemap(limits = c(bbox_arctic[1], bbox_arctic[2],
                   bbox_arctic[3], bbox_arctic[4]), 
        glaciers = TRUE, bathymetry = TRUE, rotate = TRUE) +
  geom_spatial_point(data = study_sites, crs = "+init=epsg:4326",
                     aes(x = lon, y = lat, colour = Campaign), size = 4) +
  labs(x = NULL, y = NULL)
fig_1b

# Combine and save
fig_1 <- ggpubr::ggarrange(fig_1a, fig_1b, nrow = 1, labels = c("A)", "B)"))
ggsave("figures/fig_1.png", fig_1, height = 6, width = 16)


# Table 1 -----------------------------------------------------------------
# Importance of variables by kelps and models


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
    filter(land_distance <= 100 | depth <= 100) %>% 
    ggplot(aes(x = lon, y = lat)) +
    geom_tile(aes_string(fill = paste0("change_",year_label))) +
    borders(fill = "grey60", colour = "black") +
    scale_y_continuous(breaks = c(60, 70), labels = c("60°N", "70°N")) +
    scale_x_continuous(breaks = c(-80, -60), labels = c("80°W", "60°W")) +
    coord_quickmap(xlim = c(bbox_arctic[1], bbox_arctic[2]),
                   ylim = c(bbox_arctic[3], bbox_arctic[4]), expand = F) +
    # scale_fill_brewer(palette = "Set1", direction = -1) +
    scale_fill_manual(values = c("forestgreen", "white", "red")) +
    labs(x = NULL, y = NULL, fill = "change", title = paste0(year_label," - present")) +
    # theme_bw() +
    theme(legend.position = "bottom",
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA))
  # diff_plot
}

# Function for plotting ensemble model results
# sps_choice <- "Acla"
ensemble_plot <- function(sps_choice, add_legend = F){
  
  # Create full species name
  if(sps_choice %in% c("Ldig", "Lsol", "Slat")){
    sps_title <- "Laminariales spp."
  } else if(sps_choice == "Acla"){
    sps_title <- "Agarum clathratum"
  } else if(sps_choice == "Aesc"){
    sps_title <- "Alaria esculenta"
  } else{
    stop("*sad robot noises*")
  }
  
  # Prep data for plotting
  if(sps_choice %in% c("Ldig", "Lsol", "Slat")){
    sps_choice <- c("Ldig", "Lsol", "Slat") # This is used by 'sps_points' below to get all three species data
    df_project <- plyr::ddply(data.frame(sps_name = c("Ldig", "Lsol", "Slat")), c("sps_name"), ensemble_prep) %>% 
      pivot_wider(values_from = presence, names_from = sps_name) %>% 
      mutate(presence = ifelse(Ldig + Lsol + Slat >= 2, TRUE, FALSE)) %>% 
      dplyr::select(lon, lat, presence, land_distance, depth, projection)
  } else{
    df_project <- ensemble_prep(sps_choice)
  }
  
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
    filter(land_distance <= 100 | depth <= 100,
           presence == 1,
           projection == "proj_pres") %>%
    # dplyr::select(lon, lat) %>% 
    # distinct()
    mutate(presence = as.logical(presence)) %>% 
    ggplot(aes(x = lon, y = lat)) +
    geom_tile(aes(fill = presence)) +
    borders(fill = "grey60", colour = "black") +
    geom_point(data = sps_points, colour = "yellow", size = 0.5) +
    scale_y_continuous(breaks = c(60, 70), labels = c("60°N", "70°N")) +
    scale_x_continuous(breaks = c(-80, -60), labels = c("80°W", "60°W")) +
    coord_quickmap(xlim = c(bbox_arctic[1], bbox_arctic[2]),
                   ylim = c(bbox_arctic[3], bbox_arctic[4]), expand = F) +
    scale_fill_manual(values = c("forestgreen")) +
    labs(x = NULL, y = NULL, title = paste0(sps_title)) +
    theme_bw() +
    theme(legend.position = "bottom",
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
ensemble_Lam <- ensemble_plot("Ldig")# NB: This automagically combines all the Laminariales
ensemble_legend <- ensemble_plot("Acla", add_legend = T)

# Combine into one mecha-figure
ensemble_ALL <- ggpubr::ggarrange(ensemble_Acla, ensemble_Aesc, ensemble_Lam, ensemble_legend,
                                  ncol = 1, labels = c("A)", "B)", "C)", ""), heights = c(1, 1, 1, 0.15))
ggsave("figures/fig_3.png", ensemble_ALL, width = 7, height = 11)


# Figure 4 ----------------------------------------------------------------
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
    mutate(pred_present_round = plyr::round_any(pred_present_mean, 10),
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
                            depth <= 100 | land_distance <= 100),
              aes(x = lon, y = lat, fill = pred_present_round)) +
    # scale_fill_viridis_c(paste0("cover (%)"), limits = c(10, 70)) +
    # scale_fill_distiller(palette = "Greens", direction = 1, limits = c(10, 70)) +
    scale_fill_gradient("cover (%)", low = "white", high = "forestgreen", 
                        limits = c(0, 70), breaks = c(10, 20, 30, 40, 50, 60), guide = "legend") +
    borders(fill = "grey60", colour = "black") +
    scale_y_continuous(breaks = c(60, 70), labels = c("60°N", "70°N")) +
    scale_x_continuous(breaks = c(-80, -60), labels = c("80°W", "60°W")) +
    coord_quickmap(xlim = c(bbox_arctic[1], bbox_arctic[2]),
                   ylim = c(bbox_arctic[3], bbox_arctic[4]), expand = F) +
    labs(x = NULL, y = NULL, title = paste0(sps_title)) +
    theme_bw() +
    theme(legend.position = "bottom",
          panel.border = element_rect(colour = "black", fill = NA))
  # p_present
  
  # 2050 plot
  p_2050 <- ggplot() +
    # geom_tile(data = df, # No filter
    geom_tile(data = filter(project_diff, depth <= 100 | land_distance <= 100),
              aes(x = lon, y = lat, fill = pred_diff_2050)) +
    scale_fill_gradient2("change (%)", low = "red", high = "forestgreen",
                         limits = c(-40, 40), breaks = c(-40, -30, -20, -10, 0, 10, 20, 30, 40),
                         guide = "legend", na.value = NA) +
    # scale_fill_discrete("2050 - present (%)") +
    borders(fill = "grey60", colour = "black") +
    scale_y_continuous(breaks = c(60, 70), labels = c("60°N", "70°N")) +
    scale_x_continuous(breaks = c(-80, -60), labels = c("80°W", "60°W")) +
    coord_quickmap(xlim = c(bbox_arctic[1], bbox_arctic[2]),
                   ylim = c(bbox_arctic[3], bbox_arctic[4]), expand = F) +
    labs(x = NULL, y = NULL, title = "2050 - present") +
    # theme_bw() +
    theme(legend.position = "bottom",
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA))
  # p_2050
  
  # 2050 plot
  p_2100 <-  ggplot() +
    # geom_tile(data = df, # No filter
    geom_tile(data = filter(project_diff, depth <= 100 | land_distance <= 100),
              aes(x = lon, y = lat, fill = pred_diff_2100)) +
    scale_fill_gradient2("change (%)", low = "red", high = "forestgreen", 
                         limits = c(-40, 40), breaks = c(-40, -30, -20, -10, 0, 10, 20, 30, 40), 
                         guide = "legend", na.value = NA) +
    # scale_fill_discrete("2100 - present (%)") +
    borders(fill = "grey60", colour = "black") +
    scale_y_continuous(breaks = c(60, 70), labels = c("60°N", "70°N")) +
    scale_x_continuous(breaks = c(-80, -60), labels = c("80°W", "60°W")) +
    coord_quickmap(xlim = c(bbox_arctic[1], bbox_arctic[2]),
                   ylim = c(bbox_arctic[3], bbox_arctic[4]), expand = F) +
    labs(x = NULL, y = NULL, title = "2100 - present") +
    # theme_bw() +
    theme(legend.position = "bottom",
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
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
rf_legend <- rf_plot("alaria", add_legend = T)

# Combine into one mecha-figure
rf_ALL <- ggpubr::ggarrange(rf_agarum, rf_alaria, rf_lam, rf_legend,
                             ncol = 1, labels = c("A)", "B)", "C)", ""), heights = c(1, 1, 1, 0.15))
ggsave("figures/fig_4.png", rf_ALL, width = 7, height = 11)


# Figure 5 ----------------------------------------------------------------

# Comparison of the present modelling approaches


