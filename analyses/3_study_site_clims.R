# analyses/3_study_site_clims
# The purpose of this script is to pull out the monthly clims at each study site
# Whenever new sites are added re-run this script to get the clims for the new pixels

# The first round of data compiled come form the NAPA model created at BIO in Halifax
  # These data are only available through a specific server and so much of the code used below
  # to access them won't work outside of that environment.
  # For this reason much of the NAPA code is commented out
# The second round of data are from Bio-Oracle (http://www.bio-oracle.org/code.php)
  # They are downloaded directly from their server and so may be accessed from anywhere


# Setup -------------------------------------------------------------------

# Load study sites and base packages
source("analyses/1_study_sites.R")

# Bio-Oracle access
library(sdmpredictors)

# Additional packages
library(FNN)

# The NAPA Arctic coords
if(!exists("NAPA_arctic")){
  load("metadata/NAPA_arctic.RData")
  NAPA_arctic <- NAPA_arctic %>% 
    mutate(NAPA_index = as.integer(row.names(.)))
}

# Model variable explanations
model_info <- read_csv("metadata/model_info.csv")

# Convenience loading and coord rounding function
round_coords <- function(df){
  df_round <- data.frame(df) %>%
    mutate(nav_lon = round(nav_lon, 4),
           nav_lat = round(nav_lat, 4))
  if(!("depth" %in% colnames(df_round))) df_round$depth <- 0
  return(df_round)
}


# Load clim files ---------------------------------------------------------

# NB: This section will not run as it requires access to the clim files
  # These files are too large to host on GitHub
  # e-mail Robert (robert.schlegel@dal.ca) for him to send you the large files

# if(!exists("Arctic_surface_clim")){
#   load("data/Arctic_surface_clim.RData")
#   Arctic_surface_clim <- round_coords(Arctic_surface_clim)
# }
# if(!exists("Arctic_ice_clim")){
#   load("data/Arctic_ice_clim.RData")
#   Arctic_ice_clim <- round_coords(Arctic_ice_clim)
# }
# if(!exists("Arctic_depth_T_clim")){
#   load("data/Arctic_depth_T_clim.RData")
#   Arctic_depth_T_clim <- round_coords(Arctic_depth_T_clim)
#   }
# if(!exists("Arctic_depth_U_clim")){
#   load("data/Arctic_depth_U_clim.RData")
#   Arctic_depth_U_clim <- round_coords(Arctic_depth_U_clim)
# }
# if(!exists("Arctic_depth_V_clim")){
#   load("data/Arctic_depth_V_clim.RData")
#   Arctic_depth_V_clim <- round_coords(Arctic_depth_V_clim)
# }
# if(!exists("Arctic_depth_W_clim")){
#   load("data/Arctic_depth_W_clim.RData")
#   Arctic_depth_W_clim <- round_coords(Arctic_depth_W_clim)
# }
# # Join everything
# if(!exists("Arctic_clim")){
#   Arctic_clim <- left_join(Arctic_depth_T_clim, Arctic_depth_U_clim, by = c("x", "y", "nav_lon", "nav_lat", "depth", "bathy", "month")) %>%
#     left_join(Arctic_depth_V_clim, by = c("x", "y", "nav_lon", "nav_lat", "depth", "bathy", "month")) %>%
#     left_join(Arctic_depth_W_clim, by = c("x", "y", "nav_lon", "nav_lat", "depth", "bathy", "month")) %>%
#     left_join(Arctic_surface_clim, by = c("x", "y", "nav_lon", "nav_lat", "depth", "bathy", "month")) %>%
#     left_join(Arctic_ice_clim, by = c("x", "y", "nav_lon", "nav_lat", "depth", "bathy", "month"))
# }


# Load overall mean files -------------------------------------------------

if(!exists("Arctic_surface_mean")){
  load("data/Arctic_surface_mean.RData")
  Arctic_surface_mean <- round_coords(Arctic_surface_mean)
}
if(!exists("Arctic_ice_mean")){
  load("data/Arctic_ice_mean.RData")
  Arctic_ice_mean <- round_coords(Arctic_ice_mean)
}
if(!exists("Arctic_depth_T_mean")){
  load("data/Arctic_depth_T_mean.RData")
  Arctic_depth_T_mean <- round_coords(Arctic_depth_T_mean)
}
if(!exists("Arctic_depth_U_mean")){
  load("data/Arctic_depth_U_mean.RData")
  Arctic_depth_U_mean <- round_coords(Arctic_depth_U_mean)
}
if(!exists("Arctic_depth_V_mean")){
  load("data/Arctic_depth_V_mean.RData")
  Arctic_depth_V_mean <- round_coords(Arctic_depth_V_mean)
}
if(!exists("Arctic_depth_W_mean")){
  load("data/Arctic_depth_W_mean.RData")
  Arctic_depth_W_mean <- round_coords(Arctic_depth_W_mean)
}

# Join everything
if(!exists("Arctic_mean")){
  Arctic_mean <- left_join(Arctic_depth_T_mean, Arctic_depth_U_mean, by = c("x", "y", "nav_lon", "nav_lat", "depth", "bathy")) %>%
    left_join(Arctic_depth_V_mean, by = c("x", "y", "nav_lon", "nav_lat", "depth", "bathy")) %>% 
    left_join(Arctic_depth_W_mean, by = c("x", "y", "nav_lon", "nav_lat", "depth", "bathy")) %>% 
    left_join(Arctic_surface_mean, by = c("x", "y", "nav_lon", "nav_lat", "depth", "bathy")) %>% 
    left_join(Arctic_ice_mean, by = c("x", "y", "nav_lon", "nav_lat", "depth", "bathy"))#, "emp_ice", "emp_oce", "qemp_oce"))
}


# Find nearest neighbours -------------------------------------------------

# Find the row numbers and subset less directly
study_sites_index <- study_sites %>% 
  mutate(NAPA_index = as.vector(knnx.index(as.matrix(NAPA_arctic[,c("nav_lon", "nav_lat")]),
                                 as.matrix(study_sites[,c("lon", "lat")]), k = 1))) %>% 
  left_join(NAPA_arctic, by = "NAPA_index") %>% 
  dplyr::select(-NAPA_index)


# Match monthly clims to each site ----------------------------------------

# NB: This section will not run as it requires access to the clim files
  # These files are too large to host on GitHub
  # e-mail Robert (robert.schlegel@dal.ca) for him to send you the large files

# Variable explanations
# view(model_info)

# Extract clims for each site
# study_site_clims <- right_join(Arctic_clim, study_sites_index, by = c("nav_lon", "nav_lat", "x", "y", "bathy")) %>%
#   dplyr::select(site:Campaign, nav_lon:bathy, everything(), -Date, -Notes)#,
#   #               eken, soce, toce, # Depth variables
#   #               mldr10_1, runoffs, # surface variables
#   #               iceconc_cat, icethic_cat)  # Ice variables
# save(study_site_clims, file = "data/study_site_clims.RData")

# Melt for plotting purposes
# study_site_clims_long <- study_site_clims %>%
#   gather("var", "val", -c(site:bathy)) %>%
#   na.omit() #%>%
#   # mutate(var = factor(var, levels = c("iceconc_cat", "icethic_cat", #Ice
#   #                                     "mldr10_1", "runoffs", # Surface
#   #                                     "eken", "soce", "toce"))) # Depth
# save(study_site_clims_long, file = "data/study_site_clims_long.RData")


# Study site overall means ------------------------------------------------

study_site_means <- right_join(Arctic_mean, study_sites_index, by = c("nav_lon", "nav_lat", "x", "y", "bathy")) %>%
  dplyr::select(site:Campaign, nav_lon:bathy, everything(), -Date, -Notes)
  #               eken, soce, toce, # Depth variables
  #               mldr10_1, runoffs, # surface variables
  #               iceconc_cat, icethic_cat)  # Ice variables
save(study_site_means, file = "data/study_site_means.RData")

# Melt for plotting purposes
study_site_means_long <- study_site_means %>%
  gather("var", "val", -c(site:bathy)) %>%
  na.omit() #%>%
  # mutate(var = factor(var, levels = c("iceconc_cat", "icethic_cat", #Ice
  #                                     "mldr10_1", "runoffs", # Surface
  #                                     "eken", "soce", "toce"))) # Depth
save(study_site_means_long, file = "data/study_site_means_long.RData")


# Visualise ---------------------------------------------------------------

# ggplot(study_site_clims_long, aes(x = month, y = val, colour = depth)) +
#   geom_line(aes(group = depth), size = 1.5) +
#   scale_colour_gradient(low = "steelblue", high = "black", breaks = seq(0, 30, 5),
#                         guide = guide_legend(title.position = "left", nrow = 1)) +
#   facet_grid(var~site, scales = "free", switch = "y") +
#   scale_x_discrete(expand = c(0, 0)) +
#   labs(x = "Monthly climatology", y = NULL, colour = "Depth (m)") +
#   theme(axis.text.x = element_text(angle = 45),
#         legend.position = "bottom")
# ggsave("graph/study_site_clims.pdf", height = 12, width = 44)

