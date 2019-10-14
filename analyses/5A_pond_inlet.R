# analyses/pond_inlet.R
# The purpose of this script is to extract NAPA data for only the Pond Inlet area
# It then creates visualisations at depth for monthly clims


# Libraries ---------------------------------------------------------------

# The study sites and bounding box
source("2_study_sites.R")

# The Pond Inlet bounding box
bbox_PI <- c(-81, -76, 71.5, 73) 

# The arctic lon/lat/depth coords
load("metadata/NAPA_arctic.RData")


# Subset for Pond Inlet ---------------------------------------------------

if(!exists("Arctic_depth_T_clim")){
  load("data/Arctic_depth_T_clim.RData")
}

if(!exists("Arctic_ice_clim")){
  load("data/Arctic_ice_clim.RData")
}

PI_sites <- study_sites %>% 
  filter(lon >= bbox_PI[1], lon <= bbox_PI[2],
         lat >= bbox_PI[3], lat <= bbox_PI[4])

PI_depth_T_clim <- Arctic_depth_T_clim %>%
  filter(nav_lon >= bbox_PI[1], nav_lon <= bbox_PI[2],
         nav_lat >= bbox_PI[3], nav_lat <= bbox_PI[4])

PI_ice_clim <- Arctic_ice_clim %>%
  filter(nav_lon >= bbox_PI[1], nav_lon <= bbox_PI[2],
         nav_lat >= bbox_PI[3], nav_lat <= bbox_PI[4])

# Look at ice only for April to September
PI_ice_clim_sub <- PI_ice_clim %>% 
  filter(month %in% c("Apr", "May", "Jun", "Jul", "Aug", "Sep"))

# Pull out the weird low salinity data
low_sss <- PI_depth_T_clim %>%
  filter(soce < 10)


# Visualise ---------------------------------------------------------------

# # Single file of data
# ggplot(filter(res, depth == 5), 
#        aes(x = nav_lon, y = nav_lat, colour = soce)) +
#   geom_point(size = 0.001) +
#   scale_colour_viridis_c()
# 
# # Daily data
# ggplot(filter(Arctic_depth_T, t == "1998-01-05"), 
#        aes(x = nav_lon, y = nav_lat, colour = soce)) +
#   geom_point(size = 0.001) +
#   scale_colour_viridis_c()
# 
#  # Monthly clims
# ggplot(filter(Arctic_depth_T_clim, month == "Jan"), 
#        aes(x = nav_lon, y = nav_lat, colour = soce)) +
#   geom_point(size = 0.001) +
#   scale_colour_viridis_c()

# Plot the salinity in the area just around Pond Inlet
sss_month_depth <- ggplot(PI_depth_T_clim, aes(x = nav_lon, y = nav_lat)) +
  geom_point(size = 3, shape = 19, aes(colour = soce)) +
  borders(fill = "grey70", colour = "black") +
  geom_point(data = PI_sites, colour = "red", aes(x = lon, y = lat)) +
  geom_label_repel(data = PI_sites, aes(x = lon, y = lat, label = site), nudge_y = -1) +
  scale_colour_viridis_c() +
  coord_equal(xlim = c(-81, -76), ylim = c(71.5, 73)) +
  labs(x = "Longitude", y = "Latitude", colour = "salinity") +
  facet_grid(month~depth)
ggsave("graph/sss_month_depth.pdf", sss_month_depth, height = 12, width = 25)

# Plot the kinetic energy in the area just around Pond Inlet
kin_month_depth <- ggplot(PI_depth_T_clim, aes(x = nav_lon, y = nav_lat)) +
  geom_point(size = 3, shape = 19, aes(colour = eken)) +
  borders(fill = "grey70", colour = "black") +
  geom_point(data = PI_sites, colour = "red", aes(x = lon, y = lat)) +
  geom_label_repel(data = PI_sites, aes(x = lon, y = lat, label = site), nudge_y = -1) +
  scale_colour_viridis_c(option = "B") +
  coord_equal(xlim = c(-81, -76), ylim = c(71.5, 73)) +
  labs(x = "Longitude", y = "Latitude", colour = "kinetic\nenergy\n(m2/s2)") +
  facet_grid(month~depth)
# kin_month_depth
ggsave("graph/kin_month_depth.pdf", kin_month_depth, height = 12, width = 25)

# Plot the ice concentration just around Pond Inlet
# A rainbow colour palette was explicitly requested
# ice_month <- ggplot(PI_ice_clim, aes(x = nav_lon, y = nav_lat)) +
ice_month <- ggplot(PI_ice_clim_sub, aes(x = nav_lon, y = nav_lat)) +
  geom_point(size = 5, shape = 15, aes(colour = iceconc_cat*100)) +
  borders(fill = "grey70", colour = "black") +
  geom_point(data = PI_sites, colour = "red", aes(x = lon, y = lat)) +
  geom_label_repel(data = PI_sites, aes(x = lon, y = lat, label = site), nudge_y = -1) +
  scale_colour_gradientn(colours = rainbow_palette, 
                         breaks = c(0, 5, 10, 15, 20),
                         labels = c(0, 5, 10, 15, 20), 
                         limits = c(0, 20)) +
  coord_equal(xlim = c(-81, -76), ylim = c(71.5, 73)) +
  labs(x = "Longitude", y = "Latitude", colour = "Ice cover (%)") +
  facet_wrap(~month)
# ice_month
# ggsave("graph/ice_month.pdf", ice_month, height = 5, width = 15)
ggsave("graph/ice_month_sub.pdf", ice_month, height = 5, width = 15)
