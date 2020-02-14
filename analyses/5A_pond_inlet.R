# analyses/pond_inlet.R
# The purpose of this script is to extract NAPA data for only the Pond Inlet area
# It then creates visualisations at depth for monthly clims


# Libraries ---------------------------------------------------------------

# The study sites and bounding box
source("analyses/4_kelp_cover.R")

# The Pond Inlet bounding box
bbox_PI <- c(-81, -76, 71.5, 73) 


# Subset for Pond Inlet ---------------------------------------------------

study_site_means_PI <- study_site_means %>% 
  filter(lon >= bbox_PI[1], lon <= bbox_PI[2],
         lat >= bbox_PI[3], lat <= bbox_PI[4])

adf_summary_PI <- adf_summary %>% 
  filter(site %in% study_site_means_PI$site)

# Pull out the weird low salinity data
# low_sss <- PI_depth_T_mean %>%
#   filter(soce < 10)

kelp_means_PI <- left_join(adf_summary_PI, study_site_means_PI, by = c("Campaign", "site")) %>% 
  # filter(depth > -1) %>%  # Removes sites with names that don't match up
  dplyr::select(Campaign:mean_cover, depth.y, eken:icethic_cat) %>% 
  mutate(eken = ifelse(depth.x == depth.y, eken, NA),  # A funny way of getting rid of non-target depth data
         soce = ifelse(depth.x == depth.y, soce, NA), 
         toce = ifelse(depth.x == depth.y, toce, NA)) %>% 
  dplyr::select(-depth.y) %>% 
  dplyr::rename(depth = depth.x) %>% 
  gather(key = "model_var", value = "val", -c(Campaign:mean_cover)) %>% 
  na.omit()


# Multivariate analyses ---------------------------------------------------

# Need to consider percent sand and rock when looking at patterns
sand_rock_PI <- adf %>%
  filter(site %in% study_site_means_PI$site) %>% 
  select(Campaign, site, depth, Bedrock..:Shell) %>% 
  replace_na(list(Boulders.. = 0, Cobbles.. = 0, Pebbles.. = 0, sand = 0, Shell = 0)) %>% 
  group_by(Campaign, site, depth) %>% 
  summarise_all(mean, na.rm = T) 

# Create data.frame that makes vegan happy
kelp_wide_PI <- kelp_means_PI %>% 
  pivot_wider(names_from = model_var, values_from = val, values_fn = list(val = mean)) %>% 
  pivot_wider(names_from = family, values_from = mean_cover, values_fn = list(val = mean)) %>% 
  ungroup() %>% 
  left_join(sand_rock_PI, by = c("Campaign", "site", "depth"))

# The reduced version that doesn't know about depth etc.
kelp_wide_blind_PI <- kelp_wide_PI %>% 
  # select(eken:toce, Bedrock..:Shell)
  select(eken:toce)

# The "environmental" variables
kelp_wide_env_PI <- kelp_wide_PI %>% 
  select(Campaign, depth, kelp.cover, Bedrock..:Shell)

# Run the MDS
kelp_MDS_PI <- metaMDS(decostand(kelp_wide_blind_PI, method = "standardize"),
                       distance = "euclidean", try = 100, autotransform = F)
# kelp_MDS <- metaMDS(kelp_wide_blind, distance = "euclidean", try = 100)
# kelp_MDS$species <- kelp_MDS$points

# Fit environmental variables
ord_fit <- envfit(kelp_MDS_PI ~ depth + kelp.cover + Bedrock.. + Boulders.. + Cobbles.. + Pebbles.. + sand, data = kelp_wide_env_PI)
ord_fit_df <- data.frame(ord_fit$vectors$arrows) %>% 
  mutate(var = row.names(.))

# Create a data.frame for ggplot
mds_df_PI <- data.frame(site = kelp_wide_PI$site, kelp_wide_env_PI, kelp_MDS_PI$points)

# Test correlations
cor(x = kelp_wide_PI$kelp.cover, y = kelp_wide_PI$sand)
cor(x = kelp_wide_PI$kelp.cover, y = kelp_wide_PI$Boulders..)
cor(x = kelp_wide_PI$kelp.cover, y = kelp_wide_PI$Bedrock..)

# The ordiplot
ggplot(data = mds_df_PI, aes(x = MDS1, y = MDS2)) +
  # ggforce::geom_voronoi_tile(aes(group = Campaign)) +
  # stat_ellipse(aes(group = Campaign)) +
  geom_label(data = ord_fit_df, aes(label = var, x = NMDS1, y = NMDS2), size = 8) +
  # geom_label(aes(label = site))
  geom_segment(data = ord_fit_df, aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2),
               arrow = arrow(angle = 40, length = unit(0.2, "cm"), type = "open"), 
               alpha = 1, colour = "black", size = 0.5) +
  geom_point(aes(size = kelp.cover, shape = Campaign, colour = as.factor(depth))) +
  scale_colour_brewer(name = "Depth (m)", palette = "Dark2") +
  # scale_fill_brewer(name = "Campaign", palette = "Dark2") +
  scale_size_continuous(name = "Kelp cover (total %)", breaks = c(0, 25, 50, 75, 100)) +
  guides(colour = guide_legend(order = 1),
         shape = guide_legend(order = 2),
         size = guide_legend(override.aes = list(shape = 16), order = 3)) +
  # labs(size = "Duration") +
  theme_grey() +
  theme(strip.background = element_rect(fill = NA),
        panel.border = element_rect(fill = NA, colour = "black", size = 1),
        axis.text = element_text(size = 12, colour = "black"),
        axis.ticks = element_line(colour = "black"))
ggsave("graph/MDS_plot_PI.pdf", width = 10, height = 8)
ggsave("graph/MDS_plot_PI.png", width = 10, height = 8)


# Visualise ---------------------------------------------------------------

# NB: The data.frames the following visuals rely on may have changed slightly

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
# sss_month_depth <- ggplot(PI_depth_T_clim, aes(x = nav_lon, y = nav_lat)) +
#   geom_point(size = 3, shape = 19, aes(colour = soce)) +
#   borders(fill = "grey70", colour = "black") +
#   geom_point(data = PI_sites, colour = "red", aes(x = lon, y = lat)) +
#   geom_label_repel(data = PI_sites, aes(x = lon, y = lat, label = site), nudge_y = -1) +
#   scale_colour_viridis_c() +
#   coord_equal(xlim = c(-81, -76), ylim = c(71.5, 73)) +
#   labs(x = "Longitude", y = "Latitude", colour = "salinity") +
#   facet_grid(month~depth)
# ggsave("graph/sss_month_depth.pdf", sss_month_depth, height = 12, width = 25)

# Plot the kinetic energy in the area just around Pond Inlet
# kin_month_depth <- ggplot(PI_depth_T_clim, aes(x = nav_lon, y = nav_lat)) +
#   geom_point(size = 3, shape = 19, aes(colour = eken)) +
#   borders(fill = "grey70", colour = "black") +
#   geom_point(data = PI_sites, colour = "red", aes(x = lon, y = lat)) +
#   geom_label_repel(data = PI_sites, aes(x = lon, y = lat, label = site), nudge_y = -1) +
#   scale_colour_viridis_c(option = "B") +
#   coord_equal(xlim = c(-81, -76), ylim = c(71.5, 73)) +
#   labs(x = "Longitude", y = "Latitude", colour = "kinetic\nenergy\n(m2/s2)") +
#   facet_grid(month~depth)
# # kin_month_depth
# ggsave("graph/kin_month_depth.pdf", kin_month_depth, height = 12, width = 25)

# Plot the ice concentration just around Pond Inlet
# A rainbow colour palette was explicitly requested
# ice_month <- ggplot(PI_ice_clim, aes(x = nav_lon, y = nav_lat)) +
# ice_month <- ggplot(PI_ice_clim_sub, aes(x = nav_lon, y = nav_lat)) +
#   geom_point(size = 5, shape = 15, aes(colour = iceconc_cat*100)) +
#   borders(fill = "grey70", colour = "black") +
#   geom_point(data = PI_sites, colour = "red", aes(x = lon, y = lat)) +
#   geom_label_repel(data = PI_sites, aes(x = lon, y = lat, label = site), nudge_y = -1) +
#   scale_colour_gradientn(colours = rainbow_palette, 
#                          breaks = c(0, 5, 10, 15, 20),
#                          labels = c(0, 5, 10, 15, 20), 
#                          limits = c(0, 20)) +
#   coord_equal(xlim = c(-81, -76), ylim = c(71.5, 73)) +
#   labs(x = "Longitude", y = "Latitude", colour = "Ice cover (%)") +
#   facet_wrap(~month)
# # ice_month
# # ggsave("graph/ice_month.pdf", ice_month, height = 5, width = 15)
# ggsave("graph/ice_month_sub.pdf", ice_month, height = 5, width = 15)
