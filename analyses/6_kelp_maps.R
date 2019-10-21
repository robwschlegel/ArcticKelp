# 6_kelp_maps.R
# The purpose of this script is to provide a spatial visualisation of kelp coverage
# It also shows how the kelp ecosystems abiotic properties stack up against the broader arctic

# Source scripts ----------------------------------------------------------

source("analyses/4_kelp_cover.R")


# Maps --------------------------------------------------------------------

# Convenience wrapper to show kelp cover of a species laid over an abiotic variable
# Cover options are:
unique(as.character(adf_summary$family))
# Abiotic variable options are:
# unique(as.character(study_site_means_long$var))
colnames(Arctic_surface_mean)
# The function
map_cover_abiotic <- function(cover = "kelp.cover", abiotic = "sst"){
  # Subset cover
  kelp_sub <- adf_summary %>% 
    filter(family == cover) %>% 
    left_join(study_sites, by = c("Campaign", "site"))
  # Subset abiotic background variable
  abiotic_sub <- Arctic_surface_mean[,c("nav_lon", "nav_lat", abiotic)]
  colnames(abiotic_sub) <- c("lon", "lat", "var")
  # Plot it
  ggplot(data = abiotic_sub, aes(x = lon, y = lat)) +
    geom_point(aes(colour = var)) +
    borders(fill = "grey70", colour = "black") +
    geom_point(data = kelp_sub, colour = "red", 
               aes(size = mean_cover, shape = as.factor(depth))) +
    scale_colour_viridis_c(option = "D") +
    # geom_label_repel(aes(label = site)) +
    coord_cartesian(xlim = c(bbox_arctic[1], bbox_arctic[2]),
                    ylim = c(bbox_arctic[3], bbox_arctic[4])) +
    labs(x = NULL, y = NULL, colour = abiotic, shape = "Depth (m)", size = paste(cover," (%)"))
}

# Visualise all kelp cover against sst
map_cover_abiotic()
ggsave("graph/kelp_cover_vs_sst.png", width = 9, height = 7)

# Laminariales vs sss
map_cover_abiotic(cover = "Laminariales", abiotic = "sss")


# Distribution figures ----------------------------------------------------

# Convenience wrapper to show where within the distribution of an abiotic variable throughout the Arctic
# that the kelp study sites fall, and what the kelp cover at those sites are
# Cover options are:
unique(as.character(adf_summary$family))
# Abiotic variable options are:
# unique(as.character(study_site_means_long$var))
colnames(Arctic_surface_mean)
# The function
distribution_cover_abiotic <- function(cover = "kelp.cover", abiotic = "sst"){
  # Subset abiotic background variable
  abiotic_sub <- Arctic_surface_mean[,c("nav_lon", "nav_lat", abiotic)]
  colnames(abiotic_sub) <- c("lon", "lat", "var")
  abiotic_sub <- arrange(abiotic_sub, var) %>% 
    mutate(dist_index = 1:n())
  # Subset cover
  kelp_sub <- adf_summary %>% 
    filter(family == cover) %>% 
    left_join(study_sites_index, by = c("Campaign", "site")) %>% 
    left_join(abiotic_sub, by = c("nav_lon" = "lon", "nav_lat" = "lat"))
  # Find nearest 
  # Plot it
  ggplot(data = abiotic_sub, aes(x = dist_index, y = var)) +
    geom_line() +
    geom_point(data = kelp_sub, colour = "red",
               aes(size = mean_cover, shape = as.factor(depth))) +
    # scale_colour_viridis_c(option = "D") +
    geom_label_repel(data = kelp_sub, aes(label = site)) +
    # coord_cartesian(xlim = c(bbox_arctic[1], bbox_arctic[2]),
                    # ylim = c(bbox_arctic[3], bbox_arctic[4])) +
    scale_x_continuous(expand = c(0, 0)) +
    labs(x = "Rank order of all pixels", y = abiotic, colour = abiotic, shape = "Depth (m)", size = paste(cover," (%)"))
}

# Visualise all kelp cover against sst
distribution_cover_abiotic()
ggsave("graph/kelp_cover_vs_sst_distribution.png", width = 9, height = 5)

# Laminariales vs sss
distribution_cover_abiotic(cover = "Laminariales", abiotic = "sss")
