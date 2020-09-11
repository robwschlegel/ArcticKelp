# analyses/6_kelp_maps.R
# The purpose of this script is to provide a spatial visualisation of kelp coverage
# It also shows how the kelp ecosystems abiotic properties stack up against the broader Arctic


# Source scripts ----------------------------------------------------------

source("analyses/4_kelp_cover.R")

load("data/Arctic_BO.RData")
load("data/Arctic_AM.RData")

Arctic_env <- right_join(Arctic_BO, Arctic_AM)


# Maps --------------------------------------------------------------------

# Convenience wrapper to show kelp cover of a species laid over an abiotic variable

# Cover options are:
unique(as.character(adf_summary$family))

# Abiotic variable options are:
colnames(Arctic_env)

# The function
map_cover_abiotic <- function(cover = "kelp.cover", 
                              abiotic = "BO2_tempmean_bdmax"){
  # Subset cover
  kelp_sub <- adf_summary %>% 
    filter(family == cover) %>% 
    left_join(study_sites, by = c("Campaign", "site"))
  
  # Subset abiotic background variable
  abiotic_sub <- Arctic_env[,c("lon", "lat", abiotic)]
  
  # Plot it
  ggplot(data = abiotic_sub, aes(x = lon, y = lat)) +
    geom_tile(aes_string(fill = abiotic)) +
    borders(fill = "grey70", colour = "black") +
    geom_point(data = kelp_sub, colour = "red", 
               aes(size = mean_cover)) +
    scale_fill_viridis_c(option = "D") +
    coord_cartesian(xlim = c(bbox_arctic[1], bbox_arctic[2]),
                    ylim = c(bbox_arctic[3], bbox_arctic[4]),
                    expand = F) +
    labs(x = NULL, y = NULL, fill = abiotic, size = paste(cover," (%)"))
}

# Visualise some different things
map_cover_abiotic()
map_cover_abiotic(cover = "Laminariales", abiotic = "BO2_curvelmean_bdmax")
map_cover_abiotic(cover = "Agarum", abiotic = "BO2_icethickmean_ss")
map_cover_abiotic(cover = "Alaria", abiotic = "BO_damean")


# Distribution figures ----------------------------------------------------

# Convenience wrapper to show where within the distribution of an abiotic variable
# the kelp study sites fall, and what the kelp cover at those sites are

# Cover options are:
unique(as.character(adf_summary$family))

# Abiotic variable options are:
colnames(Arctic_env)

# The function
distribution_cover_abiotic <- function(cover = "kelp.cover", 
                                       abiotic = "BO2_tempmean_bdmax"){
  # Subset abiotic background variable
  abiotic_sub <- Arctic_env[,c("lon", "lat", abiotic)]
  colnames(abiotic_sub) <- c("lon_env", "lat_env", "var")
  abiotic_sub <- arrange(abiotic_sub, var) %>% 
    mutate(dist_index = 1:n()) %>% 
    na.omit()
  
  # Subset cover
  kelp_sub <- adf_summary %>% 
    filter(family == cover) %>% 
    left_join(study_site_env, by = c("Campaign", "site")) %>% 
    left_join(abiotic_sub, by = c("lon_env", "lat_env"))
  
  # Plot it
  ggplot(data = abiotic_sub, aes(x = dist_index, y = var)) +
    geom_line() +
    geom_point(data = kelp_sub, size = 4,
               aes(colour = mean_cover, shape = as.factor(depth))) +
    scale_colour_gradient(low = "yellow", high = "red") +
    scale_x_continuous(expand = c(0, 0)) +
    labs(x = "Rank order of all pixels", y = abiotic, colour = cover, 
         shape = "Depth (m)", size = paste(cover," (%)"))
}

# Visualise all kelp cover against sst
distribution_cover_abiotic()
distribution_cover_abiotic(cover = "Alaria", abiotic = "BO2_curvelltmax_bdmax")
distribution_cover_abiotic(cover = "Agarum", abiotic = "BO2_salinitymean_bdmax")
distribution_cover_abiotic(cover = "Laminariales", abiotic = "BO2_icethickmean_ss")

