# analyses/4_kelp_cover.R
# The purpose of this script is to load the kelp cover data
# It then goes about making some summaries and visualisations


# Libraries ---------------------------------------------------------------

# If new sites have been added it is necessary to run this script
# source("analyses/3_study_site_clims.R")

# Load study sites and base packages
source("analyses/1_study_sites.R")

library(ggridges)

# Convenience loading and coord rounding function
round_coords <- function(df){
  df_round <- data.frame(df) %>%
    mutate(nav_lon = round(nav_lon, 4),
           nav_lat = round(nav_lat, 4))
  if(!("depth" %in% colnames(df_round))) df_round$depth <- 0
  return(df_round)
}


# Load data ---------------------------------------------------------------

# If the following excel file is opened in Excel on a Windows computer it will 
# change the comma (,) separation into semicolon (;) separation and the commented out 
# bit in the following line needs to be turned back on
adf <- read.csv("data/Kelp cover photograph quadrats 2019.csv", sep=';', dec=',') %>%
  dplyr::rename(site = Site, depth = Depth.m, 
                sand = Sand.., Agarum = Agarum.., Alaria = Alaria.., Alaria = Alaria.., Sacharrina = S..latissima.., 
                L.solidungula = L..solid.., L.digitata = L..digitata..) %>% 
  filter(depth >= 3) %>%
  mutate(depth = ifelse(depth == 3, 5, depth), # Match abiotic 5 m depth increments
         Agarum = replace_na(Agarum, 0),
         L.solidungula = replace_na(L.solidungula, 0),
         L.digitata = replace_na(L.digitata, 0),
         Alaria = replace_na(Alaria, 0),
         Laminariales.unspecified = replace_na(Laminariales.unspecified, 0),
         Sacharrina = replace_na(Sacharrina, 0),
         Saccorhiza.. = replace_na(Saccorhiza.., 0),
         kelp.density = S.latissima.No. + Agarum.No. + Alaria.No.,
         kelp.cover = Agarum + Alaria + L.solidungula + L.digitata + Laminariales.unspecified + Sacharrina,
         Laminariales = L.solidungula + L.digitata + Laminariales.unspecified,
         Bedrock.. = replace_na(Bedrock.., 0),
         Boulders.. = replace_na(Boulders.., 0),
         Cobbles.. = replace_na(Cobbles.., 0),
         Pebbles.. = replace_na(Pebbles.., 0),
         rock = Bedrock.. + Boulders.. + Cobbles.. + Pebbles..,
         site = str_replace(site, "Durban Habour", "Durban Harbour"),
         site = as.character(site),
         too_little_cover = ifelse(sand + rock < 100, TRUE, FALSE),
         too_much_cover = ifelse(sand + rock > 100, TRUE, FALSE),
         Campaign = as.character(Campaign))

load("data/study_site_means.RData")

load_Arctic_mean <- function(){
  load("data/Arctic_surface_mean.RData")
  Arctic_surface_mean <- round_coords(Arctic_surface_mean) %>% 
    dplyr::select(-emp_ice) %>%  # Rather use the emp_ice from the ice data
    dplyr::select(-qla_oce, -qsb_oce) # These two variables are empty
  load("data/Arctic_ice_mean.RData")
  Arctic_ice_mean <- round_coords(Arctic_ice_mean) %>% 
    dplyr::select(-emp_oce, -qemp_oce) # Rather use the emp_oce and qemp_oce from the sea data
  load("data/Arctic_depth_T_mean.RData")
  Arctic_depth_T_mean <- round_coords(Arctic_depth_T_mean)
  load("data/Arctic_depth_U_mean.RData")
  Arctic_depth_U_mean <- round_coords(Arctic_depth_U_mean)
  load("data/Arctic_depth_V_mean.RData")
  Arctic_depth_V_mean <- round_coords(Arctic_depth_V_mean)
  load("data/Arctic_depth_W_mean.RData")
  Arctic_depth_W_mean <- round_coords(Arctic_depth_W_mean)
  
  # Join everything
  Arctic_mean <- left_join(Arctic_depth_T_mean, Arctic_depth_U_mean, by = c("x", "y", "nav_lon", "nav_lat", "depth", "bathy")) %>%
    left_join(Arctic_depth_V_mean, by = c("x", "y", "nav_lon", "nav_lat", "depth", "bathy")) %>% 
    left_join(Arctic_depth_W_mean, by = c("x", "y", "nav_lon", "nav_lat", "depth", "bathy")) %>% 
    left_join(Arctic_surface_mean, by = c("x", "y", "nav_lon", "nav_lat", "depth", "bathy")) %>% 
    left_join(Arctic_ice_mean, by = c("x", "y", "nav_lon", "nav_lat", "depth", "bathy"))#, "emp_ice", "emp_oce", "qemp_oce"))
  
  # Clean up
  rm(Arctic_depth_T_mean, Arctic_depth_U_mean, Arctic_depth_V_mean, 
     Arctic_depth_W_mean, Arctic_surface_mean, Arctic_ice_mean)
  return(Arctic_mean)
}
Arctic_mean <- load_Arctic_mean()


# Quadrat data ------------------------------------------------------------

adf_quadrat <- adf %>% 
  dplyr::select(Campaign, site, depth, Quadrat, kelp.cover, Laminariales, Agarum, Alaria) %>% 
  gather(key = "family", value = "cover", -Campaign, -site, -depth, -Quadrat) %>% 
  mutate(family = factor(family, levels = c("Agarum", "Alaria", 
                                            "Laminariales", "Sacharrina", "kelp.cover")))

# Summarise data ----------------------------------------------------------

adf_summary <- adf %>% 
  dplyr::select(Campaign, site, depth, kelp.cover, Laminariales, Agarum, Alaria) %>% 
  gather(key = "family", value = "cover", -Campaign, -site, -depth) %>% 
  mutate(family = factor(family, levels = c("Agarum", "Alaria", 
                                            "Laminariales", "Sacharrina", "kelp.cover"))) %>% 
  group_by(Campaign, site, depth, family) %>% 
  summarise(mean_cover = round(mean(cover, na.rm = T), 2),
            sd_cover = round(sd(cover, na.rm = T), 2),
            count = n()) %>% 
  mutate(mean_cover = ifelse(is.na(mean_cover), 0, mean_cover),
         sd_cover = ifelse(is.na(sd_cover), 0, sd_cover))


# Create figures ----------------------------------------------------------

# Ridge plot showing percent kelp cover in quadrats by site and depth
fig1 <- ggplot(adf, aes(y = site, x = kelp.cover, fill = as.factor(depth))) +
  stat_density_ridges(alpha = 0.7) +
  labs(y = "Site", x = expression(Kelp~cover~('%')), fill = "depth") +
  theme_bw(base_size = 14) +
  theme(legend.position = 'top')
# fig1
# ggsave("graph/kelp_cover_vs_sites.png", fig1, width = 8, height = 10, units = "in", dpi = 300)

# Ridge plots showing coverage of Laminariales by site and depth
fig2 <- ggplot(adf, aes(y = site, x = Laminariales, fill = as.factor(depth))) +
  stat_density_ridges(alpha = 0.7) +
  labs(y = "Site", x = expression(Laminariales~cover~('%')), fill = "depth") +
  theme_bw(base_size = 14) +
  theme(legend.position = 'top')
# fig2
# ggsave("graph/Laminariales_cover_vs_sites.png", fig2, width = 8, height = 10, units = "in", dpi = 300)

# Ridge plots showing coverage of Agarum by site and depth
# fig3 <- ggplot(adf, aes(y = site, x = Agarum, fill = as.factor(depth))) +
#   stat_density_ridges(alpha = 0.7) +
#   labs(y = "Site", x = expression(Agarum~cover~('%')), fill = "depth") +
#   theme_bw(base_size = 14) +
#   theme(legend.position = 'top')
# fig3
# ggsave("graph/Agarum_cover_vs_sites.png", fig3, width = 8, height = 10, units = "in", dpi = 300)

