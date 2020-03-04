# analyses/3_study_site_clims
# The purpose of this script is to pull out the monthly clims and mean states at each study site
# Whenever new sites are added re-run this script to get the clims for the new pixels


# Setup -------------------------------------------------------------------

# Load study sites and base packages
source("analyses/1_study_sites.R")

# Additional packages
library(FNN)

load("data/Arctic_bathy.RData")


# Study site BO data ------------------------------------------------------

# NB: The following code chunks require Arctic_BO.RData
# This is not hosted on GitHub as it is 27.7 MB
# E-mail robert.schlegel@dal.ca for the file
# Or download it from '2_monthly_clims.R' ln 78-89
load("data/Arctic_BO.RData")
Arctic_BO$BO_index <- 1:nrow(Arctic_BO)

# Find the nearest BO points to each site and add bathy data
study_site_BO <- study_sites %>%
  mutate(BO_index = as.vector(knnx.index(as.matrix(Arctic_BO[,c("lon", "lat")]),
                                         as.matrix(study_sites[,c("lon", "lat")]), k = 1))) %>%
  left_join(Arctic_BO, by = "BO_index") %>%
  dplyr::select(-BO_index, -Date, -Notes) %>%
  dplyr::rename(lon = lon.x, lat = lat.x, lon_BO = lon.y, lat_BO = lat.y) %>% 
  left_join(Arctic_bathy, by = c("lon_BO" = "lon", "lat_BO" = "lat"))
save(study_site_BO, file = "data/study_site_BO.RData")

# Melt for plotting purposes
study_site_BO_long <- study_site_BO %>%
  gather("var", "val", -c(site:lat_BO)) %>%
  na.omit()
save(study_site_BO_long, file = "data/study_site_BO_long.RData")


# Study site index --------------------------------------------------------

# Find the nearest NAPA points to each site
study_site_index <- study_sites %>% 
  mutate(BO_index = as.vector(knnx.index(as.matrix(Arctic_BO[,c("lon", "lat")]),
                                         as.matrix(study_sites[,c("lon", "lat")]), k = 1))) %>% 
  left_join(Arctic_BO[,c("lon", "lat", "BO_index")], by = "BO_index") %>% 
  dplyr::rename(lon = lon.x, lat = lat.x, lon_BO = lon.y, lat_BO = lat.y)
save(study_site_index, file = "data/study_site_index.RData")


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

