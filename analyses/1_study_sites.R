# analyses/1_study_sites.R
# The purpose of this script is to load the study site coordinate data,
# clean it up and run a test visual on the site coordinates

# NB: The workflow really begins at "analyses/3_study_site_clims.R"


# Libraries ---------------------------------------------------------------

.libPaths(c("~/R-packages", .libPaths()))
# library(rlang, lib.loc = "~/R-packages")
library(tidyverse)
library(ggrepel)
library(measurements)


# Load data ---------------------------------------------------------------

# The study sites
study_sites <- readxl::read_xlsx("metadata/Kelp cover photograph 2019 GPS locations.xlsx") %>% 
  dplyr::rename(lon = Longitude, lat = Latitude, site = Site) %>% 
  mutate(lon = round(as.numeric(lon), 4),
         lat = round(as.numeric(lat), 4),
         site = str_replace(site, "Durban Habour", "Durban Harbour")) %>% 
  filter(Notes != "not entered") %>% 
  dplyr::select(site, lon, lat, everything())

# Set bounding box for study
bbox_arctic <- c(-90, -50, 50, 80)


# Visualise ---------------------------------------------------------------

# ggplot(data = study_sites, aes(x = lon, y = lat)) +
#   borders(fill = "grey70", colour = "black") +
#   geom_point(colour = "red", size=4) +
#  # geom_label_repel(aes(label = site)) +
#   coord_cartesian(xlim = c(bbox_arctic[1], bbox_arctic[2]),
#                   ylim = c(bbox_arctic[3], bbox_arctic[4])) +
#   labs(x = NULL, y = NULL)
# ggsave("graph/study_area_points.png", width = 9, height = 7)

# ggplot(data = study_sites, aes(x = lon, y = lat)) +
#   borders(fill = "grey70", colour = "black") +
#   geom_point(colour = "red", size=4) +
#   geom_label_repel(aes(label = site)) +
#   coord_cartesian(xlim = c(bbox_arctic[1], bbox_arctic[2]),
#                   ylim = c(bbox_arctic[3], bbox_arctic[4])) +
#   labs(x = NULL, y = NULL)
# ggsave("graph/study_area_pointswlabels.png", width = 9, height = 7)

 