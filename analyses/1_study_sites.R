# analyses/1_study_sites.R
# The purpose of this script is to load the study site coordinate data,
# clean it up and run a test visual on the site coordinates

# NB: The workflow really begins at "analyses/3_study_site_layers.R"


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
bbox_arctic <- c(-95, -50, 50, 80)

# Load the CANA data
  # These data are from the archives of the Canadian Museum of Nature
  # Each line represents one collection
  # RWS: I originally manually edited the column 3 name, but then 
  # everything went pear shaped, so I'm leaving the file untouched
CANA_kelp <- read.csv("data/CANA_Arctic_Kelp.csv")
colnames(CANA_kelp)[3] <- "Year"


# Visualise ---------------------------------------------------------------

# Unlabelled ArcticKelp points
ggplot(data = study_sites, aes(x = lon, y = lat)) +
  borders(fill = "grey70", colour = "black") +
  geom_point(colour = "red", size = 4) +
  # geom_label_repel(aes(label = site)) +
  coord_cartesian(xlim = c(bbox_arctic[1], bbox_arctic[2]),
                  ylim = c(bbox_arctic[3], bbox_arctic[4]),
                  expand = F) +
  labs(x = NULL, y = NULL)
# ggsave("graph/study_area_points.png", width = 9, height = 7)

# Labelled ArcticKelp points
ggplot(data = study_sites, aes(x = lon, y = lat)) +
  borders(fill = "grey70", colour = "black") +
  geom_point(colour = "red", size = 4) +
  geom_label_repel(aes(label = site)) +
  coord_cartesian(xlim = c(bbox_arctic[1], bbox_arctic[2]),
                  ylim = c(bbox_arctic[3], bbox_arctic[4])) +
  labs(x = NULL, y = NULL)
# ggsave("graph/study_area_pointswlabels.png", width = 9, height = 7)

# CANA collection points
ggplot(data = CANA_kelp, aes(x = Longitude, y = Latitude)) +
  borders(fill = "grey70", colour = "black") +
  geom_point(colour = "red", size = 3) +
  # geom_label_repel(aes(label = site)) +
  coord_cartesian(xlim = c(-180, 180),
                  ylim = c(25, 85)) +
  labs(x = NULL, y = NULL)
# ggsave("graph/CANA_points.png", width = 16, height = 5)

