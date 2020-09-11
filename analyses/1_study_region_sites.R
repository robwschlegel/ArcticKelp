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
# library(rgdal) # Used for processing the Arctic shape file


# Load site data ----------------------------------------------------------

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
# ggplot(data = study_sites, aes(x = lon, y = lat)) +
#   borders(fill = "grey70", colour = "black") +
#   geom_point(colour = "red", size = 4) +
#   # geom_label_repel(aes(label = site)) +
#   coord_cartesian(xlim = c(bbox_arctic[1], bbox_arctic[2]),
#                   ylim = c(bbox_arctic[3], bbox_arctic[4]),
#                   expand = F) +
#   labs(x = NULL, y = NULL)
# ggsave("graph/study_area_points.png", width = 9, height = 7)

# Labelled ArcticKelp points
# ggplot(data = study_sites, aes(x = lon, y = lat)) +
#   borders(fill = "grey70", colour = "black") +
#   geom_point(colour = "red", size = 4) +
#   geom_label_repel(aes(label = site)) +
#   coord_cartesian(xlim = c(bbox_arctic[1], bbox_arctic[2]),
#                   ylim = c(bbox_arctic[3], bbox_arctic[4])) +
#   labs(x = NULL, y = NULL)
# ggsave("graph/study_area_pointswlabels.png", width = 9, height = 7)

# CANA collection points
# ggplot(data = CANA_kelp, aes(x = Longitude, y = Latitude)) +
#   borders(fill = "grey70", colour = "black") +
#   geom_point(colour = "red", size = 3) +
#   # geom_label_repel(aes(label = site)) +
#   coord_cartesian(xlim = c(-180, 180),
#                   ylim = c(25, 85)) +
#   labs(x = NULL, y = NULL)
# ggsave("graph/CANA_points.png", width = 16, height = 5)


# Load Arctic shape file --------------------------------------------------

# Load the Arctic study region shape file
# Arctic_poly <- readOGR(dsn = "metadata/", layer = "amaplim_lam_poly")
# plot(Arctic_poly)

# Convert to a more useful coordinate system
# Arctic_flat <- spTransform(Arctic_poly, "+init=epsg:4326 +proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0")
# plot(Arctic_flat)

# Convert that to a data.frame
Arctic_boundary <- as.data.frame(Arctic_flat@polygons[[1]]@Polygons[[1]]@coords) %>%
  `colnames<-`(c("lon", "lat")) %>%
  arrange(lon)
Arctic_boundary <- rbind(Arctic_boundary,
                         data.frame(lon = rev(Arctic_boundary$lon),
                                    lat = rep(90, nrow(Arctic_boundary))))

# Manually remove points that are too close to the coast for the filtering of the data layers in "2_study_region_layers.R"
# The contours of St James Bay confuse the machine when it is looking for points within the polygon
Arctic_boundary$lat[340:444] <- Arctic_boundary$lat[340:444]-1
Arctic_boundary$lon[405:444] <- Arctic_boundary$lon[405:444]-1
Arctic_boundary$lon[445:683] <- NA
Arctic_boundary <- na.omit(Arctic_boundary)

# Visualise
ggplot(Arctic_boundary, aes(x = lon, y = lat)) +
  geom_point() + geom_polygon() + borders() +
  geom_point(data = Arctic_boundary[405,], colour = "red", size = 5) +
  coord_quickmap(ylim = c(50, 90), expand = F)

# Save
save(Arctic_boundary, file = "metadata/Arctic_boundary.RData")

# Load
load("metadata/Arctic_boundary.RData")

