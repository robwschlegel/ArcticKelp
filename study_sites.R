# study_sites.R
# The purpose of this script is to load the study site coordinate data,
# clean it up and run a test visual on the site coordinates


# Libraries ---------------------------------------------------------------

library(tidyverse)
library(ggrepel)
library(measurements)


# Load data ---------------------------------------------------------------

# The study sites
study_sites <- readxl::read_xlsx("data/ArcticKelp study sites.xlsx") %>% 
  dplyr::rename(lon = Longitude, lat = Latitude, site = Site_Name)

# Manually look at the spreadsheet to see which rows need to be converted to decimal lon/lat
con_rows <- c(11:20, 23:25)

# convert from decimal minutes to decimal degrees
study_sites$lon[con_rows] <- conv_unit(study_sites$lon[con_rows] , from = 'deg_dec_min', to = 'dec_deg')
study_sites$lat[con_rows] <- conv_unit(study_sites$lat[con_rows] , from = 'deg_dec_min', to = 'dec_deg')

# Convert coordinates to numeric values
study_sites <- study_sites %>% 
  mutate(lon = round(as.numeric(lon), 4),
         lat = round(as.numeric(lat), 4),
         # It seems that there are several coordinates that should be negative
         lon = -abs(lon))

# Set bounding box for study
bbox_arctic <- c(-100, -40, 50, 80)


# Visualise ---------------------------------------------------------------

ggplot(data = study_sites, aes(x = lon, y = lat)) +
  borders(fill = "grey70", colour = "black") +
  geom_point(colour = "red") +
  geom_label_repel(aes(label = site)) +
  coord_cartesian(xlim = c(bbox_arctic[1], bbox_arctic[2]),
                  ylim = c(bbox_arctic[3], bbox_arctic[4])) +
  labs(x = NULL, y = NULL)
ggsave("study_area.png", width = 9, height = 7)
