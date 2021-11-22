# analyses/1_study_sites.R
# The purpose of this script is to load the study site coordinate data,
# clean it up and run a test visual on the site coordinates
#note there is no 02 script. 
# NB: The workflow really begins at "analyses/3_study_site_clims.R"


# Libraries ---------------------------------------------------------------

#.libPaths(c("~/R-packages", .libPaths()))
#C:\Program Files\R

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
  mutate(Region=site,
         Region = str_replace(Region, "Black Point","North Baffin Island"),
         Region = str_replace(Region, "Iceberg Point","North Baffin Island"),
         Region = str_replace(Region,"Pooh Corner","North Baffin Island"),
         Region = str_replace(Region, "Trevor's Cove","North Baffin Island"),
         Region = str_replace(Region, "Sheaties Cabin","North Baffin Island"),
         Region = str_replace(Region, "Narwhal Cove","North Baffin Island"),
         Region = str_replace(Region, "Pangnirtung 1","Davis Strait"),
         Region = str_replace(Region, "Pangnirtung 2","Davis Strait"),
         Region = str_replace(Region, "Iqualuit T1","Davis Strait"),
         Region = str_replace(Region, "Iqualuit T2","Davis Strait"),
         Region = str_replace(Region, "Iqualuit T3","Davis Strait"),
         Region = str_replace(Region, "Iqualuit T4","Davis Strait"),
         Region = str_replace(Region, "Iqualuit T5","Davis Strait"),
         Region = str_replace(Region, "Station 14","Foxe Basin"),
         Region = str_replace(Region, "Station 17","Foxe Basin"),
         Region = str_replace(Region, "East Bay","Foxe Basin"),
         Region = str_replace(Region, "Station 20","Hudson Strait"),
         Region = str_replace(Region, "Station 21","Hudson Strait"),
         Region = str_replace(Region, "Station 22","Hudson Strait"),
         Region = str_replace(Region, "Station 25","Hudson Strait"),
         Region = str_replace(Region, "Station 3","Hudson Strait"),
         Region = str_replace(Region, "Bear Island","Hudson Strait"),
         Region = str_replace(Region, "Station 6","Roes Welcome Sound"),
         Region = str_replace(Region, "Station 8","Roes Welcome Sound"),
         Region = str_replace(Region, "Station 10","Roes Welcome Sound"),
         Region = str_replace(Region, "Station 12","Roes Welcome Sound"),
         Region = str_replace(Region, "Station 13","Roes Welcome Sound"),
         Region = str_replace(Region, "Steensby Inlet T1","Foxe Basin"),
         Region = str_replace(Region, "Deception Bay T1","Hudson Strait"),
         Region = str_replace(Region, "Deception Bay T2","Hudson Strait"),
         Region = str_replace(Region, "Deception Bay T3","Hudson Strait"),
         Region = str_replace(Region, "Deception Bay T4","Hudson Strait"),
         Region = str_replace(Region, "Deception Bay T5","Hudson Strait"),
         Region = str_replace(Region, "Deception Bay T6","Hudson Strait"),
         Region = str_replace(Region, "Steensby Inlet T2","Foxe Basin"),
         Region = str_replace(Region, "Steensby Inlet T3","Foxe Basin"),
         Region = str_replace(Region, "Steensby Inlet T4","Foxe Basin"),
         Region = str_replace(Region, "Steensby Inlet T5","Foxe Basin"),
         Region = str_replace(Region, "Steensby Inlet T6","Foxe Basin"),
         Region = str_replace(Region, "Makkovik","Labrador Sea"),
         Region = str_replace(Region, "Mussel Point","Labrador Sea"),
         Region = str_replace(Region, "Turnagain","Davis Strait"),
         Region = str_replace(Region, "Hogg Island","Davis Strait"),
         Region = str_replace(Region, "Evan's Bight","Davis Strait"),
         Region = str_replace(Region, "Duck Island","Davis Strait"),
         Region = str_replace(Region, "S. of Qik","Baffin Bay"),
         Region = str_replace(Region, "Qik","Baffin Bay"),
         Region = str_replace(Region, "Durban Harbour","Baffin Bay"),
         Region = str_replace(Region, "Rocks near Durban","Baffin Bay"),
         Region = str_replace(Region, "Churchill T2","Hudson Bay"),
         Region = str_replace(Region, "Churchill T3","Hudson Bay"),
         Region = str_replace(Region, "Churchill T4","Hudson Bay"), 
         Region = str_replace(Region, "Grise fiord","Ellesmere I"),
         Region = str_replace(Region, "Harbour fiord","Ellesmere I"),
         Region = str_replace(Region, "Vaga20.23a","Ellesmere I"),
         Region = str_replace(Region, "Starnes Fiord","Ellesmere I")
  )%>% 
  dplyr::select(site, lon, lat, everything())

# Set bounding box for study
bbox_arctic <- c(-95, -55, 50, 80)


# Visualise ---------------------------------------------------------------

# ArcticKelp points
#reorder the regions names
study_sites$Region=factor(study_sites$Region, 
    levels = c( "Ellesmere I", "North Baffin Island", "Baffin Bay","Davis Strait",                 
            "Foxe Basin"  , "Roes Welcome Sound",  "Hudson Strait"  
            ,  "Labrador Sea","Hudson Bay"  ))

study_sites$ID=1:nrow(study_sites)

#map by regions
mapfig=ggplot(data = study_sites, aes(x = lon, y = lat, col=Region)) +
  borders(fill = "grey70", colour = "black") +
  geom_point( size = 2.5) +
 # geom_label_repel(aes(label = ID), size =2) +
  coord_cartesian(xlim = c(bbox_arctic[1], bbox_arctic[2]),
                  ylim = c(bbox_arctic[3], bbox_arctic[4])) +
  labs(x = NULL, y = NULL)

#save 
#tiff("figures/study_area_pointswlabels.tiff", width = 16, height = 12, units = 'cm', res = 300)
#mapfig # Make plot
#dev.off()

#with names
mapnames= ggplot(data = study_sites, aes(x = lon, y = lat, col=Campaign)) +
  borders(fill = "grey70", colour = "black") +
   geom_label_repel(aes(label = site), size=3, max.overlaps =25) +
  geom_point(size = 3, aes(col=Campaign)) +
  coord_cartesian(xlim = c(bbox_arctic[1], bbox_arctic[2]),
                  ylim = c(bbox_arctic[3], bbox_arctic[4])) +
  theme(text = element_text(size = 20), legend.title = element_blank())+
  labs(x = expression(paste('latitude (', degree, ')')), y = expression(paste('longitude (', degree,')')))
mapnames

#tiff("figures/study_area_pointswlabels.tiff", width = 28, height = 22,units = 'cm', res = 300)
#mapnames # Make plot
#dev.off()


#PICKER sites
ggplot(data = study_sites%>%filter(Campaign=='PICKeR'), aes(x = lon, y = lat)) +
  borders(fill = "grey70", colour = "black") +
  geom_point(colour = "red", size = 3) +
  geom_label_repel(aes(label = site), size=3) +
  coord_cartesian(xlim = c(-81,-76),
                  ylim = c(71.5, 74)) +
  labs(x = NULL, y = NULL)


#Vaga sites
ggplot(data = study_sites%>%filter(Campaign=='Vagabond'), aes(x = lon, y = lat)) +
  borders(fill = "grey70", colour = "black") +
  geom_point(colour = "red", size = 3) +
  geom_label_repel(aes(label = site), size=3) +
  coord_cartesian(xlim = c(-85,-81),
                  ylim = c(75.5, 77.5)) +
  labs(x = NULL, y = NULL)

#SIMEP sites
ggplot(data = study_sites%>%filter(Campaign=='SIMEP'), aes(x = lon, y = lat)) +
  borders(fill = "grey70", colour = "black") +
  geom_point(colour = "red", size = 3) +
  geom_label_repel(aes(label = site), size=3) +
  coord_cartesian(xlim = c(-89,-78),
                  ylim = c(61.5, 67)) +
  labs(x = NULL, y = NULL)

