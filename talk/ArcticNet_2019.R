# talk/ArcticNet_2019.R
# The purpose of this script is to house all of the code used to 
# create the figures etc. used in the talk with the same name


# Setup -------------------------------------------------------------------

source("analyses/4_kelp_cover.R")

# The base map to use for everything else
Arctic_map <- ggplot() +
  borders(fill = "grey70", colour = "black") +
  coord_cartesian(xlim = c(bbox_arctic[1], bbox_arctic[2]),
                  ylim = c(bbox_arctic[3], bbox_arctic[4])) +
  labs(x = NULL, y = NULL)


# Data --------------------------------------------------------------------

alaria_trading_points <- data.frame(lon = c(-67.491047),
                                   lat = c(66.552445))

adf_summary_coords <- left_join(adf_summary, study_sites, by = c("Campaign", "site"))

# What is known of Arctic Kelp? -------------------------------------------

kelp_question <- Arctic_map +
  geom_label(aes(x = -70, y = 65, label = "?"), size = 60, alpha = 0.6)
kelp_question
ggsave(kelp_question, filename = "talk/figure/kelp_question.png", height = 6, width = 6)

alaria_trading <- Arctic_map +
  geom_point(data = alaria_trading_points, aes(x = lon, y = lat), colour = "brown", size = 10)
alaria_trading
ggsave(alaria_trading, filename = "talk/figure/alaria_trading.png", height = 6, width = 6)

where_kelp <- Arctic_map +
  geom_label(aes(x = -80, y = 60, label = "?"), size = 40, alpha = 0.6, colour = "green") +
  geom_label(aes(x = -78, y = 70, label = "?"), size = 40, alpha = 0.6, colour = "brown") +
  geom_label(aes(x = -60, y = 73, label = "?"), size = 40, alpha = 0.6, colour = "blue")
where_kelp
ggsave(where_kelp, filename = "talk/figure/where_kelp.png", height = 6, width = 6)


# ArcticKelp --------------------------------------------------------------

study_site_labels <- Arctic_map +
  geom_point(data = study_sites, aes(x = lon, y = lat), colour = "red", size = 4) +
  geom_label_repel(data = study_sites, aes(x = lon, y = lat, label = site), alpha = 0.8)
study_site_labels
ggsave(study_site_labels, filename = "talk/figure/study_site_labels.png", width = 6, height = 6)

study_site_mean_cover <- Arctic_map +
  geom_point(data = filter(adf_summary_coords, family == "kelp.cover"), 
             aes(x = lon, y = lat, colour = mean_cover, shape = as.factor(depth)), size = 4) +
  scale_colour_viridis_c() +
  scale_shape_manual(values = c(17, 16, 15, 18)) +
  labs(colour = "Mean cover (%)", shape = "Depth (m)") +
  guides(color = guide_colourbar(order = 1)) +
  theme(strip.text.x = element_blank(),
        strip.background = element_rect(colour = "white", fill = "white"),
        legend.position = c(0.289, 0.08),
        legend.direction = "horizontal",
        legend.spacing.y = unit(0, "mm"))
study_site_mean_cover
ggsave(study_site_mean_cover, filename = "talk/figure/study_site_mean_cover.png", width = 6, height = 6)

