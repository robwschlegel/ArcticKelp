# kelp_cover.R
# The purpose of this script is to load the kelp cover data
# It then goes about making some summaries and visualisations


# Libraries ---------------------------------------------------------------

library(tidyverse)
library(ggridges)


# Load data ---------------------------------------------------------------

adf <- read.csv("data/Kelp cover photograph quadrats 2019.csv", dec = ',', sep = ';') %>%
  dplyr::rename(site = Site, Depth = Depth.m, 
                sand = Sand.., Agarum = Agarum.., Alaria = Alaria..) %>% 
  filter(sand <= 75, Depth >= 3) %>%
  mutate(kelp.density = S.latissima.No. + Agarum.No. + Alaria.No.,
         kelp.cover = Laminariales + Agarum + Alaria,
         site = str_replace(site, "Durban Habour", "Durban Harbour"),
         site = str_replace(site, "Qik", "Qikitarjuak"))


# Summarise data ----------------------------------------------------------

adf_summary <- adf %>% 
  dplyr::select(Campaign, site, Depth, kelp.cover, Laminariales, Agarum, Alaria) %>% 
  mutate(Depth = ifelse(Depth == 3, 5, Depth), # Match abiotic 5 m Depth increments
         site = as.character(site)) %>% 
  gather(key = "family", value = "cover", -Campaign, -site, -Depth) %>% 
  mutate(family = factor(family, levels = c("Agarum", "Alaria", 
                                            "Laminariales", "kelp.cover"))) %>% 
  group_by(Campaign, site, Depth, family) %>% 
  summarise(mean_cover = round(mean(cover, na.rm = T), 2),
            sd_cover = round(sd(cover, na.rm = T), 2),
            count = n())


# Create figures ----------------------------------------------------------

# Ridge plot showing percent kelp cover in quadrats by site and Depth
# fig1 <- ggplot(adf, aes(y = site, x = kelp.cover, fill = as.factor(Depth))) +
#   stat_density_ridges(alpha = 0.7) +
#   labs(y = "Site", x = expression(Kelp~cover~('%')), fill = "Depth") +
#   theme_bw(base_size = 14) +
#   theme(legend.position = 'top')
# fig1
# ggsave("graph/kelp_cover_vs_sites.png", fig1, width = 8, height = 10, units = "in", dpi = 300)

# Ridge plots showing coverage of Laminariales by site and Depth
# fig2 <- ggplot(adf, aes(y = site, x = Laminariales, fill = as.factor(Depth))) +
#   stat_density_ridges(alpha = 0.7) +
#   labs(y = "Site", x = expression(Laminariales~cover~('%')), fill = "Depth") +
#   theme_bw(base_size = 14) +
#   theme(legend.position = 'top')
# fig2
# ggsave("graph/Laminariales_cover_vs_sites.png", fig2, width = 8, height = 10, units = "in", dpi = 300)

# Ridge plots showing coverage of Agarum by site and Depth
# fig3 <- ggplot(adf, aes(y = site, x = Agarum, fill = as.factor(Depth))) +
#   stat_density_ridges(alpha = 0.7) +
#   labs(y = "Site", x = expression(Agarum~cover~('%')), fill = "Depth") +
#   theme_bw(base_size = 14) +
#   theme(legend.position = 'top')
# fig3
# ggsave("graph/Agarum_cover_vs_sites.png", fig3, width = 8, height = 10, units = "in", dpi = 300)

