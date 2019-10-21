# analyses/kelp_cover.R
# The purpose of this script is to load the kelp cover data
# It then goes about making some summaries and visualisations


# Libraries ---------------------------------------------------------------

source("analyses/3_study_site_clims.R")

library(ggridges)


# Load data ---------------------------------------------------------------

adf <- read.csv("data/Kelp cover photograph quadrats 2019.csv", dec = ',', sep = ';') %>%
  dplyr::rename(site = Site, depth = Depth.m, 
                sand = Sand.., Agarum = Agarum.., Alaria = Alaria..) %>% 
  filter(sand <= 75, depth >= 3) %>%
  mutate(depth = ifelse(depth == 3, 5, depth), # Match abiotic 5 m depth increments
         kelp.density = S.latissima.No. + Agarum.No. + Alaria.No.,
         kelp.cover = Laminariales + Agarum + Alaria,
         Bedrock.. = replace_na(Bedrock.., 0),
         Boulders.. = replace_na(Boulders.., 0),
         Cobbles.. = replace_na(Cobbles.., 0),
         Pebbles.. = replace_na(Pebbles.., 0),
         rock = Bedrock.. + Boulders.. + Cobbles.. + Pebbles..,
         site = str_replace(site, "Durban Habour", "Durban Harbour"),
         site = as.character(site),
         too_little_cover = ifelse(sand + rock < 100, TRUE, FALSE),
         too_much_cover = ifelse(sand + rock > 100, TRUE, FALSE))


# Summarise data ----------------------------------------------------------

adf_summary <- adf %>% 
  dplyr::select(Campaign, site, depth, kelp.cover, Laminariales, Agarum, Alaria) %>% 
  gather(key = "family", value = "cover", -Campaign, -site, -depth) %>% 
  mutate(family = factor(family, levels = c("Agarum", "Alaria", 
                                            "Laminariales", "kelp.cover"))) %>% 
  group_by(Campaign, site, depth, family) %>% 
  summarise(mean_cover = round(mean(cover, na.rm = T), 2),
            sd_cover = round(sd(cover, na.rm = T), 2),
            count = n())


# Create figures ----------------------------------------------------------

# Ridge plot showing percent kelp cover in quadrats by site and depth
# fig1 <- ggplot(adf, aes(y = site, x = kelp.cover, fill = as.factor(depth))) +
#   stat_density_ridges(alpha = 0.7) +
#   labs(y = "Site", x = expression(Kelp~cover~('%')), fill = "depth") +
#   theme_bw(base_size = 14) +
#   theme(legend.position = 'top')
# fig1
# ggsave("graph/kelp_cover_vs_sites.png", fig1, width = 8, height = 10, units = "in", dpi = 300)

# Ridge plots showing coverage of Laminariales by site and depth
# fig2 <- ggplot(adf, aes(y = site, x = Laminariales, fill = as.factor(depth))) +
#   stat_density_ridges(alpha = 0.7) +
#   labs(y = "Site", x = expression(Laminariales~cover~('%')), fill = "depth") +
#   theme_bw(base_size = 14) +
#   theme(legend.position = 'top')
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

