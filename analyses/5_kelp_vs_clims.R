# analyses/5_kelp_vs_clims.R
# The purpose of this script is to load the kelp and monthly clim data
# and then squish them together to see if anything interesting falls out


# Source scripts ----------------------------------------------------------

source("analyses/4_kelp_cover.R")

library(vegan)


# Check for site mismatches -----------------------------------------------

# NB: The following two lines of code should return 'character(0)'
# This means there are no differences in the sites for the ArctcKelp and env data
# If there are, re-run '3_study_site_clims.R'

# Sites from Kelp data not in study site names
unique(adf_summary$site)[!(unique(adf_summary$site) %in% unique(study_site_env$site))]

# Study site names not in kelp data
unique(study_site_env$site)[!(unique(study_site_env$site) %in% unique(adf_summary$site))]


# Join kelp to env variables ----------------------------------------------

# Mean kelp cover and env variables
kelp_means <- left_join(adf_summary, study_site_env, by = c("Campaign", "site")) %>% 
  ungroup()


# Multivariate analyses ---------------------------------------------------

# Need to consider percent sand and rock when looking at patterns
sand_rock <- adf %>% 
  dplyr::select(Campaign, site, depth, sand, rock) %>% 
  group_by(Campaign, site, depth) %>% 
  summarise_all(mean, na.rm = T)

# Create data.frame that makes vegan happy
kelp_wide <- kelp_means %>% 
  left_join(sand_rock, by = c("Campaign", "site", "depth"))

# Visualise
ggplot(data = kelp_wide, aes(x = BO2_icecovermean_ss)) +
  geom_histogram(bins = 6) +
  geom_point(aes(colour = mean_cover, shape = family), y = 0, size = 4) +
  scale_colour_viridis_c() +
  labs(x = "Mean ice cover")

# The env variables
   # NB: Ignoring the GMED variables for now as there are missing values
kelp_wide_var <- kelp_wide %>% 
  dplyr::select(BO2_templtmin_bdmax:BO2_curvelltmax_bdmax)

# The response variables
kelp_wide_env <- kelp_wide %>% 
  dplyr::select(Campaign, depth, family, mean_cover, sd_cover, rock, sand)

# Run the MDS
kelp_MDS <- metaMDS(decostand(kelp_wide_var, method = "standardize"),
                    distance = "euclidean", try = 100, autotransform = F)

# Fit environmental variables
ord_fit <- envfit(formula = kelp_MDS ~ depth + mean_cover + rock + sand, data = kelp_wide_env)
# plot(kelp_MDS)
# plot(ord_fit)
# ordiArrowMul(ord_fit)
ord_fit_df <- data.frame(ord_fit$vectors$arrows) %>% 
  mutate(NMDS1 = NMDS1 * ord_fit$vectors$r * ordiArrowMul(ord_fit),
         NMDS2 = NMDS2 * ord_fit$vectors$r * ordiArrowMul(ord_fit),
         var = row.names(.))

# Create a data.frame for ggplot
mds_df <- data.frame(site = kelp_wide$site, kelp_wide_env, kelp_MDS$points)

# Test correlations
cor(x = kelp_wide$mean_cover, y = kelp_wide$sand)
cor(x = kelp_wide$mean_cover, y = kelp_wide$rock)

# The ordiplot
ggplot(data = mds_df, aes(x = MDS1, y = MDS2)) +
  # ggforce::geom_voronoi_tile(aes(group = Campaign)) +
  # stat_ellipse(aes(group = Campaign)) +
  geom_label(data = ord_fit_df, aes(label = var, x = NMDS1, y = NMDS2), size = 8) +
  # geom_label(aes(label = site))
  geom_segment(data = ord_fit_df, aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2),
               arrow = arrow(angle = 40, length = unit(0.2, "cm"), type = "open"), 
               alpha = 1, colour = "black", size = 0.5) +
  geom_point(aes(size = mean_cover, colour = as.factor(Campaign))) +
  scale_colour_brewer(name = "Campaign", palette = "Dark2") +
  # scale_fill_brewer(name = "Campaign", palette = "Dark2") +
  scale_size_continuous(name = "Kelp cover (total %)", breaks = c(0, 25, 50, 75, 100)) +
  guides(colour = guide_legend(order = 1),
         shape = guide_legend(order = 2),
         size = guide_legend(override.aes = list(shape = 16), order = 3)) +
  theme_grey() +
  theme(strip.background = element_rect(fill = NA),
        panel.border = element_rect(fill = NA, colour = "black", size = 1),
        axis.text = element_text(size = 12, colour = "black"),
        axis.ticks = element_line(colour = "black"))
ggsave("graph/MDS_plot.pdf", width = 10, height = 8)
ggsave("graph/MDS_plot.png", width = 10, height = 8)

