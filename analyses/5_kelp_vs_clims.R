# analyses/5_kelp_vs_clims.R
# The purpose of this script is to load the kelp and monthly clim data
# and then squoosh them together to see if anything interesting falls out


# Source scripts ----------------------------------------------------------

source("analyses/4_kelp_cover.R")

library(vegan)
library(randomForest)

model_info <- read_csv("metadata/model_info.csv")


# Check for site mismatches -----------------------------------------------

# Sites from Kelp data not in study site names
unique(adf_summary$site)[!(unique(adf_summary$site) %in% unique(study_site_means$site))]

# Study site names not in kelp data
unique(study_site_means$site)[!(unique(study_site_means$site) %in% unique(adf_summary$site))]

# Change site names to match kelp data
# study_site_clims <- study_site_clims %>% 
#   mutate(site = str_replace(site, "Durban Habour", "Durban Harbour"),
#          site = str_replace(site, "S. of Qik", "S. of Qikitarjuak"),
#          site = str_replace(site, "Qik", "Qikitarjuak"))

# Re-run the above lines to check if everything has been fixed


# Join kelp and means -----------------------------------------------------

kelp_means <- left_join(adf_summary, study_site_means, by = c("Campaign", "site")) %>% 
  # filter(depth > -1) %>%  # Removes sites with names that don't match up
  dplyr::select(Campaign:mean_cover, depth.y, eken:icethic_cat) %>% 
  mutate(eken = ifelse(depth.x == depth.y, eken, NA),  # A funny way of getting rid of non-target depth data
         soce = ifelse(depth.x == depth.y, soce, NA), 
         toce = ifelse(depth.x == depth.y, toce, NA)) %>% 
  dplyr::select(-depth.y) %>% 
  dplyr::rename(depth = depth.x) %>% 
  gather(key = "model_var", value = "val", -c(Campaign:mean_cover)) %>% 
  na.omit()

# Kelp cover vs. overall mean for abiotic variables
scatter_means <- kelp_means %>% 
  filter(family == "kelp.cover") %>%
  ggplot(aes(x = val, y = mean_cover)) +
  geom_point(aes(group = site, colour = as.factor(depth))) +
  geom_smooth(method = "lm", se = F, aes(colour = as.factor(depth), linetype = Campaign)) +
  geom_label(aes(label = site, fill = as.factor(depth)), size = 2) +
  facet_grid(~model_var, scales = "free_x", switch = "both") +
  labs(x = "Model variable value", y = "Cover (%)", 
       colour = "Depth (m)", fill = "Depth (m)", linetype = "Campaign") +
  theme(legend.position = "top")
# scatter_means
# ggsave("graph/scatter_clim_mean.pdf", scatter_clim_mean, height = 4, width = 16)
# ggsave("graph/scatter_clim_mean.png", scatter_clim_mean, height = 4, width = 16)


# Multivariate analyses ---------------------------------------------------

# Need to consider percent sand and rock when looking at patterns
sand_rock <- adf %>% 
  select(Campaign, site, depth, sand, rock, Shell) %>% 
  mutate(Shell = replace_na(Shell, 0)) %>% 
  group_by(Campaign, site, depth) %>% 
  summarise_all(mean, na.rm = T) 

# Create data.frame that makes vegan happy
kelp_wide <- kelp_means %>% 
  spread(key = model_var, value = val) %>% 
  spread(key = family, value = mean_cover) %>% 
  ungroup() %>% 
  left_join(sand_rock, by = c("Campaign", "site", "depth"))

# The reduced version that doesn't know about depth etc.
kelp_wide_blind <- kelp_wide %>% 
  # select(eken:toce, Bedrock..:Shell)
  select(eken:toce)

# The "environmental" variables
kelp_wide_env <- kelp_wide %>% 
  select(Campaign, depth, kelp.cover:Shell)

# Run the MDS
kelp_MDS <- metaMDS(decostand(kelp_wide_blind, method = "standardize"),
                    distance = "euclidean", try = 100, autotransform = F)
# kelp_MDS <- metaMDS(kelp_wide_blind, distance = "euclidean", try = 100)
# kelp_MDS$species <- kelp_MDS$points

# Fit environmental variables
ord_fit <- envfit(kelp_MDS ~ depth + kelp.cover + rock + sand, data = kelp_wide_env)
ord_fit_df <- data.frame(ord_fit$vectors$arrows) %>% 
  mutate(var = row.names(.))

# Create a data.frame for ggplot
mds_df <- data.frame(site = kelp_wide$site, kelp_wide_env, kelp_MDS$points)

# Test correlations
cor(x = kelp_wide$kelp.cover, y = kelp_wide$sand)
cor(x = kelp_wide$kelp.cover, y = kelp_wide$rock)
cor(x = kelp_wide$kelp.cover, y = kelp_wide$Shell)

# The ordiplot
ggplot(data = mds_df, aes(x = MDS1, y = MDS2)) +
  # ggforce::geom_voronoi_tile(aes(group = Campaign)) +
  # stat_ellipse(aes(group = Campaign)) +
  geom_label(data = ord_fit_df, aes(label = var, x = NMDS1, y = NMDS2), size = 8) +
  # geom_label(aes(label = site))
  geom_segment(data = ord_fit_df, aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2),
               arrow = arrow(angle = 40, length = unit(0.2, "cm"), type = "open"), 
               alpha = 1, colour = "black", size = 0.5) +
  geom_point(aes(size = kelp.cover, shape = Campaign, colour = as.factor(depth))) +
  scale_colour_brewer(name = "Depth (m)", palette = "Dark2") +
  # scale_fill_brewer(name = "Campaign", palette = "Dark2") +
  scale_size_continuous(name = "Kelp cover (total %)", breaks = c(0, 25, 50, 75, 100)) +
  guides(colour = guide_legend(order = 1),
         shape = guide_legend(order = 2),
         size = guide_legend(override.aes = list(shape = 16), order = 3)) +
  # labs(size = "Duration") +
  theme_grey() +
  theme(strip.background = element_rect(fill = NA),
        panel.border = element_rect(fill = NA, colour = "black", size = 1),
        axis.text = element_text(size = 12, colour = "black"),
        axis.ticks = element_line(colour = "black"))
ggsave("graph/MDS_plot.pdf", width = 10, height = 8)
ggsave("graph/MDS_plot.png", width = 10, height = 8)


# Random forest -----------------------------------------------------------

# Test example
# set.seed(17)
# iris.urf <- randomForest(iris[, -5])
# MDSplot(iris.urf, iris$Species)

# Proposed formulae:
# Kelp cover ~ Ice + % rock + depth + Kinetic energy + Lat + Exposure + (1|Region/Site)
# Laminariales cover ~ Ice + % rock + depth + Kinetic energy + Lat + Exposure + (1|Region/Site)  

# Prep the data for a random forest
random_kelp_prep <- kelp_wide %>% 
  left_join(study_sites, by = c("Campaign", "site")) %>% 
  select(kelp.cover, iceconc_cat, rock, depth, eken, toce, lon, lat)

random_kelp <- randomForest(random_kelp_prep)
MDSplot(random_kelp, factor(kelp_wide$Campaign), k = 1)
ggsave("graph/random_forest_test_plot.png", height = 5, width = 5)
