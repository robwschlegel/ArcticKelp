# kelp_vs_clims.R
# The purpose of this script is to load the kelp and monthly clim data
# and then squoosh them together to see if anything interesting falls out


# Source scripts ----------------------------------------------------------

source("study_sites.R")
source("kelp_cover.R")


# Need to consider percent sand and rock when looking at patterns


# Load monthly clims ------------------------------------------------------

# The wide format data
load("data/study_site_clims.RData")

# The data in a long format
# load("data/study_site_clims_long.RData")


# Check for site mismatches -----------------------------------------------

# Sites from Kelp data not in study site names
unique(adf_summary$site)[!(unique(adf_summary$site) %in% unique(study_site_clims$site))]

# Study site names not in kelp data
unique(study_site_clims$site)[!(unique(study_site_clims$site) %in% unique(adf_summary$site))]

# Change site names to match kelp data
study_site_clims <- study_site_clims %>% 
  mutate(site = str_replace(site, "STN", "Station "),
         site = str_replace(site, "Durban Harbor", "Durban Harbour"),
         site = str_replace(site, "Bear I", "Bear Island"),
         site = str_replace(site, "S. of Qik", "S. of Qikitarjuak"))

# Re-run the above lines to check if everything has been fixed


# Join kelp and clims -----------------------------------------------------

kelp_clims <- left_join(adf_summary, study_site_clims, by = c("site", "Depth")) %>% 
  filter(depth > -1) %>%  # Removes sites with names that don't match up
  dplyr::select(Campaign:mean_cover, month, depth, eken:icethic_cat) %>% 
  mutate(eken = ifelse(Depth == depth, eken, NA),  # A funny way of getting rid of non-target depth data
         soce = ifelse(Depth == depth, soce, NA), 
         toce = ifelse(Depth == depth, toce, NA)) %>% 
  dplyr::select(-depth) %>% 
  gather(key = "model_var", value = "val", -c(Campaign:month)) %>% 
  na.omit()

kelp_clim_mean <- kelp_clims %>% 
  data.frame() %>% 
  dplyr::select(-month) %>% 
  group_by(Campaign, site, Depth, family, mean_cover, model_var) %>% 
  summarise(val = mean(val, na.rm = T))


# Multivariate analyses ---------------------------------------------------



# Visualise ---------------------------------------------------------------

# Going with scatterplots a.t.m.

# Kelp cover vs. abiotic variables per monthly climatology
scatter_clims <- kelp_clims %>% 
  # filter(family != "kelp.cover") %>% 
  ggplot(aes(x = val, y = mean_cover)) +
  geom_point(aes(group = site, colour = family)) +
  geom_smooth(method = "lm", se = F, aes(colour = family)) +
  facet_grid(month~model_var, scales = "free_x", switch = "both") +
  labs(x = "Model variable value", y = "Cover (%)", 
       colour = "Family", shape = "Depth (m)") +
  theme(legend.position = "top")
scatter_clims

# Kelp cover vs. overall mean for abiotic variables
scatter_clim_mean <- kelp_clim_mean %>% 
  filter(family == "kelp.cover") %>%
  ggplot(aes(x = val, y = mean_cover)) +
  geom_point(aes(group = site, colour = as.factor(Depth))) +
  geom_smooth(method = "lm", se = F, aes(colour = as.factor(Depth), linetype = Campaign)) +
  geom_label(aes(label = site, fill = as.factor(Depth)), size = 2) +
  facet_grid(~model_var, scales = "free_x", switch = "both") +
  labs(x = "Model variable value", y = "Cover (%)", 
       colour = "Depth (m)", fill = "Depth (m)", linetype = "Campaign") +
  theme(legend.position = "top")
scatter_clim_mean
ggsave("graph/scatter_clim_mean.pdf", scatter_clim_mean, height = 4, width = 16)
ggsave("graph/scatter_clim_mean.png", scatter_clim_mean, height = 4, width = 16)
