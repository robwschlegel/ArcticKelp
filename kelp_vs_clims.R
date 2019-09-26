# kelp_vs_clims.R
# The purpose of this script is to load the kelp and monthly clim data
# and then squoosh them together to see if anything interesting falls out


# Source scripts ----------------------------------------------------------

source("study_sites.R")
source("kelp_cover.R")


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
  dplyr::select(site:mean_cover, month, depth, eken:icethic_cat) %>% 
  mutate(eken = ifelse(Depth == depth, eken, NA),  # A funny way of getting rid of non-target depth data
         soce = ifelse(Depth == depth, soce, NA), 
         toce = ifelse(Depth == depth, toce, NA)) %>% 
  dplyr::select(-depth) %>% 
  gather(key = "model_var", value = "val", -c(site:month)) %>% 
  na.omit()


# Visualise ---------------------------------------------------------------

# Going with scatterplots a.t.m.

scatter_all <- kelp_clims %>% 
  filter(family != "kelp.cover") %>% 
  ggplot(aes(x = val, y = mean_cover)) +
  geom_point(aes(group = site, colour = month, shape = as.factor(Depth))) +
  facet_grid(family~model_var, scales = "free", switch = "both") +
  labs(x = "Model variable value", y = "Cover (%)", 
       colour = "Monthly clim", shape = "Depth (m)")
scatter_all

