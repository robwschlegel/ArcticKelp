# kelp_vs_clims.R
# The purpose of this script is to load the kelp and monthly clim data
# and then squoosh them together to see if anything interesting falls out


# Source scripts ----------------------------------------------------------

source("study_sites.R")
source("kelp_cover.R")


# Load monthly clims ------------------------------------------------------

load("data/study_site_clims.RData")
load("data/study_site_clims_long.RData")


# Check for site mismatches -----------------------------------------------

unique(adf_summary$site)[!(unique(adf_summary$site) %in% unique(study_site_clims$site))]


# Join kelp and clims -----------------------------------------------------

kelp_clims <- left_join(adf_summary, study_site_clims, by = c("site", "Depth"))
