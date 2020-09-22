# tests/java_test.R
# The purpose of this script is to provide the minimum working example to run a MAXENT model via BIOMOD
# This runs as expected on a local machine, nut tikoraluk has a Java related error

# NB: MAXENT will run correctly if one first initialises -> source activate openjdk
# in a console and then opens R directly in that console and runs this script via source("tests/java_test.R")

# Libraries ---------------------------------------------------------------

# Tell tikoraluk where the personal R libraries are
.libPaths(c("~/R-packages", .libPaths()))

# Set Java home directory for tikoraluk
# Sys.setenv(PATH = "/software/miniconda3/envs/openjdk/bin:/software/miniconda3/condabin:/software/bin:/software/pgi/linux86-64/2018/bin:/software/pgi/linux86-64/18.4/mpi/openmpi/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/usr/games:/usr/local/games:/snap/bin:/software/miniconda3/bin")
Sys.getenv("PATH")

# Load needed libraries
library(tidyverse)
library(biomod2)
# install.packages("raster") # Checking to see if the CRAN version doesn't cause projection errors. It does. Better to still use development versinn.
# remotes::install_github("rspatial/raster")
library(raster)
library(FNN)
library(usdm)


# Data --------------------------------------------------------------------

# Set bounding box for study
bbox_arctic <- c(-95, -50, 50, 80)

# The present data
load("data/Arctic_BO.RData")

# Coordinates only
global_coords <- dplyr::select(Arctic_BO, lon, lat) %>% 
  mutate(env_index = 1:nrow(Arctic_BO))

# Identify collinear variables that should be excluded (correlation > 0.7)
v2 <- vifcor(Arctic_BO[, 3:34], th = 0.7)

# Exclude the collinear variables that were identified previously
Arctic_excl <- Arctic_BO %>% 
  dplyr::select(lon, lat, v2@results$Variables) %>% 
  mutate(BO_parmean = replace_na(BO_parmean, 0),
         BO21_curvelltmin_bdmax = replace_na(BO21_curvelltmin_bdmax, 0)) %>% 
  arrange(lon, lat)
Arctic_excl_stack <- stack(rasterFromXYZ(Arctic_excl))

# Subset of present data used for projections
Arctic_excl_sub <- Arctic_excl %>% 
  filter(lon >= bbox_arctic[1], lon <= bbox_arctic[2],
         lat >= bbox_arctic[3], lat <= bbox_arctic[4])
Arctic_excl_sub_stack <- stack(rasterFromXYZ(Arctic_excl_sub))

# Load the species
sps <- read_csv("metadata/Acla_Arct_rarefied_points.csv") %>% 
  mutate(env_index = as.vector(knnx.index(as.matrix(global_coords[,c("lon", "lat")]),
                                          as.matrix(.[,2:3]), k = 1))) %>%
  left_join(global_coords, by = "env_index") %>% 
  dplyr::select(Sp, lon.y, lat.y) %>%
  dplyr::rename(lon = lon.y, lat = lat.y)

# The species abbreviation
sps_name <- sps$Sp[1]

# Set temp folder save locations
# http://www.r-forge.r-project.org/forum/forum.php?thread_id=30946&forum_id=995&group_id=302
dir.create(file.path(sps_name), showWarnings = FALSE)
dir.create(file.path(sps_name,"/Temp"), showWarnings = FALSE)
rasterOptions(tmpdir = paste0(sps_name,"/Temp"))


# Model -------------------------------------------------------------------

# Prep data for modeling
biomod_data <- BIOMOD_FormatingData(
  resp.var = rep(1, nrow(sps)),
  resp.xy = as.matrix(sps[,2:3]),
  resp.name = sps_name,
  expl.var = Arctic_excl_stack,
  PA.nb.rep = 1,
  PA.nb.absences = 1000
)

# Change the NA absences values to 0, this makes them true absences... but then CV will run as expected
biomod_data@data.species <- replace_na(biomod_data@data.species, 0)

# Model options
biomod_option <- BIOMOD_ModelingOptions()
biomod_option@MAXENT.Phillips$memory_allocated = 2048

# Run model
biomod_model <- BIOMOD_Modeling(
  biomod_data,
  # models = c('RF', 'GLM'), # Testing without MAXENT as it doesn't run on my work server
  models = c('GLM', 'MAXENT.Phillips'), # Testing with MAXENT
  models.options = biomod_option,
  NbRunEval = 1, 
  DataSplit = 70,
  VarImport = 3,
  models.eval.meth = c('TSS', 'ROC', 'FAR', 'ACCURACY', 'SR'),
  rescal.all.models = FALSE,
  do.full.models = FALSE,
  SaveObj = FALSE,
  modeling.id = sps_name)

# print summary
biomod_model

# get evaluation scores
get_evaluations(biomod_model)

# Build the ensemble models
biomod_ensemble <- BIOMOD_EnsembleModeling(
  modeling.output = biomod_model,
  chosen.models = 'all',  # defines models kept (useful for removing non-preferred models)
  em.by = 'all',
  eval.metric = c('TSS'),
  eval.metric.quality.threshold = c(0.7),
  models.eval.meth = c('TSS', 'ROC', 'FAR', 'ACCURACY', 'SR'),
  prob.mean = TRUE, # Mean probabilities across predictions
  prob.cv = TRUE, # Coefficient of variation across predictions
  prob.ci = TRUE, # confidence interval around prob.mean
  prob.ci.alpha = 0.05,
)

# Tetst rasters
raster(x = "Acla/models/Acla/Acla_PA1_RUN1_GLM")

# Create projections
biomod_projection <- BIOMOD_Projection(
  modeling.output = biomod_model,
  new.env = Arctic_excl_sub_stack,
  proj.name = 'present',
  binary.meth = 'TSS',
  compress = "xz",
  build.clamping.mask = FALSE)

