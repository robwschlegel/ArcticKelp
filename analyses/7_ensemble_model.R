# analyses/7_ensemble_model.R
# This script contains the code used to run an ensemble model with the presence data
# The order of operations is:
# 1: Setup the environment
# 2: Load data
# 3: Prep data for modelling
# 4: Run model ensemble
# 5: Present projections
# 6: Future projections
# 7: Run the pipeline


# 1: Setup ----------------------------------------------------------------

# Load study sites and base packages
source("analyses/1_study_region_sites.R")

# Load additional packages
# devtools::install_github("biomodhub/biomod2", dependencies = TRUE) # Developmental version. Wasn't needed.
# remotes::install_github("rspatial/raster") # Uncomment and run this line of code to install the development version of raster
library(biomod2)
library(raster)
library(FNN)
library(doParallel)
library(usdm)
library(corrplot)

# Function for re-loading .RData files as necessary
loadRData <- function(fileName){
  load(fileName)
  get(ls()[ls() != "fileName"])
}

# The species occurrence data
sps_files <- dir("metadata", full.names = T, pattern = "rarefied")
sps_names <- str_remove(dir("metadata", full.names = F, pattern = "rarefied"), pattern = "_Arct_rarefied_points.csv")

# The coastal data
# load("data/Arctic_coast.RData")

# The coastal coordinates
# load("metadata/coastal_coords.RData")

# The present data
load("data/Arctic_BO.RData")

# Coordinates only
global_coords <- dplyr::select(Arctic_BO, lon, lat) %>% 
  mutate(env_index = 1:nrow(Arctic_BO))

# The usdm package can do stepwise elimination of highly inflating variables
# Calculate vif for the variables in Arctic_BO_stack
# v0 <- vif(Arctic_BO[, 3:34])

# Identify collinear variables that should be excluded (VIF > 10)
# v1 <- vifstep(Arctic_BO[, 3:34])

# Identify collinear variables that should be excluded (correlation > 0.7)
# Select only mean values and max SST
# v2 <- Arctic_BO %>% 
#   dplyr::select(BO21_templtmax_ss,
#                 BO21_salinitymean_ss, 
#                 BO21_icethickmean_ss, 
#                 BO21_curvelmean_bdmax, 
#                 BO_parmean, 
#                 BO21_dissoxmean_bdmax,
#                 BO21_ironmean_bdmax,
#                 BO21_nitratemean_bdmax,
#                 BO21_phosphatemean_bdmax) %>% 
#   na.omit() %>%
#   vifcor(th = 0.8)

# Save column names for use with Random Forest variable screening
# BO_vars <- c(v2@results$Variables)
# BO_vars <- BO_vars[-5] # Remove PAR because it correlates with SST and Ice
# BO_vars
# save(BO_vars, file = "metadata/BO_vars.RData")
load("metadata/BO_vars.RData")
BO_vars

# Exclude the collinear variables that were identified previously
Arctic_excl <- Arctic_BO %>% 
  dplyr::select(lon, lat, all_of(BO_vars)) %>%
  replace(is.na(.), 0) %>% 
  arrange(lon, lat)

# One more layer of correlation screening
# Correlation plots
# excl_VIF_df <- na.omit(as.data.frame(Arctic_excl[,-c(1,2)]))
# dataVIF.cor <- cor(excl_VIF_df, method = c('pearson'))
# corrplot(dataVIF.cor, type = "lower", method = "number")
# heatmap(x = dataVIF.cor, symm = T)
# Pearson_cor <- cor(excl_VIF_df)


# 2: Load data ------------------------------------------------------------

# Create raster stack
Arctic_excl_stack <- stack(rasterFromXYZ(Arctic_excl, crs = "+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs"))
# plot(Arctic_excl_stack) # Visualise raster stack

# Subset of present data used for projections
Arctic_excl_sub <- Arctic_excl %>% 
  filter(lon >= bbox_arctic[1], lon <= bbox_arctic[2],
         lat >= bbox_arctic[3], lat <= bbox_arctic[4])
Arctic_excl_sub_stack <- stack(rasterFromXYZ(Arctic_excl_sub, crs = "+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs"))
# plot(Arctic_excl_sub_stack)
rm(Arctic_BO, Arctic_excl_sub); gc()

# Very small area for testing
# Arctic_excl_baby <- Arctic_excl %>% 
#   filter(lon >= -90, lon <= -89,
#          lat >= 60, lat <= 61)
# Arctic_excl_baby_stack <- stack(rasterFromXYZ(Arctic_excl_baby, crs = "+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs"))

# The 2050 data
load("data/Arctic_BO_2050.RData")
# Note that we do not need the full Arctic data for projections, just the study area
Arctic_excl_2050_sub <- Arctic_BO_2050 %>% 
  right_join(dplyr::select(Arctic_excl, lon, lat), by = c("lon", "lat")) %>% 
  dplyr::select(colnames(Arctic_excl)) %>% 
  filter(lon >= bbox_arctic[1], lon <= bbox_arctic[2],
         lat >= bbox_arctic[3], lat <= bbox_arctic[4])
Arctic_excl_2050_sub_stack <- stack(rasterFromXYZ(Arctic_excl_2050_sub, crs = "+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs"))
# plot(Arctic_excl_2050_sub_stack)
rm(Arctic_BO_2050, Arctic_excl_2050_sub); gc()

# The 2100 data
load("data/Arctic_BO_2100.RData")
# Note that we do not need the full Arctic data for projections, just the study area
Arctic_excl_2100_sub <- Arctic_BO_2100 %>% 
  right_join(dplyr::select(Arctic_excl, lon, lat), by = c("lon", "lat")) %>% 
  dplyr::select(colnames(Arctic_excl)) %>% 
  filter(lon >= bbox_arctic[1], lon <= bbox_arctic[2],
         lat >= bbox_arctic[3], lat <= bbox_arctic[4])
Arctic_excl_2100_sub_stack <- stack(rasterFromXYZ(Arctic_excl_2100_sub, crs = "+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs"))
# plot(Arctic_excl_2100_sub_stack)
rm(Arctic_BO_2100, Arctic_excl_2100_sub); gc()

# Choose a species for testing the code
# sps_choice <- sps_files[1]

# The full pipeline wrapped into a function
biomod_pipeline <- function(sps_choice){
  
  print(paste0("Began run on ",sps_choice))
  
  # Load the species
  sps <- read_csv(sps_choice) %>% 
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
  
  
  # 3: Prep data ------------------------------------------------------------
  
  # Prep data for modeling
  # system.time(
  biomod_data <- BIOMOD_FormatingData(
    resp.var = rep(1, nrow(sps)),
    resp.xy = as.matrix(sps[,2:3]),
    resp.name = sps_name,
    expl.var = Arctic_excl_stack,
    #eval.resp.var, eval.expl.var, eval.resp.xy is for sp data to evaluate models. But we are doing the DataSplitTable
    PA.strategy = 'random', # leave random (tried 'disk' but models were not robust enough)
    # PA.dist.min = 0, # units = meters
    # PA.dist.max = NULL, # units = meters
    PA.nb.rep = 5, # several runs to prevent sampling bias since moderate number of pseudo-absence
    PA.nb.absences = 1000)
  # ) # 2 seconds
  
  # Change the NA absences values to 0, this makes them true absences... but then CV will run as expected
  # biomod_data@data.species <- replace_na(biomod_data@data.species, 0)
  
  # biomod_data # object summary
  # plot(biomod_data) # plot selected pseudo-absences
  
  # Save the pre-model data for possible later use
  saveRDS(biomod_data, file = paste0(sps_name,"/",sps_name,".base.Rds"))
  
  # The PA points can be visualised using the function in section 8
  # NB: This requires that the section 8 function be loaded into the environment first
  # presence_absence_fig(sps_choice)
  
  # biomod_data <- readRDS(paste0(sps_name,"/",sps_name,".base.Rds"))
  
  # Model options
  biomod_option <- BIOMOD_ModelingOptions()
  
  # Setting up Maxent run
  # See here for more: https://gist.github.com/hannahlowens/974066848f8f85554ff7
  # biomod_option@MAXENT.Phillips$path_to_maxent.jar = paste(system.file(package = "dismo"), "/java", sep = '')
  
  ## Creating DataSplitTable
  # DataSplitTable <- BIOMOD_cv(biomod_data)
  # DataSplitTable <- BIOMOD_cv(biomod_data, k = 5, repetition = 2, do.full.models = F, 
  #do.full.models=T models calibrated and evaluated with the whole dataset are done
  #                             stratified.cv = F, stratify = "both", balance = "pres")
  # DataSplitTable.y <- BIOMOD_cv(biomod_data, stratified.cv = T, stratify = "y", k = 2)
  # colnames(DataSplitTable.y)[1:2] <- c("RUN11","RUN12")
  # DataSplitTable <- cbind(DataSplitTable,DataSplitTable.y)
  # head(DataSplitTable)
  
  
  # 4: Model ----------------------------------------------------------------
  
  # Run the model
  # system.time(
  biomod_model <- BIOMOD_Modeling(
    biomod_data,
    # models = c('RF', 'GLM'), #, 'MAXENT.Phillips', 'ANN', 'GAM'), # Testing without MAXENT as it doesn't run on my work server
    #models = c('RF', 'ANN'), # changed to ANN for testing (GLM Warning: 'glm.fit: fitted probabilities numerically 0 or 1 occurred') 
    models = c('MAXENT.Phillips', 'GLM', 'ANN', 'RF', 'GAM'),
    models.options = biomod_option,
    NbRunEval = 5, # 5 reps for final models
    DataSplit = 70, # Either chose a 70/30 split
    # DataSplitTable = DataSplitTable, # Or the cross-validation method. This takes much longer, but does run.
    ## JG: although it runs now here, it messes sections below and get the GLM warning
    VarImport = 5, # Number of permutations to estimate variable importance
    models.eval.meth = c('TSS', 'ROC', 'FAR', 'ACCURACY', 'SR'),
    # models.eval.meth = c('KAPPA', 'TSS', 'ROC', 'FAR', 'SR', 'ACCURACY', 'BIAS', 'POD', 'CSI', 'ETS'),
    rescal.all.models = TRUE,
    do.full.models = FALSE,
    modeling.id = sps_name)
  # ) # 894 seconds for 1 rep
  # biomod_model <- loadRData(paste0(sps_name,"/",sps_name,".",sps_name,".models.out"))
  
  # Build the ensemble models
  # system.time(
  biomod_ensemble <- BIOMOD_EnsembleModeling(
    modeling.output = biomod_model,
    chosen.models = 'all',  # defines models kept (useful for removing non-preferred models)
    em.by = 'all',
    eval.metric = c('TSS'),
    eval.metric.quality.threshold = c(0.7), # Turn this off during testing if the ensemble won't run...
    models.eval.meth = c('TSS', 'ROC'), #, 'FAR', 'ACCURACY', 'SR'),
    prob.mean = TRUE, # Mean probabilities across predictions
    prob.cv = TRUE, # Coefficient of variation across predictions
    prob.ci = TRUE, # Confidence interval around prob.mean
    prob.ci.alpha = 0.05,
    VarImport = 10)
  # ) # 792 seconds for 1 rep
  
  # Load if the model has already been run
  # biomod_ensemble <- loadRData(paste0(sps_name,"/",sps_name,".",sps_name,"ensemble.models.out"))
  
  
  # 5. Present projections --------------------------------------------------
  
  # Create projections
  # system.time(
  biomod_projection <- BIOMOD_Projection(
    modeling.output = biomod_model,
    new.env = Arctic_excl_sub_stack,
    # new.env = Arctic_excl_baby_stack, # A very small region for faster testing
    # new.env = Arctic_excl_stack, # The full Arctic region. Not necessary.
    proj.name = 'present',
    selected.models = 'all',
    binary.meth = 'TSS',
    output.format = '.RData',
    compress = "xz",
    build.clamping.mask = FALSE,
    do.stack = TRUE)
  # ) # 76 seconds for 1 rep
  # plot(biomod_projection)
  
  # Create ensemble projections
  # system.time(
  biomod_ensemble_projection <- BIOMOD_EnsembleForecasting(
    EM.output = biomod_ensemble,
    projection.output = biomod_projection,
    binary.meth = 'TSS',
    output.format = '.RData',
    do.stack = TRUE)
  # ) # 7 seconds for 1 rep
  
  # Visualise
  # plot(biomod_ensemble_projection)
  
  # Clean out some space
  rm(biomod_projection, biomod_ensemble_projection); gc()
  
  # Flush local tmp drive. Better not to do this if running on multiple cores
  # unlink(paste0(normalizePath(tempdir()), "/", dir(tempdir())), recursive = TRUE)
  # dir(tempdir())
  
  
  # 6: Future projections ---------------------------------------------------
  
  # Run 2050 projections
  # system.time(
  biomod_projection_2050 <- BIOMOD_Projection(
    modeling.output = biomod_model,
    new.env = Arctic_excl_2050_sub_stack,
    proj.name = '2050',
    binary.meth = 'TSS',
    output.format = '.RData',
    compress = 'xz',
    build.clamping.mask = FALSE)
  # ) # 85 seconds for 1 rep
  # plot(biomod_projection_2050)
  
  # Create 2050 ensemble projections
  # system.time(
  biomod_ensemble_projection_2050 <- BIOMOD_EnsembleForecasting(
    EM.output = biomod_ensemble,
    projection.output = biomod_projection_2050,
    binary.meth = 'TSS',
    output.format = '.RData',
    do.stack = TRUE)
  # ) # 8 seconds for 1 rep
  # plot(biomod_ensemble_projection_2050)
  
  # Clean out 2050
  rm(biomod_projection_2050, biomod_ensemble_projection_2050); gc()
  
  # Flush local tmp drive. Better not to do this if running on multiple cores
  # unlink(paste0(normalizePath(tempdir()), "/", dir(tempdir())), recursive = TRUE)
  
  # Run 2100 projections
  # system.time(
  biomod_projection_2100 <- BIOMOD_Projection(
    modeling.output = biomod_model,
    new.env = Arctic_excl_2100_sub_stack,
    proj.name = '2100',
    binary.meth = 'TSS',
    output.format = '.RData',
    compress = 'xz',
    build.clamping.mask = FALSE)
  # ) # 84 seconds
  
  # Create 2100 ensemble projections
  # system.time(
  biomod_ensemble_projection_2100 <- BIOMOD_EnsembleForecasting(
    EM.output = biomod_ensemble,
    projection.output = biomod_projection_2100,
    binary.meth = 'TSS',
    output.format = '.RData',
    do.stack = TRUE)
  # ) # 8 seconds
  # plot(biomod_ensemble_projection_2100)
  
  # Clean out 2100
  rm(biomod_projection_2100, biomod_ensemble_projection_2100); gc()
  
  # Flush local tmp drive. Better not to do this if running on mulitple cores
  # unlink(paste0(normalizePath(tempdir()), "/", dir(tempdir())), recursive = TRUE)
  
  # Unlink Temp folder
  unlink(paste0(sps_name,"/Temp"), recursive = TRUE)
}


# 7: Run the pipeline -----------------------------------------------------

# Detect available cores at set accordingly
# registerDoParallel(cores = detectCores()-1)

# Run one
# registerDoParallel(cores = 1)
# biomod_pipeline(sps_files[1])

# Run them all
registerDoParallel(cores = 5)
plyr::l_ply(sps_files, biomod_pipeline, .parallel = TRUE)


# 8: Extract pseudo-absence points ----------------------------------------

# There is a need to see where the pseudo-absence points were selected for each species
# This function saves figures in the "graph/comparison_PA" folder

presence_absence_fig <- function(sps_choice){
  
  # Load the presence data
  sps_presence <- read_csv(sps_choice) %>% 
    mutate(env_index = as.vector(knnx.index(as.matrix(global_coords[,c("lon", "lat")]),
                                            as.matrix(.[,2:3]), k = 1))) %>%
    left_join(global_coords, by = "env_index") %>% 
    dplyr::select(Sp, lon.y, lat.y) %>%
    dplyr::rename(lon = lon.y, lat = lat.y) %>% 
    mutate(presence = 1)
  
  # Load the absence data
  sps <- str_remove(sps_choice, pattern = "_Arct_rarefied_points.csv")
  sps <- str_remove(sps, pattern = "metadata/")
  biomod_data <- readRDS(paste0(sps,"/",sps,".base.Rds"))
  sps_absence <- data.frame(presence = 0, biomod_data@coord)
  
  # Plot
  PA_fig <- ggplot(data = sps_absence, aes(x = lon, y = lat)) +
    borders(fill = "grey20", colour = "black") +
    geom_point(aes(colour = as.factor(presence)), size = 0.01) +
    geom_point(data = sps_presence, aes(colour = as.factor(presence)), size = 0.01) +
    labs(x = NULL, y = NULL, colour = "Presence") +
    scale_colour_manual(values = c("darkred", "forestgreen")) +
    coord_quickmap(expand = F, ylim = c(50, 85)) +
    theme(legend.position = "bottom")
  # PA_fig
  ggsave(plot = PA_fig, filename = paste0("graph/comparison_PA/",sps,"_PA.png"), width = 5, height = 2)
}

# Run for one species
# presence_absence_fig(sps_files[1])

# Run for all species
# plyr::l_ply(sps_files, presence_absence_fig, .parallel = TRUE)


# 9: Full model analysis --------------------------------------------------

# Choose a species for the following code
sps_choice <- sps_names[1]

# Load chosen biomod_model and print evaluation scores
biomod_model <- loadRData(paste0(sps_choice,"/",sps_choice,".",sps_choice,".models.out"))
(Model_scores <- get_evaluations(biomod_model))
apply(Model_scores, c(1,2,3), mean, na.rm = T)
# dim(Model_scores)
# dimnames(Model_scores)

# Model evaluation by algorithm
models_scores_graph(biomod_model, by = "models", metrics = c('ROC','TSS'),
                    xlim = c(0.5,1), ylim = c(0.5,1)) + ggtitle("Algorithm") +
  geom_hline(aes(yintercept = 0.7), colour = "red", size = 2)

# Model evaluation by cross-validation
models_scores_graph(biomod_model, by = "cv_run", metrics = c('ROC','TSS'),
                    xlim = c(0.5,1), ylim = c(0.5,1)) + ggtitle("Run") +
  geom_hline(aes(yintercept = 0.7), colour = "red", size = 2)

# Model evaluation by dataset
models_scores_graph(biomod_model, by = "data_set", metrics = c('ROC','TSS'),
                    xlim = c(0.5,1), ylim = c(0.5,1)) + ggtitle("PA") +
  geom_hline(aes(yintercept = 0.7), colour = "red", size = 2)

## Calculate mean of variable importance by algorithm
# JG: I have read that scores reported are raw in the table (to be easier to interpret, 
# it should be normalized on our own - sum to 1 across algorithms)
(models_var_import <- get_variables_importance(biomod_model))
apply(models_var_import, c(1,2), mean, na.rm = T)
apply(apply(models_var_import, c(1,2), mean, na.rm = T), 1, mean) # Overall mean per variable

# Visualize species' modeled response to the given variables
# These necessary files are not on GitHub as they are too large
sp_name_Maxent <- BIOMOD_LoadModels(biomod_model, models = 'MAXENT.Phillips')
sp_name_GLM <- BIOMOD_LoadModels(biomod_model, models = 'GLM')
sp_name_ANN <- BIOMOD_LoadModels(biomod_model, models = 'ANN')
sp_name_RF <- BIOMOD_LoadModels(biomod_model, models = 'RF')
sp_name_GAM <- BIOMOD_LoadModels(biomod_model, models = 'GAM')
sp_name_ALL <- BIOMOD_LoadModels(biomod_model, models = c('MAXENT.Phillips', 'GLM', 'ANN', 'RF', 'GAM'))

# Evaluate individual models
sp_curve_data <- response.plot2(
  models  = sp_name_ALL,
  Data = get_formal_data(biomod_model, 'expl.var'),
  show.variables = get_formal_data(biomod_model, 'expl.var.names'),
  data_species = get_formal_data(biomod_model, 'resp.var'),
  do.bivariate = FALSE,
  fixed.var.metric = 'median',
  # col = c("blue", "red"),
  legend = FALSE,
  plot = FALSE
)

# ggplot2 response curves
sp_curve_data %>%
  ## transform the pred.name to extract model, cv run and data info
  mutate(
    species = pred.name %>% strsplit('_') %>% sapply(function(x) x[1]),
    pa.dat = pred.name %>% strsplit('_') %>% sapply(function(x) x[2]),
    cv.rep = pred.name %>% strsplit('_') %>% sapply(function(x) x[3]),
    model = pred.name %>% strsplit('_') %>% sapply(function(x) x[4])
  ) %>%
  ggplot(
    aes(
      x = expl.val,
      y = pred.val,
      colour = model,
      group = pred.name
    )
  ) +
  geom_line(size = 1) +
  facet_wrap(~ expl.name, scales = 'free_x') + 
  labs(
    x = '',
    y = 'probability of occurence',
    colour = 'model type'
  ) + 
  scale_color_brewer(type = 'qual', palette = 4) +
  theme_minimal() +
  theme(
    legend.position = 'bottom'
  )


Maxent_eval_strip <- biomod2::response.plot2(
  models = sp_name_Maxent,
  Data = get_formal_data(biomod_model, 'expl.var'),
  show.variables = get_formal_data(biomod_model, 'expl.var.names'),
  do.bivariate = F,
  fixed.var.metric = 'mean',
  legend = F,
  display_title = F,
  data_species = get_formal_data(biomod_model, 'resp.var')
)
ggsave("", Maxent_eval_strip)
GLM_eval_strip <- biomod2::response.plot2(
  models = sp_name_GLM,
  Data = get_formal_data(biomod_model, 'expl.var'),
  show.variables = get_formal_data(biomod_model, 'expl.var.names'),
  do.bivariate = F,
  fixed.var.metric = 'mean',
  legend = F,
  display_title = F,
  data_species = get_formal_data(biomod_model, 'resp.var')
)
ANN_eval_strip <- biomod2::response.plot2(
  models = sp_name_ANN,
  Data = get_formal_data(biomod_model, 'expl.var'),
  show.variables = get_formal_data(biomod_model, 'expl.var.names'),
  do.bivariate = F,
  fixed.var.metric = 'mean',
  legend = F,
  display_title = F,
  data_species = get_formal_data(biomod_model, 'resp.var')
)
RF_eval_strip <- biomod2::response.plot2(
  models = sp_name_RF,
  Data = get_formal_data(biomod_model, 'expl.var'),
  show.variables = get_formal_data(biomod_model, 'expl.var.names'),
  do.bivariate = F,
  fixed.var.metric = 'mean',
  legend = F,
  display_title = F,
  data_species = get_formal_data(biomod_model, 'resp.var')
)
GAM_eval_strip <- biomod2::response.plot2(
  models = sp_name_GAM,
  Data = get_formal_data(biomod_model, 'expl.var'),
  show.variables = get_formal_data(biomod_model, 'expl.var.names'),
  do.bivariate = F,
  fixed.var.metric = 'mean',
  legend = F,
  display_title = F,
  data_species = get_formal_data(biomod_model, 'resp.var')
)

# Load and print ensemble model results
biomod_ensemble <- loadRData(paste0(sps_choice,"/",sps_choice,".",sps_choice,"ensemble.models.out"))
(models_scores_biomod_ensemble <- get_evaluations(biomod_ensemble))
(models_var_import <- get_variables_importance(biomod_ensemble))
apply(models_var_import, c(1,2), mean, na.rm = T)
apply(apply(models_var_import, c(1,2), mean, na.rm = T), 1, mean) # Overall mean per variable


# 10: Visualise ensemble models -------------------------------------------

# Load data used for maps etc.
load("data/Arctic_AM.RData")
colnames(Arctic_AM)[4] <- "depth"
Arctic_AM <- Arctic_AM %>%
  mutate(lon = round(lon, 4), lat = round(lat, 4))

# Choose a species
# sps_choice <- sps_names[1]

# Function that outputs BIOMOD projection comparison figures
plot_biomod <- function(sps_choice){
  # Load the species points
  sps_points <- read_csv(sps_files[str_which(sps_files,sps_choice)]) %>% 
    mutate(env_index = as.vector(knnx.index(as.matrix(global_coords[,c("lon", "lat")]),
                                            as.matrix(.[,2:3]), k = 1))) %>%
    left_join(global_coords, by = "env_index") %>% 
    dplyr::select(Sp, lon.y, lat.y) %>%
    dplyr::rename(lon = lon.y, lat = lat.y)
  
  # Load the ensemble projections
  biomod_project_present <- loadRData(paste0(sps_choice,"/proj_present/proj_present_",sps_choice,"_ensemble_TSSbin.RData"))
  biomod_project_2050 <- loadRData(paste0(sps_choice,"/proj_2050/proj_2050_",sps_choice,"_ensemble_TSSbin.RData"))
  biomod_project_2100 <- loadRData(paste0(sps_choice,"/proj_2100/proj_2100_",sps_choice,"_ensemble_TSSbin.RData"))
  
  # Convert to data.frames
  rast_df <- function(rast){
    df_out <- as.data.frame(rast[[1]], xy = T) %>% 
      `colnames<-`(c("lon", "lat", "presence")) %>% 
      mutate(lon = round(lon, 4), lat = round(lat, 4)) %>% 
      left_join(Arctic_AM, by = c("lon", "lat")) %>% 
      na.omit() 
  }
  df_project_present <- rast_df(biomod_project_present[[1]])
  df_project_2050 <- rast_df(biomod_project_2050[[1]])
  df_project_2100 <- rast_df(biomod_project_2100[[1]])
  
  # Visualise present data
  plot_present <- df_project_present %>% 
    filter(land_distance <= 50 | depth <= 100) %>% 
    filter(presence == TRUE) %>% 
    ggplot(aes(x = lon, y = lat)) +
    geom_tile(aes(fill = presence)) +
    borders(fill = "grey90", colour = "black") +
    geom_point(data = sps_points, colour = "yellow", size = 0.5) +
    scale_y_continuous(breaks = c(60, 70), labels = c("60°N", "70°N")) +
    scale_x_continuous(breaks = c(-80, -60), labels = c("80°W", "60°W")) +
    coord_quickmap(xlim = c(bbox_arctic[1], bbox_arctic[2]),
                   ylim = c(bbox_arctic[3], bbox_arctic[4]), expand = F) +
    scale_fill_manual(values = c("forestgreen")) +
    labs(x = NULL, y = NULL, title = paste0(sps_choice,": Present")) +
    theme_bw() +
    theme(legend.position = "bottom")
  
  # Function for visualising changes over time
  plot_diff <- function(df_future, year_label){
    plot_out <- left_join(df_project_present, df_future, 
                          by = c("lon", "lat", "land_distance", "depth")) %>% 
      mutate(change = presence.x - presence.y,
             change = ifelse(presence.x == F & presence.y == F, NA, change),
             change =  factor(change, 
                              levels = c("-1", "0", "1"),
                              labels = c("gain", "no change", "loss"))) %>% 
      na.omit() %>% 
      filter(land_distance <= 100 | depth <= 100) %>% 
      ggplot(aes(x = lon, y = lat)) +
      geom_tile(aes(fill = change)) +
      borders(fill = "grey90", colour = "black") +
      scale_y_continuous(breaks = c(60, 70), labels = c("60°N", "70°N")) +
      scale_x_continuous(breaks = c(-80, -60), labels = c("80°W", "60°W")) +
      coord_quickmap(xlim = c(bbox_arctic[1], bbox_arctic[2]),
                     ylim = c(bbox_arctic[3], bbox_arctic[4]), expand = F) +
      scale_fill_brewer(palette = "Set1", direction = -1) +
      labs(x = NULL, y = NULL, title = paste0(sps_choice,": Present - ",year_label)) +
      theme_bw() +
      theme(legend.position = "bottom",
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank())
  }
  
  # Visualise present - 2050
  plot_2050 <- plot_diff(df_project_2050, "2050")
  
  # Visualise present - 2100
  plot_2100 <- plot_diff(df_project_2100, "2100")
  
  # Combine and save
  plot_ALL <- cowplot::plot_grid(
    cowplot::plot_grid(
      plot_present + theme(legend.position = "none"),
      plot_2050 + theme(legend.position = "none"),
      plot_2100 + theme(legend.position = "none"),
      ncol = 3,
      align = "v"),
    cowplot::plot_grid(
      cowplot::get_legend(plot_present),
      ggplot() + theme_void(),
      cowplot::get_legend(plot_2100),
      ncol = 3, rel_widths = c(1, 0, 1.5)),
    nrow = 2, rel_heights = c(10,1)
  )
  ggsave(paste0("graph/biomod_diff_",sps_choice,".png"), plot_ALL, width = 8, height = 4.7)
}

# Create all visuals
registerDoParallel(cores = 5)
# plyr::l_ply(sps_names, plot_biomod, .parallel = T)


# 11: Save ensemble models as .grd files ----------------------------------

# Write function to load these .RData files and save them as .grd files

# Save as a raster file # NB: Still testing
# testers...
# sps_choice = "Acla"
# proj_choice = "2050"
RData_to_grd <- function(df){
  sps_choice <- df$sps_choice; proj_choice <- df$proj_choice
  proj_data <- loadRData(paste0(sps_choice,"/proj_",proj_choice,"/proj_",proj_choice,
                                "_",sps_choice,"_ensemble_TSSbin.RData"))
  proj_names <- sapply(strsplit(names(proj_data), "_"), "[[", 2)
  writeRaster(x = proj_data, bylayer = TRUE, suffix = proj_names, overwrite = TRUE,
              filename = paste0("data/ascii_results/proj_",proj_choice,"_",sps_choice,"_ensemble_TSSbin.asc"))
  # outfile <- writeRaster(proj_data, format = "ascii", overwrite = TRUE, 
  #                        options = c("INTERLEAVE=BAND", "COMPRESS=LZW"),
  #                        # suffix = c("EMmeanByTSS", "EMcvByTSS", "EMciInfByTSS", "EMciSupByTSS"),
  #                        filename = paste0(sps_choice,"/proj_",proj_choice,"/proj_",proj_choice,
  #                                          "_",sps_choice,"_ensemble_TSSbin.asc"))
}

# Run them all
quick_grid <- expand.grid(sps_choice = sps_names, 
                          proj_choice = c("present", "2050", "2100"), stringsAsFactors = F) %>% 
  mutate(plyr_id = 1:n()) %>% 
  data.frame()
plyr::d_ply(quick_grid, c("plyr_id"), RData_to_grd, .parallel = T)

# Check the output
test_rast <- raster("data/ascii_results/proj_2050_Acla_ensemble_TSSbin_EMmeanByTSS.asc")
plot(test_rast)

