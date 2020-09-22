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
library(biomod2)
library(sp)
# remotes::install_github("rspatial/raster") # Uncomment and run this line of code to install the development version of raster
library(raster)
library(FNN)
library(doParallel)
library(usdm)
library(corrplot)

# The species occurrence data
sps_files <- dir("metadata", full.names = T, pattern = "rarefied")
sps_names <- str_remove(dir("metadata", full.names = F, pattern = "rarefied"), pattern = "_Arct_rarefied_points.csv")

# The present data
load("data/Arctic_BO.RData")

# Coordinates only
global_coords <- dplyr::select(Arctic_BO, lon, lat) %>% 
  mutate(env_index = 1:nrow(Arctic_BO))

# The usdm package can do stepwise elimination of highly inflating variables
# Calculate vif for the variables in Arctic_BO_stack
# vif(Arctic_BO[, 3:34])

#identify collinear variables that should be excluded (VIF > 10)
# v1 <- vifstep(Arctic_BO[, 3:34])

# Identify collinear variables that should be excluded (correlation > 0.7)
v2 <- vifcor(Arctic_BO[, 3:34], th = 0.7)

# Exclude the collinear variables that were identified previously
Arctic_excl <- Arctic_BO %>% 
  dplyr::select(lon, lat, v2@results$Variables) %>% 
  mutate(BO_parmean = replace_na(BO_parmean, 0),
         BO21_curvelltmin_bdmax = replace_na(BO21_curvelltmin_bdmax, 0)) %>% 
  arrange(lon, lat)
Arctic_excl_stack <- stack(rasterFromXYZ(Arctic_excl))
# plot(Arctic_excl_stack) # Visualise raster stack

# Correlation plots
# excl_VIF_df <- na.omit(as.data.frame(Arctic_excl))
# dataVIF.cor <- cor(excl_VIF_df, method = c('pearson')) ##WHERE DO I SET 0.7? You don't. This just calculates the correlation values, it doesn't filter by them.
# corrplot(dataVIF.cor)
# heatmap(x = dataVIF.cor, symm = T)

# Subset of present data used for projections
Arctic_excl_sub <- Arctic_excl %>% 
  filter(lon >= bbox_arctic[1], lon <= bbox_arctic[2],
         lat >= bbox_arctic[3], lat <= bbox_arctic[4])
Arctic_excl_sub_stack <- stack(rasterFromXYZ(Arctic_excl_sub))
# plot(Arctic_excl_sub_stack)
#rm(Arctic_BO, Arctic_BO_sub); gc()

# The 2050 data
load("data/Arctic_BO_2050.RData")
# Note that we do not need the full Arctic data for projections, just the study area
Arctic_excl_2050_sub <- Arctic_BO_2050 %>% 
  dplyr::select(colnames(Arctic_excl)) %>% 
  filter(lon >= bbox_arctic[1], lon <= bbox_arctic[2],
         lat >= bbox_arctic[3], lat <= bbox_arctic[4])
Arctic_excl_2050_sub_stack <- stack(rasterFromXYZ(Arctic_excl_2050_sub))
# plot(Arctic_excl_2050_sub_stack)
#rm(Arctic_BO_2050, Arctic_BO_2050_sub); gc()

# The 2100 data
load("data/Arctic_BO_2100.RData")
# Note that we do not need the full Arctic data for projections, just the study area
Arctic_excl_2100_sub <- Arctic_BO_2100 %>% 
  dplyr::select(colnames(Arctic_excl)) %>% 
  filter(lon >= bbox_arctic[1], lon <= bbox_arctic[2],
         lat >= bbox_arctic[3], lat <= bbox_arctic[4])
Arctic_excl_2100_sub_stack <- stack(rasterFromXYZ(Arctic_excl_2100_sub))
# plot(Arctic_excl_2100_sub_stack)
#rm(Arctic_BO_2100, Arctic_BO_2100_sub); gc()

# Function for re-loading .RData files as necessary
loadRData <- function(fileName){
  load(fileName)
  get(ls()[ls() != "fileName"])
}

# Choose a species for testing the code
# sps_choice <- sps_files[1]

# The full pipeline wrapped into a function
biomod_pipeline <- function(sps_choice){
  
  print(paste0("Began run on ",sps_choice))
  
  
  # 2: Load data ------------------------------------------------------------
  
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
  biomod_data <- BIOMOD_FormatingData(
    resp.var = rep(1, nrow(sps)),
    resp.xy = as.matrix(sps[,2:3]),
    resp.name = sps_name,
    expl.var = Arctic_excl_stack,
    # PA.strategy = 'disk', # leave random for now, but might be 'disk' the one to use
    # PA.dist.min = 10000, # units = meters
    # PA.dist.max = 50000, # units = meters
    PA.nb.rep = 1,#5, # several runs to prevent sampling bias since moderate number of pseudo-absence
    PA.nb.absences = 1000
    )
  
  # Change the NA absences values to 0, this makes them true absences... but then CV will run as expected
  biomod_data@data.species <- replace_na(biomod_data@data.species, 0)
  
  # biomod_data # object summary
  # plot(biomod_data) # plot selected pseudo-absences
  
  # The PA points can be visualised using the function in section 8
    # NB: This requires that the section 8 function be loaded into the environment first
  # presence_absence_fig(sps_choice)
  
  # To see where PA points are placed 
  # from http://rstudio-pubs-static.s3.amazonaws.com/416446_3ef37751ae1e4e569964dabc09a75b56.html
  # function to get PA dataset
  # get_PAtab <- function(bfd) {
  #   dplyr::bind_cols(
  #     x = bfd@coord[, 1],
  #     y = bfd@coord[, 2],
  #     status = bfd@data.species,
  #     bfd@PA
  #   )
  # }
  
  # function to get background mask
  # get_mask <- function(bfd){
  #   bfd@data.mask
  # }
  
  # get the coordinates of presences
  # (pres.xy <- get_PAtab(biomod_data) %>% 
  #     filter(status == 1) %>%
  #     dplyr::select(x, y))
  
  # get the coordiantes of pseudo - absences
  # all repetition of pseudo absences sampling merged 
  # (pa.all.xy <- get_PAtab(biomod_data) %>% 
  #     filter(is.na(status)) %>%
  #     dplyr::select(x, y) %>%
  #     distinct())
  
  # pseudo absences sampling for the first repetition only 
  # (pa.1.xy <- get_PAtab(biomod_data) %>% 
  #     filter(is.na(status) & PA1 == TRUE) %>%
  #     dplyr::select(x, y) %>%
  #     distinct())
  
  # plot the first PA selection and add the presences on top
  # plot(get_mask(biomod_data)[['PA1']])
  # points(pres.xy, pch = 11) 
  
  # biomod_data <- readRDS(paste0(sps_name,"/",sps_name,".base.Rds"))
  
  # Save the pre-model data for possible later use
  saveRDS(biomod_data, file = paste0(sps_name,"/",sps_name,".base.Rds"))
  
  # Model options
  biomod_option <- BIOMOD_ModelingOptions()
  
  # Setting up Maxent run
  # See here for more: https://gist.github.com/hannahlowens/974066848f8f85554ff7
  # biomod_option@MAXENT.Phillips$path_to_maxent.jar = paste(system.file(package = "dismo"), "/java", sep = '')
  biomod_option@MAXENT.Phillips$memory_allocated = 2048 # Allocates 2048 MB/2 GB of memory to modeling. Be careful not to give too much.
  # biomod_option@MAXENT.Phillips$maximumiterations = 10000
  # biomod_option@MAXENT.Phillips$threshold = F
  # biomod_option@MAXENT.Phillips$hinge = F
  # biomod_option@MAXENT.Phillips$visible = F
  # biomod_option@MAXENT.Phillips$beta_lqp = .95
  
  ## Creating DataSplitTable (check the code because it gives error of missing values)
  # DataSplitTable <- BIOMOD_cv(biomod_data)
  DataSplitTable <- BIOMOD_cv(biomod_data, k = 5, repetition = 2, do.full.models = F,
                              stratified.cv = F, stratify = "both", balance = "pres")
  DataSplitTable.y <- BIOMOD_cv(biomod_data, stratified.cv = T, stratify = "y", k = 2)
  colnames(DataSplitTable.y)[1:2] <- c("RUN11","RUN12")
  DataSplitTable <- cbind(DataSplitTable,DataSplitTable.y)
  # head(DataSplitTable)
  
  
  # 4: Model ----------------------------------------------------------------
  
  # Run the model
  biomod_model <- BIOMOD_Modeling(
    biomod_data,
    models = c('RF', 'GLM'), #, 'MAXENT.Phillips', 'ANN', 'GAM'), # Testing without MAXENT as it doesn't run on my work server
    # models = c('MAXENT.Phillips', 'GLM'), #, 'ANN', 'RF', 'GAM'),
    models.options = biomod_option,
    NbRunEval = 1, 
    DataSplit = 70,
    VarImport = 3, #number of permutations to estimate variable importance
    models.eval.meth = c('TSS', 'ROC', 'FAR', 'ACCURACY', 'SR'), # The fewer evaluation methods used the faster it runs
    # models.eval.meth = c('KAPPA', 'TSS', 'ROC', 'FAR', 'SR', 'ACCURACY', 'BIAS', 'POD', 'CSI', 'ETS'),
    rescal.all.models = FALSE,
    do.full.models = FALSE,
    modeling.id = sps_name)
   # biomod_model <- loadRData(paste0(sps_name,"/",sps_name,".",sps_name,".models.out"))
  
  biomod_model # print summary
  get_evaluations(biomod_model) # get evaluation scores
  
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
  
  # biomod_ensemble <- loadRData(paste0(sps_name,"/",sps_name,".",sps_name,"ensemble.models.out"))
  
  #variable importance (to complete with details but the code is as follows)
  variables_importance(
    model = biomod_ensemble, # ensemble models are also supported
    data = Arctic_excl_stack, # RWS: Not sure what is supposed to go here... 
    method = 'full_rand', 
    nb_rand = 3
  )
  
  # 5. Present projections --------------------------------------------------
  
  # Create projections
  biomod_projection <- BIOMOD_Projection(
    modeling.output = biomod_model,
    new.env = Arctic_excl_sub_stack,
    proj.name = 'present',
    binary.meth = 'TSS',
    compress = "xz",
    build.clamping.mask = FALSE)
  
  # Create ensemble projections
    # RWS: Can't run this until the above error has been addressed
  biomod_ensemble_projection <- BIOMOD_EnsembleForecasting(
    EM.output = biomod_ensemble,
    projection.output = biomod_projection)
  
  # Clean out some space
  rm(biomod_projection, biomod_ensemble_projection); gc()
  
  # Flush local tmp drive. Better not to do this if running on multiple cores
  # unlink(paste0(normalizePath(tempdir()), "/", dir(tempdir())), recursive = TRUE)
  # dir(tempdir())
  
  
  # 6: Future projections ---------------------------------------------------
  
  # Run 2050 projections
  biomod_projection_2050 <- BIOMOD_Projection(
    modeling.output = biomod_model,
    new.env = Arctic_excl_2050_sub_stack,
    proj.name = '2050',
    binary.meth = 'TSS',
    compress = 'xz',
    build.clamping.mask = FALSE)
  
  # Create 2050 ensemble projections
  biomod_ensemble_projection_2050 <- BIOMOD_EnsembleForecasting(
    EM.output = biomod_ensemble,
    projection.output = biomod_projection_2050)
  
  # Clean out 2050
  rm(biomod_projection_2050, biomod_ensemble_projection_2050); gc()
  
  # Flush local tmp drive. Better not to do this if running on mulitple cores
  # unlink(paste0(normalizePath(tempdir()), "/", dir(tempdir())), recursive = TRUE)
  
  # Run 2100 projections
  biomod_projection_2100 <- BIOMOD_Projection(
    modeling.output = biomod_model,
    new.env = Arctic_excl_2100_sub_stack,
    proj.name = '2100',
    binary.meth = 'TSS',
    compress = 'xz',
    build.clamping.mask = FALSE)
  
  # Create 2100 ensemble projections
  biomod_ensemble_projection_2100 <- BIOMOD_EnsembleForecasting(
    EM.output = biomod_ensemble,
    projection.output = biomod_projection_2100)
  
  # Clean out 2100
  rm(biomod_projection_2100); gc()
  
  # Flush local tmp drive. Better not to do this if running on mulitple cores
  # unlink(paste0(normalizePath(tempdir()), "/", dir(tempdir())), recursive = TRUE)
  
  # Unlink Temp folder
  unlink(paste0(sps_name,"/Temp"), recursive = TRUE)
}


# 7: Run the pipeline -----------------------------------------------------

# Detect available cores at set accordingly
registerDoParallel(cores = detectCores()-1)

# Run one
registerDoParallel(cores = 1)
biomod_pipeline(sps_files[5])

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
presence_absence_fig(sps_files[1])

# Run for all species
  # NB: This currently won't run as we haven't modeled all of the species yet
plyr::l_ply(sps_files, presence_absence_fig, .parallel = TRUE)

