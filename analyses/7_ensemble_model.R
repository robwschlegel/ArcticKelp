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

.libPaths(c("~/R-packages", .libPaths()))
library(tidyverse)
library(biomod2)
library(sp)
library(raster)
library(FNN)
library(doParallel) 

# The species occurrence data
sps_files <- dir("metadata", full.names = T, pattern = "rarefied")
sps_names <- str_remove(dir("metadata", full.names = F, pattern = "rarefied"), pattern = "_Arct_rarefied_points.csv")

# The present data
load("data/Arctic_BO.RData")
Arctic_BO_stack <- stack(rasterFromXYZ(Arctic_BO)) # This warning is caused by some layers not having the same number of pixels
# plot(Arctic_BO_stack)
# Arctic_BO_stack

# Subset of present data used for projections
Arctic_BO_sub <- Arctic_BO %>% 
  filter(lon >= bbox_arctic[1], lon <= bbox_arctic[2],
         lat >= bbox_arctic[3], lat <= bbox_arctic[4])
Arctic_BO_sub_stack <- stack(rasterFromXYZ(Arctic_BO_sub))
rm(Arctic_BO, Arctic_BO_sub); gc()

# The 2050 data
load("data/Arctic_BO_2050.RData")
# Note that we do not need the full Arctic data for projections, just the study area
Arctic_BO_2050_sub <- Arctic_BO_2050 %>% 
  filter(lon >= bbox_arctic[1], lon <= bbox_arctic[2],
         lat >= bbox_arctic[3], lat <= bbox_arctic[4])
Arctic_BO_2050_sub_stack <- stack(rasterFromXYZ(Arctic_BO_2050_sub))
rm(Arctic_BO_2050, Arctic_BO_2050_sub); gc()

# The 2100 data
load("data/Arctic_BO_2100.RData")
# Note that we do not need the full Arctic data for projections, just the study area
Arctic_BO_2100_sub <- Arctic_BO_2100 %>% 
  filter(lon >= bbox_arctic[1], lon <= bbox_arctic[2],
         lat >= bbox_arctic[3], lat <= bbox_arctic[4])
Arctic_BO_2100_sub_stack <- stack(rasterFromXYZ(Arctic_BO_2100_sub))
rm(Arctic_BO_2100, Arctic_BO_2100_sub); gc()

# Coordinates only
global_coords <- dplyr::select(Arctic_BO, lon, lat) %>% 
  mutate(env_index = 1:nrow(Arctic_BO))

# The best variables per species
  # Not currently known
# top_var <- read_csv("metadata/top_var.csv") %>% 
#   dplyr::select(Code:var6) %>% 
#   pivot_longer(cols = var1:var6) %>% 
#   dplyr::select(-name) %>% 
#   na.omit()

# Function for re-loading .RData files as necessary
loadRData <- function(fileName){
  load(fileName)
  get(ls()[ls() != "fileName"])
}

# Choose a species
# sps_choice <- sps_files[1]
# sps <- sps_names[1]

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
  
  # Filter out the top variables
    # We don't currently know what these will be before running the models
  # top_var_sub <- top_var %>% 
  #   filter(Code == sps_name) %>% 
  #   mutate(value = paste0(value,".asc"))
  
  # Load the top variables for the species
  # expl <- raster::stack(var_files[which(sapply(str_split(var_files, "/"), "[[", 3) %in% top_var_sub$value)])
  # expl_2050 <- raster::stack(var_2050_files[which(sapply(str_split(var_2050_files, "/"), "[[", 3) %in% top_var_sub$value)])
  # expl_2100 <- raster::stack(var_2100_files[which(sapply(str_split(var_2100_files, "/"), "[[", 3) %in% top_var_sub$value)])
  
  # Set temp folder save locations
  # http://www.r-forge.r-project.org/forum/forum.php?thread_id=30946&forum_id=995&group_id=302
  dir.create(file.path(sps_name), showWarnings = FALSE)
  dir.create(file.path(sps_name,"/Temp"), showWarnings = FALSE)
  rasterOptions(tmpdir = paste0(sps_name,"/Temp"))
  
  
  # 3: Prep data ------------------------------------------------------------
  
  # Prep data for modelling
  biomod_data <- BIOMOD_FormatingData(
    resp.var = rep(1, nrow(sps)),
    resp.xy = as.matrix(sps[,2:3]),
    resp.name = sps_name,
    expl.var = Arctic_BO_stack,
    PA.nb.rep = 1, #5, # It seems like 5 runs is unnecessary if 10,000 points are used
    PA.nb.absences = 10000)
  # biomod_data <- readRDS(paste0(sps_name,"/",sps_name,".base.Rds"))
  
  # Save the pre-model data for possible later use
  saveRDS(biomod_data, file = paste0(sps_name,"/",sps_name,".base.Rds"))
  
  # Model options
  biomod_option <- BIOMOD_ModelingOptions()
  
  # Setting up Maxent run
  # See here for more: https://gist.github.com/hannahlowens/974066848f8f85554ff7
  # biomod_option@MAXENT.Phillips$path_to_maxent.jar = paste(system.file(package="dismo"), "/java", sep='')
  biomod_option@MAXENT.Phillips$memory_allocated = 2048 # Allocates 2048 MB/2 GB of memory to modeling. Be careful not to give too much.
  biomod_option@MAXENT.Phillips$maximumiterations = 10000
  # biomod_option@MAXENT.Phillips$threshold = F
  # biomod_option@MAXENT.Phillips$hinge = F
  # biomod_option@MAXENT.Phillips$visible = F
  # biomod_option@MAXENT.Phillips$beta_lqp = .95
  
  
  # 4: Model ----------------------------------------------------------------
  
  # Run the model
  biomod_model <- BIOMOD_Modeling(
    biomod_data,
    models = c('MAXENT.Phillips', 'GLM', 'ANN', 'RF', 'GAM'),
    models.options = biomod_option,
    NbRunEval = 5,
    DataSplit = 70,
    VarImport = 0,
    # models.eval.meth = 'TSS', # The fewer evaluation methods used the faster it runs
    models.eval.meth = c('KAPPA', 'TSS', 'ROC', 'FAR', 'SR', 'ACCURACY', 'BIAS', 'POD', 'CSI', 'ETS'),
    rescal.all.models = TRUE,
    do.full.models = FALSE,
    modeling.id = sps_name)
  # biomod_model <- loadRData(paste0(sps_name,"/",sps_name,".",sps_name,".models.out"))
  
  # Build the ensemble models
     # RWS: This is currently giving an error but I haven't looked into why
  biomod_ensemble <- BIOMOD_EnsembleModeling(
    modeling.output = biomod_model,
    eval.metric = 'TSS',
    eval.metric.quality.threshold = 0.7,
    models.eval.meth = 'TSS'
  )
  # biomod_ensemble <- loadRData(paste0(sps_name,"/",sps_name,".",sps_name,"ensemble.models.out"))
  
  
  # 5. Present projections --------------------------------------------------
  
  # Create projections
  biomod_projection <- BIOMOD_Projection(
    modeling.output = biomod_model,
    new.env = Arctic_BO_sub_stack,
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
    new.env = Arctic_BO_2050_sub_stack,
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
    new.env = Arctic_BO_2100_sub_stack,
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
biomod_pipeline(sps_files[23])

# Run them all
registerDoParallel(cores = 5)
plyr::l_ply(sps_files, biomod_pipeline, .parallel = TRUE)

