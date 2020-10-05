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
library(raster)
# library(sp)
# remotes::install_github("rspatial/raster") # Uncomment and run this line of code to install the development version of raster
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
v2 <- vifcor(na.omit(Arctic_BO)[, 3:34], th = 0.7)

# Exclude the collinear variables that were identified previously
Arctic_excl_pre <- Arctic_BO %>% 
  dplyr::select(lon, lat, v2@results$Variables) %>% ###JG: should we also remove SSTmax and PAR? (see lines 64-66)
  mutate(BO_parmean = replace_na(BO_parmean, 0),
         BO21_curvelltmin_bdmax = replace_na(BO21_curvelltmin_bdmax, 0)) %>% 
  arrange(lon, lat)

# One more layer of correlation screening
# Correlation plots
excl_VIF_df <- na.omit(as.data.frame(Arctic_excl_pre))
# dataVIF.cor <- cor(excl_VIF_df, method = c('pearson')) ##WHERE DO I SET 0.7? You don't. This just calculates the correlation values, it doesn't filter by them.
# corrplot(dataVIF.cor)
# heatmap(x = dataVIF.cor, symm = T)
Pearson_cor <- cor(excl_VIF_df) ## JG: NOT ALL CORRELATED VARIABLES HAVE BEEN EXCLUDED IN LATER STEP,
##WITH PEARSON CORRELATION, THERE ARE STILL 3 CORRELATIONS: 
##1) SSTmax/SSTmin, 2) PAR and SSTmax, 3)PAR and ice

# Remove SST long-term min and PAR
Arctic_excl <- Arctic_excl_pre %>% 
  dplyr::select(-BO2_templtmin_ss, -BO_parmean)
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
  dplyr::select(colnames(Arctic_excl)) %>% 
  filter(lon >= bbox_arctic[1], lon <= bbox_arctic[2],
         lat >= bbox_arctic[3], lat <= bbox_arctic[4])
Arctic_excl_2100_sub_stack <- stack(rasterFromXYZ(Arctic_excl_2100_sub, crs = "+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs"))
# plot(Arctic_excl_2100_sub_stack)
rm(Arctic_BO_2100, Arctic_excl_2100_sub); gc()

# Function for re-loading .RData files as necessary
# loadRData <- function(fileName){
#   load(fileName)
#   get(ls()[ls() != "fileName"])
# }

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
    # expl.var = Arctic_excl_sub_stack, # The MAXENT raster errors may be due to how large the raster files are
    #eval.resp.var, eval.expl.var, eval.resp.xy is for sp data to evaluate models. But we are doing 
    #the DataSplitTable ##JG= IS THAT ENOUGH?
    PA.strategy = 'random', # leave random (tried 'disk' but models were not robust enough)
    PA.dist.min = 0, # units = meters
    PA.dist.max = NULL, # units = meters
    PA.nb.rep = 5,#5, # several runs to prevent sampling bias since moderate number of pseudo-absence
    PA.nb.absences = 1000
    )
  
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
  # biomod_option@MAXENT.Phillips$memory_allocated = 4096 # Allocates 2048 MB/2 GB of memory to modeling. Be careful not to give too much.
  # biomod_option@MAXENT.Phillips$maximumiterations = 10000
  # biomod_option@MAXENT.Phillips$threshold = F
  # biomod_option@MAXENT.Phillips$hinge = F
  # biomod_option@MAXENT.Phillips$visible = F
  # biomod_option@MAXENT.Phillips$beta_lqp = .95
  
  ## Creating DataSplitTable
  # DataSplitTable <- BIOMOD_cv(biomod_data)
  # DataSplitTable <- BIOMOD_cv(biomod_data, k = 5, repetition = 2, do.full.models = F, #do.full.models=T models calibrated and evaluated with the whole dataset are done
  #                             stratified.cv = F, stratify = "both", balance = "pres")
  # DataSplitTable.y <- BIOMOD_cv(biomod_data, stratified.cv = T, stratify = "y", k = 2)
  # colnames(DataSplitTable.y)[1:2] <- c("RUN11","RUN12")
  # DataSplitTable <- cbind(DataSplitTable,DataSplitTable.y)
  # head(DataSplitTable)
  
  
  # 4: Model ----------------------------------------------------------------
  
  # Run the model
  biomod_model <- BIOMOD_Modeling(
    biomod_data,
    # models = c('RF', 'GLM'), #, 'MAXENT.Phillips', 'ANN', 'GAM'), # Testing without MAXENT as it doesn't run on my work server
    #models = c('RF', 'ANN'), #changed to ANN for testing (GLM Warning: 'glm.fit: fitted probabilities numerically 0 or 1 occurred') 
    models = c('MAXENT.Phillips', 'GLM', 'ANN', 'RF', 'GAM'),
    models.options = biomod_option,
    NbRunEval = 5, # 5 reps for final models
    DataSplit = 70, # Either chose a 70/30 split
    # DataSplitTable = DataSplitTable, # Or the cross-validation method. This takes much longer, but does run.
    ##JG= although it runs now here, it messes sections below and get the GLM warning
    VarImport = 5, # Number of permutations to estimate variable importance
    models.eval.meth = c('TSS', 'ROC', 'FAR', 'ACCURACY', 'SR'),
    # models.eval.meth = c('KAPPA', 'TSS', 'ROC', 'FAR', 'SR', 'ACCURACY', 'BIAS', 'POD', 'CSI', 'ETS'),
    rescal.all.models = TRUE,
    do.full.models = FALSE,
    modeling.id = sps_name)
  #biomod_model <- loadRData(paste0(sps_name,"/",sps_name,".",sps_name,".models.out"))
  
  # biomod_model # print summary
  # Model_scores <- get_evaluations(biomod_model) # get evaluation scores
  # dim(Model_scores)
  # dimnames(Model_scores)
  
  # Model evaluation by algorithm
  # models_scores_graph(biomod_model, by = "models", metrics = c('ROC','TSS'), 
  #                     xlim = c(0.5,1), ylim = c(0.5,1))
  
  # Model evaluation by cross-validation
  # models_scores_graph(biomod_model, by = "cv_run", metrics = c('ROC','TSS'), 
  #                     xlim = c(0.5,1), ylim = c(0.5,1))
  
  # Model evaluation by dataset
  # models_scores_graph(biomod_model, by = "data_set", metrics = c('ROC','TSS'), 
  #                     xlim = c(0.5,1), ylim = c(0.5,1))
  
  ## Calculate mean of variable importance by algorithm
      # JG: I have read that scores reported are raw in the table (to be easier to interpret, 
      # it should be normalized on our own - sum to 1 across algorithms)
  # (models_var_import <- get_variables_importance(biomod_model))
  # apply(models_var_import, c(1,2), mean, na.rm = T)
  # apply(apply(models_var_import, c(1,2), mean, na.rm = T), 1, mean) # Overall mean per variable
  
 # To visualize species' modeled response to the given variables
    ## JG: MAYBE THIS CAN BE MOVED TO ANOTHER PLACE LATER IF NEED BE AND/OR A WAY TO SUMMARIZE LINES 223-283 IF POSSIBLE
  # sp_name_Maxent <- BIOMOD_LoadModels(biomod_model, models = 'MAXENT.Phillips') 
  # sp_name_GLM <- BIOMOD_LoadModels(biomod_model, models = 'GLM')
  # sp_name_ANN <- BIOMOD_LoadModels(biomod_model, models = 'ANN')
  # sp_name_RF <- BIOMOD_LoadModels(biomod_model, models = 'RF')
  # sp_name_GAM <- BIOMOD_LoadModels(biomod_model, models = 'GAM')
  
  # Evaluate individual models
  # Maxent_eval_strip <- biomod2::response.plot2(
  #   models = sp_name_Maxent,
  #   Data = get_formal_data(biomod_model, 'expl.var'),
  #   show.variables = get_formal_data(biomod_model, 'expl.var.names'),
  #   do.bivariate = F,
  #   fixed.var.metric = 'mean',
  #   legend = F,
  #   display_title = F,
  #   data_species = get_formal_data(biomod_model, 'resp.var')
  # )
  
  # GLM_eval_strip <- biomod2::response.plot2(
  #   models = sp_name_GLM,
  #   Data = get_formal_data(biomod_model, 'expl.var'),
  #   show.variables = get_formal_data(biomod_model, 'expl.var.names'),
  #   do.bivariate = F,
  #   fixed.var.metric = 'mean',
  #   legend = F,
  #   display_title = F,
  #   data_species = get_formal_data(biomod_model, 'resp.var')
  # )
  
  # ANN_eval_strip <- biomod2::response.plot2(
  #   models = sp_name_ANN,
  #   Data = get_formal_data(biomod_model, 'expl.var'),
  #   show.variables = get_formal_data(biomod_model, 'expl.var.names'),
  #   do.bivariate = F,
  #   fixed.var.metric = 'mean',
  #   legend = F,
  #   display_title = F,
  #   data_species = get_formal_data(biomod_model, 'resp.var')
  # )
  
  # RF_eval_strip <- biomod2::response.plot2(
  #   models = sp_name_RF,
  #   Data = get_formal_data(biomod_model, 'expl.var'),
  #   show.variables = get_formal_data(biomod_model, 'expl.var.names'),
  #   do.bivariate = F,
  #   fixed.var.metric = 'mean',
  #   legend = F,
  #   display_title = F,
  #   data_species = get_formal_data(biomod_model, 'resp.var')
  # )
  
  # GAM_eval_strip <- biomod2::response.plot2(
  #   models = sp_name_GAM,
  #   Data = get_formal_data(biomod_model, 'expl.var'),
  #   show.variables = get_formal_data(biomod_model, 'expl.var.names'),
  #   do.bivariate = F,
  #   fixed.var.metric = 'mean',
  #   legend = F,
  #   display_title = F,
  #   data_species = get_formal_data(biomod_model, 'resp.var')
  # )
  
  # Build the ensemble models
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
    VarImport = 10 #10 for final models
    )
  
  # (models_scores_biomod_ensemble <- get_evaluations(biomod_ensemble))
  
  # Load if the model has already been run
  # biomod_ensemble <- loadRData(paste0(sps_name,"/",sps_name,".",sps_name,"ensemble.models.out"))
  
 # get_variables_importance(biomod_ensemble)
 # (models_var_import <- get_variables_importance(biomod_ensemble))
 # apply(models_var_import, c(1,2), mean, na.rm = T)
 # apply(apply(models_var_import, c(1,2), mean, na.rm = T), 1, mean) # Overall mean per variable
 ## JG: THESE VALUES SHOULD ALSO BE NORMALIZED? # RWS: I think they are?
 
 ### The other hint of how this could be than is similar to this, but still struggling to get the list
 # https://r-forge.r-project.org/forum/forum.php?thread_id=31877&forum_id=4342&group_id=302
 ## get BIOMOD_Modeling output object
 # bm.mod <- get(load(biomod_ensemble@models.out.obj@link))
 # 
 # ## load ensemble models
 # em.mods.names <- BIOMOD_LoadModels(biomod_ensemble)
 # em.mods.names
 # 
 # ## by default variable importance is not computed with ensemble models
 # get_variables_importance(biomod_ensemble)
 # 
 # ## but you can do it a posteriori
 # em.vi.list <- lapply(em.mods.names,
 #                      function(emn) {
 #                        variables_importance(get(emn), data = get_formal_data(bm.mod,'expl.var'))
 #                      })
 # names(em.vi.list) <- em.mods.names
 # str(em.vi.list )
 # 
 # myBiomodModelEval <- getModelsEvaluations(biomod_ensemble)
 # dimnames(myBiomodModelEval)
 # myBiomodModelEval["TSS"]
 # getModelsVarImport(myBiomodEM)
 

  # 5. Present projections --------------------------------------------------
  
  # Create projections
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
    do.stack = TRUE
  )
  
  # plot(biomod_projection)
  
  # Create ensemble projections  
  ##25-09 Error in dimnames(x) <- dn : length of 'dimnames' [2] not equal to array extent when usign DataSpiltTable
  ##When using 70/30 it runs but warnings about projection and WS84 ellipsoid # RWS: I no longer receive this warning.
  biomod_ensemble_projection <- BIOMOD_EnsembleForecasting(
    EM.output = biomod_ensemble,
    projection.output = biomod_projection, # Should be biomod_projection when problem fixed
    binary.meth = 'TSS',
    output.format = '.RData',
    do.stack = TRUE)
  
  # Visualise
  # plot(biomod_ensemble_projection)
  
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
    output.format = '.RData',
    compress = 'xz',
    build.clamping.mask = FALSE)
  
  # plot(biomod_projection_2050)
  
  # Create 2050 ensemble projections
  biomod_ensemble_projection_2050 <- BIOMOD_EnsembleForecasting(
    EM.output = biomod_ensemble,
    projection.output = biomod_projection_2050,
    binary.meth = 'TSS',
    output.format = '.RData',
    do.stack = TRUE)
  
  plot(biomod_ensemble_projection_2050)
  
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
    output.format = '.RData',
    compress = 'xz',
    build.clamping.mask = FALSE)
  
  # Create 2100 ensemble projections
  biomod_ensemble_projection_2100 <- BIOMOD_EnsembleForecasting(
    EM.output = biomod_ensemble,
    projection.output = biomod_projection_2100,
    binary.meth = 'TSS',
    output.format = '.RData',
    do.stack = TRUE)
  
  # plot(biomod_ensemble_projection_2100)
  
  # Clean out 2100
  rm(biomod_projection_2100); gc()
  
  ##Tests to plot ensemble models current and future conditions but I could not make it work.
  ##Delete this section if necessary
  # stk_biomod_ensemble_projection_2100 <- get_predictions(biomod_ensemble_projection_2100)
  # stk_biomod_ensemble_projection_2100 <- subset(stk_biomod_ensemble_projection_2100,
  #  grep("EMca/EMwmean", names(stk_biomod_ensemble_projection_2100)))
  # names(stk_biomod_ensemble_projection_2100) <- sapply(strsplit(names(stk_biomod_ensemble_projection_2100),
  #           "_"),
  #  getElement, 2)
  #  
  #  levelplot(biomod_ensemble_projection_2100,
  #           main = "Future 2100",
  #           col.regions = colorRampPalette(c("grey90", "yellow4", "green4"))(100))
  # 
  
  ##Species Range change
  # binary_2050 <- stack("Acla/proj_2050/proj_2050_Acla_TSSbin.grd")
  # binary_2100 <- raster::stack("Acla/proj_2100/proj_2100_Acla_TSSbin.grd")
  #       ##Did not use present, couldn't find the TSSbin.grd file
  #       ##There is another way to do this in Guisan book using .img files but was not able to do it
  # 
  #   RangeSize <- BIOMOD_RangeSize(
  #   CurrentPred = binary_2050,
  #   FutureProj = binary_2100
  # )
  # 
  # RangeSize$Compt.By.Models
  # plot(RangeSize$Diff.By.Pixel)
  ##Don't know what each layer is
  
  
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
  # NB: This currently won't run as we haven't modeled all of the species yet
# plyr::l_ply(sps_files, presence_absence_fig, .parallel = TRUE)

