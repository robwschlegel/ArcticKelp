# tests/MAXENT_test.R
# MAXENT refuses to allow model projections, but all other models work
# So this script runs an out-of-the-box MAXENT demo to see what is going wrong


# Libraries ---------------------------------------------------------------

library(tidyverse)
library(biomod2)
library(raster)
library(usdm)


# Data --------------------------------------------------------------------

# the name
# myRespName <-'Gulo'
myRespName <-'Acla' # Our data

# load a species raster
# we consider only the presences of Gulo Gulo species
myResp.ras <- raster(system.file("external/species/GuloGulo.grd", package = "biomod2"))

# A kelp sps
sps <- read_csv("metadata/Acla_Arct_rarefied_points.csv") %>% 
  dplyr::rename(x = lon, y = lat) %>% 
  mutate(z = 1,
         x = plyr::round_any(x, 0.25),
         y = plyr::round_any(y, 0.25)) %>% 
  dplyr::select(x, y, z)

# Convert CSV points to raster
x <- raster(xmn=min(sps$x), xmx=max(sps$x), ymn=min(sps$y), ymx=max(sps$y), res=3, crs="+proj=longlat +datum=WGS84")
sps.ras <- rasterize(sps[, c('x', 'y')], x, sps[, 'z'], fun=mean)
plot(sps.ras)

# extract the presences data
# the XY coordinates of the presence
# myRespXY <- xyFromCell(object = myResp.ras, cell = which(myResp.ras[]>0))
myRespXY <- xyFromCell(object = sps.ras, cell = which(sps.ras[]>0)) # Our data

# and the presence data
# myResp <- extract(x = myResp.ras, y = myRespXY)
myResp <- extract(x = sps.ras, y = myRespXY) # Our data

# load the environmental raster layers (could be .img, ArcGIS rasters or any supported format by the raster package)
myExpl <- stack(raster(system.file("external/bioclim/current/bio3.grd", package = "biomod2")),
                raster(system.file("external/bioclim/current/bio4.grd", package = "biomod2")),
                raster(system.file("external/bioclim/current/bio7.grd", package = "biomod2")),
                raster(system.file("external/bioclim/current/bio11.grd", package = "biomod2")),
                raster(system.file("external/bioclim/current/bio12.grd", package = "biomod2")))

# Our explanatory data
# The present data
load("data/Arctic_BO.RData")

# Identify collinear variables that should be excluded (correlation > 0.7)
v2 <- vifcor(Arctic_BO[, 3:34], th = 0.7)

# Exclude the collinear variables that were identified previously
Arctic_excl <- Arctic_BO %>% 
  dplyr::select(lon, lat, v2@results$Variables) %>% ###JG: should we also remove SSTmax and PAR? (see lines 64-66)
  mutate(BO_parmean = replace_na(BO_parmean, 0),
         BO21_curvelltmin_bdmax = replace_na(BO21_curvelltmin_bdmax, 0)) %>% 
  arrange(lon, lat) %>% 
  mutate(lon = plyr::round_any(lon, 3),
         lat = plyr::round_any(lat, 3)) %>% 
  group_by(lon, lat) %>% 
  summarise_all(mean, 3, na.rm = T, .groups = "drop")
Arctic_excl_stack <- stack(rasterFromXYZ(Arctic_excl, crs = "+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs"))


# Model -------------------------------------------------------------------

# Format data
myBiomodData <- BIOMOD_FormatingData(resp.var = myResp,
                                     # expl.var = myExpl,
                                     expl.var = Arctic_excl_stack, # Our data
                                     resp.xy = myRespXY,
                                     resp.name = myRespName,
                                     PA.nb.rep = 2,
                                     PA.nb.absences = 200,
                                     PA.strategy ='random')

# 2Defining Models Options using default options.
myBiomodOption <- BIOMOD_ModelingOptions()

# Computing the models
myBiomodModelOut <- BIOMOD_Modeling(myBiomodData,
                                    models = c('SRE','RF','MAXENT.Phillips'),
                                    models.options = myBiomodOption,
                                    NbRunEval = 1,
                                    DataSplit = 80,
                                    Yweights = NULL,
                                    VarImport = 3,
                                    models.eval.meth = c('TSS','ROC'),
                                    SaveObj = TRUE,
                                    rescal.all.models = TRUE)

# let's have a look at different models scores
getModelsEvaluations(myBiomodModelOut)

# Project our models over studied area
myBiomomodProj <- BIOMOD_Projection(modeling.output = myBiomodModelOut,
                                    # new.env = myExpl,
                                    new.env = Arctic_excl_stack, # Our data
                                    proj.name = 'current',
                                    # xy.new.env = myRespXY, # I don't think this is necessary because the data are a raster
                                    selected.models = 'all',
                                    binary.meth = 'ROC',
                                    filtered.meth = 'TSS',
                                    compress = 'xz',
                                    clamping.mask = T,
                                    do.stack = T)

# make some plots sub-selected by str.grep argument
plot(myBiomomodProj, str.grep = 'MAXENT')

