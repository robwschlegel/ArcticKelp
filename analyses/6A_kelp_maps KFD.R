library(randomForest)

#Karen playing with things
  # fmmflx - Water flux due to freezing/melting
  # mldr10_1 - Mixed Layer Depth (dsigma = 0.01 wrt 10m)
  # qt - Net Downward Heat Flux
  # runoffs - Water Flux into Sea Water From Rivers
  # ssh - sea surface height
  # sss - Sea Surface Salinity
  # sst - Sea Surface Temperature
  # taum - wind stress module
  ## Ice:
  # iceconc_cat - Ice concentration for categories
  # icethic_cat - Ice thickness for categories
  ## Depth:
  # eken - kinetic energy
  # soce - Sea Water Salinity
  # toce - Sea Water Potential Temperature
colnames(Arctic_surface_mean)
#Agarum       #Alaria       #Laminariales #kelp.cover 

#interesting relationships:
#taum .

map_cover_abiotic(cover = "kelp.cover", abiotic = "emp_ice")
distribution_cover_abiotic(cover = "Laminariales", abiotic = "qns")
distribution_cover_abiotic(cover = "kelp.cover", abiotic = "bathy")
distribution_cover_abiotic(cover = "Agarum", abiotic = "emp_ice")
distribution_cover_abiotic(cover = "Laminariales", abiotic = "emp_ice")


abiotic_var <- Arctic_surface_mean
colnames(abiotic_var)[1:2]<-c("lon", "lat")
#abiotic_sub <- arrange(abiotic_sub, var) %>% 
 # mutate(dist_index = 1:n())

kelp_var <- adf_summary %>% 
  filter(family == 'kelp.cover') %>% filter(depth==10|depth==15) %>% 
  left_join(study_sites_index, by = c("Campaign", "site")) %>% 
  left_join(abiotic_var, by = c("nav_lon" = "lon", "nav_lat" = "lat"))
data1=kelp_var[,c(5,3,20)]
train <- sample(nrow(data1), 0.7*nrow(data1), replace = FALSE)
TrainSet <- data1[train,]

kelp.rf <- randomForest(kelp.cover ~ ., data=kelp_var, mtry=3,
                         importance=TRUE, na.action=na.omit)

set.seed(100)
train <- sample(nrow(data1), 0.7*nrow(data1), replace = FALSE)
TrainSet <- data1[train,]
ValidSet <- data1[-train,]
summary(TrainSet)
summary(ValidSet)