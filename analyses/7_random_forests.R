# analyses/7_random_forests.R
# The purpose of this script is to house the code used for the random forest analyses
# There is an analysis in script 5 and 6A as well, but my (RWS) thinking was to put everything
# in one place moving forward in the sake of tidiness.
# I (RWS) also found an additional method of performing a random forest in R that I
# will test out here as well



# Setup -------------------------------------------------------------------

# Load all data nad previous libraries
source("analyses/4_kelp_cover.R")

# Libraries for this script specifically
library(randomForest)
# library(forestRK)

# Site clims
load("data/study_site_clims.RData")

# BO data
load("data/study_site_BO.RData")

# Prep data.frames for combining
study_site_means$month <- "Annual"
study_site_BO$month <- "Annual"
study_site_BO$depth <- 0

# Combine the data
study_site_ALL <- rbind(study_site_means, study_site_clims) %>% 
  left_join(study_site_BO, by = c("site", "lon", "lat", "Campaign", "depth", "month")) %>% 
  dplyr::select(Campaign, site, month, lon, lat, nav_lon, nav_lat, lon_BO, lat_BO, x, y,
                bathy, depth, everything())
rm(study_site_means, study_site_clims, study_site_BO)


# Data --------------------------------------------------------------------

# Create wide data.frame with all of the abiotic variables matched to depth of site
# Also includes substrate percentage variables
kelp_all <- adf %>% 
  dplyr::select(Campaign, site, depth, Quadrat, Quadrat.size.m2, Bedrock..:sand, 
                kelp.cover, Laminariales, Agarum, Alaria) %>% 
  left_join(study_site_ALL, by = c("Campaign", "site")) %>% 
  # dplyr::select(Campaign:Alaria, depth.y, eken:icethic_cat, lon, lat) %>%
  mutate(eken = ifelse(depth.x == depth.y, eken, NA),  # A funny way of getting rid of non-target depth data
         soce = ifelse(depth.x == depth.y, soce, NA), 
         toce = ifelse(depth.x == depth.y, toce, NA),
         uoce = ifelse(depth.x == depth.y, uoce, NA),
         voce = ifelse(depth.x == depth.y, voce, NA),
         avt = ifelse(depth.x == depth.y, avt, NA),
         wo = ifelse(depth.x == depth.y, wo, NA)) %>% 
  dplyr::select(-c(nav_lon:depth.y)) %>% 
  dplyr::rename(depth = depth.x) %>% 
  gather(key = "model_var", value = "val", -c(Campaign:month, lon, lat)) %>%
  na.omit() %>% 
  pivot_wider(names_from = c(model_var, month), values_from = val)

# Filter down to only total kelp cover 
  #RWS: This is rather done below as it's own step
# kelp_var <- kelp_all %>% 
  # dplyr::select(-Laminariales, -Agarum, -Alaria) #%>% 
  # filter(depth %in% c(10, 15)) #%>% 
  # mutate(kelp.cover = round(kelp.cover, -1)) # Round kelp cover to the nearest 10% step


# Random Forest: Total kelp cover -----------------------------------------

# The data.frame that will be fed to the model
# The Quadrat information is fairly random, as it should be, and so isn't used in the model TRUE
# The substrate and depth values also score consistently in the bottom of importance so aren't used
data1 <- select(kelp_all, kelp.cover, everything(),
                -c(Campaign:Quadrat.size.m2), -Laminariales, -Agarum, -Alaria) 

test <- data1$vfxthin_Jul

# Which variables are highly correlated?
# RWS: It could be useful here to screen out any variables with correlations over +-0.5
tmp <- cor(data1) %>% 
  replace_na(0)
# tmp <- tmp[complete.cases(tmp)]
# tmp[upper.tri(tmp)] <- 0
# diag(tmp) <- 0
# Above two commands can be replaced with 
tmp[!lower.tri(tmp)] <- 0


data1 <- data1[,!apply(tmp,2,function(x) any(abs(x) >= 0.50))]


cor_df <- round(cor(data1), 2) %>% 
  data.frame() %>% 
  mutate(var1 = row.names(.)) %>% 
  pivot_longer(cols = -var1, names_to = "var2") %>% 
  # filter(value == 1 | value <= 0.5, value >= -0.5) %>%  # Find low correlation results
  # # filter(!(var1 %in% var2))
  # select(-value) %>% 
  # unite(col = "var3", var1, var2, sep = "_") %>% 
  # unique() %>% 
  # pivot_wider(names_from = var1, values_from = var2)
  filter(value != 1, value >= 0.5 | value <= -0.5) %>% # Find high correlation results
  mutate()
  # filter(var2 %in% var1)
  # select(var2) %>%
  # unique()

# Filter by columns that correlate poorly with something else
data1 <- data1[,colnames(data1) %in% cor_df$var2] # RWS: This does nothing at the moment

# Split data up for training and testing
set.seed(666)
train <- sample(nrow(data1), 0.7*nrow(data1), replace = FALSE)
train_set <- data1[train,]
valid_set <- data1[-train,]

# Test function to see what the best `mtry` value is
rf_mtry_test <- function(mtry_num){
  set.seed(666)
  test_rf <- randomForest(kelp.cover ~ ., data = train_set, mtry = mtry_num, ntree = 1000,
                          importance = TRUE)
  pred_test <- predict(test_rf, train_set)
  pred_accuracy <- data.frame(mtry = mtry_num,
                              acc = round(mean(abs(pred_test - train_set$kelp.cover))))
  return(pred_accuracy)
}
lapply(1:10, rf_mtry_test) # It looks like an mtry = 1 is best

# Random forest model based on all quadrat data
kelp_rf <- randomForest(kelp.cover ~ ., data = train_set, mtry = 1, ntree = 1000,
                        importance = TRUE, na.action = na.omit)
summary(kelp_rf)
print(kelp_rf)
# explains 42% of variance

round(importance(kelp_rf), 1)
varImpPlot(kelp_rf)
partialPlot(kelp_rf, train_set, lat) # RWS This doesn't run for me
partialPlot(kelp_rf, train_set, iceconc_cat) # RWS This doesn't run for me

# Predicting on training set
pred_train <- predict(kelp_rf, train_set)
table(pred_train, train_set$kelp.cover)  
mean(abs(pred_train - train_set$kelp.cover)) # within 21% accuracy...

# Predicting on Validation set
pred_valid <- predict(kelp_rf, valid_set)
table(pred_valid, valid_set$kelp.cover)  
mean(abs(pred_valid - valid_set$kelp.cover)) # within 18% accuracy...


# Random Forest: Laminariales ---------------------------------------------

data2 <- data1 <- select(kelp_all, Laminariales, everything(),
                         -c(Campaign:Quadrat.size.m2), -kelp.cover, -Agarum, -Alaria)
train <- sample(nrow(data2), 0.7*nrow(data2), replace = FALSE)
train_set <- data2[train,]
valid_set <- data2[-train,]
kelp_rf <- randomForest(Laminariales ~ ., data = train_set, mtry = 1, ntree = 1000,
                        importance = TRUE, na.action = na.omit)
summary(kelp_rf)
print(kelp_rf) # 25% variance explained
round(importance(kelp_rf), 1)
varImpPlot(kelp_rf)


# Predicting on training set
pred_train <- predict(kelp_rf, train_set)
table(pred_train, train_set$Laminariales)  
mean(abs(pred_train - train_set$Laminariales)) # within 12% accuracy...

# Predicting on Validation set
pred_valid <- predict(kelp_rf, valid_set)
table(pred_valid, valid_set$Laminariales)  
mean(abs(pred_valid - valid_set$Laminariales)) # within 13% accuracy...


# A different approach ----------------------------------------------------

# The results from the forestRK method above are worse than the classic random forest

# Going through everything I am starting to suspect that the difference in sampling 
# methodology between the campaigns may be an important factor
# So in this section I use the standard random forest method and see
# how well we can predict which campaign a data point comes from
# We'll do this both with and without lon lat

data3 <- select(kelp_all, Campaign, lon, eken_Annual, soce_Annual, toce_Annual) # With lon/lat
# data3 <- select(kelp_all, Campaign, eken_Annual, soce_Annual, toce_Annual) # Without lon/lat
data3$Campaign <- as.factor(data3$Campaign)
set.seed(666)
train <- sample(nrow(data1), 0.7*nrow(data3), replace = FALSE)
train_set <- data3[train,]
valid_set <- data3[-train,]

# Test function to see what the best `mtry` value is
rf_mtry_test <- function(mtry_num){
  set.seed(666)
  test_rf <- randomForest(Campaign ~ ., data = train_set, mtry = mtry_num, 
                          importance = TRUE, proximity = TRUE)
  pred_test <- predict(test_rf, train_set, type = "class")
  pred_accuracy <- data.frame(mtry = mtry_num,
                              acc = round(mean(pred_test == train_set$Campaign)*100))
  return(pred_accuracy)
}
lapply(1:9, rf_mtry_test) # It's always 100 percent accurate...

# Random forest model based on all quadrat data
kelp_rf <- randomForest(Campaign ~ ., data = train_set, importance = TRUE,
                        proximity = TRUE)
summary(kelp_rf)
round(importance(kelp_rf), 2)
varImpPlot(kelp_rf)

# Predicting on training set
pred_train <- predict(kelp_rf, train_set, type = "class") # type = 'response' also works but returns the same results
table(pred_train, train_set$Campaign)  
mean(pred_train == train_set$Campaign)*100 # 100% accuracy...

# Predicting on Validation set
pred_valid <- predict(kelp_rf, valid_set, type = "class")
table(pred_valid, valid_set$Campaign)  
mean(pred_valid == valid_set$Campaign)*100 # 100% accuracy...

