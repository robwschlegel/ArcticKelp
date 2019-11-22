# analyses/7_random_forests.R
# The purpose of this script is to house the code used for the random forest analyses
# There is an analysis in script 5 and 6A as well, but my (RWS) thinking was to put everything
# in one place moving forward in the sake of tidiness.


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
  # Create mean kelp cover values per site
  # group_by(Campaign, site, depth) %>% 
  # summarise_all(mean) %>% 
  # ungroup() %>% 
  # End mean site creation
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


# Which variables are highly correlated? ----------------------------------

# Identify variables that correlate with 50 or more other variables at abs(0.75) or more
cor_df <- round(cor(select(kelp_all, -c(Campaign:Quadrat.size.m2))), 2) %>% 
  data.frame() %>% 
  mutate(var1 = row.names(.)) %>% 
  pivot_longer(cols = -var1, names_to = "var2") %>% 
  na.omit() %>% 
  filter(value != 1, abs(value) >= 0.75) %>% # Find high correlation results
  group_by(var2) %>% 
  summarise(var2_count = n()) %>% 
  filter(var2_count >= 50) %>% 
  ungroup()


# Random Forest function --------------------------------------------------

# kelp_choice <- "Laminariales"
random_kelp_forest <- function(kelp_choice){
  
  # Trim down data.frame
  # The Quadrat information is fairly random, as it should be, and so isn't used in the model TRUE
  df_1 <- select(kelp_all, -c(Campaign:Quadrat.size.m2)) %>% 
    data.frame()
  df_1 <- df_1[,!(colnames(df_1) %in% cor_df$var2)]
  
  # Create double of chosen kelp cover column
  # df_1$chosen_kelp <- df_1[ ,kelp_choice] 
  names(df_1)[names(df_1) == kelp_choice] <- "chosen_kelp"
  
  other_kelps <- c("kelp.cover", "Laminariales", "Agarum", "Alaria")
  other_kelps <- other_kelps[other_kelps != kelp_choice]
  
  # The data.frame that will be fed to the model
  # The substrate and depth values also score consistently in the bottom of importance so aren't used
  # RWS: I put them back in for the moment...
  df_2 <- select(df_1, chosen_kelp, everything(), -c(other_kelps)) 
  
  # Split data up for training and testing
  set.seed(666)
  train <- sample(nrow(df_2), 0.7*nrow(df_2), replace = FALSE)
  train_set <- df_2[train,]
  colnames(train_set)[1] <- "chosen_kelp"
  valid_set <- df_2[-train,]
  colnames(valid_set)[1] <- "chosen_kelp"
  
  # Random forest model based on all quadrat data
  kelp_rf <- randomForest(chosen_kelp ~ ., data = train_set, mtry = 6, ntree = 1000,
                          importance = TRUE, na.action = na.omit)
  print(kelp_rf) # percent of variance explained
  varImpPlot(kelp_rf)
  
  # Predicting on training set
  pred_train <- predict(kelp_rf, train_set)
  print(paste0("Average accuracy of prediction on test data: ", 
         round(mean(abs(pred_train - train_set$chosen_kelp)), 2),"%"))
  
  # Predicting on Validation set
  pred_valid <- predict(kelp_rf, valid_set)
  print(paste0("Average accuracy of prediction on validation data: ", 
         round(mean(abs(pred_valid - valid_set$chosen_kelp)), 2),"%"))
}


# Random Forest per family ------------------------------------------------

random_kelp_forest("kelp.cover")
random_kelp_forest("Laminariales")
random_kelp_forest("Agarum")
random_kelp_forest("Alaria")


# More thorough Random Forest ---------------------------------------------

# The data.frame that will be fed to the model
# The Quadrat information is fairly random, as it should be, and so isn't used in the model TRUE
# The substrate and depth values also score consistently in the bottom of importance so aren't used
data1 <- select(kelp_all, kelp.cover, everything(),
                -c(Campaign:Quadrat.size.m2), -Laminariales, -Agarum, -Alaria) 

# Filter by columns that correlate highly with many others
data1 <- data1[,!(colnames(data1) %in% cor_df$var2)]

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
lapply(1:10, rf_mtry_test) # It looks like an mtry >=6 1 is best

# Random forest model based on all quadrat data
kelp_rf <- randomForest(kelp.cover ~ ., data = train_set, mtry = 6, ntree = 1000,
                        importance = TRUE, na.action = na.omit)
summary(kelp_rf)
print(kelp_rf) # explains 48% of variance

round(importance(kelp_rf), 1)
varImpPlot(kelp_rf)
partialPlot(kelp_rf, train_set, lat) # RWS This doesn't run for me
partialPlot(kelp_rf, train_set, iceconc_cat) # RWS This doesn't run for me

# Predicting on training set
pred_train <- predict(kelp_rf, train_set)
table(pred_train, train_set$kelp.cover)  
mean(abs(pred_train - train_set$kelp.cover)) # within 18% accuracy, 10% with meaned sites

# Predicting on Validation set
pred_valid <- predict(kelp_rf, valid_set)
table(pred_valid, valid_set$kelp.cover)  
mean(abs(pred_valid - valid_set$kelp.cover)) # within 17% accuracy, 21% with meaned sites


# Random Forest: Campaign -------------------------------------------------

# The results from the forestRK method above are worse than the classic random forest

# Going through everything I am starting to suspect that the difference in sampling 
# methodology between the campaigns may be an important factor
# So in this section I use the standard random forest method and see
# how well we can predict which campaign a data point comes from
# We'll do this both with and without lon lat

dataC <- select(kelp_all, Campaign, lon, eken_Annual, soce_Annual, toce_Annual) # With lon/lat
# dataC <- select(kelp_all, Campaign, eken_Annual, soce_Annual, toce_Annual) # Without lon/lat
dataC$Campaign <- as.factor(dataC$Campaign)
set.seed(666)
train <- sample(nrow(data1), 0.7*nrow(dataC), replace = FALSE)
train_set <- dataC[train,]
valid_set <- dataC[-train,]

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
lapply(1:4, rf_mtry_test) # It's always 100 percent accurate...

# Random forest model based on all quadrat data
kelp_rf <- randomForest(Campaign ~ ., data = train_set, importance = TRUE,
                        proximity = TRUE, mtry = 2)
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

