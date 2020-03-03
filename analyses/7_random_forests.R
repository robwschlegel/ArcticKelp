# analyses/7_random_forests.R
# The purpose of this script is to house the code used for the random forest analyses
# There is an analysis in script 5 and 6A as well, but my (RWS) thinking was to put everything
# in one place moving forward in the sake of tidiness.


# An open water period would be a useful measurement to consider


# Setup -------------------------------------------------------------------

# Load all data nad previous libraries
source("analyses/4_kelp_cover.R")

# Libraries for this script specifically
library(tidymodels)  # Loads parsnip, rsample, recipes, yardstick
library(skimr)       # Quickly get a sense of data
library(randomForest)
library(knitr)
library(OneR) # For single rule machine learning
library(doParallel); doParallel::registerDoParallel(cores = 50) # This will be between 4 - 8 on a laptop

# Site clims
load("data/study_site_clims.RData")

# BO data
load("data/study_site_BO.RData")

# Prep data.frames for combining
# study_site_means$month <- "mean" # Use these two lines if using clim values
# study_site_BO$month <- "mean"
# study_site_BO$depth <- 0
# study_site_BO$Campaign <- NULL

# Combine the data
# study_site_ALL <- rbind(study_site_means, study_site_clims) %>% # Include monthly climatologies in modelling
# study_site_ALL <- study_site_means %>% # Use only overall means in modelling
#   # left_join(study_site_BO, by = c("site", "lon", "lat", "Campaign", "depth", "month")) %>% # For clim use 
#   left_join(study_site_BO, by = c("site", "lon", "lat", "Campaign", "depth")) %>% # For mean only use
#   # dplyr::select(Campaign, site, month, lon, lat, # Select month if using clim values 
#   dplyr::select(Campaign, site, lon, lat, 
#                 nav_lon, nav_lat, lon_BO, lat_BO, x, y,
#                 bathy, depth, everything()) %>% 
#   filter(depth == 0) # Filter out all NAPA model depth data, this keeps the BIO depth data
# rm(study_site_means, study_site_clims, study_site_BO)

# Remove scientific notation from data.frame displays in RStudio
options(scipen = 9999)

# The base map to use for everything else
Arctic_map <- ggplot() +
  borders(fill = "grey70", colour = "black") +
  coord_cartesian(expand = F,
                  xlim = c(bbox_arctic[1], bbox_arctic[2]),
                  ylim = c(bbox_arctic[3], bbox_arctic[4])) +
  labs(x = NULL, y = NULL)


# Data --------------------------------------------------------------------

# Create wide data.frame with all of the BO variables matched to sites
kelp_all <- adf %>% 
  # The substrate variables need to be removed as they can't be used in the final data
  # to predict kelp cover presence because we don't know what they are everywhere
  dplyr::select(Campaign, site, depth, -c(Bedrock..:sand), kelp.cover, Laminariales, Agarum, Alaria) %>% 
  # Create mean kelp cover values per site and depths
  # group_by(Campaign, site, depth) %>%
  # summarise_all(mean) %>%
  # ungroup() %>%
  # End mean site creation
  left_join(study_site_BO, by = c("Campaign", "site")) %>%
  mutate(kelp.cover = ifelse(kelp.cover > 100, 100, kelp.cover)) %>% # Correct values over 100
  dplyr::select(-lon_BO, -lat_BO) %>% 
  na.omit()


# Quick visuals -----------------------------------------------------------

# ggplot(kelp_all, aes(x = lon, y = kelp.cover)) +
#   geom_point()


# Which variables are highly correlated? ----------------------------------

# Identify variables that correlate with 5 or more other variables at abs(0.75) or more
cor_df <- round(cor(dplyr::select(kelp_all, -c(Campaign:site))), 2) %>% 
  data.frame() %>% 
  mutate(var1 = row.names(.)) %>% 
  pivot_longer(cols = -var1, names_to = "var2") %>% 
  na.omit() %>% 
  filter(value != 1, abs(value) >= 0.75) %>% # Find high correlation results
  group_by(var2) %>% 
  summarise(var2_count = n()) %>% 
  filter(var2_count >= 5) %>% 
  ungroup()


# Data prep function ------------------------------------------------------

# Convenience function for prepping dataframe for use in the random forest
# This removes all other kelp cover values
rf_data_prep <- function(kelp_choice, df = kelp_all, exclude_cor = T){
  # Trim down data.frame
  # The Quadrat information is fairly random, as it should be, and so isn't used in the model TRUE
  df_1 <- data.frame(dplyr::select(df, -c(Campaign:site))) 
  
  # Rmove highly correlated values
  if(exclude_cor) df_1 <- df_1[,!(colnames(df_1) %in% cor_df$var2)]
  
  # Create double of chosen kelp cover column
  names(df_1)[names(df_1) == kelp_choice] <- "chosen_kelp"
  
  # Determine which kelp cover is being used
  other_kelps <- c("kelp.cover", "Laminariales", "Agarum", "Alaria")
  other_kelps <- other_kelps[other_kelps != kelp_choice]
  
  # The data.frame that will be fed to the model
  df_2 <- dplyr::select(df_1, -chosen_kelp, everything(), chosen_kelp, -c(other_kelps)) #, depth) %>% 
  # mutate(depth = as.numeric(depth))
}


# Single rule model -------------------------------------------------------

# Prep the data for one rule modelling
kelp_cut <- kelp_all %>% 
  group_by(Campaign, site, depth) %>%
  summarise_all(mean) %>%
  ungroup() %>%
  # dplyr::select(-Campaign, -site) %>% 
  mutate(depth = as.numeric(depth), # Must be numeric for upcoming model training
         kelp.cover = ifelse(kelp.cover == 0, kelp.cover + 1, kelp.cover),  # Need to bump 0's for the binning process
         Laminariales = ifelse(Laminariales == 0, Laminariales + 1, Laminariales),
         Agarum = ifelse(Agarum == 0, Agarum + 1, Agarum),
         Alaria = ifelse(Alaria == 0, Alaria + 1, Alaria),
         kelp.cover = base::cut(kelp.cover, breaks = seq(0, 100, 20), ordered_result = F),
         Laminariales = base::cut(Laminariales, breaks = seq(0, 100, 20), ordered_result = F),
         Agarum = base::cut(Agarum, breaks = seq(0, 100, 20), ordered_result = F),
         Alaria = base::cut(Alaria, breaks = seq(0, 100, 20), ordered_result = F))

# Pull out training index
set.seed(666) # for reproducibility
one_random <- sample(1:nrow(kelp_all), 0.8*nrow(kelp_all))

# Function for OneR model run
OneR_model <- function(kelp_choice, df = kelp_cut){
  kelp_train <- rf_data_prep(kelp_choice = kelp_choice, df = df, exclude_cor = F) %>% 
    slice(one_random) %>%
    data.frame() %>% 
    optbin(., method = "infogain") # The choice of method does have an effect on the results
  kelp_test <- rf_data_prep(kelp_choice = kelp_choice, df = df, exclude_cor = F) %>% 
    slice(-one_random) %>% 
    data.frame()
  kelp_model <- OneR(kelp_train, verbose = TRUE)
  summary(kelp_model)
  plot(kelp_model)
  prediction <- predict(kelp_model, kelp_test)
  eval_model(prediction, kelp_test)
}

# Run models for each type of kelp cover
OneR_model("kelp.cover")
OneR_model("Laminariales")
OneR_model("Agarum")
OneR_model("Alaria")


# Test the best mtry ------------------------------------------------------

# The data.frame that will be fed to the model
# The Quadrat information is fairly random, as it should be, and so isn't used in the model TRUE
# The substrate and depth values also score consistently in the bottom of importance so aren't used
data1 <- select(kelp_all, kelp.cover, everything(),
                -c(Campaign:site), -Laminariales, -Agarum, -Alaria) 

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
lapply(1:10, rf_mtry_test) # It looks like an mtry >=6 is best


# Which variables are the most important? ---------------------------------

# This function runs many random forest models to determine which variables
# are consistantly the most important
top_variables <- function(lplyr_bit, kelp_choice){
  
  # Extract only the kelp cover of choice
  df_kelp_choice <- rf_data_prep(kelp_choice)
  
  # Split data up for training and testing
  train <- sample(1:nrow(df_kelp_choice), 0.7*nrow(df_kelp_choice), replace = FALSE)
  train_set <- df_kelp_choice[train,]
  valid_set <- df_kelp_choice[-train,]
  
  # Random forest model based on all quadrat data
  kelp_rf <- randomForest(chosen_kelp ~ ., data = train_set, mtry = 6, ntree = 1000,
                          importance = TRUE, na.action = na.omit)
  res <- data.frame(var = row.names(kelp_rf$importance), 
                    kelp_rf$importance, 
                    mean_MSE = mean(kelp_rf$mse, na.rm = T)) %>% 
    arrange(-X.IncMSE) %>% 
    # slice(1:30) %>% 
    mutate_if(is.numeric, round, 0)
}
# top_variables(kelp_choice = "kelp.cover")

# We then run this 1000 times to increase our certainty in the findings
top_variables_multi <- function(kelp_choice){
  multi_kelp <- plyr::ldply(.data = 1:1000, .fun = top_variables, .parallel = T, kelp_choice = kelp_choice)
  multi_kelp_importance <- multi_kelp %>% 
    group_by(var) %>% 
    summarise(mean_MSE = round(mean(mean_MSE, na.rm = T)),
              mean_IncMSE = mean(X.IncMSE),
              count = n()) %>% 
    ungroup() %>% 
    arrange(-mean_IncMSE, -count)
}

# Find the top variables for the different kelp covers
# system.time(top_var_kelpcover <- top_variables_multi("kelp.cover")) # ~19 seconds on 50 cores, ~543 on 3
# save(top_var_kelpcover, file = "data/top_var_kelpcover.RData")
# top_var_laminariales <- top_variables_multi("Laminariales")
# save(top_var_laminariales, file = "data/top_var_laminariales.RData")
# top_var_agarum <- top_variables_multi("Agarum")
# save(top_var_agarum, file = "data/top_var_agarum.RData")
# top_var_alaria <- top_variables_multi("Alaria")
# save(top_var_alaria, file = "data/top_var_alaria.RData")


# Random Forest function --------------------------------------------------

# Load top variable dataframes
load("data/top_var_kelpcover.RData")
load("data/top_var_laminariales.RData")
load("data/top_var_agarum.RData")
load("data/top_var_alaria.RData")

# kelp_choice <- "Laminariales"
# kelp_choice <- Laminariales
# column_choice <- top_var_laminariales
random_kelp_forest_check <- function(kelp_choice, column_choice){
  
  # Extract only the kelp cover of choice
  df_kelp_choice <- rf_data_prep(kelp_choice)
  df_kelp_choice$chosen_kelp <- base::cut(df_kelp_choice$chosen_kelp)
  
  # Chose only desired columns
  df_var_choice <- select(df_kelp_choice, chosen_kelp, as.character(column_choice$var))#[1:30])
  
  # Split data up for training and testing
  set.seed(666)
  train <- sample(nrow(df_var_choice), 0.7*nrow(df_var_choice), replace = FALSE)
  train_set <- df_var_choice[train,]
  # colnames(train_set)[1] <- "chosen_kelp"
  valid_set <- df_var_choice[-train,]
  # colnames(valid_set)[1] <- "chosen_kelp"
  
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
  # print(arrange(data.frame(var = row.names(kelp_rf$importance), 
  #                          kelp_rf$importance), -X.IncMSE)[1:30,])
}

# Check the random forests
random_kelp_forest_check("kelp.cover", top_var_kelpcover)
random_kelp_forest_check("Laminariales", top_var_laminariales)
random_kelp_forest_check("Agarum", top_var_agarum)
random_kelp_forest_check("Alaria", top_var_alaria)


# Many random forests -----------------------------------------------------

# This function is designed to output the model created and it's predictive accuracy
random_kelp_forest_test <- function(lplyr_bit, kelp_choice, column_choice){
  
  # Extract only the kelp cover of choice
  df_kelp_choice <- rf_data_prep(kelp_choice)
  
  # Chose only desired columns
  df_var_choice <- select(df_kelp_choice, chosen_kelp, as.character(column_choice$var))#[1:30])
  
  # Split data up for training and testing
  train <- sample(nrow(df_var_choice), 0.7*nrow(df_var_choice), replace = FALSE)
  train_set <- df_var_choice[train,]
  # colnames(train_set)[1] <- "chosen_kelp"
  valid_set <- df_var_choice[-train,]
  # colnames(valid_set)[1] <- "chosen_kelp"
  
  # Random forest model based on all quadrat data
  kelp_rf <- randomForest(chosen_kelp ~ ., data = train_set, mtry = 6, ntree = 1000,
                          importance = TRUE, na.action = na.omit)
  
  # Predicting on training and validation sets
  pred_train <- predict(kelp_rf, train_set)
  pred_valid <- predict(kelp_rf, valid_set)
  
  # Create data frame of accuracy results
  train_accuracy <- data.frame(portion = "train",
                               pred = round(pred_train, 2), 
                               original = train_set$chosen_kelp) %>% 
    mutate(accuracy = pred-original)
  validate_accuracy <- data.frame(portion = "validate",
                                  pred = round(pred_valid, 2), 
                                  original = valid_set$chosen_kelp) %>% 
    mutate(accuracy = pred-original)
  res_accuracy <- rbind(train_accuracy, validate_accuracy)
  
  # Package the model up with the accuracy results and exit
  res <- list(model = kelp_rf, accuracy = res_accuracy)
}


# Select best random forest -----------------------------------------------

# Now we run the test on each kelp cover 1000 times to see what the spread is
# in the accuracy of the random forests
# This is caused by different random splitting of test/validation sets
# as well as the many possible routes that the random forest may then take
random_kelp_forest_select <- function(kelp_choice, column_choice){
  # system.time(
    multi_test <- plyr::llply(.data = 1:1000, .fun = random_kelp_forest_test, .parallel = T, 
                              kelp_choice = kelp_choice, column_choice = column_choice)
    # ) # 56 seconds
  
  # Extract the model accuracies
  model_accuracy <- lapply(multi_test, function(x) x$accuracy) %>% 
    do.call(rbind.data.frame, .) %>% 
    mutate(model_id = rep(1:1000, each = nrow(kelp_all)))
  # RWS: A whole range of further analyses could be done with these values
  
  # Find which model had the best validation scores
  accuracy_check <- model_accuracy %>% 
    filter(portion == "validate") %>% 
    group_by(model_id) %>% 
    summarise(mean_acc = mean(abs(accuracy)),
              r_acc = round(cor(x = original, y = pred), 2))
  
  # Extract that model
  best_model <- arrange(accuracy_check, -r_acc, mean_acc) %>% 
    slice(1)
  choice_model <- multi_test[[best_model$model_id]]$model
  res <- list(choice_model = choice_model,
              model_accuracy = model_accuracy)
}

# doParallel::registerDoParallel(cores = 50)
# system.time(best_rf_kelpcover <- random_kelp_forest_select("kelp.cover", top_var_kelpcover)) # ~38 seconds with 50 cores
# save(best_rf_kelpcover, file = "data/best_rf_kelpcover.RData", compress = T)
# best_rf_laminariales <- random_kelp_forest_select("Laminariales", top_var_laminariales)
# save(best_rf_laminariales, file = "data/best_rf_laminariales.RData", compress = T)
# best_rf_agarum <- random_kelp_forest_select("Agarum", top_var_agarum)
# save(best_rf_agarum, file = "data/best_rf_agarum.RData", compress = T)
# best_rf_alaria <- random_kelp_forest_select("Alaria", top_var_alaria)
# save(best_rf_alaria, file = "data/best_rf_alaria.RData", compress = T)


# Analyse model accuracy --------------------------------------------------

# First load the best random forest models produced above
load("data/best_rf_kelpcover.RData")
load("data/best_rf_laminariales.RData")
load("data/best_rf_agarum.RData")
load("data/best_rf_alaria.RData")

# Find the distributions of accuracy from 0 - 100%
test <- best_rf_kelpcover$model_accuracy #%>%
  # filter(model_id == 1000)

# Quick visuals
ggplot(filter(test, portion == "validate"), aes(x = accuracy)) +
  geom_histogram()

ggplot(filter(test, portion == "validate"), aes(x = original, y = pred)) +
  geom_point() +
  geom_smooth(method = "lm")

ggplot(filter(test, portion == "validate"), aes(x = as.factor(original), y = pred)) +
  geom_boxplot()

# Function for creating figure showing confidence intervals of prediction accuracy
conf_plot <- function(df, plot_title){
  
  # 90 CI around predictions per step
  conf_acc <- df %>% 
    filter(portion == "validate") %>% 
    group_by(original) %>% 
    mutate(accuracy = round(accuracy)) %>% 
    summarise(q05 = quantile(accuracy, 0.05),
              q25 = quantile(accuracy, 0.25),
              q50 = median(accuracy),
              mean = mean(accuracy),
              q75 = quantile(accuracy, 0.75),
              q95 = quantile(accuracy, 0.95)) %>% 
    ungroup()
  
  conf_mean <- df %>% 
    filter(portion == "validate") %>% 
    group_by(model_id) %>% 
    mutate(accuracy = round(accuracy)) %>% 
    summarise(mean_acc = mean(abs(accuracy)),
              sd_acc = sd(abs(accuracy)),
              r_acc = cor(x = original, y = pred)) %>% 
    ungroup()
  
  conf_mean_label <- conf_mean %>% 
    summarise(mean_acc = round(mean(abs(mean_acc))),
              sd_acc = round(mean(abs(sd_acc))),
              r_acc = round(mean(r_acc), 2))
  
  conf_best_label <- conf_mean %>% 
    filter(mean_acc == min(mean_acc)) %>% 
    mutate(mean_acc = round(abs(mean_acc)),
           sd_acc = round(sd_acc),
           r_acc = round(r_acc, 2))
  
  conf_best <- df %>% 
    filter(model_id == conf_best_label$model_id,
           portion == "validate") %>% 
    group_by(original) %>% 
    summarise(mean_acc = mean(accuracy),
              sd_acc = sd(accuracy)) %>% 
    ungroup()
  
  # Visualise
  ggplot(conf_acc, aes(x = original, y = mean)) +
    geom_hline(yintercept = 0, size = 2, colour = "red") +
    geom_crossbar(aes(y = 0, ymin = q05, ymax = q95),
                  fatten = 0, fill = "grey70", colour = NA, width = 1) +
    geom_crossbar(aes(ymin = q25, ymax = q75),
                  fatten = 0, fill = "grey50", width = 1) +
    geom_crossbar(aes(ymin = q50, ymax = q50),
                  fatten = 0, fill = NA, colour = "black", width = 1) +
    geom_segment(data = conf_best, aes(xend = original, y = mean_acc, yend = 0), 
                 colour = "purple", size = 1.2, alpha = 0.8) +
    geom_point(data = conf_best, aes(y = mean_acc), colour = "purple", size = 3, alpha = 0.8) +
    geom_label(data = conf_mean_label, 
               aes(x = 75, y = 75, label = paste0("Mean accuracy: ",mean_acc,"±",sd_acc,"; r = ",r_acc))) +
    geom_label(data = conf_best_label, colour = "purple",
               aes(x = 75, y = 60, label = paste0("Best accuracy: ",mean_acc,"±",sd_acc,"; r = ",r_acc))) +
    scale_y_continuous(limits = c(-100, 100)) +
    labs(y = "Range in accuracy of predictions", x = "Original value (% cover)", title = plot_title)
}

# Create the plots
conf_plot(best_rf_kelpcover$model_accuracy, "Total cover confidence")
conf_plot(best_rf_laminariales$model_accuracy, "Laminariales cover confidence")
conf_plot(best_rf_agarum$model_accuracy, "Agarum cover confidence")
conf_plot(best_rf_alaria$model_accuracy, "Alaria cover confidence")


# Project kelp cover in the Arctic ----------------------------------------

# First load the best random forest models produced above
load("data/best_rf_kelpcover.RData")
load("data/best_rf_laminariales.RData")
load("data/best_rf_agarum.RData")
load("data/best_rf_alaria.RData")

# Load top variable choices
load("data/top_var_kelpcover.RData")
load("data/top_var_laminariales.RData")
load("data/top_var_agarum.RData")
load("data/top_var_alaria.RData")

# Load the Arctic data
load("data/Arctic_BO.RData")

# Prep the data for lazy joining
Arctic_mean_prep <- Arctic_mean %>% 
  # dplyr::select(-qla_oce, -qsb_oce) %>%  # These two columns have no values
  na.omit() %>% 
  dplyr::rename(lon = nav_lon, lat = nav_lat) %>% 
  mutate(lon = plyr::round_any(lon, 0.25),
         lat = plyr::round_any(lat, 0.25)) %>% 
  group_by(lon, lat) %>% 
  summarise_if(is.numeric, mean, na.rm = T) %>% 
  ungroup()

Arctic_BO_prep <- na.omit(Arctic_BO) %>% 
  mutate(lon = plyr::round_any(lon, 0.25),
         lat = plyr::round_any(lat, 0.25)) %>% 
  group_by(lon, lat) %>% 
  summarise_if(is.numeric, mean, na.rm = T) %>% 
  ungroup()

# Join the data for being fed to the model
Arctic_data <- left_join(Arctic_BO_prep, Arctic_mean_prep, by = c("lon", "lat")) %>%
  # na.omit() %>% 
  dplyr::select(-depth, -x, -y) %>% 
  dplyr::rename(depth = bathy) %>% 
  ungroup() %>% 
  tidyr::fill(lon:vtau_ice, .direction = "downup")

# Convenience function for final step before prediction
Arctic_cover_predict <- function(model_choice){
  pred_df <- data.frame(lon = Arctic_data$lon, lat = Arctic_data$lat,
                        depth = Arctic_data$depth,
                        pred_val = predict(model_choice, Arctic_data))
}

# Visualise a family of cover
cover_squiz <- function(df, legend_title, x_nudge){
  Arctic_map +
    geom_tile(data = filter(df, depth <= 50), 
              aes(x = lon, y = lat, fill = pred_val)) +
    scale_fill_viridis_c(legend_title) +
    theme(strip.background = element_rect(colour = "white", fill = "white"),
          legend.position = c(x_nudge, 0.96),
          legend.direction = "horizontal",
          legend.spacing.y = unit(0, "mm"))
}

# Predictions
pred_kelpcover <- Arctic_cover_predict(best_rf_kelpcover$choice_model)
cover_squiz(pred_kelpcover, "Total cover (%)", 0.785)
pred_laminariales <- Arctic_cover_predict(best_rf_laminariales$choice_model)
cover_squiz(pred_laminariales, "Laminariales cover (%)", 0.745)
pred_agarum <- Arctic_cover_predict(best_rf_agarum$choice_model)
cover_squiz(pred_agarum, "Agarum cover (%)", 0.77)
pred_alaria <- Arctic_cover_predict(best_rf_alaria$choice_model)
cover_squiz(pred_alaria, "Alaria cover (%)", 0.78)


# Substrate random forest -------------------------------------------------

# This would be useful to do, but would require substrate data for the entire Arctic


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


# A different prep method -------------------------------------------------

# Create a category version of the data for use with other models
kelp_all_factor <- kelp_all %>% 
  mutate(kelp.cover = as.factor(as.character(round(kelp.cover, -1))))

# Train/Test Split
train_test_split <-
  rsample::initial_split(
    data = kelp_all, # For continous response     
    # data = kelp_all_factor, # For categorical response
    prop = 0.70   
  ) 

train_test_split

train_tbl <- train_test_split %>% training() 
test_tbl  <- train_test_split %>% testing() 

# Prepare
recipe_simple <- function(dataset) {
  recipe(kelp.cover ~ ., data = dataset) %>%
    step_string2factor(all_nominal(), -all_outcomes()) %>%
    prep(data = dataset)
}

# NB: In order to avoid Data Leakage 
# (e.g: transferring information from the train set into the test set), 
# data should be “prepped” using the train_tbl only.
recipe_prepped <- recipe_simple(dataset = train_tbl)

# Bake the recipe
train_baked <- bake(recipe_prepped, new_data = train_tbl)
test_baked  <- bake(recipe_prepped, new_data = test_tbl)


# GLM ---------------------------------------------------------------------

# RWS-NB: The GLM method is white hot garbage
  # I've just left it in here for now to demonstrate how poorly it works

logistic_glm <- logistic_reg(mode = "classification") %>%
  # Change this to any number of things to change the model:
  # glm, glmnet, stan, spark, and keras
  set_engine("glm") %>%
  fit(kelp.cover ~ ., data = train_baked)

predictions_glm <- logistic_glm %>%
  predict(new_data = test_baked) %>%
  bind_cols(test_baked %>% select(kelp.cover))

predictions_glm %>% head() %>% kable()

predictions_glm %>%
  conf_mat(kelp.cover, .pred_class) %>%
  pluck(1) %>%
  as_tibble() %>%
  # Visualize with ggplot
  ggplot(aes(Prediction, Truth, alpha = n)) +
  geom_tile(show.legend = FALSE) +
  geom_text(aes(label = n), colour = "white", alpha = 1, size = 8)

# Accuracy
predictions_glm %>%
  metrics(kelp.cover, .pred_class) %>%
  select(-.estimator) %>%
  filter(.metric == "accuracy") %>%
  kable()

# Precision and Recall
# Precision shows how sensitive models are to False Positives 
# (i.e. predicting a customer is leaving when he-she is actually staying) 
# whereas Recall looks at how sensitive models are to False Negatives 
# (i.e. forecasting that a customer is staying whilst he-she is in fact leaving).
tibble(
  "precision" = 
    precision(predictions_glm, kelp.cover, .pred_class) %>%
    select(.estimate),
  "recall" = 
    recall(predictions_glm, kelp.cover, .pred_class) %>%
    select(.estimate)
) %>%
  unnest(cols = c(precision, recall)) %>%
  kable()

# F1 Score
# Another popular performance assessment metric is the F1 Score, 
# which is the harmonic average of the precision and recall. 
# An F1 score reaches its best value at 1 with perfect precision and recall.
predictions_glm %>%
  f_meas(kelp.cover, .pred_class) %>%
  select(-.estimator) %>%
  kable()


# A different random forest approach --------------------------------------

# NB: This method is not better than the original

rf_fun <- function(split, id, try, tree) {
  
  analysis_set <- split %>% analysis()
  analysis_prepped <- analysis_set %>% recipe_simple()
  analysis_baked <- analysis_prepped %>% bake(new_data = analysis_set)
  
  model_rf <-
    rand_forest(
      mode = "regression",
      mtry = try,
      trees = tree
    ) %>%
    set_engine("ranger",
               importance = "impurity"
    ) %>%
    fit(kelp.cover ~ ., data = analysis_baked)
  
  assessment_set     <- split %>% assessment()
  assessment_prepped <- assessment_set %>% recipe_simple()
  assessment_baked   <- assessment_prepped %>% bake(new_data = assessment_set)
  
  tibble(
    "id" = id,
    "truth" = assessment_baked$kelp.cover,
    "prediction" = model_rf %>%
      predict(new_data = assessment_baked) %>%
      unlist()
  )
}

cross_val_tbl <- vfold_cv(train_tbl, v = 10)
cross_val_tbl
cross_val_tbl %>% pluck("splits", 1)

# Accuracy
pred_rf <- map2_df(
  .x = cross_val_tbl$splits,
  .y = cross_val_tbl$id,
  ~ rf_fun(split = .x, id = .y, try = 3, tree = 200)) %>% 
  mutate(accuracy = abs(truth-prediction)) %>% 
  group_by(id) %>% 
  summarise(accuracy = mean(accuracy))
pred_rf


