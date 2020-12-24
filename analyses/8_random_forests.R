# analyses/8_random_forests.R
# The purpose of this script is to house the code used for the random forest analyses


# Setup -------------------------------------------------------------------

# Load all data nad previous libraries
source("analyses/4_kelp_cover.R")

# Libraries for this script specifically
library(randomForest)
library(OneR) # For single rule machine learning
# library(caret) # For cross validation option
library(doParallel); doParallel::registerDoParallel(cores = 50) # This will be between 4 - 8 on a laptop

# Load BO layer names used for the ensemble models
  # NB: Somehow the ensembles were given a different salinity variable
  # This will need to be corrected later
load("metadata/BO_vars.RData")
BO_vars[4] <- "BO2_salinityltmax_ss" 

# Environmental data per site
load("data/study_site_env.RData")
study_site_env <- study_site_env %>% 
  dplyr::select(site:land_distance, all_of(BO_vars))

# Load Arctic data for testing variable correlations and for making model projections
load("data/Arctic_BO.RData")
Arctic_BO <- Arctic_BO %>% 
  mutate(lon = round(lon, 4), lat = round(lat, 4)) %>% 
  dplyr::select(lon, lat, all_of(BO_vars))
load("data/Arctic_AM.RData")
Arctic_AM <- Arctic_AM %>% 
  mutate(lon = round(lon, 4), lat = round(lat, 4))
Arctic_env <- right_join(Arctic_BO, Arctic_AM, by = c("lon", "lat"))
rm(Arctic_BO); gc()

# Load future layers
load("data/Arctic_BO_2050.RData")
Arctic_env_2050 <- Arctic_BO_2050 %>% 
  mutate(lon = round(lon, 4), lat = round(lat, 4)) %>% 
  dplyr::select(lon, lat, all_of(BO_vars)) %>% 
  right_join(Arctic_AM, by = c("lon", "lat")) 
load("data/Arctic_BO_2100.RData")
Arctic_env_2100 <- Arctic_BO_2100 %>% 
  mutate(lon = round(lon, 4), lat = round(lat, 4)) %>% 
  dplyr::select(lon, lat, all_of(BO_vars)) %>% 
  right_join(Arctic_AM, by = c("lon", "lat")) 
rm(Arctic_AM, Arctic_BO_2050, Arctic_BO_2100); gc()

# Load the BO correlation matrix
# load("data/BO_cor_matrix.RData")

# Remove scientific notation from data.frame displays in RStudio
options(scipen = 9999)

# The base map to use for everything else
Arctic_map <- ggplot() +
  borders(fill = "grey70", colour = "black") +
  coord_quickmap(expand = F,
                  xlim = c(bbox_arctic[1], bbox_arctic[2]),
                  ylim = c(bbox_arctic[3], bbox_arctic[4])) +
  labs(x = NULL, y = NULL)

# See how well the models perform given the restrictions in range between the different 
# bounding boxes that can be used


# Data --------------------------------------------------------------------

# NB: The substrate variables need to be removed as they can't be used in the final data
# to predict kelp cover presence because we don't know what they are everywhere

# All of the BO variables matched to sites
kelp_all <- adf %>% 
  dplyr::select(Campaign, site, depth, -c(Bedrock..:sand), kelp.cover, Laminariales, Agarum, Alaria) %>% 
  left_join(study_site_env, by = c("Campaign", "site")) %>%
  mutate(kelp.cover = ifelse(kelp.cover > 100, 100, kelp.cover)) %>% # Correct values over 100
  dplyr::select(-lon_env, -lat_env, -lon, -lat) %>%
  dplyr::select(-depth) %>% # This is not used in Jesi's model
  dplyr::select(-bathy, -land_distance) %>% # Decided against these variables
  na.omit() # No missing data

# The mean kelp covers per site/depth
kelp_all_mean <- kelp_all %>% 
  group_by(Campaign, site) %>%
  summarise_all(mean) %>%
  ungroup()

# The highest kelp covers per site/depth
kelp_all_max <- kelp_all %>% 
  group_by(Campaign, site) %>%
  summarise_all(max) %>%
  ungroup()


# Data prep function ------------------------------------------------------

# Convenience function for prepping dataframe for use in the random forest
# This removes all other kelp cover values
rf_data_prep <- function(kelp_choice, df = kelp_all){
  df_1 <- data.frame(dplyr::select(df, -c(Campaign:site))) %>% 
    pivot_longer(cols = kelp.cover:Alaria, names_to = "chosen_kelp", values_to = "cover") %>% 
    filter(chosen_kelp == kelp_choice)
}


# Single rule model -------------------------------------------------------

# Function for OneR model run
OneR_model <- function(kelp_choice, df = kelp_all){

  # Pull out training index
  one_random <- sample(1:nrow(df), 0.7*nrow(df))
  
  # Prep the data
  kelp_cat <- rf_data_prep(kelp_choice = kelp_choice, df = df)
           
  # The training data
  kelp_train <- kelp_cat %>% 
    slice(one_random) %>%
    data.frame() %>% 
    optbin(., method = "infogain") # The choice of method does have an effect on the results
  
  # The test data
  kelp_test <- kelp_cat %>% 
    slice(-one_random) %>% 
    data.frame()
  
  # Various tests and results
  kelp_model <- OneR(kelp_train, verbose = TRUE)
  summary(kelp_model)
  prediction <- predict(kelp_model, kelp_test)
  eval_model(prediction, kelp_test)
  # plot(kelp_model) # Can't plot more than 20 variables
}

# Run models for each type of kelp cover
# OneR_model("kelp.cover")
# OneR_model("Laminariales")
# OneR_model("Agarum")
# OneR_model("Alaria")


# Which variables are the most important? ---------------------------------

# Convenience wrapper for extracting variable importance from an RF
extract_var_imp <- function(kelp_rf){
    res <- data.frame(var = row.names(kelp_rf$importance), 
                      kelp_rf$importance, 
                      importanceSD = kelp_rf$importanceSD,
                      mean_MSE = mean(kelp_rf$mse),
                      mean_rsq = mean(kelp_rf$rsq)) %>% 
      arrange(-X.IncMSE) %>% 
      mutate_if(is.numeric, round, 4)
  res$var <- as.character(res$var)
  return(res)
}

# This function runs many random forest models to determine which variables
# are consistently the most important
top_var <- function(lplyr_bit, kelp_choice, df = kelp_all){
  
  # Prep the four possibilities for a single kelp cover choice
  df_prep <- rf_data_prep(kelp_choice, df = df) %>% 
    dplyr::select(-chosen_kelp)
  
  # Random sampling to split data up for training
  train <- sample(1:nrow(df_prep), 0.7*nrow(df_prep), replace = FALSE)
  
  # Random forest models for the four possibilities
  rf_reg <- randomForest(cover ~ ., data = df_prep[train,], ntree = 200, importance = TRUE, do.trace = F)

  # Extract results
  var_reg <- extract_var_imp(rf_reg)
  return(var_reg)
}
# top_var(kelp_choice = "Agarum", df = kelp_all)

# Convenience wrapper to remove correlated variables
# cor_var_rm <- function(df_multi){
#   
#   # Order the dataframe based on the Increase in MSE
#   df_ordered <- arrange(df_multi, -X.IncMSE)
#   
#   # Remove variables that correlate with better predictors
#   row_i <- 2
#   cor_df <- df_ordered
#   while(row_i < nrow(cor_df)){
#     cor_cols <- cor_df[1:row_i-1, "var"]
#     cor_check <- cor_df[row_i, "var"]
#     BO_cor_check <- BO_cor_matrix %>% 
#       dplyr::select(Parameter1, cor_cols$var) %>% 
#       filter(Parameter1 == cor_check$var) %>% 
#       pivot_longer(cols = -Parameter1) %>% 
#       filter(abs(value) >= 0.7)
#     if(nrow(BO_cor_check) > 0){
#       cor_df <- cor_df[-row_i,]
#     } else{
#       row_i <- row_i+1
#     }
#   }
#   return(cor_df)
# }

# We then run this 1000 times to increase our certainty in the findings
top_var_multi <- function(kelp_choice, df = kelp_all){
  
  # Run 100 models
  multi_kelp <- plyr::ldply(.data = 1:1000, .fun = top_var, .parallel = T, 
                            kelp_choice = kelp_choice, df = df)
  
  # Clean up the results
  multi_kelp_mean <- multi_kelp %>% 
    group_by(var) %>% 
    summarise_all(mean) %>% 
    ungroup() %>% 
    arrange(-X.IncMSE)
  
  # Remove correlated variables and exit
  # res <- cor_var_rm(multi_kelp_mean)
  return(multi_kelp_mean)
}

## Find the top variables for the different kelp covers
# registerDoParallel(cores = 50)

# kelp.cover
# system.time(top_var_kelpcover <- top_var_multi("kelp.cover")) # ~16 seconds on 50 cores
# save(top_var_kelpcover, file = "data/top_var_kelpcover.RData")

# Laminariales
# top_var_laminariales <- top_var_multi("Laminariales")
# save(top_var_laminariales, file = "data/top_var_laminariales.RData")

# # Agarum
# top_var_agarum <- top_var_multi("Agarum")
# save(top_var_agarum, file = "data/top_var_agarum.RData")

# Alaria
# top_var_alaria <- top_var_multi("Alaria")
# save(top_var_alaria, file = "data/top_var_alaria.RData")


# Random Forest function --------------------------------------------------

# The base function used for the random forest
# It is informed by the top variables from the previous section

# Load top variable dataframes
load("data/top_var_kelpcover.RData")
load("data/top_var_laminariales.RData")
load("data/top_var_agarum.RData")
load("data/top_var_alaria.RData")

# testers...
# kelp_choice <- "kelp.cover"
# column_choice <- top_var_kelpcover
# kelp_choice <- "Agarum"
# column_choice <- top_var_agarum
# kelp_choice <- "Laminariales"
# column_choice <- top_var_laminariales
random_kelp_forest <- function(lply_bit, kelp_choice, column_choice,
                               df = kelp_all, print_res = F){
  
  # Extract only the kelp cover of choice
  df_prep <- rf_data_prep(kelp_choice, df)
  
  # Prep the selection columns
  top_var_reg <- as.character(column_choice$var)
  
  # Chose only desired columns
  df_var_reg <- dplyr::select(df_prep, cover, all_of(top_var_reg))
  
  # Split data up for training and testing
  train <- sample(nrow(df), 0.7*nrow(df), replace = FALSE)
  
  # Random forest model regression
  rf_reg <- randomForest(cover ~ ., data = df_var_reg[train, ], 
                         ntree = 200, importance = TRUE, do.trace = F)
  
  # Predicting on training set
  pred_reg_train <- round(predict(rf_reg, df_var_reg[train, ]), 2)
  
  # Predicting on Validation set
  pred_reg_valid <- round(predict(rf_reg, df_var_reg[-train, ]), 2)
  
  # Create dataframe of accuracy results
  # mean(rf_reg_base$rsq)
  # rf_reg_base$predicted
  accuracy_reg <- data.frame(portion = c(rep("train", length(train)), 
                                         rep("validate", nrow(df)-length(train))),
                             pred = c(pred_reg_train, pred_reg_valid), 
                             original = c(df_var_reg[train, ]$cover,
                                          df_var_reg[-train, ]$cover)) %>% 
    mutate(accuracy = as.numeric(pred)-as.numeric(original))
  
  # Either print the results or return them as a list
  if(print_res){
    # Regression base
    print(rf_reg); varImpPlot(rf_reg, main = kelp_choice)
  } else {
    # Package the model up with the accuracy results and exit
    res <- list(rf_reg = rf_reg, accuracy_reg = accuracy_reg)
  }
}

# Check the random forests
# random_kelp_forest(kelp_choice = "kelp.cover", column_choice = top_var_kelpcover, print_res = T)
# random_kelp_forest(kelp_choice = "Laminariales", column_choice = top_var_laminariales, print_res = T)
# random_kelp_forest(kelp_choice = "Agarum", column_choice = top_var_agarum, print_res = T)
# random_kelp_forest(kelp_choice = "Alaria", column_choice = top_var_alaria, print_res = T)

# For grabbing a single tree
# getTree(rf_reg, k = 1, labelVar = FALSE)


# Calculate and merge 1000 models -----------------------------------------

# Function for projecting a model
# model_choice <- multi_test[[1]]
project_cover <- function(model_choice){
  pred_df <- data.frame(lon = Arctic_env$lon, lat = Arctic_env$lat,
                        depth = Arctic_env$bathy, land_distance = Arctic_env$land_distance,
                        pred_present = predict(model_choice$rf_reg, Arctic_env),
                        pred_2050 = predict(model_choice$rf_reg, Arctic_env_2050),
                        pred_2100 = predict(model_choice$rf_reg, Arctic_env_2100)) %>% 
    mutate(pred_present = case_when(pred_present < 0 ~ 0, TRUE ~ pred_present),
           pred_2050 = case_when(pred_present < 0 ~ 0, TRUE ~ pred_2050),
           pred_2100 = case_when(pred_present < 0 ~ 0, TRUE ~ pred_2100))
  gc(); return(pred_df)
}

# Now we run the test on each kelp cover 1000 times to see what the spread is
# in the accuracy of the random forests
# This is caused by different random splitting of test/validation sets
# as well as the many possible routes that the random forest may then take
# kelp_choice <- "Laminariales"
# column_choice <- top_var_laminariales
random_kelp_forest_select <- function(kelp_choice, column_choice, df = kelp_all){
  # system.time(
  multi_test <- plyr::llply(.data = 1:1000, .fun = random_kelp_forest, .parallel = T, 
                            kelp_choice = kelp_choice, column_choice = column_choice, df = df)
  # ) # ~18 seconds
  gc()
    
  # Create mean projection from all models for present data
  # system.time(
  project_multi <- plyr::ldply(multi_test, project_cover, .parallel = T) %>%
    group_by(lon, lat, depth, land_distance) %>% 
    summarise(pred_present_mean = round(mean(pred_present, na.rm = T)),
              pred_2050_mean = round(mean(pred_2050, na.rm = T)),
              pred_2100_mean = round(mean(pred_2100, na.rm = T)),
              pred_present_sd = round(sd(pred_present, na.rm = T), 2),
              pred_2050_sd = round(sd(pred_2050, na.rm = T), 2), 
              pred_2100_sd = round(sd(pred_2100, na.rm = T), 2), 
              .groups = "drop")
  # ) # 169 seconds
  gc()
  
  # Extract the accuracy of the models 
  accuracy_reg <- lapply(multi_test, function(x) x$accuracy_reg) %>% 
    do.call(rbind.data.frame, .) %>% 
    mutate(model_id = rep(1:1000, each = nrow(df)))
  
  # Find which model had the best validation scores
  accuracy_reg_check <- accuracy_reg %>%
    filter(portion == "validate") %>% 
    group_by(model_id) %>% 
    summarise(mean_acc = round(mean(abs(accuracy)), 2),
              r_acc = round(cor(x = original, y = pred), 2),
              .groups = "drop")
  
  # Find the best model ID
  best_reg <- arrange(accuracy_reg_check, -r_acc, mean_acc)[1,]
  
  # Extract the best model
  choice_reg <- multi_test[[best_reg$model_id]]$rf_reg
  
  # Combine and exit
  res <- list(choice_reg = choice_reg, 
              accuracy_reg = accuracy_reg, 
              project_multi = project_multi)
  return(res)
}

# Set cores
doParallel::registerDoParallel(cores = 50)

# Kelp.cover
# system.time(best_rf_kelpcover <- random_kelp_forest_select("kelp.cover", top_var_kelpcover)) # 230 seconds with 50 cores
# save(best_rf_kelpcover, file = "data/best_rf_kelpcover.RData", compress = T)

# Laminariales
# best_rf_laminariales <- random_kelp_forest_select("Laminariales", top_var_laminariales)
# save(best_rf_laminariales, file = "data/best_rf_laminariales.RData", compress = T)

# Agarum
# best_rf_agarum <- random_kelp_forest_select("Agarum", top_var_agarum)
# save(best_rf_agarum, file = "data/best_rf_agarum.RData", compress = T)

# Alaria
# best_rf_alaria <- random_kelp_forest_select("Alaria", top_var_alaria)
# save(best_rf_alaria, file = "data/best_rf_alaria.RData", compress = T)


# Analyse model accuracy --------------------------------------------------

# First load the best random forest models produced above
load("data/best_rf_kelpcover.RData")
load("data/best_rf_laminariales.RData")
load("data/best_rf_agarum.RData")
load("data/best_rf_alaria.RData")

# Find the distributions of accuracy from 0 - 100%
test_acc <- best_rf_kelpcover$accuracy_reg #%>%
  # filter(model_id == 1000)

# Quick visuals
ggplot(filter(test_acc, portion == "validate"), aes(x = accuracy)) +
  geom_histogram()

ggplot(filter(test_acc, portion == "validate"), aes(x = original, y = pred)) +
  geom_point() +
  geom_smooth(method = "lm")

ggplot(filter(test_acc, portion == "validate"), aes(x = as.factor(original), y = pred)) +
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
  
  # The mean CI per model run
  conf_mean <- df %>% 
    filter(portion == "validate") %>% 
    group_by(model_id) %>% 
    mutate(accuracy = round(accuracy)) %>% 
    summarise(mean_acc = mean(abs(accuracy)),
              sd_acc = sd(abs(accuracy)),
              r_acc = cor(x = original, y = pred)) %>% 
    ungroup()
  
  # Easy labels for plotting
  conf_mean_label <- conf_mean %>% 
    summarise(mean_acc = round(mean(abs(mean_acc))),
              sd_acc = round(mean(abs(sd_acc))),
              r_acc = round(mean(r_acc), 2))
  
  # Label for best model
  conf_best_label <- conf_mean %>% 
    filter(mean_acc == min(mean_acc)) %>% 
    mutate(mean_acc = round(abs(mean_acc)),
           sd_acc = round(sd_acc),
           r_acc = round(r_acc, 2))
  
  # CI of best model
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
# conf_plot(best_rf_kelpcover$accuracy_reg, "Total cover confidence")
# conf_plot(best_rf_laminariales$accuracy_reg, "Laminariales cover confidence")
# conf_plot(best_rf_agarum$accuracy_reg, "Agarum cover confidence")
# conf_plot(best_rf_alaria$accuracy_reg, "Alaria cover confidence")


# Visualise kelp cover projections ----------------------------------------

# First load the best random forest models produced above
load("data/best_rf_kelpcover.RData")
load("data/best_rf_laminariales.RData")
load("data/best_rf_agarum.RData")
load("data/best_rf_alaria.RData")

# Visualise a family of cover
cover_squiz <- function(best_rf, legend_title, x_nudge, kelp_choice, pred_choice){
  
  # Prep point data
  data_point <- adf %>% 
    dplyr::select(Campaign, site, depth, -c(Bedrock..:sand), kelp.cover, Laminariales, Agarum, Alaria) %>% 
    left_join(study_site_env, by = c("Campaign", "site")) %>%
    mutate(kelp.cover = ifelse(kelp.cover > 100, 100, kelp.cover)) %>%
    dplyr::select(lon, lat, {{kelp_choice}}) %>% 
    group_by(lon, lat) %>%
    summarise_all(mean) %>%
    ungroup()
  
  # Prep projection data
  project_single <- best_rf$project_multi
  
  # Create plot
  Arctic_map +
    # geom_tile(data = df, # No filter
    geom_tile(data = filter(project_single, depth <= 100 | land_distance <= 100),
              aes_string(x = "lon", y = "lat", fill = pred_choice)) +
    geom_point(data = data_point, colour = "red", shape = 21,
               aes_string(x = "lon", y = "lat", size = kelp_choice)) +
    scale_fill_viridis_c(legend_title) +
    theme(strip.background = element_rect(colour = "white", fill = "white"),
          legend.position = c(x_nudge, 0.96),
          legend.direction = "horizontal",
          legend.spacing.y = unit(0, "mm")) +
    labs(size = legend_title)
}

## Visualisations
# Total cover
# cover_squiz(best_rf_kelpcover, "Total cover (%)", 0.785, "kelp.cover", "pred_present_mean")
# cover_squiz(best_rf_kelpcover, "Total cover (%)", 0.785, "kelp.cover", "pred_2050_mean")
# cover_squiz(best_rf_kelpcover, "Total cover (%)", 0.785, "kelp.cover", "pred_2100_mean")
# Laminariales
# cover_squiz(best_rf_laminariales, "Laminariales cover (%)", 0.745, "Laminariales", "pred_present_mean")
# cover_squiz(best_rf_laminariales, "Laminariales cover (%)", 0.745, "Laminariales", "pred_2050_mean")
# cover_squiz(best_rf_laminariales, "Laminariales cover (%)", 0.745, "Laminariales", "pred_2100_mean")
# Alaria
# cover_squiz(best_rf_alaria, "Alaria cover (%)", 0.78, "Alaria", "pred_present_mean")
# cover_squiz(best_rf_alaria, "Alaria cover (%)", 0.78, "Alaria", "pred_2050_mean")
# cover_squiz(best_rf_alaria, "Alaria cover (%)", 0.78, "Alaria", "pred_2100_mean")
# Agarum
# cover_squiz(best_rf_agarum, "Agarum cover (%)", 0.77, "Agarum", "pred_present_mean")
# cover_squiz(best_rf_agarum, "Agarum cover (%)", 0.77, "Agarum", "pred_2050_mean")
# cover_squiz(best_rf_agarum, "Agarum cover (%)", 0.77, "Agarum", "pred_2100_mean")


# Compare future projections to present -----------------------------------

project_compare <- function(best_rf, kelp_choice){
  
  # Calculate differences
  project_diff <- best_rf$project_multi %>% 
    mutate(pred_diff_2050 = plyr::round_any(pred_2050_mean - pred_present_mean, 20),
           pred_diff_2100 = pred_2100_mean - pred_present_mean) %>% 
    mutate(pred_diff_2050 = ifelse(pred_diff_2050 == 0, NA, pred_diff_2050))
    # pivot_longer(cols = pred_present_mean:pred_diff_2100) %>% 
    # filter(name == "pred_diff_2050" & value != 0) %>% 
    # pivot_wider()
    # mutate(pred_diff_2050 = base::cut(pred_2050_mean - pred_present_mean, breaks = c(-0.1, 0, 0.1, 0.2)),
    #        pred_diff_2100 = base::cut(pred_2100_mean - pred_present_mean, breaks = c(-0.1, 0, 0.1, 0.2)))
  
  # Scale for difference plots
  diff_range <- range(c(project_diff$pred_diff_2050, project_diff$pred_diff_2100), na.rm = T)
  
  # Present plot
  p_present <-  Arctic_map +
    # geom_tile(data = df, # No filter
    geom_tile(data = filter(project_diff, depth <= 100 | land_distance <= 100),
              aes(x = lon, y = lat, fill = pred_present_mean)) +
    scale_fill_viridis_c(paste0(kelp_choice," cover (%)")) +
    theme(legend.position = "bottom")
  
  # 2050 plot
  p_2050 <-  Arctic_map +
    # geom_tile(data = df, # No filter
    geom_tile(data = filter(project_diff, depth <= 100 | land_distance <= 100),
              aes(x = lon, y = lat, fill = pred_diff_2050)) +
    scale_fill_gradient2("2050 - present (%)", low = "blue", high = "red", limits = diff_range, guide = "legend", na.value = NA) +
    # scale_fill_discrete("2050 - present (%)") +
    theme(legend.position = "bottom",
          axis.ticks.y = element_blank(),
          axis.text.y = element_blank())
  
  # 2050 plot
  p_2100 <-  Arctic_map +
    # geom_tile(data = df, # No filter
    geom_tile(data = filter(project_diff, depth <= 100 | land_distance <= 100),
              aes(x = lon, y = lat, fill = pred_diff_2100)) +
    scale_fill_gradient2("2100 - present (%)", low = "blue", high = "red", limits = diff_range) +
    # scale_fill_discrete("2100 - present (%)") +
    theme(legend.position = "bottom",
          axis.ticks.y = element_blank(),
          axis.text.y = element_blank())
  
  # Save and Print
  p_all <- ggpubr::ggarrange(p_present, p_2050, p_2100, ncol = 3, nrow = 1, align = "hv")
  ggsave(paste0("graph/future_diff_",kelp_choice,".png"), p_all, width = 9, height = 5)
  p_all
}

# Visualise the comparisons
project_compare(best_rf_kelpcover, "Total_cover")
# project_compare(best_rf_laminariales, "Laminariales")
# project_compare(best_rf_agarum, "Agarum")
project_compare(best_rf_alaria, "Alaria")


# Relationship between cover and physical variables -----------------------

# Run a regression to see in which direction the relationships with percent cover 
# are with the top variables. e.g. more cover with more iron

