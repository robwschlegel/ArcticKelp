# analyses/7_random_forests.R
# The purpose of this script is to house the code used for the random forest analyses


# Setup -------------------------------------------------------------------

# Load all data nad previous libraries
source("analyses/4_kelp_cover.R")

# Libraries for this script specifically
library(randomForest)
library(OneR) # For single rule machine learning
# library(caret) # For cross validation option
library(doParallel); doParallel::registerDoParallel(cores = 50) # This will be between 4 - 8 on a laptop

# Environmental data per site
load("data/study_site_env.RData")

# Load Arctic data for testing variable correlations
load("data/Arctic_env.RData")

# Load the BO correlation matrix
load("data/BO_cor_matrix.RData")

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

# NB: The substrate variables need to be removed as they can't be used in the final data
# to predict kelp cover presence because we don't know what they are everywhere

# All of the BO variables matched to sites
kelp_all <- adf %>% 
  dplyr::select(Campaign, site, depth, -c(Bedrock..:sand), kelp.cover, Laminariales, Agarum, Alaria) %>% 
  left_join(study_site_env, by = c("Campaign", "site")) %>%
  mutate(kelp.cover = ifelse(kelp.cover > 100, 100, kelp.cover)) %>% # Correct values over 100
  dplyr::select(-lon_env, -lat_env, -env_index, -lon, -lat) %>%
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

# Scenarios
  # Rather use the same number of variables for each scenario
  # This means using the variables that don't have projections
  # with the projected variables for the two different time tests
base <- colnames(dplyr::select(kelp_all, 
                               BO2_templtmin_bdmax:BO2_curvelltmax_bdmax))
future_2050 <- colnames(dplyr::select(kelp_all, 
                                      BO_parmean:BO2_phosphateltmax_bdmax,
                                      BO2_RCP85_2050_curvelltmax_bdmax:BO2_RCP85_2050_tempmean_ss))
future_2100 <- colnames(dplyr::select(kelp_all, 
                                      BO_parmean:BO2_phosphateltmax_bdmax,
                                      BO2_RCP85_2100_curvelltmax_bdmax:BO2_RCP85_2100_tempmean_ss))


# Data prep function ------------------------------------------------------

# Convenience function for prepping dataframe for use in the random forest
# This removes all other kelp cover values
rf_data_prep <- function(kelp_choice, df = kelp_all, cat_cover = F, scenario = base){
  
  # Trim down data.frame
  df_1 <- data.frame(dplyr::select(df, -c(Campaign:site))) %>% 
    pivot_longer(cols = kelp.cover:Alaria, names_to = "chosen_kelp", values_to = "cover") %>% 
    filter(chosen_kelp == kelp_choice) %>% 
    dplyr::select(scenario, cover)
  
  # Cut cover into categories if desired
  if(cat_cover){
    df_1$cover <- cut(df_1$cover, breaks = c(-Inf, 10, 50, 100), ordered_result = T)
    levels(df_1$cover)[1] <- "[0,10]"
  }
  
  return(df_1)
}


# Single rule model -------------------------------------------------------

# Function for OneR model run
OneR_model <- function(kelp_choice, df = kelp_all){

  # Pull out training index
  one_random <- sample(1:nrow(df), 0.7*nrow(df))
  
  # Prep the data
  kelp_cat <- rf_data_prep(kelp_choice = kelp_choice, df = df, cat_cover = T)
           
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
  plot(kelp_model)
  prediction <- predict(kelp_model, kelp_test)
  eval_model(prediction, kelp_test)
}

# Run models for each type of kelp cover
# OneR_model("kelp.cover")
# OneR_model("Laminariales")
# OneR_model("Agarum")
# OneR_model("Alaria")


# Which variables are the most important? ---------------------------------

# Convenience wrapper for extracting variable importance from an RF
extract_var_imp <- function(kelp_rf){
  if(ncol(kelp_rf$importance) == 5){
    res <- data.frame(var = row.names(kelp_rf$importance), 
                      kelp_rf$importance[,4:5]) %>% 
      arrange(-MeanDecreaseAccuracy) %>% 
      mutate_if(is.numeric, round, 4)
  } else {
    res <- data.frame(var = row.names(kelp_rf$importance), 
                      kelp_rf$importance, 
                      importanceSD = kelp_rf$importanceSD,
                      mean_MSE = mean(kelp_rf$mse),
                      mean_rsq = mean(kelp_rf$rsq)) %>% 
      arrange(-X.IncMSE) %>% 
      mutate_if(is.numeric, round, 4)
  }
  res$var <- as.character(res$var)
  return(res)
}

# This function runs many random forest models to determine which variables
# are consistently the most important
top_var <- function(lplyr_bit, kelp_choice, df = kelp_all){
  
  # Prep the four possibilities for a single kelp cover choice
  df_reg <- rf_data_prep(kelp_choice, df = df, cat_cover = F, scenario = base)
  df_cat <- rf_data_prep(kelp_choice, df = df, cat_cover = T, scenario = base)
  
  # Random sampling to split data up for training
  train <- sample(1:nrow(df_reg), 0.7*nrow(df_reg), replace = FALSE)
  
  # Random forest models for the four possibilities
  rf_reg <- randomForest(cover ~ ., data = df_reg[train,], ntree = 200, importance = TRUE, do.trace = F)
  rf_cat <- randomForest(cover ~ ., data = df_cat[train,], ntree = 200, importance = TRUE, do.trace = F)

  # Extract results
  var_reg <- extract_var_imp(rf_reg)
  var_cat <-  extract_var_imp(rf_cat)
  
  # Combine and exit
  res <- left_join(var_reg, var_cat, by = "var")
  return(res)
}
# top_var(kelp_choice = "Agarum", df = kelp_all)

# Convenience wrapper to remove correlated variables
cor_var_rm <- function(df){
  
  # Order the dataframe based on the second column
  df_order <- order(as.vector(as.data.frame(df[,2])), decreasing = T)
  df_ordered <- df[df_order,]
  
  # Remove variables that correlate with better predictors
  row_i <- 2
  cor_df <- df_ordered
  while(row_i < nrow(cor_df)){
    cor_cols <- cor_df[1:row_i-1, "var"]
    cor_check <- cor_df[row_i, "var"]
    BO_cor_check <- BO_cor_matrix %>% 
      dplyr::select(var, cor_cols$var) %>% 
      filter(var == cor_check$var) %>% 
      pivot_longer(cols = -var) %>% 
      filter(value >= 0.7)
    if(nrow(BO_cor_check)  > 0){
      cor_df <- cor_df[-row_i,]
    } else{
      row_i <- row_i+1
    }
  }
  return(cor_df)
}

# We then run this 100 times to increase our certainty in the findings
top_var_multi <- function(kelp_choice, df = kelp_all){
  
  # Run 100 models
  multi_kelp <- plyr::ldply(.data = 1:100, .fun = top_var, .parallel = T, 
                            kelp_choice = kelp_choice, df = df)
  
  # Clean up the results
  multi_kelp_mean <- multi_kelp %>% 
    group_by(var) %>% 
    summarise_all(mean) %>% 
    ungroup()
  
  # Remove correlated variables and exit
  res <- list(reg_var = cor_var_rm(multi_kelp_mean[,1:6]),
              cat_var = cor_var_rm(multi_kelp_mean[,c(1,7:8)]))
  return(res)
}

## Find the top variables for the different kelp covers
# kelp.cover
# system.time(top_var_kelpcover <- top_var_multi("kelp.cover")) # ~5 seconds on 50 cores
# save(top_var_kelpcover, file = "data/top_var_kelpcover.RData")

# Laminariales
# top_var_laminariales <- top_var_multi("Laminariales")
# save(top_var_laminariales, file = "data/top_var_laminariales.RData")

# Agarum
# top_var_agarum <- top_var_multi("Agarum")
# save(top_var_agarum, file = "data/top_var_agarum.RData")

# Alaria
# top_var_alaria <- top_var_multi("Alaria")
# save(top_var_alaria, file = "data/top_var_alaria.RData")


# Random Forest function --------------------------------------------------

# This section is largely redundant to the following section
# The difference here is that this function spits out the results
# without saving anything to the environment

# Load top variable dataframes
load("data/top_var_kelpcover.RData")
load("data/top_var_laminariales.RData")
load("data/top_var_agarum.RData")
load("data/top_var_alaria.RData")

# testers...
# kelp_choice <- "Agarum"
# column_choice <- top_var_agarum
random_kelp_forest <- function(lply_bit, kelp_choice, column_choice, 
                               df = kelp_all, scenario = base, print_res = F){
  
  # Extract only the kelp cover of choice
  df_reg <- rf_data_prep(kelp_choice, df, cat_cover = F, scenario = scenario)
  df_cat <- rf_data_prep(kelp_choice, df, cat_cover = T, scenario = scenario)
  
  # Prep the selection columns
  top_var_reg <- as.character(column_choice$reg_var$var)
  top_var_cat <- as.character(column_choice$cat_var$var)
  
  # Chose only desired columns
  df_var_reg <- dplyr::select(df_reg, cover, all_of(top_var_reg))
  df_var_cat <- dplyr::select(df_cat, cover, all_of(top_var_cat))
  
  # Split data up for training and testing
  train <- sample(nrow(df), 0.7*nrow(df), replace = FALSE)
  
  # Random forest model regression
  rf_reg <- randomForest(cover ~ ., data = df_var_reg[train, ], 
                         ntree = 200, importance = TRUE, do.trace = F)
  
  # Random forest model category
  rf_cat <- randomForest(cover ~ ., data = df_var_cat[train, ],
                         ntree = 200, importance = TRUE, do.trace = F)
  
  # Predicting on training set
  pred_reg_train <- round(predict(rf_reg, df_var_reg[train, ]), 2)
  pred_cat_train <- predict(rf_cat, df_var_cat[train, ])
  
  # Predicting on Validation set
  pred_reg_valid <- round(predict(rf_reg, df_var_reg[-train, ]), 2)
  pred_cat_valid <- predict(rf_cat, df_var_cat[-train, ])
  
  # Create dataframe of accuracy results
  # mean(rf_reg_base$rsq)
  # rf_reg_base$predicted
  accuracy_reg <- data.frame(portion = c(rep("train", length(train)), 
                                         rep("validate", nrow(df)-length(train))),
                             pred = c(pred_reg_train, pred_reg_valid), 
                             original = c(df_var_reg[train, ]$cover,
                                          df_var_reg[-train, ]$cover)) %>% 
    mutate(accuracy = as.numeric(pred)-as.numeric(original))
  accuracy_cat <- data.frame(portion = c(rep("train", length(train)), 
                                         rep("validate", nrow(df)-length(train))),
                             pred = c(pred_cat_train, pred_cat_valid), 
                             original = c(df_var_cat[train, ]$cover,
                                          df_var_cat[-train, ]$cover)) %>% 
    mutate(accuracy = ifelse(pred == original, 1, 0)) 
  
  # Either print the results or return them as a list
  if(print_res){
    # Regression base
    print(rf_reg); varImpPlot(rf_reg, main = kelp_choice)
    # Category base
    print(rf_cat); varImpPlot(rf_cat, main = kelp_choice)
  } else {
    # Package the model up with the accuracy results and exit
    res <- list(rf_reg = rf_reg, accuracy_reg = accuracy_reg,
                rf_cat = rf_cat, accuracy_cat = accuracy_cat)
  }
}

# Check the random forests
# random_kelp_forest(kelp_choice = "kelp.cover", column_choice = top_var_kelpcover, print_res = T)
# random_kelp_forest(kelp_choice = "Laminariales", column_choice = top_var_laminariales, print_res = T)
# random_kelp_forest(kelp_choice = "Agarum", column_choice = top_var_agarum, print_res = T)
# random_kelp_forest(kelp_choice = "Alaria", column_choice = top_var_alaria, print_res = T)

# For grabbing a single tree
# getTree(rfobj, k=1, labelVar=FALSE)


# Select best random forest -----------------------------------------------

# Now we run the test on each kelp cover 1000 times to see what the spread is
# in the accuracy of the random forests
# This is caused by different random splitting of test/validation sets
# as well as the many possible routes that the random forest may then take
random_kelp_forest_select <- function(kelp_choice, column_choice, 
                                      df = kelp_all, scenario = base){
  # system.time(
    multi_test <- plyr::llply(.data = 1:100, .fun = random_kelp_forest, .parallel = T, 
                              kelp_choice = kelp_choice, column_choice = column_choice,
                              df = df, scenario = scenario)
    # ) # ~50 seconds
  
  # Extract the model accuracies
  accuracy_reg <- lapply(multi_test, function(x) x$accuracy_reg) %>% 
    do.call(rbind.data.frame, .) %>% 
    mutate(model_id = rep(1:100, each = nrow(df)))
  accuracy_cat <- lapply(multi_test, function(x) x$accuracy_cat) %>% 
    do.call(rbind.data.frame, .) %>% 
    mutate(model_id = rep(1:100, each = nrow(df)))
  
  # Find which model had the best validation scores
  accuracy_reg_check <- accuracy_reg %>%
    filter(portion == "validate") %>% 
    group_by(model_id) %>% 
    summarise(mean_acc = round(mean(abs(accuracy)), 2),
              r_acc = round(cor(x = original, y = pred), 2)) %>% 
    ungroup()
  accuracy_cat_check <- accuracy_cat %>% 
    filter(portion == "validate") %>% 
    group_by(model_id) %>% 
    mutate(count_all = n()) %>% 
    filter(accuracy == 0) %>% 
    mutate(count_0 = n(),
           mean_acc = round(count_0/count_all, 4)*100) %>% 
    ungroup() %>% 
    dplyr::select(model_id, mean_acc) %>% 
    unique()
  
  # Find the best model ID
  best_reg <- arrange(accuracy_reg_check, -r_acc, mean_acc)[1,]
  best_cat <- arrange(accuracy_cat_check, -mean_acc)[1,]
  
  # Extract the best model
  choice_reg <- multi_test[[best_reg$model_id]]$rf_reg
  choice_cat <- multi_test[[best_cat$model_id]]$rf_cat
  
  # Combine and exit
  res <- list(choice_reg = choice_reg, accuracy_reg = accuracy_reg, 
              choice_cat = choice_cat, accuracy_cat = accuracy_cat)
  return(res)
}

# doParallel::registerDoParallel(cores = 50)
# Kelp.cover
# system.time(best_rf_kelpcover <- random_kelp_forest_select("kelp.cover", top_var_kelpcover)) # 3 seconds with 50 cores
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

## NB: It appears that the category models are massively overfitting


# Analyse model accuracy --------------------------------------------------

# First load the best random forest models produced above
load("data/best_rf_kelpcover.RData")
load("data/best_rf_laminariales.RData")
load("data/best_rf_agarum.RData")
load("data/best_rf_alaria.RData")

# Find the distributions of accuracy from 0 - 100%
# test_acc <- best_rf_kelpcover$accuracy_reg #%>%
  # filter(model_id == 1000)

# Quick visuals
# ggplot(filter(test_acc, portion == "validate"), aes(x = accuracy)) +
  # geom_histogram()

# ggplot(filter(test_acc, portion == "validate"), aes(x = original, y = pred)) +
  # geom_point() +
  # geom_smooth(method = "lm")

# ggplot(filter(test_acc, portion == "validate"), aes(x = as.factor(original), y = pred)) +
  # geom_boxplot()

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


# Project kelp cover in the Arctic ----------------------------------------

# First load the best random forest models produced above
# load("data/best_rf_kelpcover.RData")
# load("data/best_rf_laminariales.RData")
# load("data/best_rf_agarum.RData")
# load("data/best_rf_alaria.RData")

# Load top variable choices
# load("data/top_var_kelpcover.RData")
# load("data/top_var_laminariales.RData")
# load("data/top_var_agarum.RData")
# load("data/top_var_alaria.RData")

# Load the Arctic data
# load("data/Arctic_env.RData")

# This function changes the variable names in the future layers to match the expected names
predict_future <- function(model_choice, scenario){
  
  # Filter out only columns in chosen scenario
  Arctic_env_sub <- dplyr::select(Arctic_env, all_of(scenario))
  
  # Translate best columns choice from future projections as necessary
  if(tail(scenario, 1) == "BO2_RCP85_2050_tempmean_ss"){
    colnames(Arctic_env_sub) <- str_replace(colnames(Arctic_env_sub), "BO2_RCP85_2050_", "BO2_")
  } else if(tail(scenario, 1) == "BO2_RCP85_2100_tempmean_ss"){
    colnames(Arctic_env_sub) <- str_replace(colnames(Arctic_env_sub), "BO2_RCP85_2100_", "BO2_")
  } else {
  }
  
  # Run prediction and exit
  pred_val <-  predict(model_choice, Arctic_env_sub)
  return(pred_val)
}


# Convenience function for final step before prediction
Arctic_cover_predict <- function(model_choice, scenario){
  pred_df <- data.frame(lon = Arctic_env$lon, lat = Arctic_env$lat,
                        depth = Arctic_env$bathy,
                        land_distance = Arctic_env$land_distance,
                        pred_val = predict_future(model_choice, scenario))
}

# Predict the covers
# pred_kelpcover <- Arctic_cover_predict(best_rf_kelpcover$choice_reg, base)
# pred_laminariales <- Arctic_cover_predict(best_rf_laminariales$choice_reg, base)
# pred_agarum <- Arctic_cover_predict(best_rf_agarum$choice_reg, base)
# pred_agarum_2050 <- Arctic_cover_predict(best_rf_agarum$choice_reg, future_2050)
# pred_agarum_2100 <- Arctic_cover_predict(best_rf_agarum$choice_reg, future_2100)
# pred_alaria <- Arctic_cover_predict(best_rf_alaria$choice_reg, base)

# Visualise a family of cover
cover_squiz <- function(df, legend_title, x_nudge, kelp_choice){
  
  # Prep point data
  data_point <- adf %>% 
    dplyr::select(Campaign, site, depth, -c(Bedrock..:sand), kelp.cover, Laminariales, Agarum, Alaria) %>% 
    left_join(study_site_env, by = c("Campaign", "site")) %>%
    mutate(kelp.cover = ifelse(kelp.cover > 100, 100, kelp.cover)) %>%
    dplyr::select(lon, lat, {{kelp_choice}}) %>% 
    group_by(lon, lat) %>%
    summarise_all(mean) %>%
    ungroup()
  # colnames(data_point)[3] <- "mean_cover"
  
  # Create plot
  Arctic_map +
    # geom_tile(data = df, # No filter
    geom_tile(data = filter(df, depth <= 100 | land_distance <= 100),
              aes(x = lon, y = lat, fill = pred_val)) +
    geom_point(data = data_point, colour = "red", shape = 21,
               aes_string(x = "lon", y = "lat", size = kelp_choice)) +
    scale_fill_viridis_c(legend_title) +
    theme(strip.background = element_rect(colour = "white", fill = "white"),
          legend.position = c(x_nudge, 0.96),
          legend.direction = "horizontal",
          legend.spacing.y = unit(0, "mm")) +
    labs(size = legend_title)
}

# Visualisations
# cover_squiz(pred_kelpcover, "Total cover (%)", 0.785, "kelp.cover")
# cover_squiz(pred_laminariales, "Laminariales cover (%)", 0.745, "Laminariales")
# cover_squiz(pred_agarum, "Agarum cover (%)", 0.77, "Agarum")
# cover_squiz(pred_agarum_2050, "Agarum cover (%)", 0.77, "Agarum")
# cover_squiz(pred_agarum_2100, "Agarum cover (%)", 0.77, "Agarum")
# cover_squiz(pred_alaria, "Alaria cover (%)", 0.78, "Alaria")

