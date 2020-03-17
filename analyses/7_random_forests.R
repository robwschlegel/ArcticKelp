# analyses/7_random_forests.R
# The purpose of this script is to house the code used for the random forest analyses


# Setup -------------------------------------------------------------------

# Load all data nad previous libraries
source("analyses/4_kelp_cover.R")

# Libraries for this script specifically
library(randomForest)
library(OneR) # For single rule machine learning
library(tidymodels) # For tdy modelling conventions
library(caret) # More modelling code
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
  dplyr::select(-BO_ph) %>% # There are large errors in pH throughout the entire study region
  dplyr::select(-BO2_silicateltmin_bdmax, 
                -BO2_silicatemean_bdmax, 
                -BO2_silicateltmax_bdmax) %>% # There are large errors in silicate in Hudson Bay region
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
base_sub <- colnames(dplyr::select(kelp_all, 
                                   BO2_templtmin_bdmax:BO2_salinityltmax_ss,
                                   BO2_icethickltmin_ss:BO2_icethickltmax_ss, 
                                   BO2_curvelltmin_bdmax:BO2_curvelltmax_bdmax))
future_2050 <- colnames(dplyr::select(kelp_all, 
                                      BO2_icecoverltmin_ss:BO2_icecoverltmax_ss,
                                      BO_parmean:BO2_phosphateltmax_bdmax,
                                      BO2_RCP85_2050_curvelltmax_bdmax:BO2_RCP85_2050_tempmean_ss))
future_2100 <- colnames(dplyr::select(kelp_all, 
                                      BO2_icecoverltmin_ss:BO2_icecoverltmax_ss,
                                      BO_parmean:BO2_phosphateltmax_bdmax,
                                      BO2_RCP85_2100_curvelltmax_bdmax:BO2_RCP85_2100_tempmean_ss))


# Data prep function ------------------------------------------------------

# Convenience function for prepping dataframe for use in the random forest
# This removes all other kelp cover values
rf_data_prep <- function(kelp_choice, df = kelp_all, cut_cover = F, scenario = base){
  
  # Trim down data.frame
  df_1 <- data.frame(dplyr::select(df, -c(Campaign:site))) %>% 
    pivot_longer(cols = kelp.cover:Alaria, names_to = "chosen_kelp", values_to = "cover") %>% 
    filter(chosen_kelp == kelp_choice) %>% 
    dplyr::select(scenario, cover)
  
  # Cut cover into categories if desired
  if(cut_cover){
    df_1$cover <- cut(df_1$cover, breaks = c(-Inf, 0, 10, 50, 100), ordered_result = T)
    levels(df_1$cover)[1] <- "(0]"
  }
  
  return(df_1)
}


# Single rule model -------------------------------------------------------

# Function for OneR model run
OneR_model <- function(kelp_choice, df = kelp_all){

  # Pull out training index
  one_random <- sample(1:nrow(df), 0.7*nrow(df))
  
  # Prep the data
  kelp_cut <- rf_data_prep(kelp_choice = kelp_choice, df = df, cut_cover = T)
           
  # The training data
  kelp_train <- kelp_cut %>% 
    slice(one_random) %>%
    data.frame() %>% 
    optbin(., method = "infogain") # The choice of method does have an effect on the results
  
  # The test data
  kelp_test <- kelp_cut %>% 
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
#

# Establish RF conventions ------------------------------------------------

# Prep data
df_kelp.cover <- rf_data_prep("kelp.cover")

# Define training control: 10 fold cross-validation. 
train_control <- trainControl(method = "cv", number = 10)

# Train the model using randomForest (rf).
model <- train(cover ~ ., data = df_kelp.cover, trControl = train_control, 
               method = "rf", importance = T)
print(model)

plot(model$finalModel)

# Make predictions and produce Mean Absolute Error (MAE) scores to evaluate
# model performance
predictions <- predict(model, df_kelp.cover)
result <- data.frame(Actual = df_kelp.cover$cover, Predicted = predictions)
result$Difference <- abs(result$Actual - result$Predicted)
summary(result$Difference) 

# Plot predicted vs observed abundance to inspect model performance
df_kelp.cover$predictions <- predictions

predictions_fig <- ggplot(df_kelp.cover, aes(predictions, cover)) +
  geom_point(colour = "black") +
  ggtitle(bquote(~italic("All macroalgae"))) +
  geom_abline(slope = 1, intercept = 0) +
  scale_x_continuous(limits = c(0, 100)) +
  labs(x = "Predicted abundance (% cover)", 
       y = "Observed abundance (% cover)") +
  theme_bw() +
  theme(aspect.ratio = 1) +
  theme(plot.title = element_text(size = 20, vjust = 0),
        axis.text.x = element_text(angle = 50, size = 10, vjust = 1, hjust = 1),
        axis.text.y = element_text(size = 10))
predictions_fig

# View variables by relative importance in the model to see what factors 
# may be most ecologically important for this species.
importance(model$finalModel)
varImp(model)
plot(varImp(model))


# Which variables are the most important? ---------------------------------

# Covenience wrapper for extracting variable importance from an RF
extract_var_imp <- function(kelp_rf){
  if(ncol(kelp_rf$importance) == 6){
    res <- data.frame(var = row.names(kelp_rf$importance), 
                      kelp_rf$importance[,5:6]) %>% 
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
# are consistantly the most important
top_var <- function(lplyr_bit, kelp_choice, df = kelp_all){
  
  # Prep the four possibilities for a single kelp cover choice
  df_reg_base <- rf_data_prep(kelp_choice, df = df, cut_cover = F, scenario = base)
  df_cut_base <- rf_data_prep(kelp_choice, df = df, cut_cover = T, scenario = base)
  df_reg_sub <- rf_data_prep(kelp_choice, df = df, cut_cover = F, scenario = base_sub)
  df_cut_sub <- rf_data_prep(kelp_choice, df = df, cut_cover = T, scenario = base_sub)
  
  # Random sampling to split data up for training
  train <- sample(1:nrow(df_reg_base), 0.7*nrow(df_reg_base), replace = FALSE)
  
  # Random forest models for the four posibilities
  rf_reg_base <- randomForest(cover ~ ., data = df_reg_base[train,], ntree = 200, importance = TRUE, do.trace = F)
  rf_cut_base <- randomForest(cover ~ ., data = df_cut_base[train,], ntree = 200, importance = TRUE, do.trace = F)
  rf_reg_sub <- randomForest(cover ~ ., data = df_reg_sub[train,], ntree = 200, importance = TRUE, do.trace = F)
  rf_cut_sub <- randomForest(cover ~ ., data = df_cut_sub[train,], ntree = 200, importance = TRUE, do.trace = F)

  # Extract results
  reg_base <- extract_var_imp(rf_reg_base)
  cut_base <-  extract_var_imp(rf_cut_base)
  reg_sub <- extract_var_imp(rf_reg_sub)
  cut_sub <- extract_var_imp(rf_cut_sub)
  
  # Combine and exit
  res <- left_join(reg_base, cut_base, by = "var") %>% 
    left_join(reg_sub, by = "var") %>% 
    left_join(cut_sub, by = "var")
  return(res)
}
# top_variables(kelp_choice = "Agarum", df = kelp_all)

# Convenience wrapper to remove corrrelated variables
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
  res <- list(reg_base = cor_var_rm(multi_kelp_mean[,1:6]),
              cut_base = cor_var_rm(multi_kelp_mean[,c(1,7:8)]),
              reg_sub = cor_var_rm(multi_kelp_mean[,c(1,9:13)]),
              cut_sub = cor_var_rm(multi_kelp_mean[,c(1,14:15)]))
  return(res)
}

## Find the top variables for the different kelp covers
# kelp.cover
system.time(top_var_kelpcover <- top_var_multi("kelp.cover")) # ~5 seconds on 50 cores
save(top_var_kelpcover, file = "data/top_var_kelpcover.RData")

# Laminariales
top_var_laminariales <- top_var_multi("Laminariales")
save(top_var_laminariales, file = "data/top_var_laminariales.RData")

# Agarum
top_var_agarum <- top_var_multi("Agarum")
save(top_var_agarum, file = "data/top_var_agarum.RData")

# Alaria
top_var_alaria <- top_var_multi("Alaria")
save(top_var_alaria, file = "data/top_var_alaria.RData")


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
random_kelp_forest_check <- function(kelp_choice, column_choice){
  
  # Extract only the kelp cover of choice
  df_kelp_reg <- rf_data_prep(kelp_choice)
  df_kelp_cat <- rf_data_prep(kelp_choice, cut_cover = T)
  
  # Chose only desired columns
  df_var_reg_base <- dplyr::select(df_kelp_reg, cover, as.character(column_choice$reg_base$var))
  df_var_reg_sub <- dplyr::select(df_kelp_reg, cover, as.character(column_choice$reg_sub$var))
  df_var_cut_base <- dplyr::select(df_kelp_cat, cover, as.character(column_choice$cut_base$var)) 
  df_var_cut_sub <- dplyr::select(df_kelp_cat, cover, as.character(column_choice$cut_sub$var)) 
  
  # Split data up for training and testing
  train <- sample(nrow(df_kelp_reg), 0.7*nrow(df_kelp_reg), replace = FALSE)
  
  # Random forest model regression
  rf_reg_base <- randomForest(cover ~ ., data = df_var_reg_base[train, ], 
                              ntree = 200, importance = TRUE, do.trace = F)
  # print(rf_reg_base) # percent of variance explained
  # varImpPlot(rf_reg_base)
  rf_reg_sub <- randomForest(cover ~ ., data = df_var_reg_sub[train, ], 
                             ntree = 200, importance = TRUE, do.trace = F)
  
  # Random forest model category
  rf_cut_base <- randomForest(cover ~ ., data = df_var_cut_base[train, ],
                              ntree = 200, importance = TRUE, do.trace = F)
  # print(rf_cut_base) # percent of variance explained
  # varImpPlot(rf_cut_base)
  rf_cut_sub <- randomForest(cover ~ ., data = df_var_cut_sub[train, ],
                             ntree = 200, importance = TRUE, do.trace = F)
  
  # Predicting on training set
  pred_reg_base_test <- predict(rf_reg_base, df_var_reg_base[train, ])
  print(rf_cut_sub)
  print(paste0("Average inaccuracy of prediction on base test data: ", 
               round(mean(abs(pred_reg_base_test - df_var_reg_base[train, ]$cover)), 2),"%"))
  pred_reg_sub_test <- predict(rf_reg_sub, df_var_reg_sub[train, ])
  print(paste0("Average inaccuracy of prediction on sub test data: ", 
               round(mean(abs(pred_reg_sub_test - df_var_reg_sub[train, ]$cover)), 2),"%"))
  
  # Predicting on Validation set
  pred_reg_base_valid <- predict(rf_reg_base, df_var_reg_base[-train, ])
  print(paste0("Average inaccuracy of prediction on validation base data: ", 
               round(mean(abs(pred_reg_base_valid - df_var_reg_base[-train, ]$cover)), 2),"%"))
  pred_reg_sub_valid <- predict(rf_reg_sub, df_var_reg_sub[-train, ])
  print(paste0("Average inaccuracy of prediction on validation sub data: ", 
               round(mean(abs(pred_reg_sub_valid - df_var_reg_sub[-train, ]$cover)), 2),"%"))
  
  # Accuracy of categorical prediction
  pred_cut_base_train <- predict(rf_cut_base, df_var_cut_base[train, ])
  print(pred_cut_base_train)
  print(paste0("Average accuracy of category prediction on test data: ", 
               cbind(pred_cut_base_train, df_var_cut_base[train, ]$cover) %>% 
                 data.frame() %>% 
                 mutate(acc = ifelse(pred_cut_base_train == V2, 1, 0)) %>% 
                 summarise(acc = round(sum(acc)/n(), 4)*100),"%"))
  pred_valid_cut <- predict(rf_cut_base, df_var_cut_base[-train, ])
  print(paste0("Average accuracy of category prediction on validation data: ", 
               cbind(pred_valid_cat, valid_set$cover_cut) %>% 
                 data.frame() %>% 
                 mutate(acc = ifelse(pred_valid_cat == V2, 1, 0)) %>% 
                 summarise(acc = round(sum(acc)/n(), 4)*100),"%"))
}

# Check the random forests
random_kelp_forest_check("kelp.cover", top_var_kelpcover)
random_kelp_forest_check("Laminariales", top_var_laminariales)
random_kelp_forest_check("Agarum", top_var_agarum)
random_kelp_forest_check("Alaria", top_var_alaria)

# For grabbing a single tree
# getTree(rfobj, k=1, labelVar=FALSE)


# Many random forests -----------------------------------------------------

# This function is designed to output the model created and it's predictive accuracy
random_kelp_forest_test <- function(lplyr_bit, kelp_choice, column_choice,
                                    df = kelp_all, cut_cover = F){
  
  # Extract only the kelp cover of choice
  df_kelp_choice <- rf_data_prep(kelp_choice, df = df, cut_cover = cut_cover)
  
  # Chose only desired columns
  df_var_choice <- dplyr::select(df_kelp_choice, cover, as.character(column_choice$var))
  
  # Split data up for training and testing
  train <- sample(nrow(df_var_choice), 0.7*nrow(df_var_choice), replace = FALSE)
  train_set <- df_var_choice[train,]
  valid_set <- df_var_choice[-train,]
  
  # Random forest model based on all quadrat data
  kelp_rf <- randomForest(cover ~ ., data = train_set, ntree = 200, importance = TRUE)
  
  # Predicting on training and validation sets
  pred_train <- predict(kelp_rf, train_set)
  pred_valid <- predict(kelp_rf, valid_set)
  
  # Create data frame of accuracy results
    train_accuracy <- data.frame(portion = "train",
                                 pred = pred_train, 
                                 original = train_set$cover) %>% 
      mutate(accuracy = as.numeric(pred)-as.numeric(original))
    validate_accuracy <- data.frame(portion = "validate",
                                    pred = pred_valid, 
                                    original = valid_set$cover) %>% 
      mutate(accuracy = as.numeric(pred)-as.numeric(original))
  
  res_accuracy <- rbind(train_accuracy, validate_accuracy)
  
  # Package the model up with the accuracy results and exit
  res <- list(model = kelp_rf, accuracy = res_accuracy)
}


# Select best random forest -----------------------------------------------

# Now we run the test on each kelp cover 1000 times to see what the spread is
# in the accuracy of the random forests
# This is caused by different random splitting of test/validation sets
# as well as the many possible routes that the random forest may then take
random_kelp_forest_select <- function(kelp_choice, column_choice, 
                                      df = kelp_all, cut_cover = F){
  # system.time(
    multi_test <- plyr::llply(.data = 1:1000, .fun = random_kelp_forest_test, .parallel = T, 
                              kelp_choice = kelp_choice, column_choice = column_choice,
                              df = df, cut_cover = cut_cover)
    # ) # ~50 seconds
  
  # Extract the model accuracies
  model_accuracy <- lapply(multi_test, function(x) x$accuracy) %>% 
    do.call(rbind.data.frame, .) %>% 
    mutate(model_id = rep(1:1000, each = nrow(kelp_all)))
  
  # Find which model had the best validation scores
  accuracy_check <- model_accuracy %>% 
    filter(portion == "validate") %>% 
    group_by(model_id)
  if(cut_cover){
    accuracy_check <- accuracy_check %>% 
      mutate(count_all = n()) %>% 
      filter(accuracy == 0) %>% 
      mutate(count_0 = n(),
             mean_acc = count_0/count_all) %>% 
      ungroup() %>% 
      dplyr::select(portion, model_id, mean_acc)
  } else {
    accuracy_check <- accuracy_check %>% 
      summarise(mean_acc = mean(abs(accuracy)),
              r_acc = round(cor(x = original, y = pred), 2)) %>% 
      ungroup()
  }
  
  # Extract the best model
  if(cut_cover){
    best_model <- arrange(accuracy_check, -mean_acc)[1,]
  } else {
    best_model <- arrange(accuracy_check, -r_acc, mean_acc)[1,]
  }
  choice_model <- multi_test[[best_model$model_id]]$model
  res <- list(choice_model = choice_model,
              model_accuracy = model_accuracy)
}

# doParallel::registerDoParallel(cores = 50)
# Kelp.cover
system.time(best_rf_kelpcover <- random_kelp_forest_select("kelp.cover", top_var_kelpcover)) # 56 seconds with 50 cores
save(best_rf_kelpcover, file = "data/best_rf_kelpcover.RData", compress = T)

# Laminariales
best_rf_laminariales <- random_kelp_forest_select("Laminariales", top_var_laminariales)
save(best_rf_laminariales, file = "data/best_rf_laminariales.RData", compress = T)

# Agarum
best_rf_agarum <- random_kelp_forest_select("Agarum", top_var_agarum)
save(best_rf_agarum, file = "data/best_rf_agarum.RData", compress = T)
# best_rf_agarum_cut <- random_kelp_forest_select("Agarum", top_var_agarum, cut_cover = T)
# save(best_rf_agarum_cut, file = "data/best_rf_agarum_cut.RData", compress = T)

# Alaria
best_rf_alaria <- random_kelp_forest_select("Alaria", top_var_alaria)
save(best_rf_alaria, file = "data/best_rf_alaria.RData", compress = T)


# Future projections ------------------------------------------------------



# Analyse model accuracy --------------------------------------------------

# First load the best random forest models produced above
load("data/best_rf_kelpcover.RData")
load("data/best_rf_laminariales.RData")
load("data/best_rf_agarum.RData")
load("data/best_rf_alaria.RData")

# Find the distributions of accuracy from 0 - 100%
test_acc <- best_rf_kelpcover$model_accuracy #%>%
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
load("data/Arctic_env.RData")

# Convenience function for final step before prediction
Arctic_cover_predict <- function(model_choice){
  pred_df <- data.frame(lon = Arctic_env$lon, lat = Arctic_env$lat,
                        depth = Arctic_env$bathy,
                        pred_val = predict(model_choice, Arctic_env))
}

# Predict the covers
pred_kelpcover <- Arctic_cover_predict(best_rf_kelpcover$choice_model)
pred_laminariales <- Arctic_cover_predict(best_rf_laminariales$choice_model)
pred_agarum <- Arctic_cover_predict(best_rf_agarum$choice_model)
pred_alaria <- Arctic_cover_predict(best_rf_alaria$choice_model)

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
    geom_tile(data = filter(df, depth <= 100),
              aes(x = lon, y = lat, fill = pred_val)) +
    geom_point(data = data_point, colour = "red", shape = 21,
               aes_string(x = "lon", y = "lat", size = kelp_choice)) +
    scale_fill_viridis_c(legend_title) +
    theme(strip.background = element_rect(colour = "white", fill = "white"),
          legend.position = c(x_nudge, 0.96),
          legend.direction = "horizontal",
          legend.spacing.y = unit(0, "mm")) +
    labs(size = paste0(kelp_choice, "(%)"))
}

# Visualisations
cover_squiz(pred_kelpcover, "Total cover (%)", 0.785, "kelp.cover")
cover_squiz(pred_laminariales, "Laminariales cover (%)", 0.745, "Laminariales")
cover_squiz(pred_agarum, "Agarum cover (%)", 0.77, "Agarum")
cover_squiz(pred_alaria, "Alaria cover (%)", 0.78, "Alaria")

