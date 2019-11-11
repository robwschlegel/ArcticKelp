# analyses/7_random_forests.R
# The purpose of this script is to house the code used for the random forest analyses
# There is an analysis in script 5 and 6A as well, but my (RWS) thinking was to put everything
# in one place moving forward in the sake of tidiness.
# I (RWS) also found an additional method of performing a random forest in R that I
# will test out here as well


# Libraries ---------------------------------------------------------------

# Load all data nad previous libraries
source("analyses/4_kelp_cover.R")

# Libraries for this script specifically
library(randomForest)
library(forestRK)


# Data --------------------------------------------------------------------

# Create wide data.frame with all of the abiotic variables matched to depth of site
# Also includes substrate percentage variables
kelp_all <- adf %>% 
  dplyr::select(Campaign, site, depth, Quadrat, Quadrat.size.m2, Bedrock..:sand, 
                kelp.cover, Laminariales, Agarum, Alaria) %>% 
  left_join(study_site_means, by = c("Campaign", "site")) %>% 
  dplyr::select(Campaign:Alaria, depth.y, eken:icethic_cat, lon, lat) %>%
  mutate(eken = ifelse(depth.x == depth.y, eken, NA),  # A funny way of getting rid of non-target depth data
         soce = ifelse(depth.x == depth.y, soce, NA), 
         toce = ifelse(depth.x == depth.y, toce, NA)) %>% 
  dplyr::select(-depth.y) %>% 
  dplyr::rename(depth = depth.x) %>% 
  gather(key = "model_var", value = "val", -c(Campaign:Alaria, lon, lat)) %>%
  na.omit() %>% 
  spread(model_var, val)

# Filter down to only total kelp cover at depths 5, 10 and 15 m
kelp_var <- kelp_all %>% 
  dplyr::select(-Laminariales, -Agarum, -Alaria) #%>% 
#  filter(depth %in% c(10, 15)) #%>% 
 # mutate(kelp.cover = round(kelp.cover, -1)) # Round kelp cover to the nearest 10% step


# Random Forest -----------------------------------------------------------

# The data.frame that will be fed to the model
# The Quadrat information is fairly random, as it should be, and so isn't used in the model TRUE
# The substrate and depth values also score consistently in the bottom of importance so aren't used
data1 <- select(kelp_var, kelp.cover:toce, -lon) #%>% 
  #mutate(kelp.cover = as.factor(kelp.cover)) #%>% 
  # mutate(kelp.cover = factor(kelp.cover, levels = as.character(seq(0, 100, by = 10))))

#what variables are highly correlated. 
round(cor(data1),2)


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
lapply(1:9, rf_mtry_test) # It looks like an mtry >= 2 is best

# Random forest model based on all quadrat data
kelp_rf <- randomForest(kelp.cover ~ ., data = train_set, mtry = 2, ntree = 1000,
                        importance = TRUE, na.action = na.omit)
summary(kelp_rf)
print(kelp_rf)
#explains 67% pf variance

round(importance(kelp_rf), 1)
varImpPlot(kelp_rf)
partialPlot(kelp_rf, train_set, mldr10_1)
partialPlot(kelp_rf, train_set, iceconc_cat  )


# Predicting on training set
pred_train <- predict(kelp_rf, train_set) # type = 'response' also works but returns the same results
table(pred_train, train_set$kelp.cover)  
mean(abs(pred_train - train_set$kelp.cover)) # 43% accuracy...

# Predicting on Validation set
pred_valid <- predict(kelp_rf, valid_set)
table(pred_valid, valid_set$kelp.cover)  
mean(abs(pred_valid - valid_set$kelp.cover)) # 29% accuracy...


#Agarum
data2 <- select(kelp_all, Agarum:toce, -lon, -Alaria)
train <- sample(nrow(data2), 0.7*nrow(data2), replace = FALSE)
train_set <- data2[train,]
valid_set <- data2[-train,]
kelp_rf <- randomForest(Agarum ~ ., data = train_set, mtry = 3, ntree = 1000,
                        importance = TRUE, na.action = na.omit)
summary(kelp_rf)
print(kelp_rf)
round(importance(kelp_rf), 1)
varImpPlot(kelp_rf)
partialPlot(kelp_rf, train_set, iceconc_cat)
partialPlot(kelp_rf, train_set, lat)

# Predicting on training set
pred_train <- predict(kelp_rf, train_set) # type = 'response' also works but returns the same results
table(pred_train, train_set$Laminariales)  
mean(abs(pred_train - train_set$Laminariales)) # 43% accuracy...

# Predicting on Validation set
pred_valid <- predict(kelp_rf, valid_set)
table(pred_valid, valid_set$Laminariales)  
mean(abs(pred_train - valid_set$Laminariales)) # 29% accuracy...



#laminariales
data2 <- select(kelp_all, Laminariales:toce, -lon, -lat, -Alaria, -Agarum)
train <- sample(nrow(data2), 0.7*nrow(data2), replace = FALSE)
train_set <- data2[train,]
valid_set <- data2[-train,]
kelp_rf <- randomForest(Laminariales ~ ., data = train_set, mtry = 3, ntree = 1000,
                        importance = TRUE, na.action = na.omit)
summary(kelp_rf)
print(kelp_rf)
round(importance(kelp_rf), 1)
varImpPlot(kelp_rf)


# Predicting on training set
pred_train <- predict(kelp_rf, train_set) # type = 'response' also works but returns the same results
table(pred_train, train_set$Laminariales)  
mean(abs(pred_train - train_set$Laminariales)) # 43% accuracy...

# Predicting on Validation set
pred_valid <- predict(kelp_rf, valid_set)
table(pred_valid, valid_set$Laminariales)  
mean(abs(pred_train - valid_set$Laminariales)) # 29% accuracy...


# forestRK method ---------------------------------------------------------
### I do not think we want to go categorical here. It makes no sense with the data

## Prep data

# Step 1: Remove all NA's and NaN's from the original dataset
# The method prefers categorical/factor response variables
data2 <- na.omit(data1) %>%  # This is the same data.frame a the previous model method
  mutate(kelp.cover = as.factor(kelp.cover))# %>% 
  # mutate(kelp.cover = factor(kelp.cover, levels = as.character(seq(0, 100, by = 10))))
# levels(data2$kelp.cover)

# Step 2: Seperate the dataset into a chunk that stores covariates of both 
# training and test observations X and a vector y that stores class types of the 
# training observations
vec <- seq.int(1, nrow(data2), by = 3) # indices of the training observations
X <- data2[ ,2:length(data2)]
y <- data2[vec, 1]

# Step 3: Numericize the data frame X by applying the x.organizer function
# and split X into a training and a test set
X1 <- x.organizer(X, encoding = "num") # Numeric Encoding applied

## train and test set from the Numeric ncoding
x.train1 <- X1[vec,]
x.test1 <- X1[-vec,]

# Step 4: Numericize the vector y by applying the y.organizer function
y.train <- y.organizer(y)$y.new # a vector storing the numericized class type
y.factor.levels <- y.organizer(y)$y.factor.levels

## Run model

set.seed(666)

# Construct rkTree based on gini splitting
tree.gini <- construct.treeRK(x.train1, y.train, min.num.obs.end.node.tree = 2, 
                              entropy = FALSE)

# Construct rkTree based on entropy splitting
  ## default for entropy is TRUE
  ## default for min.num.obs.end.node.tree is 5
tree.entropy <- construct.treeRK(x.train1, y.train)

## Visualise model

# plotting the tree.gini rkTree obtained via construct.treeRK function
draw.treeRK(tree.gini, y.factor.levels[c(-1)], font = "Times", node.colour = "white", 
            text.colour = "dark blue", text.size = 0.6, tree.vertex.size = 30, 
            tree.title = "Decision Tree", title.colour = "dark green")
# ggsave("graph/forestRK_tree.png")

## Model predictions

# display numericized predicted class type for the first 5 observations (in the 
# order of increasing original observation index number) from the test set with
# Binary Encoding
prediction.df <- pred.treeRK(X = x.test1, tree.gini)$prediction.df
# prediction.df[1:5, dim(prediction.df)[2]]

# display the list of hierarchical flag from the test set with Numeric Encoding
pred.treeRK(X = x.test1, tree.entropy)$flag.pred

## More models

# Make a forestRK model based on the training data with Numeric Encoding and with Gini Index as the splitting criteria
# normally nbags and samp.size would be much larger than 30 and 50
forestRK.1 <- forestRK(x.train1, y.train, nbags = 30, samp.size = 50, entropy = FALSE)

# Make a forestRK model based on the training data with Numeric Encoding and with Entropy as the splitting criteria
# normally nbags and samp.size would be much larger than 30 and 50
forestRK.2 <- forestRK(x.train1, y.train, nbags = 30, samp.size = 50)

## Predictions

## make predictions on the test set x.test1 based on the forestRK model constructed from
# x.train1 and y.train (recall: these are the training data with Binary Encoding)
# after using Gini Index as splitting criteria (entropy == FALSE)
pred.forest.rk1 <- pred.forestRK(x.test = x.test1, x.training = x.train1, 
                                 y.training = y.train, nbags = 100, 
                                 samp.size = 100, y.factor.levels = y.factor.levels,
                                 entropy = FALSE)

# Make predictions on the test set x.test2 based on the forestRK model constructed from
# x.train2 and y.train (recall: these are the training data with Numeric Encoding)
# after using Entropy as splitting criteria (entropy == TRUE)
pred.forest.rk2 <- pred.forestRK(x.test = x.test1, x.training = x.train2, 
                                 y.training = y.train, nbags = 100, 
                                 samp.size = 100, y.factor.levels = y.factor.levels, 
                                 entropy = TRUE)

## Extract individual trees

# Get tree
tree.index.ex <- c(1,3,8)

# Extract information about the 1st, 3rd, and 8th tree in the forestRK.1
get.tree <- get.tree.forestRK(forestRK.1, tree.index = tree.index.ex)

# Display 8th tree in the forest (which is the third element of the vector 'tree.index.ex')
get.tree[["8"]]

# Get the list of variable used for splitting
covariate.used.for.split.tree <- var.used.forestRK(forestRK.1, tree.index = c(4,5,6))

# Get the list of names of covariates used for splitting to construct tree #6 in the forestRK.1
covariate.used.for.split.tree[["6"]]

## MDS plot 

# Generate 2-dimensional MDS plot of the test observations colour coded by the 
# predictions stored under the object pred.forest.rk1
mds.plot.forestRK(pred.forest.rk1, plot.title = "Kelp cover predictions", xlab ="1st Coordinate", 
                  ylab = "2nd Coordinate", colour.lab = "Predictions By The forestRK.1")


## Importance of variables

imp <- importance.forestRK(forestRK.1)

# Gives the list of covariate names ordered from most important to least important
covariate.names <- imp$importance.covariate.names 

# Gives the list of average decrease in (weighted) node impurity across all trees in the forest
# ordered from the highest to lowest in value.
# That is, the first element of this vector pertains to the first covariate listed under
# imp$importance.covariate.names 
ave.decrease.criteria <- imp$average.decrease.in.criteria.vec

## Generate importance plot of the forestRK.1 object
p <- importance.plot.forestRK(imp, colour.used = "dark green", fill.colour = "dark green", label.size = 8)
p
ggsave("graph/forestRK_variable_importance.png", height = 5, width = 6)

## Comparing performance

# Overall prediction accuracy (in percentage) when training data was modified with numeric Encoding
# and the Entropy was used as splitting criteria.
y.test <- data2[-vec, 1] # this is an actual class type of test observations
mean(pred.forest.rk1$pred.for.obs.forest.rk == as.vector(y.test))*100 # 6% accuracy...


# A different approach ----------------------------------------------------

# The results from the forestRK method above are worse than the classic random forest

# Going through everything I am starting to suspect that the difference in sampling 
# methodology between the campaigns may be an important factor
# So in this third method we will use the standard random forest method and see
# how well we can predict which campaign a data point comes from
# We'll do this both with and without lon lat

# data3 <- select(kelp_var, Campaign, lon:toce) # With lon/lat
data3 <- select(kelp_var, Campaign, eken:toce) # Without lon/lat
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

