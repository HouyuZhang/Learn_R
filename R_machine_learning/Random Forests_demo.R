#http://uc-r.github.io/random_forests
library(rsample)      # data splitting 
library(randomForest) # basic implementation
library(ranger)       # a faster implementation of randomForest
library(caret)        # an aggregator package for performing many machine learning models
library(h2o)          # an extremely fast java-based platform

#=========================================================================================
#0. Data arrangement
#=========================================================================================
# Create training (70%) and test (30%) sets for the AmesHousing::make_ames() data.
set.seed(123)
ames_split <- initial_split(AmesHousing::make_ames(), prop = .7)
ames_train <- training(ames_split)
ames_test  <- testing(ames_split)

#=========================================================================================
#1. Trainning using randomForest
#=========================================================================================
set.seed(123)
# default RF model
m1 <- randomForest(
  formula = Sale_Price ~ .,
  data    = ames_train
)

m1
plot(m1)

# number of trees with lowest MSE
which.min(m1$mse)

# RMSE of this optimal random forest
sqrt(m1$mse[which.min(m1$mse)])

# create training and validation data 
set.seed(123)
valid_split <- initial_split(ames_train, .8)

# training data
ames_train_v2 <- analysis(valid_split)

# validation data
ames_valid <- assessment(valid_split)
x_test <- ames_valid[setdiff(names(ames_valid), "Sale_Price")]
y_test <- ames_valid$Sale_Price

rf_oob_comp <- randomForest(
  formula = Sale_Price ~ .,
  data    = ames_train_v2,
  xtest   = x_test,
  ytest   = y_test
)

# extract OOB & validation errors
oob <- sqrt(rf_oob_comp$mse)
validation <- sqrt(rf_oob_comp$test$mse)

# compare error rates
tibble::tibble(
  `Out of Bag Error` = oob,
  `Test error` = validation,
  ntrees = 1:rf_oob_comp$ntree
) %>%
  gather(Metric, RMSE, -ntrees) %>%
  ggplot(aes(ntrees, RMSE, color = Metric)) +
  geom_line() +
  scale_y_continuous(labels = scales::dollar) +
  xlab("Number of trees")

# names of features
features <- setdiff(names(ames_train), "Sale_Price")

set.seed(123)
m2 <- tuneRF(
  x          = ames_train[features],
  y          = ames_train$Sale_Price,
  ntreeTry   = 500,
  mtryStart  = 5,
  stepFactor = 1.5,
  improve    = 0.01,
  trace      = FALSE      # to not show real-time progress 
)

# create feature names
y <- "Sale_Price"
x <- setdiff(names(ames_train), y)

# turn training set into h2o object
train.h2o <- as.h2o(ames_train)

# hyperparameter grid
hyper_grid.h2o <- list(
  ntrees      = seq(200, 500, by = 100),
  mtries      = seq(20, 30, by = 2),
  sample_rate = c(.55, .632, .70, .80)
)

# build grid search 
grid <- h2o.grid(
  algorithm = "randomForest",
  grid_id = "rf_grid",
  x = x, 
  y = y, 
  training_frame = train.h2o,
  hyper_params = hyper_grid.h2o,
  search_criteria = list(strategy = "Cartesian")
)

# collect the results and sort by our model performance metric of choice
grid_perf <- h2o.getGrid(
  grid_id = "rf_grid", 
  sort_by = "mse", 
  decreasing = FALSE
)
print(grid_perf)
























