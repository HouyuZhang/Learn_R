#http://uc-r.github.io/gbm_regression

library(rsample)      # data splitting 
library(gbm)          # basic implementation
library(xgboost)      # a faster implementation of gbm
library(caret)        # an aggregator package for performing many machine learning models
library(h2o)          # a java-based platform
library(pdp)          # model visualization
library(ggplot2)      # model visualization
library(lime)         # model visualization

#=========================================================================================
#0. Data arrangement
#=========================================================================================
# Create training (70%) and test (30%) sets for the AmesHousing::make_ames() data.
# Use set.seed for reproducibility

set.seed(123)
ames_split <- initial_split(AmesHousing::make_ames(), prop = .7)
ames_train <- training(ames_split)
ames_test  <- testing(ames_split)

#=========================================================================================
#1. Trainning using gbm
#=========================================================================================
###1.1 regular trainning ----
set.seed(123)
# train GBM model
gbm.fit <- gbm(
  formula = Sale_Price ~ .,
  distribution = "gaussian",
  data = ames_train,
  n.trees = 10000,
  interaction.depth = 1,
  shrinkage = 0.001,
  cv.folds = 5,
  n.cores = NULL, # will use all cores by default
  verbose = FALSE
)  

# print results
print(gbm.fit)

# get MSE and compute RMSE
sqrt(min(gbm.fit$cv.error))

# plot loss function as a result of n trees added to the ensemble
gbm.perf(gbm.fit, method = "cv")

###1.2 Tuning parameters----
set.seed(123)
# train GBM model
gbm.fit2 <- gbm(
  formula = Sale_Price ~ .,
  distribution = "gaussian",
  data = ames_train,
  n.trees = 5000,
  interaction.depth = 3,
  shrinkage = 0.1,
  cv.folds = 5,
  n.cores = NULL, # will use all cores by default
  verbose = FALSE
)  

# find index for n trees with minimum CV error
min_MSE <- which.min(gbm.fit2$cv.error)
# get MSE and compute RMSE
sqrt(gbm.fit2$cv.error[min_MSE])
# plot loss function as a result of n trees added to the ensemble
gbm.perf(gbm.fit2, method = "cv")

###1.3 grid of hyperparameter combinations ----
# create hyperparameter grid
hyper_grid <- expand.grid(
  shrinkage = c(.01, .1, .3),
  interaction.depth = c(1, 3, 5),
  n.minobsinnode = c(5, 10, 15),
  bag.fraction = c(.65, .8, 1), 
  optimal_trees = 0,               # a place to dump results
  min_RMSE = 0                     # a place to dump results
)

# total number of combinations
nrow(hyper_grid)

# randomize data
random_index <- sample(1:nrow(ames_train), nrow(ames_train))
random_ames_train <- ames_train[random_index, ]

# grid search 
for(i in 1:nrow(hyper_grid)) {
  
  set.seed(123)
  gbm.tune <- gbm(
    formula = Sale_Price ~ .,
    distribution = "gaussian",
    data = random_ames_train,
    n.trees = 5000,
    interaction.depth = hyper_grid$interaction.depth[i],
    shrinkage = hyper_grid$shrinkage[i],
    n.minobsinnode = hyper_grid$n.minobsinnode[i],
    bag.fraction = hyper_grid$bag.fraction[i],
    train.fraction = .75,
    n.cores = NULL, # will use all cores by default
    verbose = FALSE
  )
  
  # add min training error and trees to grid
  hyper_grid$optimal_trees[i] <- which.min(gbm.tune$valid.error)
  hyper_grid$min_RMSE[i] <- sqrt(min(gbm.tune$valid.error))
}

hyper_grid %>% 
  dplyr::arrange(min_RMSE) %>%
  head(10)

###1.4 Visualization ----
#Variable importance
par(mar = c(5, 8, 1, 1))
summary(
  gbm.fit.final, 
  cBars = 10,
  method = relative.influence, # also can use permutation.test.gbm
  las = 2
)

vip::vip(gbm.fit.final)

#Partial dependence plots
gbm.fit.final %>%
  partial(pred.var = "Gr_Liv_Area", n.trees = gbm.fit.final$n.trees, grid.resolution = 100) %>%
  autoplot(rug = TRUE, train = ames_train) +
  scale_y_continuous(labels = scales::dollar)

#ICE
ice1 <- gbm.fit.final %>%
  partial(
    pred.var = "Gr_Liv_Area", 
    n.trees = gbm.fit.final$n.trees, 
    grid.resolution = 100,
    ice = TRUE
  ) %>%
  autoplot(rug = TRUE, train = ames_train, alpha = .1) +
  ggtitle("Non-centered") +
  scale_y_continuous(labels = scales::dollar)

ice2 <- gbm.fit.final %>%
  partial(
    pred.var = "Gr_Liv_Area", 
    n.trees = gbm.fit.final$n.trees, 
    grid.resolution = 100,
    ice = TRUE
  ) %>%
  autoplot(rug = TRUE, train = ames_train, alpha = .1, center = TRUE) +
  ggtitle("Centered") +
  scale_y_continuous(labels = scales::dollar)

gridExtra::grid.arrange(ice1, ice2, nrow = 1)

#LIME
model_type.gbm <- function(x, ...) {
  return("regression")
}

predict_model.gbm <- function(x, newdata, ...) {
  pred <- predict(x, newdata, n.trees = x$n.trees)
  return(as.data.frame(pred))
}

# get a few observations to perform local interpretation on
local_obs <- ames_test[1:2, ]

# apply LIME
explainer <- lime(ames_train, gbm.fit.final)
explanation <- explain(local_obs, explainer, n_features = 5)
plot_features(explanation)

###1.5 Predicting ----
# predict values for test data
pred <- predict(gbm.fit.final, n.trees = gbm.fit.final$n.trees, ames_test)

# results
caret::RMSE(pred, ames_test$Sale_Price)

#=========================================================================================
#2. Trainning using xgboost
#=========================================================================================
###2.1 regular trainning----
# variable names
features <- setdiff(names(ames_train), "Sale_Price")

# Create the treatment plan from the training data
treatplan <- vtreat::designTreatmentsZ(ames_train, features, verbose = FALSE)

# Get the "clean" variable names from the scoreFrame
new_vars <- treatplan %>%
  magrittr::use_series(scoreFrame) %>%        
  dplyr::filter(code %in% c("clean", "lev")) %>% 
  magrittr::use_series(varName)     

# Prepare the training data
features_train <- vtreat::prepare(treatplan, ames_train, varRestriction = new_vars) %>% as.matrix()
response_train <- ames_train$Sale_Price

# Prepare the test data
features_test <- vtreat::prepare(treatplan, ames_test, varRestriction = new_vars) %>% as.matrix()
response_test <- ames_test$Sale_Price

# dimensions of one-hot encoded data
dim(features_train)
dim(features_test)

set.seed(123)
xgb.fit1 <- xgb.cv(
  data = features_train,
  label = response_train,
  nrounds = 1000,
  nfold = 5,
  objective = "reg:linear",  # for regression models
  verbose = 0                # silent,
)

# get number of trees that minimize error
xgb.fit1$evaluation_log %>%
  dplyr::summarise(
    ntrees.train = which(train_rmse_mean == min(train_rmse_mean))[1],
    rmse.train   = min(train_rmse_mean),
    ntrees.test  = which(test_rmse_mean == min(test_rmse_mean))[1],
    rmse.test   = min(test_rmse_mean),
  )

# plot error vs number trees
ggplot(xgb.fit1$evaluation_log) +
  geom_line(aes(iter, train_rmse_mean), color = "red") +
  geom_line(aes(iter, test_rmse_mean), color = "blue")

###2.2 Visualization ----
###2.3 Predicting ----
#=========================================================================================
#3. Trainning using h2o
#=========================================================================================
###3.1 regular trainning ----
h2o.no_progress()
h2o.init(max_mem_size = "15g")

# create feature names
y <- "Sale_Price"
x <- setdiff(names(ames_train), y)

# turn training set into h2o object
train.h2o <- as.h2o(ames_train)

# training basic GBM model with defaults
h2o.fit1 <- h2o.gbm(
  x = x,
  y = y,
  training_frame = train.h2o,
  nfolds = 5
)

# assess model results
h2o.fit1

###3.2 Tuning parameters by Full grid search ----
# create training & validation sets
split <- h2o.splitFrame(train.h2o, ratios = 0.75)
train <- split[[1]]
valid <- split[[2]]

# create hyperparameter grid
hyper_grid <- list(
  max_depth = c(1, 3, 5),
  min_rows = c(1, 5, 10),
  learn_rate = c(0.01, 0.05, 0.1),
  learn_rate_annealing = c(.99, 1),
  sample_rate = c(.5, .75, 1),
  col_sample_rate = c(.8, .9, 1)
)

# perform grid search 
grid <- h2o.grid(
  algorithm = "gbm",
  grid_id = "gbm_grid1",
  x = x, 
  y = y, 
  training_frame = train,
  validation_frame = valid,
  hyper_params = hyper_grid,
  ntrees = 5000,
  stopping_rounds = 10,
  stopping_tolerance = 0,
  seed = 123
)

# collect the results and sort by our model performance metric of choice
grid_perf <- h2o.getGrid(
  grid_id = "gbm_grid1", 
  sort_by = "mse", 
  decreasing = FALSE
)
grid_perf

# Grab the model_id for the top model, chosen by validation error
best_model_id <- grid_perf@model_ids[[1]]
best_model <- h2o.getModel(best_model_id)

# Now let’s get performance metrics on the best model
h2o.performance(model = best_model, valid = TRUE)

###3.2 Tuning parameters by Random discrete grid search ----
# random grid search criteria
search_criteria <- list(
  strategy = "RandomDiscrete",
  stopping_metric = "mse",
  stopping_tolerance = 0.005,
  stopping_rounds = 10,
  max_runtime_secs = 60*60
)

# perform grid search 
grid <- h2o.grid(
  algorithm = "gbm",
  grid_id = "gbm_grid2",
  x = x, 
  y = y, 
  training_frame = train,
  validation_frame = valid,
  hyper_params = hyper_grid,
  search_criteria = search_criteria, # add search criteria
  ntrees = 5000,
  stopping_rounds = 10,
  stopping_tolerance = 0,
  seed = 123
)

# collect the results and sort by our model performance metric of choice
grid_perf <- h2o.getGrid(
  grid_id = "gbm_grid2", 
  sort_by = "mse", 
  decreasing = FALSE
)
grid_perf

# Grab the model_id for the top model, chosen by validation error
best_model_id <- grid_perf@model_ids[[1]]
best_model <- h2o.getModel(best_model_id)

# Now let’s get performance metrics on the best model
h2o.performance(model = best_model, valid = TRUE)

#一旦得到好的模型,开始更精确的计算
# train final model
h2o.final <- h2o.gbm(
  x = x,
  y = y,
  training_frame = train.h2o,
  nfolds = 5,
  ntrees = 10000,
  learn_rate = 0.01,
  learn_rate_annealing = 1,
  max_depth = 3,
  min_rows = 10,
  sample_rate = 0.75,
  col_sample_rate = 1,
  stopping_rounds = 10,
  stopping_tolerance = 0,
  seed = 123
)

# model stopped after xx trees
h2o.final@parameters$ntrees
# cross validated RMSE
h2o.rmse(h2o.final, xval = TRUE)

###3.4 Visualization ----
#importance
h2o.varimp_plot(h2o.final, num_of_features = 10)

#PDP and ICE
pfun <- function(object, newdata) {
  as.data.frame(predict(object, newdata = as.h2o(newdata)))[[1L]]
}

pdp <- h2o.final %>%
  partial(
    pred.var = "Gr_Liv_Area", 
    pred.fun = pfun,
    grid.resolution = 20, 
    train = ames_train
  ) %>%
  autoplot(rug = TRUE, train = ames_train, alpha = .1) +
  scale_y_continuous(labels = scales::dollar) +
  ggtitle("PDP")

ice <- h2o.final %>%
  partial(
    pred.var = "Gr_Liv_Area", 
    pred.fun = pfun,
    grid.resolution = 20, 
    train = ames_train,
    ice = TRUE
  ) %>%
  autoplot(rug = TRUE, train = ames_train, alpha = .1, center = TRUE) +
  scale_y_continuous(labels = scales::dollar) +
  ggtitle("ICE")

gridExtra::grid.arrange(pdp, ice, nrow = 1)

h2o.partialPlot(h2o.final, data = train.h2o, cols = "Overall_Qual")

pdp <- h2o.final %>%
  partial(
    pred.var = "Overall_Qual", 
    pred.fun = pfun,
    grid.resolution = 20, 
    train = as.data.frame(ames_train)
  ) %>%
  autoplot(rug = TRUE, train = ames_train, alpha = .1) +
  scale_y_continuous(labels = scales::dollar) +
  ggtitle("PDP")

ice <- h2o.final %>%
  partial(
    pred.var = "Overall_Qual", 
    pred.fun = pfun,
    grid.resolution = 20, 
    train = as.data.frame(ames_train),
    ice = TRUE
  ) %>%
  autoplot(rug = TRUE, train = ames_train, alpha = .1, center = TRUE) +
  scale_y_continuous(labels = scales::dollar) +
  ggtitle("ICE")

gridExtra::grid.arrange(pdp, ice, nrow = 1)

#LIME
explainer <- lime(ames_train, h2o.final)
explanation <- explain(local_obs, explainer, n_features = 5)
plot_features(explanation)

###3.5 Predicting ----
# convert test set to h2o object
test.h2o <- as.h2o(ames_test)

# evaluate performance on new data
h2o.performance(model = h2o.final, newdata = test.h2o)

# predict with h2o.predict
h2o.predict(h2o.final, newdata = test.h2o)

# predict values with predict
predict(h2o.final, test.h2o)







