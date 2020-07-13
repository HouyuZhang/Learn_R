# http://uc-r.github.io/iml-pkg
library(rsample)   # data splitting
library(ggplot2)   # allows extension of visualizations
library(dplyr)     # basic data transformation
library(h2o)       # machine learning modeling
library(iml)       # ML interprtation

# initialize h2o session
h2o.no_progress()
h2o.init()

#=========================================================================================
#0. Data arrangement
#=========================================================================================
# classification data
df <- rsample::attrition %>% 
  mutate_if(is.ordered, factor, ordered = FALSE) %>%
  mutate(Attrition = recode(Attrition, "Yes" = "1", "No" = "0") %>% factor(levels = c("1", "0")))

# convert to h2o object
df.h2o <- as.h2o(df)

# create train, validation, and test splits
set.seed(123)
splits <- h2o.splitFrame(df.h2o, ratios = c(.7, .15), destination_frames = c("train","valid","test"))
names(splits) <- c("train","valid","test")

# variable names for resonse & features
y <- "Attrition"
x <- setdiff(names(df), y)

#=========================================================================================
#1. Starting trainning
#=========================================================================================
###1.1 supported models ----
# elastic net model 
glm <- h2o.glm(
  x = x, 
  y = y, 
  training_frame = splits$train,
  validation_frame = splits$valid,
  family = "binomial",
  seed = 123
)

# random forest model
rf <- h2o.randomForest(
  x = x, 
  y = y,
  training_frame = splits$train,
  validation_frame = splits$valid,
  ntrees = 1000,
  stopping_metric = "AUC",    
  stopping_rounds = 10,         
  stopping_tolerance = 0.005,
  seed = 123
)

# gradient boosting machine model
gbm <-  h2o.gbm(
  x = x, 
  y = y,
  training_frame = splits$train,
  validation_frame = splits$valid,
  ntrees = 1000,
  stopping_metric = "AUC",    
  stopping_rounds = 10,         
  stopping_tolerance = 0.005,
  seed = 123
)

# model performance
h2o.auc(glm, valid = TRUE)
h2o.auc(rf, valid = TRUE)
h2o.auc(gbm, valid = TRUE)

###1.2 unsupported models ----
# 1). create a data frame with just the features
features <- as.data.frame(splits$valid) %>% select(-Attrition)

# 2). Create a vector with the actual responses
response <- as.numeric(as.vector(splits$valid$Attrition))

# 3). Create custom predict function that returns the predicted values as a vector (probability of purchasing in our example)
pred <- function(model, newdata)  {
  results <- as.data.frame(h2o.predict(model, as.h2o(newdata)))
  return(results[[3L]])
}

# example of prediction output
pred(rf, features) %>% head()


# create predictor object to pass to explainer functions
predictor.glm <- Predictor$new(
  model = glm, 
  data = features, 
  y = response, 
  predict.fun = pred,
  class = "classification"
)

predictor.rf <- Predictor$new(
  model = rf, 
  data = features, 
  y = response, 
  predict.fun = pred,
  class = "classification"
)

predictor.gbm <- Predictor$new(
  model = gbm, 
  data = features, 
  y = response, 
  predict.fun = pred,
  class = "classification"
)

# structure of predictor object
str(predictor.gbm)

#=========================================================================================
#2. Global interpretation
#=========================================================================================
### 2.1 Feature importance ----
# compute feature importance with specified loss metric
imp.glm <- FeatureImp$new(predictor.glm, loss = "mse")
imp.rf <- FeatureImp$new(predictor.rf, loss = "mse")
imp.gbm <- FeatureImp$new(predictor.gbm, loss = "mse")

# plot output
p1 <- plot(imp.glm) + ggtitle("GLM")
p2 <- plot(imp.rf) + ggtitle("RF")
p3 <- plot(imp.gbm) + ggtitle("GBM")

gridExtra::grid.arrange(p1, p2, p3, nrow = 1)

### 2.2 Partial dependence (PDPs) ----
glm.ot <- Partial$new(predictor.glm, "OverTime") %>% plot() + ggtitle("GLM")
rf.ot <- Partial$new(predictor.rf, "OverTime") %>% plot() + ggtitle("RF") 
gbm.ot <- Partial$new(predictor.gbm, "OverTime") %>% plot() + ggtitle("GBM")

gridExtra::grid.arrange(glm.ot, rf.ot, gbm.ot, nrow = 1)

# GLM model
glm.age <- Partial$new(predictor.glm, "Age", ice = TRUE, grid.size = 50)
glm.age$center(min(features$Age))
p1 <- plot(glm.age) + ggtitle("GLM")

# RF model
rf.age <- Partial$new(predictor.rf, "Age", ice = TRUE, grid.size = 50)
rf.age$center(min(features$Age))
p2 <- plot(rf.age) + ggtitle("RF")

# GBM model
gbm.age <- Partial$new(predictor.gbm, "Age", ice = TRUE, grid.size = 50)
gbm.age$center(min(features$Age))
p3 <- plot(gbm.age) + ggtitle("GBM")

gridExtra::grid.arrange(p1, p2, p3, nrow = 1)

p1 <- Partial$new(predictor.glm, c("MonthlyIncome", "OverTime")) %>% plot() + ggtitle("GLM") + ylim(c(0, .4))
p2 <- Partial$new(predictor.rf, c("MonthlyIncome", "OverTime")) %>% plot() + ggtitle("RF") + ylim(c(0, .4))
p3 <- Partial$new(predictor.gbm, c("MonthlyIncome", "OverTime")) %>% plot() + ggtitle("GBM") + ylim(c(0, .4))

gridExtra::grid.arrange(p1, p2, p3, nrow = 1)

### 2.3 Measuring interactions (H-statistic) ----
# identify variables with largest interactions in each model
interact.glm <- Interaction$new(predictor.glm) %>% plot() + ggtitle("GLM")
interact.rf  <- Interaction$new(predictor.rf) %>% plot() + ggtitle("RF")
interact.gbm <- Interaction$new(predictor.gbm) %>% plot() + ggtitle("GBM")

gridExtra::grid.arrange(interact.glm, interact.rf, interact.gbm, nrow = 1)

# identify variables with largest interactions w/OverTime
interact.glm <- Interaction$new(predictor.glm, feature = "OverTime") %>% plot()
interact.rf  <- Interaction$new(predictor.rf, feature = "OverTime") %>% plot()
interact.gbm <- Interaction$new(predictor.gbm, feature = "OverTime") %>% plot()

gridExtra::grid.arrange(interact.glm, interact.rf, interact.gbm, nrow = 1)

### 2.4 Surrogate model ----
# fit surrogate decision tree model
tree <- TreeSurrogate$new(predictor.gbm, maxdepth = 3)

# how well does this model fit the original results
tree$r.squared
plot(tree)

#=========================================================================================
#3. Local interpretation
#=========================================================================================
# identify obs with highest and lowest probabilities
(high <- predict(rf, splits$valid) %>% .[, 3] %>% as.vector() %>% which.max())
(low  <- predict(rf, splits$valid) %>% .[, 3] %>% as.vector() %>% which.min())  

# get these observations
high_prob_ob <- features[high, ]
low_prob_ob  <- features[low, ]

# fit local model
lime.glm <- LocalModel$new(predictor.glm, k = 10, x.interest = high_prob_ob) %>% plot() + ggtitle("GLM")
lime.rf  <- LocalModel$new(predictor.rf, k = 10, x.interest = high_prob_ob) %>% plot() + ggtitle("RF")
lime.gbm <- LocalModel$new(predictor.gbm, k = 10, x.interest = high_prob_ob) %>% plot() + ggtitle("GBM")

gridExtra::grid.arrange(lime.glm, lime.rf, lime.gbm, nrow = 1)


# fit local model
lime.glm <- LocalModel$new(predictor.glm, k = 10, x.interest = low_prob_ob) %>% plot() + ggtitle("GLM")
lime.rf  <- LocalModel$new(predictor.rf, k = 10, x.interest = low_prob_ob) %>% plot() + ggtitle("RF")
lime.gbm <- LocalModel$new(predictor.gbm, k = 10, x.interest = low_prob_ob) %>% plot() + ggtitle("GBM")

gridExtra::grid.arrange(lime.glm, lime.rf, lime.gbm, nrow = 1)

# high probability observation
predict(rf, splits$valid) %>% .[high, 3] # actual probability

lime_high <- LocalModel$new(predictor.rf, k = 10, x.interest = high_prob_ob)
lime_high$predict(high_prob_ob) # predicted probability

# low probability observation
predict(rf, splits$valid) %>% .[low, 3] # actual probability

lime_low <- LocalModel$new(predictor.rf, k = 10, x.interest = low_prob_ob)
lime_low$predict(low_prob_ob) # predicted probability

### 3.2 Shapley values ----
shapley.glm <- Shapley$new(predictor.glm, x.interest = high_prob_ob) %>% plot() + ggtitle("GLM")
shapley.rf  <- Shapley$new(predictor.rf, x.interest = high_prob_ob) %>% plot() + ggtitle("RF")
shapley.gbm <- Shapley$new(predictor.gbm, x.interest = high_prob_ob) %>% plot() + ggtitle("GBM")

gridExtra::grid.arrange(shapley.glm, shapley.rf, shapley.gbm, nrow = 1)


shapley.glm <- Shapley$new(predictor.glm, x.interest = low_prob_ob) %>% plot() + ggtitle("GLM")
shapley.rf  <- Shapley$new(predictor.rf, x.interest = low_prob_ob) %>% plot() + ggtitle("RF")
shapley.gbm <- Shapley$new(predictor.gbm, x.interest = low_prob_ob) %>% plot() + ggtitle("GBM")

gridExtra::grid.arrange(shapley.glm, shapley.rf, shapley.gbm, nrow = 1)



