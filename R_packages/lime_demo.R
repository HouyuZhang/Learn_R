#http://uc-r.github.io/lime
#LIME provides a great, model-agnostic approach to assessing local interpretation of predictions
library(ggplot2)    # visualization pkg leveraged by above packages
library(caret)      # ML model building
library(h2o)        # ML model building
library(vip)        # ML global interpretation
library(pdp)        # ML global interpretation
library(lime)       # ML local interpretation

# initialize h2o
h2o.init()
h2o.no_progress()
#=========================================================================================
#0. Data arrangement
#=========================================================================================
# create data sets
df <- rsample::attrition %>% 
  dplyr::mutate_if(is.ordered, factor, ordered = FALSE) %>%
  dplyr::mutate(Attrition = factor(Attrition, levels = c("Yes", "No")))

index <- 1:5
train_obs <- df[-index, ]
local_obs <- df[index, ]

# create h2o objects for modeling
y <- "Attrition"
x <- setdiff(names(train_obs), y)
train_obs.h2o <- as.h2o(train_obs)
local_obs.h2o <- as.h2o(local_obs)

#=========================================================================================
#1. Starting trainning
#=========================================================================================
# Create Random Forest model with ranger via caret
fit.caret <- train(
  Attrition ~ ., 
  data = train_obs, 
  method = 'ranger',
  trControl = trainControl(method = "cv", number = 5, classProbs = TRUE),
  tuneLength = 1,
  importance = 'impurity'
)

# create h2o models
h2o_rf <- h2o.randomForest(x, y, training_frame = train_obs.h2o)
h2o_glm <- h2o.glm(x, y, training_frame = train_obs.h2o, family = "binomial")
h2o_gbm <- h2o.gbm(x, y, training_frame = train_obs.h2o)

# ranger model --> model type not built in to LIME
fit.ranger <- ranger::ranger(
  Attrition ~ ., 
  data = train_obs, 
  importance = 'impurity',
  probability = TRUE
)

#=========================================================================================
#2. Global Interpretation (vip,pdp)
#=========================================================================================
###2.1 variable importance measures ----
vip(fit.ranger) + ggtitle("ranger: RF")

###2.2 partial dependence plots (PDPs) ----
# built-in PDP support in H2O
h2o.partialPlot(h2o_rf, data = train_obs.h2o, cols = "MonthlyIncome")

###2.3 individual conditional expectation (ICE) curves ----
fit.ranger %>%
  partial(pred.var = "MonthlyIncome", grid.resolution = 25, ice = TRUE) %>%
  autoplot(rug = TRUE, train = train_obs, alpha = 0.1, center = TRUE)

#=========================================================================================
#3. Local Interpretation using Local Interpretable Model-agnostic Explanations (LIME)
#=========================================================================================
###3.1 lime::lime
explainer_caret <- lime(train_obs, fit.caret, n_bins = 5)
class(explainer_caret)
summary(explainer_caret)

### 3.2 lime::explain
explanation_caret <- explain(
  x = local_obs, 
  explainer = explainer_caret, 
  n_permutations = 5000,
  dist_fun = "gower",
  kernel_width = 0.75,
  n_features = 10, 
  feature_select = "highest_weights",
  labels = "Yes"
)

tibble::glimpse(explanation_caret)

plot_features(explanation_caret)

plot_explanations(explanation_caret)

###3.3 Tuning hyperparameters
# tune LIME algorithm
explanation_caret <- explain(
  x = local_obs, 
  explainer = explainer_caret, 
  n_permutations = 5000,
  dist_fun = "manhattan",
  kernel_width = 3,
  n_features = 10, 
  feature_select = "lasso_path",
  labels = "Yes"
)

plot_features(explanation_caret)

explainer_h2o_rf  <- lime(train_obs, h2o_rf, n_bins = 5)
explainer_h2o_glm <- lime(train_obs, h2o_glm, n_bins = 5)
explainer_h2o_gbm <- lime(train_obs, h2o_gbm, n_bins = 5)

explanation_rf <- explain(local_obs, explainer_h2o_rf, n_features = 5, labels = "Yes", kernel_width = .1, feature_select = "highest_weights")
explanation_glm <- explain(local_obs, explainer_h2o_glm, n_features = 5, labels = "Yes", kernel_width = .1, feature_select = "highest_weights")
explanation_gbm <- explain(local_obs, explainer_h2o_gbm, n_features = 5, labels = "Yes", kernel_width = .1, feature_select = "highest_weights")

p1 <- plot_features(explanation_rf, ncol = 1) + ggtitle("rf")
p2 <- plot_features(explanation_glm, ncol = 1) + ggtitle("glm")
p3 <- plot_features(explanation_gbm, ncol = 1) + ggtitle("gbm")
gridExtra::grid.arrange(p1, p2, p3, nrow = 1)

#=========================================================================================
#4. Work with unsupported models
#=========================================================================================
# get the model class
class(fit.ranger)

# need to create custom model_type function
model_type.ranger <- function(x, ...) {
  # Function tells lime() what model type we are dealing with
  # 'classification', 'regression', 'survival', 'clustering', 'multilabel', etc
  return("classification")
}

model_type(fit.ranger)

# need to create custom predict_model function
predict_model.ranger <- function(x, newdata, ...) {
  # Function performs prediction and returns data frame with Response
  pred <- predict(x, newdata)
  return(as.data.frame(pred$predictions))
}

predict_model(fit.ranger, newdata = local_obs)

explainer_ranger <- lime(train_obs, fit.ranger, n_bins = 5)
explanation_ranger <- explain(local_obs, explainer_ranger, n_features = 5, n_labels = 2, kernel_width = .1)
plot_features(explanation_ranger, ncol = 2) + ggtitle("ranger")


