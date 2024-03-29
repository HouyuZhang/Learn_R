#http://uc-r.github.io/naive_bayes
library(rsample)  # data splitting 
library(dplyr)    # data transformation
library(ggplot2)  # data visualization
library(caret)    # implementing with caret
library(h2o)      # implementing with h2o

#=========================================================================================
#0. Data arrangement
#=========================================================================================
# conver some numeric variables to factors
attrition <- attrition %>%
  mutate(
    JobLevel = factor(JobLevel),
    StockOptionLevel = factor(StockOptionLevel),
    TrainingTimesLastYear = factor(TrainingTimesLastYear)
  )

# Create training (70%) and test (30%) sets for the attrition data.
# Use set.seed for reproducibility
set.seed(123)
split <- initial_split(attrition, prop = 0.7, strata = "Attrition")
train <- training(split)
test  <- testing(split)

# distribution of Attrition rates across train & test set
table(train$Attrition) %>% prop.table()
table(test$Attrition) %>% prop.table()

train %>%
  filter(Attrition == "Yes") %>%
  select_if(is.numeric) %>%
  cor() %>%
  corrplot::corrplot()

train %>% 
  select(Age, DailyRate, DistanceFromHome, HourlyRate, MonthlyIncome, MonthlyRate) %>% 
  gather(metric, value) %>% 
  ggplot(aes(value, fill = metric)) + 
  geom_density(show.legend = FALSE) + 
  facet_wrap(~ metric, scales = "free")
#=========================================================================================
#1. Trainning using Caret
#=========================================================================================
###1.1 regular trainning ----
# create response and feature data
features <- setdiff(names(train), "Attrition")
x <- train[, features]
y <- train$Attrition

# set up 10-fold cross validation procedure
train_control <- trainControl(
  method = "cv", 
  number = 10
)

# train model
nb.m1 <- train(
  x = x,
  y = y,
  method = "nb",
  trControl = train_control
)

# results
confusionMatrix(nb.m1)

###1.2 Tuning parameters ----
# set up tuning grid
search_grid <- expand.grid(
  usekernel = c(TRUE, FALSE),
  fL = 0:5,
  adjust = seq(0, 5, by = 1)
)

# train model
nb.m2 <- train(
  x = x,
  y = y,
  method = "nb",
  trControl = train_control,
  tuneGrid = search_grid,
  preProc = c("BoxCox", "center", "scale", "pca")
)

# top 5 modesl
nb.m2$results %>% 
  top_n(5, wt = Accuracy) %>%
  arrange(desc(Accuracy))

# plot search grid results
plot(nb.m2)

# results for best model
confusionMatrix(nb.m2)

pred <- predict(nb.m2, newdata = test)
confusionMatrix(pred, test$Attrition)

#=========================================================================================
#2. Trainning using h2o
#=========================================================================================
###2.1 regular trainning ----
# start up h2o
h2o.no_progress()
h2o.init()

# create feature names
y <- "Attrition"
x <- setdiff(names(train), y)

# h2o cannot ingest ordered factors
train.h2o <- train %>%
  mutate_if(is.factor, factor, ordered = FALSE) %>%
  as.h2o()

# train model
nb.h2o <- h2o.naiveBayes(
  x = x,
  y = y,
  training_frame = train.h2o,
  nfolds = 10,
  laplace = 0
)

# assess results
h2o.confusionMatrix(nb.h2o)

###2.2 Tuning parameters ----
# do a little preprocessing
preprocess <- preProcess(train, method = c("BoxCox", "center", "scale", "pca"))
train_pp   <- predict(preprocess, train)
test_pp    <- predict(preprocess, test)

# convert to h2o objects
train_pp.h2o <- train_pp %>%
  mutate_if(is.factor, factor, ordered = FALSE) %>%
  as.h2o()

test_pp.h2o <- test_pp %>%
  mutate_if(is.factor, factor, ordered = FALSE) %>%
  as.h2o()

# get new feature names --> PCA preprocessing reduced and changed some features
x <- setdiff(names(train_pp), "Attrition")

# create tuning grid
hyper_params <- list(
  laplace = seq(0, 5, by = 0.5)
)

# build grid search 
grid <- h2o.grid(
  algorithm = "naivebayes",
  grid_id = "nb_grid",
  x = x, 
  y = y, 
  training_frame = train_pp.h2o, 
  nfolds = 10,
  hyper_params = hyper_params
)

# Sort the grid models by mse
sorted_grid <- h2o.getGrid("nb_grid", sort_by = "accuracy", decreasing = TRUE)
sorted_grid

# grab top model id
best_h2o_model <- sorted_grid@model_ids[[1]]
best_model <- h2o.getModel(best_h2o_model)

# confusion matrix of best model
h2o.confusionMatrix(best_model)

# ROC curve
auc <- h2o.auc(best_model, xval = TRUE)
fpr <- h2o.performance(best_model, xval = TRUE) %>% h2o.fpr() %>% .[['fpr']]
tpr <- h2o.performance(best_model, xval = TRUE) %>% h2o.tpr() %>% .[['tpr']]
data.frame(fpr = fpr, tpr = tpr) %>%
  ggplot(aes(fpr, tpr)) +
  geom_line() + 
  ggtitle(sprintf('AUC: %f', auc))

# evaluate on test set
h2o.performance(best_model, newdata = test_pp.h2o)

# predict new data
h2o.predict(nb.h2o, newdata = test_pp.h2o)

# shut down h2o
h2o.shutdown(prompt = FALSE)


