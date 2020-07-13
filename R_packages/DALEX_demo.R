#http://uc-r.github.io/dalex
#provide Descriptive mAchine Learning EXplanations ranging from global to local interpretability methods.
library(rsample)
library(dplyr)
library(h2o)
library(DALEX)
library(ggplot2)

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
#1. DALEX procedures
#=========================================================================================
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

# convert feature data to non-h2o objects
x_valid <- as.data.frame(splits$valid)[, x]

# make response variable numeric binary vector
y_valid <- as.vector(as.numeric(as.character(splits$valid$Attrition)))
head(y_valid)

# create custom predict function
pred <- function(model, newdata)  {
  results <- as.data.frame(h2o.predict(model, as.h2o(newdata)))
  return(results[[3L]])
}

pred(rf, x_valid) %>% head()

# elastic net explainer
explainer_glm <- explain(
  model = glm,
  data = x_valid,
  y = y_valid,
  predict_function = pred,
  label = "h2o glm"
)

# random forest explainer
explainer_rf <- explain(
  model = rf,
  data = x_valid,
  y = y_valid,
  predict_function = pred,
  label = "h2o rf"
)

# GBM explainer
explainer_gbm <- explain(
  model = gbm,
  data = x_valid,
  y = y_valid,
  predict_function = pred,
  label = "h2o gbm"
)

# example of explainer object
class(explainer_glm)
summary(explainer_glm)

#=========================================================================================
#2. Global interpretation
#=========================================================================================
###2.1 Residual diagnostics ----
# compute predictions & residuals
resids_glm <- model_performance(explainer_glm)
resids_rf  <- model_performance(explainer_rf)
resids_gbm <- model_performance(explainer_gbm)

# assess quantiles for residuals
resids_glm
resids_rf
resids_gbm

# create comparison plot of residuals for each model
p1 <- plot(resids_glm, resids_rf, resids_gbm)
p2 <- plot(resids_glm, resids_rf, resids_gbm, geom = "boxplot")

gridExtra::grid.arrange(p1, p2, nrow = 1)

###2.2 Variable importance ----
# compute permutation-based variable importance
vip_glm <- variable_importance(explainer_glm, n_sample = -1, loss_function = loss_root_mean_square) 
vip_rf  <- variable_importance(explainer_rf, n_sample = -1, loss_function = loss_root_mean_square)
vip_gbm <- variable_importance(explainer_gbm, n_sample = -1, loss_function = loss_root_mean_square)

plot(vip_glm, vip_rf, vip_gbm, max_vars = 10)

###2.3 Predictor-response relationship ----
#we likely want to understand how the relationship between 
#these influential variables and the predicted response differ between the models.
# compute PDP for a given variable --> uses the pdp package
pdp_glm  <- variable_response(explainer_glm, variable =  "Age", type = "pdp")
pdp_rf   <- variable_response(explainer_rf,  variable =  "Age", type = "pdp")
pdp_gbm  <- variable_response(explainer_gbm, variable =  "Age", type = "pdp")

plot(pdp_glm, pdp_rf, pdp_gbm)

cat_glm  <- variable_response(explainer_glm, variable = "EnvironmentSatisfaction", type = "factor")
cat_rf  <- variable_response(explainer_rf, variable = "EnvironmentSatisfaction", type = "factor")
cat_gbm  <- variable_response(explainer_gbm, variable = "EnvironmentSatisfaction", type = "factor")
plot(cat_glm, cat_rf, cat_gbm)

#=========================================================================================
#3. Local interpretation
#=========================================================================================
# create a single observation
new_cust <- splits$valid[1, ] %>% as.data.frame()

# compute breakdown distances
new_cust_glm <- prediction_breakdown(explainer_glm, observation = new_cust)
new_cust_rf  <- prediction_breakdown(explainer_rf, observation = new_cust)
new_cust_gbm <- prediction_breakdown(explainer_gbm, observation = new_cust)

# class of prediction_breakdown output
class(new_cust_gbm)
# check out the top 10 influential variables for this observation
new_cust_gbm[1:10, 1:5]

plot(new_cust_gbm)

# filter for top 10 influential variables for each model and plot
list(new_cust_glm, new_cust_rf, new_cust_gbm) %>%
  purrr::map(~ top_n(., 11, wt = abs(contribution))) %>%
  do.call(rbind, .) %>%
  mutate(variable = paste0(variable, " (", label, ")")) %>%
  ggplot(aes(contribution, reorder(variable, contribution))) +
  geom_point() +
  geom_vline(xintercept = 0, size = 3, color = "white") +
  facet_wrap(~ label, scales = "free_y", ncol = 1) +
  ylab(NULL)










