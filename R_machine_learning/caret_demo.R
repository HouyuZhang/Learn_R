#https://www.machinelearningplus.com/machine-learning/caret-package/
library(caret)
suppressPackageStartupMessages(library(doParallel))
cl <- makePSOCKcluster(5)
registerDoParallel(cl)
#=========================================================================================
#1. Data Preparation and Preprocessing 
#=========================================================================================
###1.1 How to split the dataset into training and validation? ----
orange <- read.csv('https://raw.githubusercontent.com/selva86/datasets/master/orange_juice_withmissing.csv')
# Structure of the dataframe
str(orange)
head(orange[, 1:10])

set.seed(100)
# Step 1: Get row numbers for the training data
trainRowNumbers <- createDataPartition(orange$Purchase, p=0.8, list=FALSE)
# Step 2: Create the training  dataset
trainData <- orange[trainRowNumbers,]
# Step 3: Create the test dataset
testData <- orange[-trainRowNumbers,]
# Store X and Y for later use.
x = trainData[, 2:18]
y = trainData$Purchase

###1.2 Descriptive statistics ----
library(skimr)
skimmed <- skim_to_wide(trainData)
skimmed[, c(1:5, 9:11, 13, 15:16)]

###1.3 How to impute missing values using preProcess()? ----
# Create the knn imputation model on the training data
preProcess_missingdata_model <- preProcess(trainData, method='knnImpute')
preProcess_missingdata_model

# Use the imputation model to predict the values of missing data points
library(RANN)  # required for knnInpute
trainData <- predict(preProcess_missingdata_model, newdata = trainData)
anyNA(trainData)

###1.4 How to create One-Hot Encoding (dummy variables)? ----
# Creating dummy variables is converting a categorical variable to as many binary variables as here are categories.
dummies_model <- dummyVars(Purchase ~ ., data=trainData)

# Create the dummy variables using predict. The Y variable (Purchase) will not be present in trainData_mat.
trainData_mat <- predict(dummies_model, newdata = trainData)

# Convert to dataframe
trainData <- data.frame(trainData_mat)

# See the structure of the new dataset
str(trainData)

###1.5 How to preprocess to transform the data? ----
preProcess_range_model <- preProcess(trainData, method='range')
trainData <- predict(preProcess_range_model, newdata = trainData)

# Append the Y variable
trainData$Purchase <- y

apply(trainData[, 1:10], 2, FUN=function(x){c('min'=min(x), 'max'=max(x))})

#=========================================================================================
#2. Visualize the importance of variables
#=========================================================================================
#Box plot
featurePlot(x = trainData[, 1:18], 
            y = trainData$Purchase, 
            plot = "box",
            strip=strip.custom(par.strip.text=list(cex=.7)),
            scales = list(x = list(relation="free"), 
                          y = list(relation="free")))

#Density plots
featurePlot(x = trainData[, 1:18], 
            y = trainData$Purchase, 
            plot = "density",
            strip=strip.custom(par.strip.text=list(cex=.7)),
            scales = list(x = list(relation="free"), 
                          y = list(relation="free")))

#Having visualised the relationships between X and Y, We can only say which variables are likely to be 
#important to predict Y. It may not be wise to conclude which variables are NOT important.
#=========================================================================================
#3. Feature Selection using recursive feature elimination (RFE)
#=========================================================================================
set.seed(100)
options(warn=-1)

subsets <- c(1:5, 10, 15, 18)

ctrl <- rfeControl(functions = rfFuncs,
                   method = "repeatedcv",
                   repeats = 5,
                   verbose = FALSE)

lmProfile <- rfe(x=trainData[, 1:18], y=trainData$Purchase,
                 sizes = subsets,
                 rfeControl = ctrl)

lmProfile
#However, it is not a mandate that only including these 3 variables will always give high accuracy over 
#larger sized models, because, the rfe() we just implemented is particular to random forest based rfFuncs

#Plus also, since the training dataset isn’t large enough, the other predictors may not have 
#had the chance to show its worth.
#=========================================================================================
#4. Training and Tuning the model
#=========================================================================================
###4.1 How to train() the model and interpret the results? ----
# See available algorithms in caret
names(getModelInfo())
modelLookup('earth')

# Set the seed for reproducibility
set.seed(100)

# Train the model using randomForest and predict on the training data itself.
model_mars <- train(Purchase ~ ., data=trainData, method='earth')
model_mars
ggplot(model_mars, main="Model Accuracies with MARS")
#train() has already done a basic cross validation and hyper parameter tuning. default behaviour.
fitted <- predict(model_mars)

###4.2 How to compute variable importance? ----
varimp_mars <- varImp(model_mars)
plot(varimp_mars, main="Variable Importance with MARS")

###4.3 Predict on testData ----
# Step 1: Impute missing values 
testData2 <- predict(preProcess_missingdata_model, testData)  
# Step 2: Create one-hot encodings (dummy variables)
testData3 <- predict(dummies_model, testData2)
# Step 3: Transform the features to range between 0 and 1
testData4 <- predict(preProcess_range_model, testData3)
head(testData4[, 1:10])

# Predict on testData
predicted <- predict(model_mars, testData4)
head(predicted)

# Compute the confusion matrix
confusionMatrix(reference = testData$Purchase, data = predicted, mode='everything', positive='MM')


###4.4 How to do hyperparameter tuning to optimize the model for better performance? ----
# Define the training control
fitControl <- trainControl(
  method = 'cv',                   # k-fold cross validation
  number = 5,                      # number of folds
  savePredictions = 'final',       # saves predictions for optimal tuning parameter
  classProbs = T,                  # should class probabilities be returned
  summaryFunction=prSummary  # results summary function
) 
#4.4.1 use tuneLength 

#Tune hyper parameters by setting tuneLength
set.seed(100)
model_mars2 <-  train(Purchase ~ ., data=trainData, method='earth', tuneLength = 5, 
                      metric='AUC', trControl = fitControl)

#Predict on testData and Compute the confusion matrix
predicted2 <- predict(model_mars2, testData4)
confusionMatrix(reference = testData$Purchase, data = predicted2, mode='everything', positive='MM')

#4.4.2 use tuneGrid

#Define the tuneGrid
marsGrid <-  expand.grid(nprune = c(2, 4, 6, 8, 10), degree = c(1, 2, 3))
#Tune hyper parameters by setting tuneGrid
set.seed(100)
model_mars3 <- train(Purchase ~ ., data=trainData, method='earth', metric='ROC', 
                    tuneGrid = marsGrid, trControl = fitControl)
#Predict on testData and Compute the confusion matrix
predicted3 <- predict(model_mars3, testData4)

#=========================================================================================
#S. Top 15 Evaluation Metrics for Classification Models
#https://www.machinelearningplus.com/machine-learning/evaluation-metrics-classification-models-r/
#=========================================================================================
#1. The Confusion Matrix
confusionMatrix(reference = testData$Purchase, data = predicted3, mode='everything', positive='MM')
ROC <- twoClassSummary(table, lev = levels(table$obs))[[1]]
AUC <- prSummary(table, lev = rev(levels(table$obs)))[[1]]
#=========================================================================================
#5. Choose model
#=========================================================================================
###5.1 Training Adaboost ----
set.seed(100)
model_adaboost = train(Purchase ~ ., data=trainData, method='adaboost', tuneLength=2, trControl = fitControl)
model_adaboost

###5.2 Training Random Forest ----
set.seed(100)
model_rf = train(Purchase ~ ., data=trainData, method='rf', tuneLength=5, trControl = fitControl)
model_rf

###5.3 Training xgBoost Dart ----
set.seed(100)
model_xgbDART = train(Purchase ~ ., data=trainData, method='xgbDART', tuneLength=5, trControl = fitControl, verbose=F)
model_xgbDART

###5.4 Training SVM ----
set.seed(100)
model_svmRadial = train(Purchase ~ ., data=trainData, method='svmRadial', tuneLength=15, trControl = fitControl)
model_svmRadial

###5.5 Run resamples() to compare the models ----
# Compare model performances using resample()
models_compare <- resamples(list(ADABOOST=model_adaboost, RF=model_rf, 
                                 XGBDART=model_xgbDART, MARS=model_mars3, SVM=model_svmRadial))

# Summary of the models performances
summary(models_compare)

# Draw box plots to compare models
scales <- list(x=list(relation="free"), y=list(relation="free"))
bwplot(models_compare, scales=scales)
#=========================================================================================
#6. Ensembling the predictions
#=========================================================================================
###6.1 How to ensemble predictions from multiple models using caretEnsemble?
library(caretEnsemble)
# Stacking Algorithms - Run multiple algos in one call.
trainControl <- trainControl(method="repeatedcv", 
                             number=10, 
                             repeats=3,
                             savePredictions=TRUE, 
                             classProbs=TRUE)

algorithmList <- c('rf', 'adaboost', 'earth', 'xgbDART', 'svmRadial')

set.seed(100)
models <- caretList(Purchase ~ ., data=trainData, trControl=trainControl, methodList=algorithmList) 
results <- resamples(models)
summary(results)

# Box plots to compare models
scales <- list(x=list(relation="free"), y=list(relation="free"))
bwplot(results, scales=scales)

### 6.2 How to combine the predictions of multiple models to form a final prediction ----
# Create the trainControl
set.seed(101)
stackControl <- trainControl(method="repeatedcv", 
                             number=10, 
                             repeats=3,
                             savePredictions=TRUE, 
                             classProbs=TRUE)

# Ensemble the predictions of `models` to form a new combined prediction based on glm
stack.glm <- caretStack(models, method="glm", metric="Accuracy", trControl=stackControl)
print(stack.glm)

# Predict on testData
stack_predicteds <- predict(stack.glm, newdata=testData4)
head(stack_predicteds)

stopCluster(cl)



