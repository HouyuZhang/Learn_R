library(glmnet)
library(doParallel)
#=========================================================================================
# 1. https://web.stanford.edu/~hastie/glmnet/glmnet_alpha.html#lin (Official)
#=========================================================================================
###1. Linear Regression
#1.1 Gaussian Family

#1.2 Multiresponse Gaussian Family
load("E:/Program Files/R/R-3.5.3/library/glmnet/data/MultiGaussianExample.RData")
mfit = glmnet(x, y, family = "mgaussian")
plot(mfit, xvar = "lambda", label = TRUE, type.coef = "2norm")
predict(mfit, newx = x[1:5,], s = c(0.1, 0.01))
cvmfit = cv.glmnet(x, y, family = "mgaussian")
plot(cvmfit)

###2. Logistic Regression
#2.1 Binomial Models
load("E:/Program Files/R/R-3.5.3/library/glmnet/data/BinomialExample.RData")
fit = glmnet(x, y, family = "binomial")

plot(fit, xvar = "dev", label = TRUE)
predict(fit, newx = x[1:5,], type = "class", s = c(0.05, 0.01))
cvfit = cv.glmnet(x, y, family = "binomial", type.measure = "class")
plot(cvfit)

#2.2 Multinomial Models

###3. Poisson Models
load("E:/Program Files/R/R-3.5.3/library/glmnet/data/PoissonExample.RData")
fit = glmnet(x, y, family = "poisson")
plot(fit)

###4. Cox Models
load("E:/Program Files/R/R-3.5.3/library/glmnet/data/CoxExample.RData")
y[1:5,]
fit = glmnet(x, y, family = "cox")
plot(fit)

#2.5 Sparse Matrices

#=========================================================================================
# 2. https://www4.stat.ncsu.edu/~post/josh/LASSO_Ridge_Elastic_Net_-_Examples.html
#=========================================================================================

library(MASS)  # Package needed to generate correlated precictors
library(glmnet)  # Package to fit ridge/lasso/elastic net models

###1. Generate data
set.seed(19873)
n <- 100    # Number of observations
p <- 50     # Number of predictors included in model
CovMatrix <- outer(1:p, 1:p, function(x,y) {.7^abs(x-y)})
x <- mvrnorm(n, rep(0,p), CovMatrix)
y <- 10 * apply(x[, 1:2], 1, sum) + 
  5 * apply(x[, 3:4], 1, sum) +
  apply(x[, 5:14], 1, sum) +
  rnorm(n)

# Split data into train and test sets
train_rows <- sample(1:n, .66*n)
x.train <- x[train_rows, ]
x.test <- x[-train_rows, ]

y.train <- y[train_rows]
y.test <- y[-train_rows]

###2. Fit models 
# (For plots on left):
fit.lasso <- glmnet(x.train, y.train, family="gaussian", alpha=1)
fit.ridge <- glmnet(x.train, y.train, family="gaussian", alpha=0)
fit.elnet <- glmnet(x.train, y.train, family="gaussian", alpha=.5)

# 10-fold Cross validation for each alpha = 0, 0.1, ... , 0.9, 1.0
# (For plots on Right)
for (i in 0:10) {
  assign(paste("fit", i, sep=""), cv.glmnet(x.train, y.train, type.measure="mse",alpha=i/10,family="gaussian"))
}

###3. Plot solution path and cross-validated MSE as function of λ:
par(mfrow=c(3,2))
# For plotting options, type '?plot.glmnet' in R console
plot(fit.lasso, xvar="lambda")
plot(fit10, main="LASSO")

plot(fit.ridge, xvar="lambda")
plot(fit0, main="Ridge")

plot(fit.elnet, xvar="lambda")
plot(fit5, main="Elastic Net")

yhat0 <- predict(fit0, s=fit0$lambda.1se, newx=x.test)
coef(fit0,s=fit0$lambda.1se)
yhat1 <- predict(fit1, s=fit1$lambda.1se, newx=x.test)
yhat2 <- predict(fit2, s=fit2$lambda.1se, newx=x.test)
yhat3 <- predict(fit3, s=fit3$lambda.1se, newx=x.test)
yhat4 <- predict(fit4, s=fit4$lambda.1se, newx=x.test)
yhat5 <- predict(fit5, s=fit5$lambda.1se, newx=x.test)
yhat6 <- predict(fit6, s=fit6$lambda.1se, newx=x.test)
yhat7 <- predict(fit7, s=fit7$lambda.1se, newx=x.test)
yhat8 <- predict(fit8, s=fit8$lambda.1se, newx=x.test)
yhat9 <- predict(fit9, s=fit9$lambda.1se, newx=x.test)
yhat10 <- predict(fit10, s=fit10$lambda.1se, newx=x.test)

mse0 <- mean((y.test - yhat0)^2)
mse1 <- mean((y.test - yhat1)^2)
mse2 <- mean((y.test - yhat2)^2)
mse3 <- mean((y.test - yhat3)^2)
mse4 <- mean((y.test - yhat4)^2)
mse5 <- mean((y.test - yhat5)^2)
mse6 <- mean((y.test - yhat6)^2)
mse7 <- mean((y.test - yhat7)^2)
mse8 <- mean((y.test - yhat8)^2)
mse9 <- mean((y.test - yhat9)^2)
mse10 <- mean((y.test - yhat10)^2)

