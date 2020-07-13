library(e1071)
#=========================================================================================
#1. Plot hyperplane in 3D plot
#=========================================================================================
n = 100
nnew = 50

# Simulate some data
set.seed(0725)
group = sample(2, n, replace=T)
dat   = data.frame(group=factor(group), matrix(rnorm(n*3, rep(group, each=3)), ncol=3, byrow=T))

fit = svm(group ~ ., data=dat)
# Plot original data
library(rgl)
plot3d(dat[,-1], col=dat$group)

# Get decision values for a new data grid
newdat.list = lapply(dat[,-1], function(x) seq(min(x), max(x), len=nnew))
newdat      = expand.grid(newdat.list)
newdat.pred = predict(fit, newdata=newdat, decision.values=T)
newdat.dv   = attr(newdat.pred, 'decision.values')
newdat.dv   = array(newdat.dv, dim=rep(nnew, 3))

# Fit/plot an isosurface to the decision boundary
library(misc3d)
contour3d(newdat.dv, level=0, x=newdat.list$X1, y=newdat.list$X2, z=newdat.list$X3, add=T)

#=========================================================================================
#2. http://uc-r.github.io/svm
#=========================================================================================
library(tidyverse)    # data manipulation and visualization
library(kernlab)      # SVM methodology
library(e1071)        # SVM methodology
library(RColorBrewer) # customized coloring of plots
###1. Maximal Margin Classifier ----
set.seed(0725)
x <- matrix(rnorm(20*2), ncol = 2)
y <- c(rep(-1,10), rep(1,10))
x[y==1,] <- x[y==1,] + 3/2
dat <- data.frame(x=x, y=as.factor(y))

ggplot(data = dat, aes(x = x.2, y = x.1, color = y, shape = y)) + 
  geom_point(size = 2) +
  scale_color_manual(values=c("#000000", "#FF0000")) +
  theme(legend.position = "none")

svmfit <- svm(y~., data = dat, kernel = "linear", scale = FALSE)
plot(svmfit, dat)

# fit model and produce plot
kernfit <- ksvm(x, y, type = "C-svc", kernel = 'vanilladot')
plot(kernfit, data = x)

###2. Support Vector Classifiers ----
x <- matrix(rnorm(20*2), ncol = 2)
y <- c(rep(-1,10), rep(1,10))
x[y==1,] <- x[y==1,] + 1
dat <- data.frame(x=x, y=as.factor(y))

ggplot(data = dat, aes(x = x.2, y = x.1, color = y, shape = y)) + 
  geom_point(size = 2) +
  scale_color_manual(values=c("#000000", "#FF0000")) +
  theme(legend.position = "none")


svmfit <- svm(y~., data = dat, kernel = "linear", cost = 10)
plot(svmfit, dat)

# Fit Support Vector Machine model to data set
kernfit <- ksvm(x,y, type = "C-svc", kernel = 'vanilladot', C = 100)
plot(kernfit, data = x)

# find optimal cost of misclassification
tune.out <- tune(svm, y~., data = dat, kernel = "linear",
                 ranges = list(cost = c(0.001, 0.01, 0.1, 1, 5, 10, 100)))
# extract the best model
(bestmod <- tune.out$best.model)

ypred <- predict(bestmod, dat)
(misclass <- table(predict = ypred, truth = dat$y))

###3. Support Vector Machines ----
x <- matrix(rnorm(200*2), ncol = 2)
x[1:100,] <- x[1:100,] + 2.5
x[101:150,] <- x[101:150,] - 2.5
y <- c(rep(1,150), rep(2,50))
dat <- data.frame(x=x,y=as.factor(y))

ggplot(data = dat, aes(x = x.2, y = x.1, color = y, shape = y)) + 
  geom_point(size = 2) +
  scale_color_manual(values=c("#000000", "#FF0000")) +
  theme(legend.position = "none")

set.seed(0725)
train <- base::sample(200,100, replace = FALSE)
svmfit <- svm(y~., data = dat[train,], kernel = "radial", gamma = 1, cost = 1)

plot(svmfit, dat)

# Fit radial-based SVM in kernlab
kernfit <- ksvm(x[train,],y[train], type = "C-svc", kernel = 'rbfdot', C = 1, scaled = c())
plot(kernfit, data = x[train,])

# tune model to find optimal cost, gamma values
tune.out <- tune(svm, y~., data = dat[train,], kernel = "radial",
                 ranges = list(cost = c(0.1,1,10,100,1000), gamma = c(0.5,1,2,3,4)))
tune.out$best.model

# validate model performance
(valid <- table(true = dat[-train,"y"], pred = predict(tune.out$best.model, newx = dat[-train,])))

###4. SVMs for Multiple Classes ----
x <- rbind(x, matrix(rnorm(50*2), ncol = 2))
y <- c(y, rep(0,50))
x[y==0,2] <- x[y==0,2] + 2.5
dat <- data.frame(x=x, y=as.factor(y))
# plot data set
ggplot(data = dat, aes(x = x.2, y = x.1, color = y, shape = y)) + 
  geom_point(size = 2) +
  scale_color_manual(values=c("#000000","#FF0000","#00BA00")) +
  theme(legend.position = "none")

svmfit <- svm(y~., data = dat, kernel = "radial", cost = 10, gamma = 1)
plot(svmfit, dat)

#construct table
ypred <- predict(svmfit, dat)
(misclass <- table(predict = ypred, truth = dat$y))


kernfit <- ksvm(as.matrix(dat[,2:1]),dat$y, type = "C-svc", kernel = 'rbfdot', C = 100, scaled = c())

# Create a fine grid of the feature space
x.1 <- seq(from = min(dat$x.1), to = max(dat$x.1), length = 100)
x.2 <- seq(from = min(dat$x.2), to = max(dat$x.2), length = 100)
x.grid <- expand.grid(x.2, x.1)

# Get class predictions over grid
pred <- predict(kernfit, newdata = x.grid)

# Plot the results
cols <- brewer.pal(3, "Set1")
plot(x.grid, pch = 19, col = adjustcolor(cols[pred], alpha.f = 0.05))

classes <- matrix(pred, nrow = 100, ncol = 100)
contour(x = x.2, y = x.1, z = classes, levels = 1:3, labels = "", add = TRUE)

points(dat[, 2:1], pch = 19, col = cols[predict(kernfit)])

