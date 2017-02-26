
# Load required libraries
library(class)  # for k-NN
library(ggplot2) #  for builing graphics

uci.data <- "http://archive.ics.uci.edu/ml/machine-learning-databases/pima-indians-diabetes/pima-indians-diabetes.data"


RawData <- read.delim(uci.data,header = F, sep = ",")
head(RawData) 

# isolate response
responseY <- as.matrix(RawData[,ncol(RawData)])
# only predictors
predictorX <- as.matrix(RawData[,1:(ncol(RawData) -1)])
# add columns names
colnames(predictorX) <- paste0("x",1:ncol(predictorX))

# Perform PCA
centr <- apply(predictorX,2,mean) # centrality measure = mean
disp <- apply(predictorX,2,sd) # dispersion measure = sd

predictorX <- scale(predictorX,center=centr,scale=disp)

pca_fit <- princomp(predictorX, cor=F)
summary(pca_fit) #  Summary of results

?biplot



####
# k-NN
####
library(class)
X <- pca_fit$scores[,1:2]
plot(X,main="First two components")

set.seed(1)  # Set random seed

# 80% training
n = nrow(X)
index <- sample(n,size = floor(n*0.80), replace = F)
train.set <- X[index,]
dim(train.set)

test.set <- X[-index,]
dim(test.set)

obs_class <- responseY[index,]


(model.knn <- knn(train=train.set,
                  test=test.set,
                  cl=obs_class,
                  k=19,
                  prob=T))



error <- table(model.knn,responseY[-index,]) #  Confusion matrix
(error[1,2] + error[2,1])/sum(error)  # Misclassification error for k = 19


knn_cv <- function(X,response,knn_max,n_folds){
  # Args:
  # X: dataset with explanatory variables
  # response : lables
  # knn_max: max value for nearest neighbors
  # n_folds: number of folds
  
  
  # Matrix for keeping track of error
  cv_error = matrix(NA,
                    nrow = knn_max,
                    ncol = n_folds)
  
  n <- nrow(X) # no. of observations
  
  # Create identifier for each fold
  id_fold = rep(1:n_folds,
                times = c(rep(floor(n/n_folds),n_folds-1),
                          n-floor(n/n_folds)*(n_folds-1))) # Remainder obs. are assigned to last fold
  
  # We then shuffle each id
  id_fold <- sample(id_fold)
  
  for (knn_value in 1:knn_max){  # m: Maximum number of values of k
    
    for(i in 1:n_folds){
      
      # All obs from the i-th fold will correspond to test data
      idx_test <- (id_fold==i)
      
      # Training data
      train.set <- X[!idx_test,]
      
      # Test data
      test.set <- X[idx_test,]
      
      # Vector of classes
      class.train <- response[!idx_test,]
      class.test <- response[idx_test,]
      
      # For this configuration fit k-NN model for a fixed value of
      # nearest-neighbors
      
      model.knn <- knn(train = train.set,
                       test = test.set,
                       cl = class.train,
                       k=knn_value,
                       prob=T)  # Fit model
      
      error <- table(model.knn,class.test)
      # Compute Error
      cv_error[knn_value,i] <- (error[1,2] + error[2,1])/sum(error)
    }
  }
  colnames(cv_error) <- paste0("fold_",1:n_folds) # Assign names
  rownames(cv_error) <- 1:knn_max # Assign names
  return(cv_error)
}

undebug(knn_cv)
CrossValid_knn = knn_cv(X,responseY,knn_max=15,n_folds=10)




# Plot error curves
# I will do the ggplot route
library(reshape2)
?melt
CrossValid_knn2 <- melt(CrossValid_knn,
                        value.name = "classifiction_error")

p = ggplot(CrossValid_knn2,aes(x = Var1,
                               y = classifiction_error,
                               colour = Var2)) +
  geom_line() +
  labs(x = "Number of nearest-neighbors",
       y = "Classification Error")
p
# Optimal values?
# Average error accross folds for a fixed value of n-n
avg_error = apply(CrossValid_knn,1,mean)
p  = p +
  geom_line(data = data.frame(k = 1:length(avg_error),
                              avg_error = avg_error),
            aes(x = k,y = avg_error),col = "black") +
  geom_hline(yintercept = min(avg_error),lty = "dashed") +
  geom_vline(xintercept = which.min(avg_error),lty = "dashed")
p


# leave-one-out crossvalidation
?knn.cv

