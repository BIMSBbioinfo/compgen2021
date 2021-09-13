---
title: 'compgen2021: Week 3 exercises'
author: 'ENTER YOUR NAME HERE'
output:
  pdf_document: default
  pdf: default
---

# Exercises for Week 3


### Classification 
For this set of exercises we will be using the gene expression and patient annotation data from the glioblastoma patient. You can read the data as shown below:
```{r,readMLdataEx,eval=FALSE}
library(compGenomRData)
# get file paths
fileLGGexp=system.file("extdata",
                      "LGGrnaseq.rds",
                      package="compGenomRData")
fileLGGann=system.file("extdata",
                      "patient2LGGsubtypes.rds",
                      package="compGenomRData")
# gene expression values
gexp=readRDS(fileLGGexp)

# patient annotation
patient=readRDS(fileLGGann)

```

1. Our first task is to not use any data transformation and do classification. Run the k-NN classifier on the data without any transformation or scaling. What is the effect on classification accuracy for k-NN predicting the CIMP and noCIMP status of the patient? [Difficulty: **Beginner**]

**solution:**
Below, we are using cross-validation to measure accuracy using the untransformed
and transformed data set. It is possible to do this with a hold-out test data
set as well. But with CV, we have a more reliable accuracy measure. Transforming
the data slightly increases the accuracy. However, note that we also select the top 1000 most variable genes. The data is pre-normalized, the transformation does not
have a large impact especially for this classifier. 
```{r}

#  transpose the matrix
tgexp <- t(gexp)


# transform the data, log and scale it
# there is a difference between scaling before or after transposing the data
# this will also affect the results
trs.g=t(scale(log10(gexp+1)))


# get most variable genes
SDs=apply(tgexp,2,sd )
topPreds=order(SDs,decreasing = TRUE)[1:1000]
tgexp=tgexp[,topPreds]

SDs2=apply(trs.g,2,sd )
topPreds2=order(SDs2,decreasing = TRUE)[1:1000]
trs.g=trs.g[,topPreds2]

# merge for class groups
tgexp=merge(patient,tgexp,by="row.names")

# push sample ids back to the row names
rownames(tgexp)=tgexp[,1]
tgexp=tgexp[,-1]


# merge for class groups
trs.g=merge(patient,trs.g,by="row.names")

# push sample ids back to the row names
rownames(trs.g)=trs.g[,1]
trs.g=trs.g[,-1]

set.seed(42) # set the random number seed for reproducibility 
require(caret)
# this method controls everything about training
# we will just set up 5-fold cross validation
trctrl <- trainControl(method = "cv",number=5)

# we will now train k-NN model
knn_fitTrs <- train(subtype~., data = trs.g, 
                 method = "knn",
                 trControl=trctrl,
                 tuneGrid = data.frame(k=2:7)) # try k between 2-7

# we will now train k-NN model with unscaled data
knn_fit<- train(subtype~., data = tgexp, 
                 method = "knn",
                 trControl=trctrl,
                 tuneGrid = data.frame(k=2:7))# try k between 2-7

#  cross-validation accuracy
knn_fitTrs$results
knn_fit$results



```
KNN model is better with scaled data.


2. Bootstrap resampling can be used to measure the variability of the prediction error. Use bootstrap resampling with k-NN for the prediction accuracy. How different is it from cross-validation for different $k$s? [Difficulty: **Intermediate**]

**solution:**
bootstrap and CV are different.
```{r}
set.seed(101)
trctrl2 <- trainControl(method = "boot",number=20,returnResamp="all")

# we will now train k-NN model

knnBoot <- train(subtype~., data = trs.g, method = "knn",
                trControl=trctrl2,
                tuneGrid = data.frame(k=1:12))
# best k value by bootsraping accuracy
kb<-knnBoot$bestTune

```

3. There are a number of ways to get variable importance for a classification problem. Run random forests on the classification problem above. Compare the variable importance metrics from random forest and the one obtained from DALEX applied on the random forests model. How many variables are the same in the top 10? [Difficulty: **Advanced**]

**solution:**
Here we are running variable importance on the the whole training set. Another
way of doing this, would be to run importance metrics only on the test set
after training. 

```{r}
library(DALEX)
#random forest model
set.seed(101)
trCtrl = trainControl(method = "cv",number = 5,classProb=TRUE)
rf = train(subtype~.,data = trs.g, method = "ranger", trControl=trCtrl,
          importance = "impurity")

#Then, let's run the explain function of DALEX with random forest model
explainer = DALEX::explain(rf, label = "random forest",
                            data =trs.g[,-1], y=trs.g[,1]=="CIMP")

#feature importance
dalrf = model_parts(explainer,B = 5,type="ratio")

#plot the results
plot(varImp(rf), top=10)
plot(dalrf, max_vars=10)
```

4. Come up with a unified importance score by normalizing importance scores from random forests and DALEX, followed by taking the average of those scores. [Difficulty: **Advanced**]

**solution:**
```{r}

rfimp=varImp(rf)$importance

uniScore= merge(data.frame(id=row.names(rfimp), rf=rfimp[,1]),
                data.frame(id=dalrf[,1], dalex=dalrf$dropout_loss),
                     by="id")

uniScore$avg= rowMeans(apply(uniScore[-1], 2, function(x) (x-min(x))/(max(x)-min(x))))
uniScore[order(uniScore$avg, decreasing = T)[1:10],]
```


### Regression
For this set of problems we will use the regression data set where we tried to predict the age of the sample from the methylation values. The data can be loaded as shown below: 
```{r, readMethAgeex,eval=FALSE}
# file path for CpG methylation and age
fileMethAge=system.file("extdata",
                      "CpGmeth2Age.rds",
                      package="compGenomRData")

# read methylation-age table
ameth=readRDS(fileMethAge)
```

1. Run random forest regression and plot the importance metrics. [Difficulty: **Beginner**]

**solution:**
```{r,echo=FALSE,eval=FALSE}
set.seed(101) # set the random number seed for reproducibility 

# removing less variable CpGs, not necessary but things run faster
# with less variables
ameth=ameth[,c(TRUE,matrixStats::colSds(as.matrix(ameth[,-1]))>0.1)]

# get indices for 80% of the data set
intrain <- createDataPartition(y = ameth[,1], p= 0.80)[[1]]


# seperate test and training sets
training <- ameth[intrain,]
testing <- ameth[-intrain,]

# we are not going to do any cross-validatin
# and rely on OOB error
trctrl <- trainControl(method = "none" )

# we will now train random forest model
rfregFit <- train(Age~., 
                  data = training, 
                  method = "ranger",
                  trControl=trctrl,
                  # calculate importance
                  importance="permutation", 
                  tuneGrid = data.frame(mtry=50,
                                        min.node.size = 5,
                                        splitrule="variance")
)


# plot variable importance for top 10 variables
plot(varImp(rfregFit),top=10)

# predict on the test set
testPred=predict(rfregFit,testing[,-1])

# R-squared for the test set
(cor(testPred,testing[,1]))^2

# R-squared from OOB training dataset
rfregFit$finalModel$r.squared


 
```



2. Split 20% of the methylation-age data as test data and run elastic net regression on the training portion to tune parameters and test it on the test portion. [Difficulty: **Intermediate**] 
**solution:**
```{r,echo=FALSE,eval=FALSE}
set.seed(101) # set the random number seed for reproducibility 

# removing less variable CpGs, not necessary but things run faster
# with less variables
ameth=ameth[,c(TRUE,matrixStats::colSds(as.matrix(ameth[,-1]))>0.1)]

# get indices for 80% of the data set
intrain <- createDataPartition(y = ameth[,1], p= 0.80)[[1]]


# seperate test and training sets
training <- ameth[intrain,]
testing <- ameth[-intrain,]

# we are not going to do any cross-validatin
# and rely on OOB error
trctrl <- trainControl(method = "none" )

# we will now train random forest model
rfregFit <- train(Age~., 
                  data = training, 
                  method = "ranger",
                  trControl=trctrl,
                  # calculate importance
                  importance="permutation", 
                  tuneGrid = data.frame(mtry=50,
                                        min.node.size = 5,
                                        splitrule="variance")
)


# plot variable importance for top 10 variables
plot(varImp(rfregFit),top=10)

# predict on the test set
testPred=predict(rfregFit,testing[,-1])

# R-squared for the test set
(cor(testPred,testing[,1]))^2

# R-squared from OOB training dataset
rfregFit$finalModel$r.squared


 
```


3. Run an ensemble model for regression using the **caretEnsemble** or **mlr** package and compare the results with the elastic net and random forest model. Did the test accuracy increase?
**HINT:** You need to install these extra packages and learn how to use them in the context of ensemble models. [Difficulty: **Advanced**] 

**solution:**
```{r}

library(caretEnsemble)
#set seed
set.seed(101)

tcontrol <- trainControl(
    method="cv",
     number=10,
     savePredictions=TRUE,
             classProbs=TRUE)

#set model list
modlist <- caretList(
Age~., data=training,
        trControl=my_control,
methodList=c("glmnet", "ranger"))


library(caretEnsemble)
set.seed(101)
ensemble <- caretEnsemble(
modlist,
trControl=trainControl(
method = "cv",
number=10,
savePredictions=TRUE))

class.res=predict(ensemble,testing[,-1])
postResample(pred = class.res, obs = testing[,1])

```
