---
title: 'compgen2021: Week 1 exercises'
author: 'ENTER YOUR NAME HERE'
output:
  pdf_document: default
  pdf: default
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

# Exercises for Week1

## Statistics for genomics 

### How to summarize collection of data points: The idea behind statistical distributions

1. Calculate the means and variances 
of the rows of the following simulated data set, and plot the distributions
of means and variances using `hist()` and `boxplot()` functions. [Difficulty: **Beginner/Intermediate**]  
```{r getDataChp3Ex,eval=FALSE}
set.seed(100)

#sample data matrix from normal distribution
gset=rnorm(600,mean=200,sd=70)
data=matrix(gset,ncol=6)
```

**solution:**
```{r,echo=FALSE,eval=FALSE}
require(matrixStats)
means=rowMeans(data)
vars=rowVars(data)
hist(means)
hist(vars) 
boxplot(means)
boxplot(vars)
```


2. Using the data generated above, calculate the standard deviation of the
distribution of the means using the `sd()` function. Compare that to the expected
standard error obtained from the central limit theorem keeping in mind the
population parameters were  $\sigma=70$ and $n=6$. How does the estimate from the random samples change if we simulate more data with
`data=matrix(rnorm(6000,mean=200,sd=70),ncol=6)`? [Difficulty: **Beginner/Intermediate**] 


**solution:**
```{r,echo=FALSE,eval=FALSE}
samples=sd(means)

clt.se=70/sqrt(6) 
```
CLT standard error and the simulated standard error are similar. 
If we simulate more data with the same parameters and with the same sample size.
It wouldn't make a difference, the SE we obtained from the simulated 
data and CLT would still be similar. 


3. Simulate 30 random variables using the `rpois()` function. Do this 1000 times and calculate the mean of each sample. Plot the sampling distributions of the means
using a histogram. Get the 2.5th and 97.5th percentiles of the
distribution. [Difficulty: **Beginner/Intermediate**] 

**solution:**
```{r,echo=FALSE,eval=FALSE}
hist((do(1000)*mean(rpois(30,lambda=5)))[,1])
```

4. Use the `t.test()` function to calculate confidence intervals
of the mean on the first random sample `pois1` simulated from the `rpois()` function below. [Difficulty: **Intermediate**] 
```{r exRpoisChp3,eval=FALSE}
#HINT
set.seed(100)

#sample 30 values from poisson dist with lamda paramater =30
pois1=rpois(30,lambda=5)

```

**solution:**
```{r,echo=FALSE,eval=FALSE}
require(mosaic)

t.test(pois1)
```


5. Use the bootstrap confidence interval for the mean on `pois1`. [Difficulty: **Intermediate/Advanced**] 

**solution:**
```{r,echo=FALSE,eval=FALSE}
quantile((do(1000)*mean(sample(pois1,30,replace=T)))[,1],probs=c(0.025,0.975))
 
```



### How to test for differences in samples
1. Test the difference of means of the following simulated genes
using the randomization, `t-test()`, and `wilcox.test()` functions.
Plot the distributions using histograms and boxplots. [Difficulty: **Intermediate/Advanced**] 
```{r exRnorm1chp3,eval=FALSE}
set.seed(101)
gene1=rnorm(30,mean=4,sd=3)
gene2=rnorm(30,mean=3,sd=3)

```

**solution:**
```{r,echo=FALSE,eval=FALSE}

t.test(gene1,gene2)
wilcox.test(gene1,gene2)

#rand test
org.diff=mean(gene1)-mean(gene2)
gene.df=data.frame(exp=c(gene1,gene2),
                  group=c( rep("test",30),rep("control",30) ) )


exp.null <- do(1000) * diff(mosaic::mean(exp ~ shuffle(group), data=gene.df))

sum(exp.null>org.diff)/1000

```



2. Test the difference of the means of the following simulated genes
using the randomization, `t-test()` and `wilcox.test()` functions.
Plot the distributions using histograms and boxplots. [Difficulty: **Intermediate/Advanced**] 
```{r exRnorm2chp3,eval=FALSE}
set.seed(100)
gene1=rnorm(30,mean=4,sd=2)
gene2=rnorm(30,mean=2,sd=2)

```

**solution:**
```{r,echo=FALSE,eval=FALSE}
t.test(gene1,gene2)
wilcox.test(gene1,gene2) 
#rand test
org.diff=mean(gene1)-mean(gene2)
gene.df=data.frame(exp=c(gene1,gene2),
                  group=c( rep("test",30),rep("control",30) ) )


exp.null <- do(1000) * diff(mosaic::mean(exp ~ shuffle(group), data=gene.df))

sum(exp.null>org.diff)/1000
```

3. We need an extra data set for this exercise. Read the gene expression data set as follows:
`gexpFile=system.file("extdata","geneExpMat.rds",package="compGenomRData") data=readRDS(gexpFile)`. The data has 100 differentially expressed genes. The first 3 columns are the test samples, and the last 3 are the control samples. Do 
a t-test for each gene (each row is a gene), and record the p-values.
Then, do a moderated t-test, as shown in section "Moderated t-tests" in this chapter, and record 
the p-values. Make a p-value histogram and compare two approaches in terms of the number of significant tests with the $0.05$ threshold.
On the p-values use FDR (BH), Bonferroni and q-value adjustment methods.
Calculate how many adjusted p-values are below 0.05 for each approach.
[Difficulty: **Intermediate/Advanced**] 

**solution:**
```{r,echo=FALSE,eval=FALSE}
gexpFile=system.file("extdata","geneExpMat.rds",package="compGenomRData") 
data=readRDS(gexpFile)
 
group1=1:3
group2=4:6
n1=3
n2=3
dx=rowMeans(data[,group1])-rowMeans(data[,group2])
  
require(matrixStats)

# get the esimate of pooled variance 
stderr = sqrt( (rowVars(data[,group1])*(n1-1) + 
       rowVars(data[,group2])*(n2-1)) / (n1+n2-2) * ( 1/n1 + 1/n2 ))

# do the shrinking towards median
mod.stderr = (stderr + median(stderr)) / 2 # moderation in variation

# esimate t statistic with moderated variance
t.mod <- dx / mod.stderr

# calculate P-value of rejecting null 
p.mod = 2*pt( -abs(t.mod), n1+n2-2 )

# esimate t statistic without moderated variance
t = dx / stderr

# calculate P-value of rejecting null 
p = 2*pt( -abs(t), n1+n2-2 )

par(mfrow=c(1,2))
hist(p,col="cornflowerblue",border="white",main="",xlab="P-values t-test")
mtext(paste("signifcant tests:",sum(p<0.05))  )
hist(p.mod,col="cornflowerblue",border="white",main="",
     xlab="P-values mod. t-test")
mtext(paste("signifcant tests:",sum(p.mod<0.05))  )

```



### Relationship between variables: Linear models and correlation

Below we are going to simulate X and Y values that are needed for the 
rest of the exercise.
```{r exLM1chp3,eval=FALSE}
# set random number seed, so that the random numbers from the text
# is the same when you run the code.
set.seed(32)

# get 50 X values between 1 and 100
x = runif(50,1,100)

# set b0,b1 and variance (sigma)
b0 = 10
b1 = 2
sigma = 20
# simulate error terms from normal distribution
eps = rnorm(50,0,sigma)
# get y values from the linear equation and addition of error terms
y = b0 + b1*x+ eps
```


1. Run the code then fit a line to predict Y based on X. [Difficulty:**Intermediate**] 

**solution:**
```{r,echo=FALSE,eval=FALSE}
lm(y~x)
 
```

2. Plot the scatter plot and the fitted line. [Difficulty:**Intermediate**] 

**solution:**
```{r,echo=FALSE,eval=FALSE}
plot(x,y)
abline(lm(y~x))
 
```

3. Calculate correlation and R^2. [Difficulty:**Intermediate**] 

**solution:**
```{r,echo=FALSE,eval=FALSE}
cor(x,y) # pearson cor
cor(x,y)^2 # R squared
summary(lm(y~x))
 
```

4. Run the `summary()` function and 
try to extract P-values for the model from the object
returned by `summary`. See `?summary.lm`. [Difficulty:**Intermediate/Advanced**] 

**solution:**
```{r,echo=FALSE,eval=FALSE}
summary(lm(y~x))
 
```


5. Plot the residuals vs. the fitted values plot, by calling the `plot()` 
function with `which=1` as the second argument. First argument
is the model returned by `lm()`. [Difficulty:**Advanced**] 

**solution:**
```{r,echo=FALSE,eval=FALSE}
mod=lm(y~x)
plot(mod,which=1)
 
```


6. For the next exercises, read the data set histone modification data set. Use the following to get the path to the file:
```
hmodFile=system.file("extdata",
                    "HistoneModeVSgeneExp.rds",
                     package="compGenomRData")
```
There are 3 columns in the dataset. These are measured levels of H3K4me3,
H3K27me3 and gene expression per gene. Once you read in the data, plot the scatter plot for H3K4me3 vs. expression. [Difficulty:**Beginner**] 

**solution:**
```{r,echo=FALSE,eval=FALSE}
hdat=readRDS(hmodFile) 
plot(hdat$H3k4me3,hdat$measured_log2)

```


7. Plot the scatter plot for H3K27me3 vs. expression. [Difficulty:**Beginner**]

**solution:**
```{r,echo=FALSE,eval=FALSE}
plot(hdat$H3k27me3,hdat$measured_log2)
 
```

8. Fit the model for prediction of expression data using: 1) Only H3K4me3 as explanatory variable, 2) Only H3K27me3 as explanatory variable, and 3) Using both H3K4me3 and H3K27me3 as explanatory variables. Inspect the `summary()` function output in each case, which terms are significant. [Difficulty:**Beginner/Intermediate**] 

**solution:**
```{r,echo=FALSE,eval=FALSE}
mod1=lm(measured_log2~H3k4me3,data=hdat)
mod2=lm(measured_log2~H3k27me3,data=hdat)
mod3=lm(measured_log2~H3k4me3+H3k27me3,data=hdat)

summary(mod1)
summary(mod2)
summary(mod3)

```
all terms are significant.

10. Is using H3K4me3 and H3K27me3 better than the model with only H3K4me3? [Difficulty:**Intermediate**] 

**solution:**
Using both histone marks is only slightly better based on R-squared values compared to just using H3K4me3. These two models can be statistically compared using `anova()` function.
