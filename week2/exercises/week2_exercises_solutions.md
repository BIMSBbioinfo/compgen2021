---
title: 'compgen2021: Week 2 exercises'
author: 'ENTER YOUR NAME HERE'
output:
  pdf_document: default
  pdf: default
---

# Exercises for Week2

For this set of exercises we will be using the expression data shown below:
```{r dataLoadClu,eval=FALSE}
expFile=system.file("extdata",
                    "leukemiaExpressionSubset.rds",
                    package="compGenomRData")
mat=readRDS(expFile)

```

### Clustering

1. We want to observe the effect of data transformation in this exercise. Scale the expression matrix with the `scale()` function. In addition, try taking the logarithm of the data with the `log2()` function prior to scaling. Make box plots of the unscaled and scaled data sets using the `boxplot()` function. [Difficulty: **Beginner/Intermediate**]

**solution:**
```{r,echo=FALSE,eval=FALSE}
boxplot(mat,outline=F)
matscl=scale(mat)
boxplot(matscl,outline=F)
matlog=log2(mat+1)
boxplot(matlog,outline=F)
matlogscl=scale(log2(mat+1))
boxplot(matlogscl,outline=F)
 
```


2. For the same problem above using the unscaled data and different data transformation strategies, use the `ward.d` distance in hierarchical clustering and plot multiple heatmaps. You can try to use the `pheatmap` library or any other library that can plot a heatmap with a dendrogram. Which data-scaling strategy provides more homogeneous clusters with respect to disease types? [Difficulty: **Beginner/Intermediate**]

**solution:**
scaling generally yields more homogeneous clusters. Just taking logs without scaling
is the worst in terms of cluster homogeneity. 
```{r}
library(pheatmap)
# set the leukemia type annotation for each sample
annotation_col = data.frame(
  LeukemiaType =substr(colnames(mat),1,3))
rownames(annotation_col)=colnames(mat)

pheatmap(mat,show_rownames=FALSE,show_colnames=FALSE,cluster_cols = TRUE,
         annotation_col=annotation_col,
         scale = "none",clustering_method="ward.D2",
         clustering_distance_cols="euclidean")

pheatmap(matscl,show_rownames=FALSE,show_colnames=FALSE,cluster_cols = TRUE,
         annotation_col=annotation_col,
         scale = "none",clustering_method="ward.D2",
         clustering_distance_cols="euclidean")

pheatmap(matlog,show_rownames=FALSE,show_colnames=FALSE,cluster_cols = TRUE,
         annotation_col=annotation_col,
         scale = "none",clustering_method="ward.D2",
         clustering_distance_cols="euclidean")

pheatmap(matlogscl,show_rownames=FALSE,show_colnames=FALSE,cluster_cols = TRUE,
         annotation_col=annotation_col,
         scale = "none",clustering_method="ward.D2",
         clustering_distance_cols="euclidean")
```



3. For the transformed and untransformed data sets used in the exercise above, use the silhouette for deciding number of clusters using hierarchical clustering. [Difficulty: **Intermediate/Advanced**]

**solution:**
k=4 seems to be the best clustering regardless of the data transformation.
However, in some cases k=3 and k=5 are close to k=4.
```{r}
library(cluster)

# clustering different transformations with hclust
hcl=hclust(dist(t(mat)))
hcllog=hclust(dist(t(matlog)))
hclscl=hclust(dist(t(matscl)))
hcllogscl=hclust(dist(t(matlogscl)))


# silhouette can be used in any clustering result
# we use it with hclust & cuttree
Ks=sapply(2:7,
          function(i){ 
            summary(silhouette(cutree(hcl,k=i),dist(t(mat))))$avg.width})

Kslog=sapply(2:7,
          function(i){ 
            summary(silhouette(cutree(hcllog,k=i),dist(t(matlog))))$avg.width})


Ksscl=sapply(2:7,
          function(i){ 
            summary(silhouette(cutree(hclscl,k=i),dist(t(matscl))))$avg.width})


Kslogscl=sapply(2:7,
          function(i){ 
            summary(silhouette(cutree(hcllogscl,k=i),dist(t(matlogscl))))$avg.width})

# plotting average silhouette scores per k and per data transformation
plot(2:7,Ks,xlab="k",ylab="av. silhouette",type="b",
     pch=19)
plot(2:7,Ksscl,xlab="k",ylab="av. silhouette",type="b",
     pch=19)
plot(2:7,Kslog,xlab="k",ylab="av. silhouette",type="b",
     pch=19)
plot(2:7,Kslogscl,xlab="k",ylab="av. silhouette",type="b",
     pch=19)



```



4. Now, use the Gap Statistic for deciding the number of clusters in hierarchical clustering. Is the same number of clusters identified by two methods? Is it similar to the number of clusters obtained using the k-means algorithm in the unsupervised learning chapter. [Difficulty: **Intermediate/Advanced**]

**solution:**
Different number of clusters are identified. 

```{r,echo=TRUE,eval=TRUE}
set.seed(101)

hfun <- function(x,k){
list(cluster = cutree(hclust(dist(x,method = "euclidean"), method="ward.D"),k = k))
}

## unscaled data:

gapmatlog= clusGap(t(matlog), FUN = hfun, K.max = 8,B=50)
plot(gapmatlog, main = "Gap statistic for the 'Unscaled Leukemia' data")


## scaled data:
gapscmat= clusGap(t(matlogscl), FUN = hfun, K.max = 8,B=50)
plot(gapscmat, main = "Gap statistic for the 'Scaled Leukemia' data")
 
```


### Dimension reduction
We will be using the leukemia expression data set again. You can use it as shown in the clustering exercises.

1. Do PCA on the expression matrix using the `princomp()` function and then use the `screeplot()` function to visualize the explained variation by eigenvectors. How many top components explain 95% of the variation? [Difficulty: **Beginner**]

**solution:**
First 25 components explain 95% of the variation
```{r,echo=TRUE,eval=TRUE}
pr = princomp(scale(mat))

# Scree plot
screeplot(pr, main = "Screeplot showing first 10 PCs.",npcs = 10)

summary(pr)
 
```


2. Our next tasks are removing the eigenvectors and reconstructing the matrix using SVD, then we need to calculate the reconstruction error as the difference between the original and the reconstructed matrix. HINT: You have to use the `svd()` function and equalize eigenvalue to $0$ for the component you want to remove. [Difficulty: **Intermediate/Advanced**]

**solution:**
I'm removing the 2nd eigenvector only. The question is not clear on which
eigenvectors to remove so I picked one and removed it.
```{r}
s=svd(matlogscl)
D <- diag(s$d)
mat.rec1 = s$u %*% D %*% t(s$v)
D[2,2] <- 0 # remove 2nd eigenvector by setting its eigenvalue to 0
mat.rec2 = s$u %*% D %*% t(s$v)
reconstruction.error = sum((mat.rec1 -  mat.rec2)^2)


# remove only 4th eigenvector
s=svd(matlogscl)
D <- diag(s$d)
mat.rec1 = s$u %*% D %*% t(s$v)
D[4,4] <- 0 # remove 2nd eigenvector by setting its eigenvalue to 0
mat.rec2 = s$u %*% D %*% t(s$v)
sum((mat.rec1 -  mat.rec2)^2)
```
removing 4th eigenvector does less damage than removing 2nd eigenvector
,which is expected


3. Produce a 10-component ICA from the expression data set. Remove each component and measure the reconstruction error without that component. Rank the components by decreasing reconstruction-error. [Difficulty: **Advanced**]

**solution:**
component 6 is the most important component, removing that creates the highest
reconstruction error
```{r,echo=FALSE,eval=FALSE}
set.seed(101)
library(fastICA)
ica.res=fastICA(t(mat),n.comp=10) # apply ICA

 
temp=ica.res # save results to a temporary variable
rec.err=c() # reconstruction errors measured by RMSE

# iterate over 10 components
for( i in 1:10){
temp$S[,i]=0 # remove the component i

rec.err[i]=sqrt(mean((temp$X- (temp$S %*% temp$A))^2)) # RMSE calculation
temp=ica.res # recreate the temp to remove the next component in the loop
}
plot(rec.err)
order(rec.err)
```



4. In this exercise we use the `Rtsne()` function on the leukemia expression data set. Try to increase and decrease perplexity t-sne, and describe the observed changes in 2D plots. [Difficulty: **Beginner**]

**solution:**
Among the values we tried, visually seperation is better for 10 and above
perplexity values. Smaller values results in less visual seperation between subtypes.
```{r,echo=FALSE,eval=FALSE}
library("Rtsne")
annotation_col = data.frame(
  LeukemiaType =substr(colnames(mat),1,3))
rownames(annotation_col)=colnames(mat)
set.seed(42) # Set a seed if you want reproducible results
tsne_out <- Rtsne(t(mat),perplexity = 1) # Run TSNE
 #image(t(as.matrix(dist(tsne_out$Y))))
# Show the objects in the 2D tsne representation
plot(tsne_out$Y,col=as.factor(annotation_col$LeukemiaType),
     pch=19)
 
tsne_out <- Rtsne(t(mat),perplexity = 3) # Run TSNE
 #image(t(as.matrix(dist(tsne_out$Y))))
# Show the objects in the 2D tsne representation
plot(tsne_out$Y,col=as.factor(annotation_col$LeukemiaType),
     pch=19)
     
tsne_out <- Rtsne(t(mat),perplexity = 5) # Run TSNE
 #image(t(as.matrix(dist(tsne_out$Y))))
# Show the objects in the 2D tsne representation
plot(tsne_out$Y,col=as.factor(annotation_col$LeukemiaType),
     pch=19)
     

tsne_out <- Rtsne(t(mat),perplexity = 10) # Run TSNE
 #image(t(as.matrix(dist(tsne_out$Y))))
# Show the objects in the 2D tsne representation
plot(tsne_out$Y,col=as.factor(annotation_col$LeukemiaType),
     pch=19)
     
tsne_out <- Rtsne(t(mat),perplexity = 100) # Run TSNE
 #image(t(as.matrix(dist(tsne_out$Y))))
# Show the objects in the 2D tsne representation
plot(tsne_out$Y,col=as.factor(annotation_col$LeukemiaType),
     pch=19)
     
```



