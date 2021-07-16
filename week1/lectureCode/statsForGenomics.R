

#--- Central Tendency and spread --- 

x=rnorm(20,mean=6,sd=0.7)

mean(x)
median(x)

var(x)
sd(x)
IQR(x)






#get the probability of P(X =< -2) | mean=0 and sd=2 
pnorm(-2, mean=0, sd=2)

#get the probability of P(X > -2) | mean=0 and sd=2 
pnorm(-2, mean=0, sd=2,lower.tail = FALSE)

#get 5 random numbers from norm. dist. | mean=0 and sd=2
rnorm(5, mean=0 , sd=2)


#--- precision of estimates ----

library(mosaic)
set.seed(21)
sample1= rnorm(50,mean=20,sd=5) # simulate a sample

# do bootstrap resampling, sampling with replacement
boot.means=do(1000) * mean(resample(sample1))

# get percentiles from the bootstrap means
q=quantile(boot.means[,1],p=c(0.025,0.975))

# plot the histogram
hist(boot.means[,1],col="cornflowerblue",border="white",
     xlab="sample means")
abline(v=c(q[1], q[2] ),col="red")
text(x=q[1],y=200,round(q[1],3),adj=c(1,0))
text(x=q[2],y=200,round(q[2],3),adj=c(0,0))



alpha=0.05 # 95% confidence interval
sd=5
n=50
mean(sample1)+qnorm(c(alpha/2,1-alpha/2))*sd/sqrt(n)

#--- hypothesis testing ----

set.seed(100)
gene1=rnorm(30,mean=4,sd=2)
gene2=rnorm(30,mean=2,sd=2)
org.diff=mean(gene1)-mean(gene2)
gene.df=data.frame(exp=c(gene1,gene2),
                   group=c( rep("test",30),rep("control",30) ) )


exp.null <- do(1000) * diff(mosaic::mean(exp ~ shuffle(group), data=gene.df))

hist(exp.null[,1],xlab="null distribution | no difference in samples",
     main=expression(paste(H[0]," :no difference in means") ),
     xlim=c(-2,2),col="cornflowerblue",border="white")
abline(v=quantile(exp.null[,1],0.95),col="red" )
abline(v=org.diff,col="blue" )
text(x=quantile(exp.null[,1],0.95),y=200,"0.05",adj=c(1,0),col="red")
text(x=org.diff,y=200,"org. diff.",adj=c(1,0),col="blue")

p.val=sum(exp.null[,1]>org.diff)/length(exp.null[,1])
p.val


# Welch's t-test
stats::t.test(gene1,gene2)


# t-test with equal variance assumption
stats::t.test(gene1,gene2,var.equal=TRUE)


#--- mult. test correction ----

library(qvalue)
data(hedenfalk)

qvalues = qvalue(hedenfalk$p)$q
bonf.pval=p.adjust(hedenfalk$p,method ="bonferroni")
fdr.adj.pval=p.adjust(hedenfalk$p,method ="fdr")

plot(hedenfalk$p,qvalues,pch=19,ylim=c(0,1),
     xlab="raw P-values",ylab="adjusted P-values")
points(hedenfalk$p,bonf.pval,pch=19,col="red")
points(hedenfalk$p,fdr.adj.pval,pch=19,col="blue")
legend("bottomright",legend=c("q-value","FDR (BH)","Bonferroni"),
       fill=c("black","blue","red"))


#--- moderated t-test ----

set.seed(100)

#sample data matrix from normal distribution

gset=rnorm(3000,mean=200,sd=70)
data=matrix(gset,ncol=6)

# set groups
group1=1:3
group2=4:6
n1=3
n2=3
dx=rowMeans(data[,group1])-rowMeans(data[,group2])

require(matrixStats)

# get the estimate of pooled variance 
stderr <- sqrt( (rowVars(data[,group1])*(n1-1) + rowVars(data[,group2])*(n2-1)) / (n1+n2-2) * ( 1/n1 + 1/n2 ))

# do the shrinking towards median
mod.stderr <- (stderr + median(stderr)) / 2 # moderation in variation

# estimate t statistic with moderated variance
t.mod = dx / mod.stderr

# calculate P-value of rejecting null 
p.mod = 2*pt( -abs(t.mod), n1+n2-2 )

# esimate t statistic without moderated variance
t = dx / stderr

# calculate P-value of rejecting null 
p = 2*pt( -abs(t), n1+n2-2 )

par(mfrow=c(1,2))
hist(p,col="cornflowerblue",border="white",main="",xlab="P-values t-test")
mtext(paste("signifcant tests:",sum(p<0.05))  )
hist(p.mod,col="cornflowerblue",border="white",main="",xlab="P-values mod. t-test")
mtext(paste("signifcant tests:",sum(p.mod<0.05))  )


#---- regression ----



# set random number seed, so that the random numbers from the text
# is the same when you run the code.
set.seed(32)

# get 50 X values between 1 and 100
x = runif(50,1,100)

# set b0,b1 and variance (sigma)
b0 = 100
b1 = 2
sigma = 20
# simulate error terms from normal distribution
eps = rnorm(50,0,sigma)
# get y values from the linear equation and addition of error terms
y = b0 + b1*x+ eps

mod1=lm(y~x)

# plot the data points
plot(x,y,pch=20,
     ylab="Gene Expression",xlab="Histone modification score")
# plot the linear fit
abline(mod1,col="blue")

mod1


summary(mod1)

# get confidence intervals 
confint(mod1)

# pull out coefficients from the model
coef(mod1)



