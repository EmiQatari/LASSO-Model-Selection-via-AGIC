This file is for a high dimensional example in model selection via AGIC

############AGIC
AGIC.Sel=function(itr,Data,n,p,l.max=4){
lambda_g=rep(0,itr)
beta_g=matrix(0,nrow=p,ncol=itr)
q=seq(0.001,l.max,by=0.001)
lq=length(q)
AGIC=matrix(0,nrow=itr,ncol=lq)
GD=list(0)
for(ia in 1:itr){
######Generate data
GD[[ia]]=Data[[ia]]
    x <- scale(GD[[ia]]$x,scale=FALSE)
p=ncol(x)
    y <- scale(GD[[ia]]$y,scale=FALSE)
#select lambda

dp=rep(0,lq)
dalpha=rep(0,lq)
AGIC=rep(0,lq)
fit.l<-glmnet(x, y,lambda=q , family="gaussian", alpha=1,intercept=FALSE)
yh= predict(fit.l,s=q, newx=x)
f=coef(fit.l)
for(ib in 1:lq){
bet=matrix(f[,ib],ncol=1)
dp[ib]=beta_r(bet,p)
}
dalpha[1:lq]=dp[lq:1]
for(ic in 1:lq){
sigma=mean((yh[,ic]-y)^2)
AGIC[ic]=sigma+(2*dalpha[ic]/n)
}
lam=q[which.min(AGIC)]
fitt=glmnet(x, y,lambda=lam , family="gaussian", alpha=1,intercept=FALSE)
print(ia)
beta_g[,ia]=as.matrix(fitt$beta)
lambda_g[ia]=q[which.min(AGIC)] # lambda
AGIC[ia,]=AGIC
}
beta_g[abs(beta_g)<0.15]=0
result=list("AGIC"=AGIC,"lambda"=lambda_g,"beta1"=beta_g)
}
############
####
#nonzero beta
beta_r=function(beta,p){
qq=0
for(i in 1:p){
if(abs(beta[i])<=0){
qq=qq+1
}
}
o=p-qq
return(o)
}

######
####Data-Generating

### input
# n : the number of sample
# p : the number of the covariates
# p_real1: the number of the related covariates
# beta_real1: the value of the related covariates
# rho : determines the correlation among covariates.

#output
# a list from observations and coefficient vector of the model

simdat=function(n,p,p_real1,beta_real1,rho){
 # n Number of observations
 # p Number of predictors included in model
real_n=n
#####
beta1=rep(beta_real1,p_real1)
beta01=rep(0,p-p_real1)
beta11=c(beta1,beta01)
beta=beta11
######
CovMatrix <- outer(1:p, 1:p, function(x,y) {rho^abs(x-y)})
x <- mvrnorm(real_n, rep(0,p), CovMatrix)
y <- x%*%beta + rnorm(real_n, 0, 1)

result <- list("y"=y,"x"=x,"beta"=beta)
  return(result)
}
#######
######
#####
###################Run
library("glmnet")
library("MASS")
#########condition-1
itr=500
rho1=0.25
p1=300
p_real1=5
beta_real1=4
n=300
l.max=1
DataSim325300=list(0)
for(i in 1:itr){
set.seed(1371+i)
DataSim325300[[i]]=simdat(n,p1,p_real1,beta_real1,rho1)
}
ntest=200
set.seed(2020)
Data.Test25300=simdat(ntest,p1,p_real1,beta_real1,rho1)
###model-fitting

AGICRUN325300=AGIC.Sel(itr,DataSim325300,n,p1,l.max)

###
######varnum and variable selection
chAGIC325300=list(0)

for(i in 1:itr){
chAGIC325300[[i]]=chosen(AGICRUN325300$beta[,i])
}
varnum(chAGIC325300)

#######MSE-Vectors
MSE.AGIC325300=MSE.AGIC.Rate325300=0
for(i in 1:itr){
MSE.AGIC325300[i]=sum((Data.Test25300$y.t-(Data.Test25300$x.t%*%AGICRUN325300$beta[,i]))^2)/(ntest-length(chAGIC325300[[i]]))
}
mean(MSE.AGIC325300)


#########condition-2
itr=250
rho1=0.25
p2=350
p_real1=5
beta_real1=4
n=300
l.max=1
DataSim325350=list(0)
for(i in 1:itr){
set.seed(1371+i)
DataSim325350[[i]]=simdat(n,p2,p_real1,beta_real1,rho1)
}
ntest=250
set.seed(2020)
Data.Test25350=simdat(ntest,p2,p_real1,beta_real1,rho1)
###model-fitting

AGICRUN325350=AGIC.Sel(itr,DataSim325350,n,p2,l.max)

###
######varnum and variable selection
chAGIC325350=list(0)

for(i in 1:itr){
chAGIC325350[[i]]=chosen(AGICRUN325350$beta[,i])
}
varnum(chAGIC325350)

#######MSE-Vectors
MSE.AGIC325350=MSE.AGIC.Rate325350=0
for(i in 1:itr){
MSE.AGIC325350[i]=sum((Data.Test25350$y.t-(Data.Test25350$x.t%*%AGICRUN325350$beta[,i]))^2)/(ntest-length(chAGIC325350[[i]]))
}
mean(MSE.AGIC325350)


#########condition-3
itr=500
rho1=0.25
p3=400
p_real1=5
beta_real1=4
n=300
l.max=1
DataSim325400=list(0)
for(i in 1:itr){
set.seed(1371+i)
DataSim325400[[i]]=simdat(n,p3,p_real1,beta_real1,rho1)
}
ntest=200
set.seed(2020)
Data.Test25400=simdat(ntest,p3,p_real1,beta_real1,rho1)
###model-fitting

AGICRUN325400=AGIC.Sel(itr,DataSim325400,n,p3,l.max)
###
######varnum and variable selection
chAGIC325400=list(0)

for(i in 1:itr){
chAGIC325400[[i]]=chosen(AGICRUN325400$beta[,i])
}
varnum(chAGIC325400)

#######MSE-Vectors
MSE.AGIC325400=MSE.AGIC.Rate325400=0
for(i in 1:itr){
MSE.AGIC325400[i]=sum((Data.Test25400$y.t-(Data.Test25400$x.t%*%AGICRUN325400$beta[,i]))^2)/(ntest-length(chAGIC325400[[i]]))
}
mean(MSE.AGIC325400)


#########condition-4
itr=500
rho2=0.50
p1=300
p_real1=5
beta_real1=4
n=300
l.max=1
DataSim350300=list(0)
for(i in 1:itr){
set.seed(1371+i)
DataSim350300[[i]]=simdat(n,p1,p_real1,beta_real1,rho2)
}
ntest=200
set.seed(2020)
Data.Test50300=simdat(ntest,p1,p_real1,beta_real1,rho2)
###model-fitting

AGICRUN350300=AGIC.Sel(itr,DataSim350300,n,p1,l.max)

###
######varnum and variable selection
chAGIC350300=list(0)

for(i in 1:itr){
chAGIC350300[[i]]=chosen(AGICRUN350300$beta[,i])
}
varnum(chAGIC350300)


#######MSE-Vectors
MSE.AGIC350300=MSE.AGIC.Rate350300=0
for(i in 1:itr){
MSE.AGIC350300[i]=sum((Data.Test50300$y.t-(Data.Test50300$x.t%*%AGICRUN350300$beta[,i]))^2)/(ntest-length(chAGIC350300[[i]]))
}
mean(MSE.AGIC350300)

#########condition-5
itr=500
rho2=0.50
p2=350
p_real1=5
beta_real1=4
n=300
l.max=1
DataSim350350=list(0)
for(i in 1:itr){
set.seed(1371+i)
DataSim350350[[i]]=simdat(n,p2,p_real1,beta_real1,rho2)
}
ntest=200
set.seed(2020)
Data.Test50350=simdat(ntest,p2,p_real1,beta_real1,rho2)
###model-fitting
AGICRUN350350=AGIC.Sel(itr,DataSim350350,n,p2,l.max)

###
######varnum and variable selection
chAGIC350350=list(0)

for(i in 1:itr){
chAGIC350350[[i]]=chosen(AGICRUN350350$beta[,i])
}
varnum(chAGIC350350)

#######MSE-Vectors
MSE.AGIC350350=MSE.AGIC.Rate350350=0
for(i in 1:itr){
MSE.AGIC350350[i]=sum((Data.Test50350$y.t-(Data.Test50350$x.t%*%AGICRUN350350$beta[,i]))^2)/(ntest-length(chAGIC350350[[i]]))
}
mean(MSE.AGIC350350)

#########condition-6
itr=500
rho2=0.50
p3=400
p_real1=5
beta_real1=4
n=300
l.max=1
DataSim350400=list(0)
for(i in 1:itr){
set.seed(1371+i)
DataSim350400[[i]]=simdat(n,p3,p_real1,beta_real1,rho2)
}
ntest=200
set.seed(1371)
Data.Test50400=simdat(ntest,p3,p_real1,beta_real1,rho2)
###model-fitting
AGICRUN350400=AGIC.Sel(itr,DataSim350400,n,p3,l.max)

###
######varnum and variable selection
chAGIC350400=list(0)
for(i in 1:itr){
chAGIC350400[[i]]=chosen(AGICRUN350400$beta[,i])
}
varnum(chAGIC350400)

#######MSE-Vectors
MSE.AGIC350400=MSE.AGIC.Rate350400=0
for(i in 1:itr){
MSE.AGIC350400[i]=sum((Data.Test50400$y.t-(Data.Test50400$x.t%*%AGICRUN350400$beta[,i]))^2)/(ntest-length(chAGIC350400[[i]]))
}
mean(MSE.AGIC350400)



