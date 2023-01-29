library("glmnet")
library("MASS")
library(Matrix)
library("mvtnorm")

### m: number of relative groups
#### noise: determines signal to noise ratio


######FUNCTIONS
####Generating simulated Data
### simdat generates simulation datasets which have grouped correlated features and m is the number of the related groups. Note that each group has 5 features. n is the number of the observations, rho is the correlation among each group and noise determines the signal to noise ration of the model

simdat=function(n,rho,p,m,noise){
k=p/5
c1=c(1,rho,rho,rho,rho)
c2=c(rho,1,rho,rho,rho)
c3=c(rho,rho,1,rho,rho)
c4=c(rho,rho,rho,1,rho)
c5=c(rho,rho,rho,rho,1)
C=cbind(c1,c2,c3,c4,c5)
M=list(0)
for(i in 1:k){
M[[i]]=C
}
	R=bdiag(M)
	mu1=0
	mu1[1:p]=0
ki=sqrt(((1/noise)-1)/(5*(1+4*rho)*(1-(0.25)^m)))
pf=c(ki,ki,ki,ki,ki)
pr=matrix(0,ncol=k,nrow=5)
for(i in 0:(k-1)){
pr[,(i+1)]=pf/(2^i)
}

	x=mvrnorm(n,mu1,R)
	beta1=as.vector(pr[,1:m])
	z1=0
	l=length(beta1)
	if(l<p){
	z1[1:(p-l)]=0
	beta=c(beta1,z1)
	beta=as.matrix(beta)
	}
	else{
	beta=beta1
	}
	m2=t(beta)%*%t(x)
	eps=rnorm(n,0,1)
	Yp=m2+eps
	Yp=as.vector(Yp)
	xtr=x[1:n,]
	ytr=Yp[1:n]
	return(list("ki"=ki,"beta"=as.vector(beta),x.train=(xtr),y.train=ytr))
}

###########nonzero beta
##### beta_r calculates the number of the selected features for each model.
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
#####
####### AGIC.Sel perform LASSO regression using AGIC. It needs a list of generated data using simdat and n & p are defined in simdat. itr is the number of iterations for simulated data and l.max is the maximum value of lambda used for glmnet function and other related calculations.
#ouputs : selected beta by AGIC, the values of AGIC for the optimum models and selected lambdas by AGIC.


AGIC.Sel<-function(itr,Data,n,p,l.max=1){
lambda_g=rep(0,itr)
q=seq(0.001,l.max,by=0.001)
lq=length(q)
beta_g=matrix(0,nrow=p,ncol=itr)
GIC=matrix(0,nrow=itr,ncol=lq)
GD=list(0)
for(ia in 1:itr){
######Generate data
GD[[ia]]=Data[[ia]]
    x <- scale(GD[[ia]]$x.train,scale=FALSE)
    y <- scale(GD[[ia]]$y.train,scale=FALSE)
	#select lambda
	dp=rep(0,lq)
	dalpha=rep(0,lq)
	gic=rep(0,lq)
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
		gic[ic]=1*(sigma)+(2*dalpha[ic]/n)
	
	}
	lam=q[which.min(gic)]
	fitt=glmnet(x, y,lambda=lam , family="gaussian", alpha=1,intercept=FALSE)
	print(ia)
	beta_g[,ia]=as.matrix(fitt$beta)
	lambda_g[ia]=q[which.min(gic)] # lambda
	GIC[ia,]=gic
}
beta_g[abs(beta_g)<0.05]=0
result=list("Agic"=GIC,"lambda"=lambda_g,"beta.AGIC"=beta_g)
}
########example

####
itr=500;rho=0.25;p=40;n=400
DataSim42540=list(0)
for(i in 1:itr){
set.seed(i+2022)
DataSim42540[[i]]=simdat(n,rho,p,m,noise)
}
GICRUN42540=AGIC.Sel(itr,DataSim42540,n,p,l.max=1)


#################MSE-Search-Path
####spath function provides a new sorting among the present features in the model. inputs are 
# beta: the estimation of coefficient vector
# Dat: a list of data involved in model selection that beta is obtained by Dat.

spath=function(beta,Dat){
y=Dat$y
x=Dat$x
m=length(beta)
o=1:m
h=1:m
a1=0
d=0
for(i in 1:m)
d[i]=sum((y-(x[,i]*beta[i]))^2)

path=which.min(d)
o=o[-which.min(d)]
while(length(path)<p){
a1=0
for(j in o){
a1[j]=sum((y-(x[,c(path,j)]%*%beta[c(path,j)]))^2)
}
a1<- a1[!is.na(a1)]
path=c(path,o[which.min(a1)])
o=o[-which.min(a1)]
}
return(path)
}

###Example:
p=40
itr=500
sPathMatrixgic42540=matrix(0,ncol=itr,nrow=p)
for(i in 1:itr){
sPathMatrixgic42540[,i]=spath(GICRUN42540$beta[,i],DataSim42540[[i]])
}

####################
####chosen function shows which features is selected by a method. The input is a vector consisting of the estimation of the coefficient vector obtained by any model selection method.

chosen=function(b){
m=length(b)
p=1:m
o=1:m
for(i in 1:m){
if(b[i]==0){
p[i]=0
}
}
chosen=p[p!=0]
return("chosen"=chosen)
}
######
####Mean of logarithm of the predictive density
MLPDsel=function(Dat,beta){
xt=Dat$x.test
yt=Dat$y.test
chosen1=chosen(beta)
x=xt[,chosen1]
p=ncol(x)
w=beta[chosen1]
m=x %*%w
mpred<-ginv(t(x)%*%x+diag(p))%*%t(x)%*%m
sdpred=diag(sd(mpred)*(t(x)%*%x+diag(p)))
pd <- dnorm(yt,m,sd=sdpred)
mlpd <- mean(log((pd)))
return((mlpd=mlpd))
}
###################################################################################
######Example for MLPD Of Selected Models
##
p=40
rho=0.25
n=400
betareal=0
betareal=simdat(n,rho,p,m,alpha)$beta

MLPDSELGic42540=0
for(i in 2:itr){
MLPDSELGic42540[i]=MLPDsel(GICRUN42540$Data[[i]],GICRUN42540$beta[,i])-MLPDsel(GICRUN42540$Data[[i]],betareal)
}
mean(MLPDSELGic42540)