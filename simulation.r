rm(list=ls())
library(locpol)
library(MASS)




#scen 1
betafun<-function(t) {
6
}

#scen 2
betafun<-function(t) {
4* t^2 - t
}


#scen 3
betafun<-function(t) {
10*t
}

# betafunction
#scen 4
betafun<-function(t) {
	C1=1
	C2=0
	sigma<-0.3
	betafun<-6*(C1*(sin(pi*t)-cos(pi*t)+sin(3*pi*t)-cos(3*pi*t)+sin(5*pi*t)/9-cos(5*pi*t)/9+sin(7*pi*t)/16-cos(7*pi*t)/16+sin(9*pi*t)/25-cos(9*pi*t)/25)+C2*(1/sqrt(2*pi)/sigma)*exp(-(t-0.5)^2/2/sigma^2))		
	
}



simu<-function(){

n<-1000

nrep<-500


measurementerrorvariance<-0.5


h<-calh<-0.01

p<-1/h+1
calp<-p

kkk=10
designpoints<-seq(0,1,calh)
tempbasematrix<-matrix(NA,2*kkk+2,length(designpoints))
for (j in 1:kkk){
	tempbasematrix[j,]<-sqrt(2)*sin(pi*(2*j-1)*designpoints)
}
for (j in (kkk+1): (2*kkk)){
	tempbasematrix[j,]<-sqrt(2)*cos(pi*(2*j-21)*designpoints)
}

tempbasematrix[2*kkk+2-1,]<-1
tempbasematrix[2*kkk+2,]<-designpoints

Sigma<-matrix(0,2*kkk+2,2*kkk+2)
Sigma[1,1]<-1
Sigma[kkk+1,kkk+1]<-1
for (j in 2:kkk){
Sigma[j,j]<-(j-1)^(-2)
}
for (j in (kkk+2):(2*kkk)){
Sigma[j,j]<-(j-(kkk+1))^(-2)	
}

Sigma[2*kkk+2-1,2*kkk+2-1]<-1
Sigma[2*kkk+2,2*kkk+2]<-1


# betafunction
betafun2<-function(t) {
	sigma<-0.3
	t
}

q<-4

gamma<-rep(0.02,q)
alpha<-rep(0.1,q)

correlationmatrix<-matrix(0,nrow(Sigma),q)
correlationmatrix[1,1:q]<-0.1
correlationmatrix[1:q,1]<-0.1

SigmaZ<-0.5^t(sapply(1:q, function(i, j) abs(i-j), 1:q))

Bigsigma<-rbind(cbind(Sigma, correlationmatrix), cbind(t(correlationmatrix), SigmaZ))

mu<-rep(0,nrow(Bigsigma))

cvalue<-seq(0.1,1,0.1)



data<-mvrnorm(n,mu,Bigsigma)

ximatrix<-data[,1:nrow(Sigma)]

rawX<-ximatrix%*%tempbasematrix

calX<-rawX

options(warn=-1)

measurementerror<-matrix(rnorm(n*p,0,sqrt(measurementerrorvariance)),n,p)
observeX<-rawX+measurementerror


Z<-data[,(nrow(Sigma)+1):(nrow(Sigma)+q)]

C1=1
beta<-rep(0,calp)
for (jj in 1:calp) {
beta[jj]<-betafun(designpoints[jj])
}

truebeta<-beta


parameterx<-calX%*%beta*calh


parameter<-(calX%*%beta*calh+Z%*%gamma)
prob = exp(parameter)/(exp(parameter)+1)

treatment = sapply(1:n,function(i){
rbinom(1,1,prob[i])
})

truetheta=1
observey=truetheta*treatment+Z%*%alpha+rnorm(n,mean=0,sd=0.01)





observedmean = sapply(1:n,function(i){
mean(observeX[i,])
})

#plot(observedmean,parameterx)


observed = lapply(1:n,function(i){
observeX[i,]
})

timepoints= lapply(1:n,function(i){
seq(0,1,by=0.01)
})


library(fdapace)
res_pace<- FPCA(observed, timepoints,list(dataType='Dense',error=TRUE, kernel='epan', verbose=TRUE,methodBwCov="GCV",methodBwMu="GCV",nRegGrid=101))

selectedK=min(20,dim(res_pace$phi)[2] )
covarmat = cbind(res_pace$xiEst[,1:selectedK],Z,rep(1,n))

covarmat_mean = cbind(observedmean,Z,rep(1,n))


library(glmnet)
lasso.model =  cv.glmnet(covarmat, treatment,family = "binomial",nfolds = 10,alpha = 1)
coef_cv=coef(lasso.model, s = "lambda.min") 
lasso.prob <- predict(lasso.model,type="response", newx = covarmat, s = 'lambda.min')


data_mean = (cbind(treatment,covarmat_mean))
colnames(data_mean)=c('yy','m','z1','z2','z3','z4','i')
data_mean=as.data.frame(data_mean)
for(i in 1:ncol(data_mean)){
data_mean[,i]=unlist(data_mean[,i])
}
glmmodel = glm(formula = yy ~m+z1+z2+z3+z4+i+0, family = binomial(logit), data = data_mean)
glm.prob=  predict(glmmodel,data_mean,type=c('response'))

sum(abs(prob-glm.prob))
sum(abs(prob-lasso.prob))

MODE="functional"
if(MODE=="mean"){
estprob= glm.prob
}
if(MODE=="functional"){
estprob= lasso.prob
}


weight = sapply(1:n,function(i){
p=estprob[i]
if(treatment[i]==1){
1/p
}else{
1/(1-p)
}
})



indexuntreat = which(treatment == 0)
indextreat = which(treatment == 1)

 
esttheta = weighted.mean(observey[indextreat],weight[indextreat]) - weighted.mean(observey[indexuntreat],weight[indexuntreat]) 
bias = abs(esttheta - truetheta)
percbias = bias / truetheta
df = cbind(observey,treatment,Z)
df=as.data.frame(df)
colnames(df)=c('y','t','z1','z2','z3','z4')

model = lm(y~t+z1+z2+z3+z4+0,data = df, weights = weight)

 
 ci = confint(model, 't', level=0.95)
cicontain = ifelse(truetheta < ci[2] && truetheta > ci[1],1,0)
 
 df$w = weight
 library(survey)
 options(survey.lonely.psu = "adjust")
des<-svydesign(ids=~1, strata=~t, weights=~w, data = df )
 fit3<-svyglm(y~t+z1+z2+z3+z4+0,des, family=gaussian)
 
 se =  SE(fit3)[1]
 
 
  diff1 = mean(df$z1[indextreat] ) - weighted.mean(df$z1[indexuntreat],weight[indexuntreat]) 
  diff2 = mean(df$z2[indextreat] ) - weighted.mean(df$z2[indexuntreat],weight[indexuntreat]) 
  diff3 = mean(df$z3[indextreat] ) - weighted.mean(df$z3[indexuntreat],weight[indexuntreat]) 
  diff4 = mean(df$z4[indextreat] ) - weighted.mean(df$z4[indexuntreat],weight[indexuntreat]) 
  
diff1 = abs(diff1)/ sd(df$z1)
diff2 = abs(diff2)/ sd(df$z2)
diff3 = abs(diff3)/ sd(df$z3)
diff4 = abs(diff4)/ sd(df$z4)

asam = mean(diff1,diff2,diff3,diff4)

list1 = (c(percbias,cicontain,se,asam))


 
 
 
 
 
 
 
 
 ###REPEAT
 
MODE="mean"
if(MODE=="mean"){
estprob= glm.prob
}
if(MODE=="functional"){
estprob= lasso.prob
}


weight = sapply(1:n,function(i){
p=estprob[i]
if(treatment[i]==1){
1/p
}else{
1/(1-p)
}
})



indexuntreat = which(treatment == 0)
indextreat = which(treatment == 1)

 #mean(observey[indextreat] )-mean(observey[indexuntreat])  
#adamodel <- ada(I_AD~., data = dataall, test.x =Test[ , ! colnames(Test) %in% c("I_AD") ], test.y = Test$I_AD, iter = 200, loss = "l", type = "discrete")
 
 #mean(observey[indextreat] ) - mean(observey[indexuntreat]) 
 
 
esttheta = weighted.mean(observey[indextreat],weight[indextreat]) - weighted.mean(observey[indexuntreat],weight[indexuntreat]) 
bias = abs(esttheta - truetheta)
percbias = bias / truetheta
df = cbind(observey,treatment,Z)
df=as.data.frame(df)
colnames(df)=c('y','t','z1','z2','z3','z4')

model = lm(y~t+z1+z2+z3+z4+0,data = df, weights = weight)

 
 ci = confint(model, 't', level=0.95)
cicontain = ifelse(truetheta < ci[2] && truetheta > ci[1],1,0)
 
 df$w = weight
 library(survey)
 options(survey.lonely.psu = "adjust")
des<-svydesign(ids=~1, strata=~t, weights=~w, data = df )
 fit3<-svyglm(y~t+z1+z2+z3+z4+0,des, family=gaussian)
 
 se =  SE(fit3)[1]
 
 
  diff1 = mean(df$z1[indextreat] ) - weighted.mean(df$z1[indexuntreat],weight[indexuntreat]) 
  diff2 = mean(df$z2[indextreat] ) - weighted.mean(df$z2[indexuntreat],weight[indexuntreat]) 
  diff3 = mean(df$z3[indextreat] ) - weighted.mean(df$z3[indexuntreat],weight[indexuntreat]) 
  diff4 = mean(df$z4[indextreat] ) - weighted.mean(df$z4[indexuntreat],weight[indexuntreat]) 
  
diff1 = abs(diff1)/ sd(df$z1)
diff2 = abs(diff2)/ sd(df$z2)
diff3 = abs(diff3)/ sd(df$z3)
diff4 = abs(diff4)/ sd(df$z4)

asam = mean(diff1,diff2,diff3,diff4)

list2 = (c(percbias,cicontain,se,asam))

 
 return(list(mean=list1,func=list2))
 
# adapred <- predict(adamodel, Test[ , ! colnames(Test) %in% c("I_AD") ], type = "prob")
 }
 
 

 
 resultmat1 = matrix(ncol=4)#mean-based
 resultmat2 = matrix(ncol=4)#functional
 
 for(i in 1:100){
 ss=simu()
 resultmat1 = rbind(resultmat1,ss$'mean')
 resultmat2 = rbind(resultmat2,ss$func)
 
 }
 
 
 xx=rbind(apply(resultmat1[2:nrow(resultmat1),],2,mean),
 apply(resultmat2[2:nrow(resultmat2),],2,mean)
 )
 
 row.names(xx)=c('mean-based','functional')
 colnames(xx)= c('Bias', 'CI',  'SD', 'ASAM')                                       
 
 