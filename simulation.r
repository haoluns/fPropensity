rm(list=ls())
library(locpol)
library(MASS)
source("fPropensity.r")
library(fdapace)
library(glmnet)
library(survey)



 
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
 
 