#setwd("~/Desktop/ranking selection/code/new code")
source("3M.funcs.R")
source("datagen.R")
source("nest.func.R")
source('rih.R')
source('rsor.R')
source('rsorgen.R')
library(CVXR)

m<-10000   #number of hypothesis
nrep<-10#number of reptitions
q<-0.1 
epsilon=10^(-5)
mu0=6
#some parameter to vary
sig.vec<-seq(from=1.5, to=2.5, by=0.2) 
np<-length(sig.vec)


dd.fdr<-rep(0, np)
rs.fdr<-rep(0, np)
or2.fdr<-rep(0, np)

dd.etp<-rep(0, np)
rs.etp<-rep(0, np)
or2.etp<-rep(0, np)


dd.etp.old<-rep(0, np)
rs.etp.old<-rep(0, np)
or2.etp.old<-rep(0, np)

dd.fdp<-matrix(rep(0, nrep*np), np, nrep)
rs.fdp<-matrix(rep(0, nrep*np), np, nrep)
or2.fdp<-matrix(rep(0, nrep*np), np, nrep)

dd.ntp<-matrix(rep(0, nrep*np), np, nrep)
rs.ntp<-matrix(rep(0, nrep*np), np, nrep)
or2.ntp<-matrix(rep(0, nrep*np), np, nrep)


dd.ntp2<-matrix(rep(0, nrep*np), np, nrep)
rs.ntp2<-matrix(rep(0, nrep*np), np, nrep)
or2.ntp2<-matrix(rep(0, nrep*np), np, nrep)

dd.ntp.old<-matrix(rep(0, nrep*np), np, nrep)
rs.ntp.old<-matrix(rep(0, nrep*np), np, nrep)
or2.ntp.old<-matrix(rep(0, nrep*np), np, nrep)
for (i in 1:np){
  

  for (j in 1:nrep) {
    mu1=rnorm(5000,5,0.5)
    mu2=rnorm(5000,7,0.5)
    mu=c(mu1,mu2)
    sig=c(rep(1,5000),rep(sig.vec[i],5000))
    oldtime=Sys.time()
    cat("j=",j,"i=",i)

    d=rdata.rs.func(m,mu,sig,mu0)
    x=d$x
    theta=d$theta
   
    fnull1<- function(z,s) integrate(function(x,z) dnorm(x,5,0.5)*dnorm(z,x,s),0,mu0,z)$value
    fmar1<- function(z,s) integrate(function(x,z) dnorm(x,5,0.5)*dnorm(z,x,s),0,12,z)$value
    
    num1=mapply(fnull1,x[1:5000],sig[1:5000])
    denom1=mapply(fmar1,x[1:5000],sig[1:5000])
    clfdr1=num1/denom1
    
    fnull2<- function(z,s) integrate(function(x,z) dnorm(x,7,0.5)*dnorm(z,x,s),0,6,z)$value
    fmar2<- function(z,s) integrate(function(x,z) dnorm(x,7,0.5)*dnorm(z,x,s),0,14,z)$value
    
    num2=mapply(fnull2,x[5001:10000],sig[5001:10000])
    denom2=mapply(fmar2,x[5001:10000],sig[5001:10000])
    clfdr2=num2/denom2

    clfdr.raw=c(clfdr1,clfdr2)
    clfdr=pmin(clfdr.raw,1)
    #########################
   
      tst.dd1=rsstathomo.func(x[1:5000],sig[1:5000],mu0,gd=50)
      tst.dd2=rsstathomo.func(x[5001:10000],sig[5001:10000],mu0,gd=50)
     
      clfdr.dd=c(tst.dd1$clfdr,tst.dd2$clfdr)
      fx.dd=x-mu0
      d.dd=rsorddgen.func(fx.dd,clfdr.dd,q)
      dd.fdp[i, j]<-sum((1-theta)*d.dd$de)/max(sum(d.dd$de), 1)
      dd.ntp[i, j]<-sum((x-mu0)*d.dd$de)
      dd.ntp.old[i,j]=sum((theta)*d.dd$de)/max(sum(theta), 1)
   
    ###########################

  
    drs=rsorddgen.func(x-mu0,clfdr,q)
    rs.fdp[i, j]<-sum((1-theta)*drs$de)/max(sum(drs$de), 1)
    rs.ntp[i, j]<-sum((x-mu0)*drs$de)
    rs.ntp.old[i,j]=sum((theta)*drs$de)/max(sum(theta), 1)
    
  
    drs2=sc.func(clfdr,q)
    or2.fdp[i, j]<-sum((1-theta)*drs2$de)/max(sum(drs2$de), 1)
    or2.ntp[i, j]<-sum((x-mu0)*drs2$de)
    or2.ntp.old[i,j]=sum((theta)*drs2$de)/max(sum(theta), 1)
    
 

    
    newtime=Sys.time()
    print(newtime-oldtime)
  }#end of j loop
  dd.fdr[i]<-mean(dd.fdp[i,])
  rs.fdr[i]<-mean(rs.fdp[i,])
  or2.fdr[i]<-mean(or2.fdp[i,])
  
  dd.etp[i]<-mean(dd.ntp[i,])
  rs.etp[i]<-mean(rs.ntp[i,])
  or2.etp[i]<-mean(or2.ntp[i,])
  
  
  dd.etp.old[i]<-mean(dd.ntp.old[i,])
  rs.etp.old[i]<-mean(rs.ntp.old[i,])
  or2.etp.old[i]<-mean(or2.ntp.old[i,])
  
  
}#end of i loop

fdr1<-cbind(dd.fdr, rs.fdr,or2.fdr)
etp1<-cbind(dd.etp, rs.etp,or2.etp)
etp.old<-cbind(dd.etp.old,rs.etp.old,or2.etp.old)



par(mfrow=c(1, 3), mgp=c(2, 0.5, 0), mar=c(3, 3, 2, 1)+0.1)
if(1==1){ 
 
  matplot(sig.vec, fdr1, type="o", pch=1:5, lwd=2,col=c('black','red','blue'), main=expression(paste("a)FDR Comparison varing ", sigma)), xlab=expression(sigma), ylim=c(-0.01, 0.20), ylab="FDR")
  
  legend("top", c("DD", "OR","Clfdr"), pch=1:5, col=c('black','red','blue'), lwd=2)
  matplot(sig.vec, etp1, type="o", pch=1:5, lwd=2, col=c('black','red','blue'),main=expression(paste("b)" ,"ETP* Comparison varing ", sigma)), xlab=expression(sigma), ylab="ETP*",ylim = c(0,10000))
 
  matplot(sig.vec, etp.old*5000, type="o", pch=1:5, lwd=2,col=c('black','red','blue'), main=expression(paste("c)ETP Comparison varing ", sigma)), xlab=expression(sigma), ylab="ETP",ylim = c(0,5000))
 
}


