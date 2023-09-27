source("3M.funcs.R")
source("datagen.R")
source("nest.func.R")
source('rih.R')
source('rsorgen.R')
library(CVXR)
library(invgamma)
#install.packages("quadprog")
#install.packages('sgof')
#install.packages('invgamma')

#library(quadprog)

set.seed(92)
m<-5000   #number of hypothesis
nrep<-100#number of reptitions
q<-0.1 
mu0=6
epsilon=10^(-5)
#some parameter to vary
sig.vec<-seq(from=2, to=4, by=0.2) 

np<-length(sig.vec)


dd.fdr<-rep(0, np)
rs.fdr<-rep(0, np)
bh.fdr=rep(0,np)
rs2.fdr=rep(0,np)

dd.etp<-rep(0, np)
rs.etp<-rep(0, np)
rs2.etp<-rep(0, np)
bh.etp<-rep(0, np)
diff1.etp=rep(0,np)
diff2.etp=rep(0,np)
diff3.etp=rep(0,np)
diff4.etp=rep(0,np)

dd.etp.old<-rep(0, np)
rs.etp.old<-rep(0, np)
rs2.etp.old<-rep(0, np)
bh.etp.old<-rep(0, np)
diff1.etp.old<-rep(0, np)
diff2.etp.old<-rep(0, np)
diff3.etp.old<-rep(0, np)
diff4.etp.old<-rep(0, np)

dd.fdp<-matrix(rep(0, nrep*np), np, nrep)
rs.fdp<-matrix(rep(0, nrep*np), np, nrep)
rs2.fdp<-matrix(rep(0, nrep*np), np, nrep)
bh.fdp<-matrix(rep(0, nrep*np), np, nrep)

dd.ntp<-matrix(rep(0, nrep*np), np, nrep)
rs.ntp<-matrix(rep(0, nrep*np), np, nrep)
rs2.ntp<-matrix(rep(0, nrep*np), np, nrep)
bh.ntp<-matrix(rep(0, nrep*np), np, nrep)
diff1.ntp<-matrix(rep(0, nrep*np), np, nrep)
diff2.ntp<-matrix(rep(0, nrep*np), np, nrep)
diff3.ntp<-matrix(rep(0, nrep*np), np, nrep)
diff4.ntp<-matrix(rep(0, nrep*np), np, nrep)

dd.ntp.old<-matrix(rep(0, nrep*np), np, nrep)
rs.ntp.old<-matrix(rep(0, nrep*np), np, nrep)
rs2.ntp.old<-matrix(rep(0, nrep*np), np, nrep)
bh.ntp.old<-matrix(rep(0, nrep*np), np, nrep)
diff1.ntp.old<-matrix(rep(0, nrep*np), np, nrep)
diff2.ntp.old<-matrix(rep(0, nrep*np), np, nrep)
diff3.ntp.old<-matrix(rep(0, nrep*np), np, nrep)
diff4.ntp.old<-matrix(rep(0, nrep*np), np, nrep)

for (i in 1:np){
  
  
  for (j in 1:nrep) {
    oldtime=Sys.time()
    cat("j=",j,"i=",i)
    sig=runif(m,0.5,sig.vec[i])
    d=rdata.unif.func4(10,0,mu0,sig,m)
    x=d$x
    theta=d$theta
    fnull<- function(z,s) integrate(function(x,z) 6*dunif(x,0,6)*dnorm(z,x,s),0,6,z)$value
    falt<- function(z,s) integrate(function(x,z) 10*dunif(x,0,10)*dnorm(z,x,s),0,10,z)$value
    
    num=mapply(fnull,x,sig)
    denom=mapply(falt,x,sig)
    clfdr=num/denom
    
    #generate data (based on i only)
    z=(x-mu0)/sig
    theta=d$theta
    pv=1-pnorm(z, 0, 1)
    
    
    
    
    #########################
 
      d.dd=rih.func(x,sig,mu0,q=q,gd=50)
      dd.fdp[i, j]<-sum((1-theta)*d.dd$de)/max(sum(d.dd$de), 1)
      dd.ntp[i, j]<-sum((x-mu0)*d.dd$de)
      dd.ntp.old[i,j]=sum((theta)*d.dd$de)/max(sum(theta), 1)
    
   
    drs=rsorddgen.func(x-mu0,clfdr,q)
    rs.fdp[i, j]<-sum((1-theta)*drs$de)/max(sum(drs$de), 1)
    rs.ntp[i, j]<-sum((x-mu0)*drs$de)
    rs.ntp.old[i,j]=sum((theta)*drs$de)/max(sum(theta), 1)
  
    drs2=sc.func(clfdr,q)
    rs2.fdp[i, j]<-sum((1-theta)*drs2$de)/max(sum(drs2$de), 1)
    rs2.ntp[i, j]<-sum((x-mu0)*drs2$de)
    rs2.ntp.old[i,j]=sum((theta)*drs2$de)/max(sum(theta), 1)
   
    bh.res<-bh.func(pv, q)
    bh.de<-bh.res$de
    bh.fdp[i, j]<-sum((1-theta)*bh.de)/max(sum(bh.de), 1)
    bh.ntp[i, j]<-sum((x-mu0)*bh.de)
    bh.ntp.old[i, j]<-sum((theta)*bh.de)/max(sum(theta), 1)
    
    ind=which(d.dd$de==1&drs2$de==0)
    ind2=which(d.dd$de==0&drs2$de==1)
    ind3=which(drs$de==1&drs2$de==0)
    ind4=which(drs$de==0&drs2$de==1)
    
    diff1.ntp[i,j]=sum(x[ind]-mu0)
    diff2.ntp[i,j]=sum(x[ind2]-mu0)
    diff3.ntp[i,j]=sum(x[ind3]-mu0)
    diff4.ntp[i,j]=sum(x[ind4]-mu0)
    
    diff1.ntp.old[i,j]=length(ind)
    diff2.ntp.old[i,j]=length(ind2)
    diff3.ntp.old[i,j]=length(ind3)
    diff4.ntp.old[i,j]=length(ind4)
    
    newtime=Sys.time()
    print(newtime-oldtime)
  }#end of j loop
  dd.fdr[i]<-mean(dd.fdp[i,])
  rs.fdr[i]<-mean(rs.fdp[i,])
  rs2.fdr[i]<-mean(rs2.fdp[i,])
  bh.fdr[i]<-mean(bh.fdp[i,])
  
  dd.etp[i]<-mean(dd.ntp[i,])
  rs.etp[i]<-mean(rs.ntp[i,])
  rs2.etp[i]<-mean(rs2.ntp[i,])
  bh.etp[i]<-mean(bh.ntp[i,])
  diff1.etp[i]=mean(diff1.ntp[i,])
  diff2.etp[i]=mean(diff2.ntp[i,])
  diff3.etp[i]=mean(diff3.ntp[i,])
  diff4.etp[i]=mean(diff4.ntp[i,])
  
  dd.etp.old[i]<-mean(dd.ntp.old[i,])
  rs.etp.old[i]<-mean(rs.ntp.old[i,])
  rs2.etp.old[i]<-mean(rs2.ntp.old[i,])
  bh.etp.old[i]<-mean(bh.ntp.old[i,])
  diff1.etp.old[i]=mean(diff1.ntp.old[i,])
  diff2.etp.old[i]=mean(diff2.ntp.old[i,])
  diff3.etp.old[i]=mean(diff3.ntp.old[i,])
  diff4.etp.old[i]=mean(diff4.ntp.old[i,])
  
  
}#end of i loop

fdr1<-cbind(dd.fdr, rs.fdr,bh.fdr,rs2.fdr)
etp1<-cbind(dd.etp, rs.etp,bh.etp,rs2.etp)
etp1.old<-cbind(dd.etp.old, rs.etp.old,bh.etp.old,rs2.etp.old)


par(mfrow=c(1, 3), mgp=c(2, 0.5, 0), mar=c(3, 3, 2, 1)+0.1)
if(T){
  
  
  matplot(sig.vec, fdr1, type="o", pch=1:5, lwd=2, main=expression(paste("a)FDR Comparison varing ", sigma[max])), xlab=expression(sigma[max]), ylim=c(-0.01, 0.20), ylab="FDR")

  legend("top", c("DD", "OR", "BH","Clfdr"), pch=1:4, col=1:4, lwd=2)
  matplot(sig.vec, etp1, type="o", pch=1:5, lwd=2, main=expression(paste("b)" ,"ETP* Comparison varing ", sigma[max])), xlab=expression(sigma[max]), ylab="ETP*",ylim = c(1000,5800))
 
  matplot(sig.vec, etp1.old*m*0.4, type="o", pch=1:5, lwd=2, main=expression(paste("c)ETP Comparison varing ", sigma[max])), xlab=expression(sigma[max]), ylab="ETP",ylim = c(0,2000))
 
}




if(T){
  ind=which(d.dd$de==1&drs2$de==0)
  ind2=which(d.dd$de==0&drs2$de==1)
  ind3=which(d.dd$de==1&drs2$de==1)
ptcolor=rep('gray',length(x))
ptcolor[ind]='red'
ptcolor[ind2]='blue'
ptcolor[ind3]='green'
par(mfrow=c(1, 2), mgp=c(2, 0.5, 0), mar=c(3, 3, 2, 1)+0.1)
ptype=rep(1,length(x))

ptype[which(ptcolor!='gray')]=3
plot(x,sig,col='gray',pch=1,cex=2,ylab=expression(sigma),main=expression(paste("a)" ,"Rejection by DD and Clfdr ", sigma[max],"=4")))

points(x[c(ind3,ind2,ind)],sig[c(ind3,ind2,ind)],col=ptcolor[c(ind3,ind2,ind)],pch=18,cex=2)

etp.diff=cbind(diff1.etp,diff2.etp)
matplot(sig.vec, etp.diff, type="o", pch=1:5,col=c('red','blue'),  lwd=2, main=expression(paste("b)" ,"ETP* Comparison varing ", sigma[max])), xlab=expression(sigma[max]), ylab="ETP*",ylim = c(0,400))

legend("top", c("DD", "Clfdr"), pch=1:4, col=c('red','blue'), lwd=2)

}



