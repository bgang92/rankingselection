rih.func=function(x,sig,mu0,q=0.1,epsilon=10^(-3),gd=100,mod=2,jacknife=T,mu=NULL,adjust=T){
  #This function implements Algorithm 1 in the paper
  #q is the target FDR level
  #gd is the grid size
  #mod=2 means use convolution to estimate the marginal density
  #mod=1 means use bivariate kernel to estimate the marginal density
  #return varaible:
  #Clfdr is the estimated CLfdr
  #t is the estimated t statistics
  #de is the decision
  
  #step 1: pilot density estimate 
  nest.res=nest.func(x,sig,jacknife=jacknife)
  dens.dd=nest.res$dens
  #step 2: Convex optimizaiton to find w
  west=west.func(dens.dd,x,sig,gdsize = gd,mu=mu,adjust=adjust) 
  w=west$w
  gd=west$gd
  #step 3: Estimate clfdr
  indm=max(which(gd<=mu0))
  gd.null=gd[1:indm]
  w.null=w[1:indm]
  fnullest=function(x,s){
    #density at x for the null function
    res=dnorm(gd.null-x,0,s)
    return(res)
  }
  A=mapply(fnullest, x,sig)
  fnull.dens=t(A)%*%w.null
  clfdr.dd=pmin(fnull.dens/nest.res$dens,1)
  clfdr.dd=pmax(clfdr.dd,0)

  fest=function(x,s){
    #density at x for the marginal distribution
    res=dnorm(gd-x,0,s)
    return(res)
  }
  B=mapply(fest, x,sig)
  fest.dens=t(B)%*%w
  
  festprime=function(x,s){
    res=-(x-gd)/(s^2)*dnorm(gd-x,0,s)
    return(res)
  }
  C=mapply(festprime, x,sig)
  fest.prime=t(C)%*%w
  posterior=x+sig^2*fest.prime/fest.dens
   if(mod==2){

    clfdr.dd=pmin(fnull.dens/fest.dens,1)
    clfdr.dd=pmax(clfdr.dd,0)
  }
  # step 4: selection

  d.dd=rsorddgen.func(x-mu0,clfdr.dd,q)

  res=list(clfdr=clfdr.dd,de=d.dd$de,w=w,dens=nest.res$dens,grid=gd,post=posterior)
}




rsstathomo.func=function(x,sig,mu0,epsilon=10^(-5),gd=200,mod=2,jacknife=T,mu=NULL,conservative=F){


  x.ds<-density(x, n=2000)
  x.dens<-lin.itp(x, x.ds$x, x.ds$y)
  west=west.func(x.dens,x,sig,gdsize = gd,mu=mu) 
  w=west$w
  gd=west$gd
  if(conservative){
  indm=max(which(gd<=mu0))
  }else{
  indm=min(which(gd>=mu0))}
  gd.null=gd[1:indm]
  w.null=w[1:indm]
  fnullest=function(x,s){
    res=dnorm(gd.null-x,0,s)
    return(res)
  }
  A=mapply(fnullest, x,sig)
  fnull.dens=t(A)%*%w.null
  if(mod==1){
  clfdr.dd=pmin(fnull.dens/x.dens,1)
  clfdr.dd=pmax(clfdr.dd,0)
  }
  #######################################
  if(mod==2){
    fest=function(x,s){
      #density at x for the null function
      res=dnorm(gd-x,0,s)
      return(res)
    }
    B=mapply(fest, x,sig)
    fest.dens=t(B)%*%w

    clfdr.dd=pmin(fnull.dens/fest.dens,1)
    clfdr.dd=pmax(clfdr.dd,0)
  }
  
  res=list(clfdr=clfdr.dd,w=w,gd=gd)
  
}



rvalue.func=function(x,sig,mu0,epsilon=10^(-3),gd=100,mod=2,jacknife=T,mu=NULL,adjust=T){
# This function computes the r-values by changing FDR level alpha (see definition 1 of the paper)
   nest.res=nest.func(x,sig,jacknife=jacknife)
  dens.dd=nest.res$dens
  west=west.func(dens.dd,x,sig,gdsize = gd,mu=mu,adjust=adjust) 

  w=west$w
  gd=west$gd
  indm=max(which(gd<=mu0))
  gd.null=gd[1:indm]
  w.null=w[1:indm]
  fnullest=function(x,s){
    res=dnorm(gd.null-x,0,s)
    return(res)
  }
  A=mapply(fnullest, x,sig)
  fnull.dens=t(A)%*%w.null
  clfdr.dd=pmin(fnull.dens/nest.res$dens,1)
  clfdr.dd=pmax(clfdr.dd,0)
 
  fest=function(x,s){

    res=dnorm(gd-x,0,s)
    return(res)
  }
  B=mapply(fest, x,sig)
  fest.dens=t(B)%*%w
  if(mod==2){
    
    clfdr.dd=pmin(fnull.dens/fest.dens,1)
    clfdr.dd=pmax(clfdr.dd,0)
  }
  rgrid=sort(clfdr.dd)
  rejectionmat=matrix(rep(0, length(rgrid)*length(x)), length(rgrid), length(x))
  for (j in 1:length(rgrid)) {
    d.dd=rsorddgen.func(x-mu0,clfdr.dd,q=rgrid[j])
    rejectionmat[j,]=d.dd$de
  }
  
  rind=rep(0,length(x))
  for (j in 1:length(x)) {
    rind[j]=min( which(rejectionmat[,j]!=0))
  }
  rval=rgrid[rind]
  return(rval)
}



rvalue2.func=function(x,sig,q=0.1,epsilon=10^(-3),gd=100,mod=2,jacknife=T,mu=NULL,adjust=T){
  # This function computes the r-values by changing reference level mu_0 (see definition 2 of the paper)
  nest.res=nest.func(x,sig,jacknife=jacknife)
  dens.dd=nest.res$dens


  west=west.func(dens.dd,x,sig,gdsize = gd,mu=mu,adjust=adjust) 

 
  w=west$w
  gd=west$gd

  rgrid=seq(gd[3],max(gd),length.out=length(x))

  rejectionmat=matrix(rep(0, length(rgrid)*length(x)), length(rgrid), length(x))
  for (j in 1:length(rgrid)) {
    mu0=rgrid[j]
    indm=max(which(gd<=mu0))
    gd.null=gd[1:indm]
    w.null=w[1:indm]
    fnullest=function(x,s){
     
      res=dnorm(gd.null-x,0,s)
      return(res)
    }
    A=mapply(fnullest, x,sig)
    fnull.dens=t(A)%*%w.null
    clfdr.dd=pmin(fnull.dens/nest.res$dens,1)
    clfdr.dd=pmax(clfdr.dd,0)
    
    fest=function(x,s){
     
      res=dnorm(gd-x,0,s)
      return(res)
    }
    B=mapply(fest, x,sig)
    fest.dens=t(B)%*%w
    if(mod==2){
      
      clfdr.dd=pmin(fnull.dens/fest.dens,1)
      clfdr.dd=pmax(clfdr.dd,0)
    }
    
    
    d.dd=rsorddgen.func(x-mu0,clfdr.dd,q)
    rejectionmat[j,]=d.dd$de
  }
  
  rind=rep(0,length(x))
  for (j in 1:length(x)) {
    rind[j]=max( which(rejectionmat[,j]!=0))
  }
  rval=rgrid[rind]
  return(rval)
}




