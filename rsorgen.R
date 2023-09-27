rsorgen.func=function(postx,clfdr,q=0.1){
  decision=rep(0,length(postx))
  
  
  pos.emuind=which((postx)>=0)
  neg.emuind=which((postx)<0)
  
  pos.clfdrind=which((clfdr-q)>0)
  neg.clfdrind=which((clfdr-q)<=0)
  
  #grouping
  gp0=intersect(pos.emuind,neg.clfdrind)
  gp2=intersect(neg.emuind,neg.clfdrind)
  gp1=intersect(pos.emuind,pos.clfdrind)
  gp3=intersect(neg.emuind,pos.clfdrind)
  
  tstat=postx/(clfdr-q)
  
  tstat2=(postx-min(postx))/clfdr
  
  #tstat.gp1=tstat2[gp1]
  tstat.gp1=tstat[gp1]
  tstat.gp2=tstat[gp2]
  clfdr.gp0=clfdr[gp0]
  clfdr.gp1=clfdr[gp1]
  clfdr.gp2=clfdr[gp2]
  
  postx.gp0=postx[gp0]
  postx.gp1=postx[gp1]
  postx.gp2=postx[gp2]
  postx.gp3=postx[gp3]
  
  #step 1:
  rejset=gp0
  

  
  #step 2:
  
  sorted.clfdrgp1=clfdr.gp1[order(tstat.gp1,decreasing = T)]
  sorted.postx.gp1=postx.gp1[order(tstat.gp1,decreasing = T)]
  sorted.gp1=gp1[order(tstat.gp1,decreasing = T)]
  
  sorted.clfdrgp2=clfdr.gp2[order(tstat.gp2,decreasing = F)]
  sorted.postx.gp2=postx.gp2[order(tstat.gp2,decreasing = F)]
  sorted.gp2=gp2[order(tstat.gp2,decreasing = F)]
  
  if (clfdr[which.min(clfdr)]>q)
  {
  
    rejset<-NULL

  }
  else{
    rejset=gp0
    
    while (mean(clfdr[rejset])<q &&length(sorted.gp1)>0) {
      rejset=c(rejset,sorted.gp1[1])
     
      
      sorted.clfdrgp1=sorted.clfdrgp1[-1]
      sorted.postx.gp1=sorted.postx.gp1[-1]
      sorted.gp1=sorted.gp1[-1]
      
    }
    
  }
  #compute ETP
  ETP.old=0
  if(length(sorted.gp1)==0){
    decision[rejset]=1
    res=list(de=decision)
    return(res)
  }
  ETP.new=sum(postx[rejset])
  
  
  #step 4 and 5
  while (ETP.new>ETP.old && length(sorted.gp2)>0 && length(sorted.gp1)>0) {
    ETP.old=ETP.new
    rejset=c(rejset,sorted.gp2[1])
 
    sorted.gp2=sorted.gp2[-1]
    while (mean(clfdr[rejset])<q &&length(sorted.gp1)>0) {
      rejset=c(rejset,sorted.gp1[1])
      
      sorted.clfdrgp1=sorted.clfdrgp1[-1]
      sorted.postx.gp1=sorted.postx.gp1[-1]
      sorted.gp1=sorted.gp1[-1]
      
    }
    if(length(sorted.gp1)==0){
      decision[rejset]=1
      res=list(de=decision)
      return(res)
    }
    ETP.new=sum(postx[rejset])
    
  }
  
  decision[rejset]=1
  
  #hypotheses from gp1 that are rejected 
  gp1.rej=intersect(rejset,gp1)
  
  
  rk=order(tstat2,decreasing = T)
  
  res=list(de=decision,rk=rk,gp2num=length(gp2),th=min(tstat[gp1.rej]))
  return(res)
  
}





rsorddgen.func=function(postx,clfdr,q=0.1){
  # this is rsorgen with computational shortcut 
  decision=rep(0,length(postx))
  
  pos.emuind=which(postx>=0)
  neg.emuind=which(postx<0)
  
  pos.clfdrind=which((clfdr-q)>0)
  neg.clfdrind=which((clfdr-q)<=0)
  
  #grouping
  gp0=intersect(pos.emuind,neg.clfdrind)
  gp2=intersect(neg.emuind,neg.clfdrind)
  gp1=intersect(pos.emuind,pos.clfdrind)
  gp3=intersect(neg.emuind,pos.clfdrind)
  

  
  Etp.dd=rep(0,(length(gp2)+1))
  
  tstat=postx/(clfdr-q)
  tstat2=(postx-min(postx))/clfdr
  rk=order(tstat2,decreasing = T)
  

  tstat.gp1=tstat[gp1]
  tstat.gp2=tstat[gp2]
  clfdr.gp0=clfdr[gp0]
  clfdr.gp1=clfdr[gp1]
  clfdr.gp2=clfdr[gp2]
  
  postx.gp0=postx[gp0]
  postx.gp1=postx[gp1]
  postx.gp2=postx[gp2]
  postx.gp3=postx[gp3]
  
  #step 1:
  rejset=gp0
  

  
  #step 2:
  #sorted.tstatgp3=sort(tstat.gp3,decreasing = T) #this is only used for ranking gp3
  sorted.clfdrgp1=clfdr.gp1[order(tstat.gp1,decreasing = T)]
  sorted.postx.gp1=postx.gp1[order(tstat.gp1,decreasing = T)]
  sorted.gp1=gp1[order(tstat.gp1,decreasing = T)]
  
  sorted.clfdrgp2=clfdr.gp2[order(tstat.gp2,decreasing = F)]
  sorted.postx.gp2=postx.gp2[order(tstat.gp2,decreasing = F)]
  sorted.gp2=gp2[order(tstat.gp2,decreasing = F)]
  
  if((length(gp0)==0))
  {
    print('gp0==0')
  
    decision=rep(0,length(x))
    res=list(de=decision,rk=rk)
    return(res)
  
  }
  else{
    rejset=gp0
 
    while (((mean(clfdr[rejset])<q) &&(length(sorted.gp1)>0))) {
      rejset=c(rejset,sorted.gp1[1])
     
      
      sorted.clfdrgp1=sorted.clfdrgp1[-1]
      sorted.postx.gp1=sorted.postx.gp1[-1]
      sorted.gp1=sorted.gp1[-1]
    
      
    }
    
  }
  #compute ETP
  if(length(sorted.gp1)==0){
    decision[rejset]=1
    res=list(de=decision,rk=rk)
    return(res)
  }
  
  if(length(gp2)==0){
    decision[rejset]=1
    gp1.rej=intersect(rejset,gp1)
    res=list(de=decision,gp2num=length(gp2),th=min(tstat[gp1.rej]))
    return(res)
  }
  ETP.new=sum(postx[rejset])
  Etp.dd[1]=ETP.new
  
  
  #step 4 and 5
  for (i in 1:length(gp2)) {
    rejset=c(rejset,sorted.gp2[1])
    sorted.gp2=sorted.gp2[-1]
   
   
    
    while (((mean(clfdr[rejset])<q) &&(length(sorted.gp1)>0))) {
      rejset=c(rejset,sorted.gp1[1])
      
      sorted.clfdrgp1=sorted.clfdrgp1[-1]
      sorted.postx.gp1=sorted.postx.gp1[-1]
      sorted.gp1=sorted.gp1[-1]
     
    }
    Etp.dd[i+1]=sum(postx[rejset])
    if(length(sorted.gp1)==0){
      break
    } 
    
  }
  
  k=which.max(Etp.dd)
  sorted.clfdrgp1=clfdr.gp1[order(tstat.gp1,decreasing = T)]
  sorted.postx.gp1=postx.gp1[order(tstat.gp1,decreasing = T)]
  sorted.gp1=gp1[order(tstat.gp1,decreasing = T)]
  
  sorted.clfdrgp2=clfdr.gp2[order(tstat.gp2,decreasing = F)]
  sorted.postx.gp2=postx.gp2[order(tstat.gp2,decreasing = F)]
  sorted.gp2=gp2[order(tstat.gp2,decreasing = F)]
  
  rejset=c(gp0,sorted.gp2[1:(k-1)])
  while (mean(clfdr[rejset])<q &&length(sorted.gp1)>0) {
    rejset=c(rejset,sorted.gp1[1])
   
    sorted.clfdrgp1=sorted.clfdrgp1[-1]
    sorted.postx.gp1=sorted.postx.gp1[-1]
    sorted.gp1=sorted.gp1[-1]
    
  }

  decision[rejset]=1
  
  gp1.rej=intersect(rejset,gp1)
  
  res=list(de=decision,gp2num=length(gp2),th=min(tstat[gp1.rej]))
  
  return(res)
  
}