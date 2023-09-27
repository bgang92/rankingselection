

bh.func<-function(pv, q)
{ 
  # the input 
    # pv: the p-values
    # q: the FDR level
  # the output 
    # nr: the number of hypothesis to be rejected
    # th: the p-value threshold
    # re: the index of rejected hypotheses
    # ac: the index of accepted hypotheses
    # de: the decision rule

  m=length(pv)
  st.pv<-sort(pv)   
  #print(length(st.pv))
  #print(m)
  pvi<-st.pv/1:m
  hps<-rep(0, m)
  if (max(pvi<=(q/m))==0)
  {
    k<-0
    pk<-1
    reject<-NULL
    accept<-1:m
  }
  else
  {
    k<-max(which(pvi<=(q/m)))
    pk<-st.pv[k]
    reject<-which(pv<=pk)
    accept<-which(pv>pk)
    hps[reject]<-1
  }
  y<-list(nr=k, th=pk, re=reject, ac=accept, de=hps)
  return (y)
}

sc.func<-function(lfdr, q)
{

## USAGE
 # mt.sc(lfdr, q)
 
## ARGUMENTS
 # lfdr: local false discovery rate sequence
 # q: the FDR level

## VALUES
 # nr: the number of rejected hypotheses
 # th: the threshold
 # re: the rejected hypotheses
 # ac: the accepted hypotheses
 # de: the decision rule

  m=length(lfdr)
  st.lfdr<-sort(lfdr)
  hps<-rep(0, m)
  #if (min(lfdr)>q)
    if (lfdr[which.min(lfdr)]>q)
  {
    k<-0
    threshold<-1
    reject<-NULL
    accept<-1:m
  }
  else
  {
    k=1
    while(k<m && (1/k)*sum(st.lfdr[1:k])<q){
      k=k+1
    }
    k<-k-1
    threshold<-st.lfdr[k]
    reject<-which(lfdr<=threshold)
    accept<-which(lfdr>threshold)
    hps[reject]<-1
  }
  y<-list(nr=k, th=threshold, re=reject, ac=accept, de=hps)
  return (y)
}








