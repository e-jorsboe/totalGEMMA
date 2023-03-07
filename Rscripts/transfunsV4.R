qtrans <- function(x){
    k<-!is.na(x)
    ran<-rank(x[k])
    y<-qnorm((1:sum(k)-0.5)/sum(k))
    x[k]<-y[ran]
    x
}

qnormtransform<-function(x){
    t1 <- x
    t2 <- x[!is.na(x)]
    r <- qnorm(rank(t2)/(length(rank(t2))+1))
    t1[!is.na(t1)] <- r
    t1
}

randomqnormtransform<-function(x,seedval=1){
    set.seed(seedval)
    t1 <- x
    t2 <- x[!is.na(x)]
    r <- qnorm(rank(t2,ties.method="random")/(length(rank(t2,ties.method="random"))+1))
    t1[!is.na(t1)] <- r
    t1
}

cctransform<-function(x){
    
    ctrls = x==0
    newx = x
    newx[ctrls]=1
    newx[!ctrls]=2
    newx
}

cctransformV2<-function(x){
    if(all(x%in%1:2)){
        
        x<-x-1
    }
    ctrls = x==0
    newx = x
    newx[ctrls]=1
    newx[!ctrls]=2
    newx
}


# tentative version of case/control function
cctransform01<-function(x){
  if(any(x==0)){
  ctrls = x==0
} else{
  ctrls = x==1
}
  newx = x
  newx[ctrls]=0
  newx[!ctrls]=1
  newx

}

