l<-commandArgs(TRUE)

assocName<-l[1]
pheno<-l[2]
out<-l[3]
env<-l[4]

if(l[5]!="NONE" & !(grepl("\\.gsm$",l[5]))){
    ## badSNPs list of SNPs to be removed basd on QC
    badSNPsFile<-commandArgs(trailingOnly=T)[5]
}


qqpemil<-function(x,ci=TRUE,add=FALSE,ylab="Observed log10(p-value)",xlab="Expected log10(p-value)",maxLogP,...){
    x<-x[!is.na(x)]
    if(!missing(maxLogP)){
        x[x<10^-maxLogP]<-10^-maxLogP
    }
    N<-length(x)
    chi1<-qchisq(1-x,1)
    x<-sort(x)
    lambda<-round(median(chi1)/qchisq(0.5,1),2)
    e<- -log((1:N-0.5)/N,10)
    if(add){
        points(e,-log(x,10),...)
    } else{
        plot(e,-log(x,10),ylab=ylab,xlab=xlab,...)
        abline(0,1,col=2,lwd=2)
    }

    ##title(paste("lambda=",lambda), cex=1.5)
    mtext(paste0("lamdbda = ",lambda))

    if(ci){
        c95<-qbeta(0.95,1:N,N-(1:N)+1)
        c05<-qbeta(0.05,1:N,N-(1:N)+1)
        lines(e,-log(c95,10))
        lines(e,-log(c05,10))
    }
}



manPlot<-function(x,chr,pass,sign,main,collar=c("darkblue","#67a9cf","orange","grey")){
  keep<-!is.na(x)
  x<-x[keep]
  chr<-chr[keep]
  x<-ifelse(x<1e-30,1e-30,x) 
  col<-chr%%2+1
  pch<-rep(16,length(x))
  if(!missing(pass)){ ## NA 
    col[!pass[keep]]<-4
    pch[!pass[keep]]<-1
  }
  par(mar=c(5.1, 5.1, 4.1, 1.1))
  
  maxy = max((-log(x,base=10)))
  if(maxy<8){
    maxy<-8
  }
  plot(-log(x,base=10),col=collar[col],ylab=expression(-log[10](italic(P))),xlab="Chromosomes",main=main,cex=1,lwd=2,pch=pch,axes=F,cex.lab=2,cex.main=2,ylim=c(0,maxy+0.2*maxy))
  box()
  axis(2,las=1,cex.axis=1.8)
  t<--log(0.05/length(x),base=10)
  print(0.05/length(x))
  abline(h=sign,lty=2,col=1)
  mtext(unique(chr),side=1,at=as.vector(tapply(1:length(chr),chr,mean)))
  
  ##          legend(0,t-0.5,"Bonferroni correction",lty=2,bty="n",col=1,cex=2)
  ## legend("topright",c("known","novel"),pch=c(4,1),bty="n",cex=2) 
  
}

qqPlot<-function(x,pass,phe,main="IHIT"){
  
  ## assoc$p_lrt,rep(T,nrow(assoc)),pheno
  x
  
  par(mar=c(5.1, 5.1, 4.1, 1.1)) ## margin of sides in plot
  keep<-!is.na(x)
  x<-x[keep]
  x<-ifelse(x<1e-30,1e-30,x)
  pass=pass[keep]
  maxy = max((-log(x,base=10)))
 
  qqpemil(x,pch=16,col=ifelse(pass,"darkblue","grey"),main=phe,las=1,cex.lab=2,cex.main=2,ylim=c(0,maxy+0.2*maxy),
      xlab=expression(Expected~~-log[10](italic(P))),ylab=expression(Observed~~-log[10](italic(P))),cex.axis=1.5)    
 
}



##################
###### FILTERING OF SNPS RELATING TO OLD METABO SNPCHIP
################

## old metabo filter not really sure if valid any longer??
##csv<-read.csv(paste("/home/ida/web/greenland/web/assoRes/run2015_v1_IHIT_addtrans_",pheno,".assoc.txt.csv.gz",sep=""),as.is=T)

## converted logical
##csv$passFilter<-csv$passFilter!="FALSE"

##ihit>0.025 & nonIhit==0
##ff<-c("chr1:119324759","chr2:165306282","chr2:226752420","chr4:81409339","chr6:72247810","chr6:118873435","chr6:160686672","chr6:160777058","chr10:114739109","chr11:2429008","chr11:47510633","chr12:48503697","chr12:110127239","chr12:110251584","chr12:110481139","chr15:59983796","chr15:89319651","chr16:28798044","chr16:66771700","chr18:56083508")
##csv[ csv$metaboName%in%ff,"passFilter"]<-FALSE

##csv[is.na(csv$passFilter),"passFilter"]<-FALSE

## remove those with too low MAF, C18_3_n6 is binary
##if(pheno %in% c("C18_3_n6")){
##   csv[csv$MAF_ALL<0.05 | is.na(csv$MAF_ALL),"passFilter"]<-FALSE
##} else{
##    csv[csv$MAF_ALL<0.01 | is.na(csv$MAF_ALL),"passFilter"]<-FALSE  
##}

## remove all from chr 0
##csv[ csv$chrHg19==0,"passFilter"]<-FALSE

##csv2<-csv[ csv$passFilter==TRUE,]

## load("/net/pontus/pontus/data/anders/shiny/dataGreenland2015c/csvInfo")
## csvInfo2<-csvInfo[ !is.na(csvInfo$metaboName),]

## ff<-c("chr1:119324759","chr2:165306282","chr2:226752420","chr4:81409339","chr6:72247810","chr6:118873435","chr6:160686672","chr6:160777058","chr10:114739109","chr11:2429008","chr11:47510633","chr12:48503697","chr12:110127239","chr12:110251584","chr12:110481139","chr15:59983796","chr15:89319651","chr16:28798044","chr16:66771700","chr18:56083508")

##################
###### FILTERING OF SNPS RELATING TO OLD METABO SNPCHIP
################

assoc<-read.table(assocName,as.is=T,h=T)

if(l[5]!="NONE" & !(grepl("\\.gsm$",l[5]))){
    ## object loaded assumed to be called csvInfo, with passFilter and snpID columns
    load(badSNPsFile)
    assoc2<-assoc[ !assoc$rs%in%csvInfo[ csvInfo$passFilter==F,"snpID"],]        
} else{
    assoc2<-assoc
}

bonf<-0.05/nrow(assoc2)

bitmap(paste(out,"_manhattan.png",sep=""),res=300)
## for gxe we can only do p_score
manPlot(assoc2[,"p_score"],assoc2$chr,rep(T,nrow(assoc2)),sign=-log10(bonf),main=pheno)
dev.off()

bitmap(paste(out,"_qqp.png",sep=""),res=300)
## for gxe we can only do p_score
qqPlot(assoc2[,"p_score"],rep(T,nrow(assoc2)),pheno)
dev.off()

dat<-read.table("data/refGeneHG19.gz",as.is=T,head=T,comment.char="")

if(nrow(assoc2[  assoc2[,"p_score"]<5e-08, ])>0){
    ## for gxe we can only do p_score
    signif<-assoc2[ assoc2[,"p_score"]<5e-08,]
    signif[,"genes"]<-apply(signif,1,function(x) paste(dat[ dat$chrom==paste0("chr",trimws(x["chr"])) & dat$txStart <= as.numeric(x["ps"]) & dat$txEnd >= as.numeric(x["ps"]),"name2"],collapse=","))    
    write.table(signif,paste0(out,"_below5e-08.assoc.txt"),col=T,row=F,quote=F)
}



