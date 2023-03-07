args<-commandArgs(trailingOnly = T)

if(length(args)==0){

    print("Arguments have to be supplied: (this version is supposed to adjust for cohort, but does not really work!")
    print("1. phenotype to analyse")
    print("2. SNP to analyse")
    print("3. if additive 0, recessive 1 or full model 2 - optional (additive is default)")
    print("4. path to genotype data to run on, - optional default is metaboDec2018 data")
    print("5. enviromental variable - optional")
    stop()
}



library(survival)
library(ISwR)
library(RColorBrewer)
library(GMMAT)

## to get proper colours, does not import from .Rprofile for some reason
colorblind6b<-function(){
  grDevices::palette(c("#b2182b","#ef8a62","#fddbc7","#d1e5f0","#67a9cf","#2166ac"))
  grDevices::dev.off()
}
colorblind6b()




pheno<-args[1]
snp<-args[2]

if(length(args)>2){    
    genoModel<-as.numeric(args[3])
}

genoData<-"/home/emil/OMNI5/CARDIO_METABOCHIP/datplus_QCed_newIDs_allcohorts_allSNPsDec2018_manuallyPlus"
if(length(args)>3){
    genoData<-args[4]
}


pheFile<-read.csv("/home/greenland/data/phenotypes/Greenland_Pheno_23maj2016.csv",as.is=T)

## B99 + IHIT raw phenotypes merged with SI and LDLR genotyping
j<-load("/home/emil/OMNI5/plotsCOE2018/B99_IHITdata_updated25thApr.Rdata")

## CVD registry data merged with data with PCs,sex,age and LDLR variant
k<-load("/home/emil/OMNI5/phenotypes/registryCVDdata270619includeDeathUntil2017withBirthyearAllAndBBHfrederikAngina.Rdata")

ifEnv<-0
if(length(args)>4){    
    env<-args[5]
    if( !(any(env%in%colnames(total)) | any(env%in%colnames(merged)) | any(env%in%colnames(pheFile))) ){
        print("env not in data!")
        stop()
    }
    ifEnv<-1
}

if(any(pheno%in%colnames(total))){
    int<-intersect(pheFile$particid,total$particid)    
    pheFile2<-pheFile[ pheFile$particid%in%int,]
    total<-total[ total$particid%in%int,]
    pheFile2<-pheFile2[ order(match(pheFile2$particid,total$particid)),]    
    table(pheFile2$particid==total$particid)
    if(!all(pheFile2$particid==total$particid)){
        print("problem with data")
        stop()
    }
    pheFile2[,pheno]<-total[,pheno]
    if(ifEnv){
        pheFile2[,env]<-total[,env]        
    }
    
} else if(any(pheno%in%colnames(merged))){
    int<-intersect(pheFile$particid,merged$particid)    
    pheFile2<-pheFile[ pheFile$particid%in%int,]
    total<-total[ total$particid%in%int,]
    pheFile2<-pheFile2[ order(match(pheFile2$particid,merged$particid)),]    
    table(pheFile2$particid==merged$particid)
    if(!all(pheFile2$particid==merged$particid)){
        print("problem with data")
        stop()
    }
    pheFile2[,pheno]<-merged[,pheno]
    if(ifEnv){
        pheFile2[,env]<-merged[,env]        
    }

} else if(any(pheno%in%colnames(pheFile))){
    pheFile2<-pheFile

} else{
    print("Phenotype not found!")
    stop()

}


table(pheFile2[,pheno],useNA='ifany')

if((all(na.omit(pheFile2[,pheno])%in%c(1,2,9)) | all(na.omit(pheFile2[,pheno])%in%c(0,1,9)))){
    print("treats 9 as missing info gives value of NA")
    pheFile2[ pheFile2[,pheno]%in%9,pheno]<-NA
}




if(!(all(na.omit(pheFile2[,pheno])%in%c(1,2)) | all(na.omit(pheFile2[,pheno])%in%c(0,1)))){
    print("pheno not binary!")
    stop()
}

if(ifEnv){    
    keep<-pheFile2[ !is.na(pheFile2[,pheno]) & !is.na(pheFile2[,"sex"])  & !is.na(pheFile2[,"age"]) & !is.na(pheFile2[,env]) ,"particid"]
} else{
    
    keep<-pheFile2[ !is.na(pheFile2[,pheno]) & !is.na(pheFile2[,"sex"])  & !is.na(pheFile2[,"age"]) ,"particid"]
}

fam<-read.table("/home/emil/OMNI5/CARDIO_METABOCHIP/datplus_QCed_newIDs_allcohorts_allSNPsDec2018_manuallyPlus.fam",as.is=T)

keep2<-keep[ keep%in%fam$V1]

pheFile3<-pheFile2[ pheFile2$particid%in%keep2 & !duplicated(pheFile2$particid),]



write.table(keep2,col=F,row=F,qu=F,"/home/emil/tmp/tmprunGLMM.keep")

system(paste0("grep -wf /home/emil/tmp/tmprunGLMM.keep ",genoData,".fam > /home/emil/tmp/runGLMM.keep"))
system(paste0("plink --bfile ",genoData," --keep /home/emil/tmp/runGLMM.keep --make-bed --out /home/emil/tmp/runGLMM"))
system("/emc/emil/riccoStuff/gemmaTotal/prog/gemma.linux -bfile /home/emil/tmp/runGLMM -gk 2 -outdir /home/emil/tmp/ -o runGLMM -maf 0.05 -miss 0.01")


pl<-plinkV2("/home/emil/tmp/runGLMM")

if(!(snp%in%colnames(pl$geno))){
    print("SNP is not in data")
    stop()
}


pheFile3<-pheFile3[ order(match(pheFile3$particid,pl$fam$V1)),]

table(pheFile3$particid==pl$fam$V1)

if(genoModel==2){
    if(ifEnv){
        snp2<-c("SNPHO",pl$bim[ pl$bim$V2%in%snp,6],pl$bim[ pl$bim$V2%in%snp,5],as.numeric(pl$geno[,snp]==0)*pheFile3[,env])
        snp3<-c("SNPHE",pl$bim[ pl$bim$V2%in%snp,6],pl$bim[ pl$bim$V2%in%snp,5],as.numeric(pl$geno[,snp]==1)*pheFile3[,env])
    } else{
        ## HO
        snp2<-c("SNPHO",pl$bim[ pl$bim$V2%in%snp,6],pl$bim[ pl$bim$V2%in%snp,5],as.numeric(pl$geno[,snp]==0))
        ## HE
        snp3<-c("SNPHE",pl$bim[ pl$bim$V2%in%snp,6],pl$bim[ pl$bim$V2%in%snp,5],as.numeric(pl$geno[,snp]==1))
    }
} else if(genoModel==1){
    if(ifEnv){
        snp2<-c("SNP1",pl$bim[ pl$bim$V2%in%snp,6],pl$bim[ pl$bim$V2%in%snp,5],as.numeric(pl$geno[,snp]==0)*pheFile3[,env])
    } else{
        snp2<-c("SNP1",pl$bim[ pl$bim$V2%in%snp,6],pl$bim[ pl$bim$V2%in%snp,5],as.numeric(pl$geno[,snp]==0))
    }
} else{
    ## first major allele and then is minor allele
    if(ifEnv){
        snp2<-c("SNP1",pl$bim[ pl$bim$V2%in%snp,6],pl$bim[ pl$bim$V2%in%snp,5],as.numeric(2-pl$geno[,snp])*pheFile3[,env])
    } else{
        snp2<-c("SNP1",pl$bim[ pl$bim$V2%in%snp,6],pl$bim[ pl$bim$V2%in%snp,5],as.numeric(2-pl$geno[,snp]))
    }
}

if(genoModel==2){

    write(paste(snp2,collapse=" "),"/home/emil/tmp/runGLMM.snpHO",sep="\t")
    system("cat /home/emil/tmp/runGLMM.snpHO | tr ' ' '\t' > /home/emil/tmp/runGLMM.snpHO.tab")

    write(paste(snp3,collapse=" "),"/home/emil/tmp/runGLMM.snpHE",sep="\t")
    system("cat /home/emil/tmp/runGLMM.snpHE | tr ' ' '\t' > /home/emil/tmp/runGLMM.snpHE.tab")


} else{

    write(paste(snp2,collapse=" "),"/home/emil/tmp/runGLMM.snp2",sep="\t")
    system("cat /home/emil/tmp/runGLMM.snp2 | tr ' ' '\t' > /home/emil/tmp/runGLMM.snp2.tab")

}

if(all(pheFile3[,pheno]%in%c(1,2))){
    pheno01<-pheFile3[,pheno]-1
} else{
    pheno01<-pheFile3[,pheno]
}

if(sum(pheno01==1) > sum(pheno01==0)){
    print("swapping around phenotypes as 0 are cases and 1 are controls")
    pheno01[pheno01==1]<-2
    pheno01[pheno01==0]<-1
    pheno01[pheno01==2]<-0
}


b99<-as.numeric(pheFile3[,"cohort"]=="B99")
b99_ihit<-as.numeric(pheFile3[,"cohort"]=="B99_IHIT")
bbh<-as.numeric(pheFile3[,"cohort"]=="BBH")

## to handle when only one cohort is present
cohorts<-data.frame(b99=b99,b99_ihit=b99_ihit,bbh=bbh,dummy=1:nrow(pheFile3),stringsAsFactors=F)
cohorts<-cohorts[,colSums(cohorts)>0]
## to handle that some cohorts have to few individuals for the GLMM to converge
##cohorts<-cohorts[,colSums(cohorts)>75]
## emil - changed from 75 to 100
cohorts<-cohorts[,colSums(cohorts)>100]

if(genoModel==2){

    if(ifEnv){
        data<-data.frame(disease=pheno01,sex=pheFile3[,"sex"],age=pheFile3[,"age"],id=pheFile3[,"particid"],genoHEcov=as.numeric(pl$geno[,snp]==1),genoHOcov=as.numeric(pl$geno[,snp]==0),genoHEcovEnv=as.numeric(pl$geno[,snp]==1)*pheFile3[,env],genoHOcovEnv=as.numeric(pl$geno[,snp]==0)*pheFile3[,env],env=pheFile3[,env],genoPrint=as.numeric(2-pl$geno[,snp]),cohorts,stringsAsFactors=F)
        data<-data[,1:(ncol(data)-1)]

    } else{
        data<-data.frame(disease=pheno01,sex=pheFile3[,"sex"],age=pheFile3[,"age"],id=pheFile3[,"particid"],genoHE=as.numeric(pl$geno[,snp]==1),genoHO=as.numeric(pl$geno[,snp]==0),genoPrint=as.numeric(2-pl$geno[,snp]),cohorts,stringsAsFactors=F)
        data<-data[,1:(ncol(data)-1)]
    }

} else if(genoModel==1){

    if(ifEnv){
        data<-data.frame(disease=pheno01,sex=pheFile3[,"sex"],age=pheFile3[,"age"],id=pheFile3[,"particid"],genoCov=as.numeric(pl$geno[,snp]==0),genoPrint=as.numeric(2-pl$geno[,snp]),env=pheFile3[,env],cohorts,stringsAsFactors=F)
        data<-data[,1:(ncol(data)-1)]

    } else{
        data<-data.frame(disease=pheno01,sex=pheFile3[,"sex"],age=pheFile3[,"age"],id=pheFile3[,"particid"],geno=as.numeric(pl$geno[,snp]==0),genoPrint=as.numeric(2-pl$geno[,snp]),cohorts,stringsAsFactors=F)
        data<-data[,1:(ncol(data)-1)]
    }

} else{

    if(ifEnv){
        data<-data.frame(disease=pheno01,sex=pheFile3[,"sex"],age=pheFile3[,"age"],id=pheFile3[,"particid"],genoCov=as.numeric(2-pl$geno[,snp]),env=pheFile3[,env],genoPrint=as.numeric(2-pl$geno[,snp]),cohorts,stringsAsFactors=F)
        data<-data[,1:(ncol(data)-1)]

    } else{
        data<-data.frame(disease=pheno01,sex=pheFile3[,"sex"],age=pheFile3[,"age"],id=pheFile3[,"particid"],geno=as.numeric(2-pl$geno[,snp]),genoPrint=as.numeric(2-pl$geno[,snp]),cohorts,stringsAsFactors=F)
        data<-data[,1:(ncol(data)-1)]

    }

}

if(ifEnv){

    if(env=="sex"){                                                                                                                                                      
        data<-data[ ,!grepl("sex",colnames(data))]
    }
    if(env=="age"){
        data<-data[ ,!grepl("age",colnames(data))]        
    }

}
    
print("this many cases and controls:")
print(table(data$disease))

if(ifEnv){

    if(genoModel==2){

        print("this many cases and controls - being analysed:")
        print(table(data[ !is.na(data[,"env"]) & !is.na(data$disease) & !is.na(data$sex) & !is.na(data$age) & !is.na(data$genoHEcov) & !is.na(data$id),"disease"]))
        
        print("controls have these genotypes:")
        print(table(data[ !is.na(data[,"env"]) & !is.na(data$disease) & !is.na(data$sex) & !is.na(data$age) & !is.na(data$genoHEcov) & !is.na(data$id) & data$disease%in%0,"genoPrint"]))
        
        print("cases have these genotypes:")
        print(table(data[ !is.na(data[,"env"]) & !is.na(data$disease) & !is.na(data$sex) & !is.na(data$age) & !is.na(data$genoHEcov) & !is.na(data$id) & data$disease%in%1,"genoPrint"]))

    } else{    
        print("this many cases and controls - being analysed:")
        print(table(data[ !is.na(data[,"env"]) & !is.na(data$disease) & !is.na(data$sex) & !is.na(data$age) & !is.na(data$genoCov) & !is.na(data$id),"disease"]))
        
        print("controls have these genotypes:")
        print(table(data[ !is.na(data[,"env"]) & !is.na(data$disease) & !is.na(data$sex) & !is.na(data$age) & !is.na(data$genoCov) & !is.na(data$id) & data$disease%in%0,"genoPrint"]))
        
        print("cases have these genotypes:")
        print(table(data[ !is.na(data[,"env"]) & !is.na(data$disease) & !is.na(data$sex) & !is.na(data$age) & !is.na(data$genoCov) & !is.na(data$id) & data$disease%in%1,"genoPrint"]))

    }
        
} else{
        
    if(genoModel==2){

        print("this many cases and controls - being analysed:")
        print(table(data[ !is.na(data$disease) & !is.na(data$sex) & !is.na(data$age) & !is.na(data$genoHE) & !is.na(data$id),"disease"]))
        
        print("controls have these genotypes:")
        print(table(data[ !is.na(data$disease) & !is.na(data$sex) & !is.na(data$age) & !is.na(data$genoHE) & !is.na(data$id) & data$disease%in%0,"genoPrint"]))
        
        print("cases have these genotypes:")
        print(table(data[ !is.na(data$disease) & !is.na(data$sex) & !is.na(data$age) & !is.na(data$genoHE) & !is.na(data$id) & data$disease%in%1,"genoPrint"]))

    } else{    
        print("this many cases and controls - being analysed:")
        print(table(data[ !is.na(data$disease) & !is.na(data$sex) & !is.na(data$age) & !is.na(data$geno) & !is.na(data$id),"disease"]))
        
        print("controls have these genotypes:")
        print(table(data[ !is.na(data$disease) & !is.na(data$sex) & !is.na(data$age) & !is.na(data$geno) & !is.na(data$id) & data$disease%in%0,"genoPrint"]))
        
        print("cases have these genotypes:")
        print(table(data[ !is.na(data$disease) & !is.na(data$sex) & !is.na(data$age) & !is.na(data$geno) & !is.na(data$id) & data$disease%in%1,"genoPrint"]))

    }
    
}

#######################################################################################


if(ifEnv){
    
    write(paste("this many cases and controls:",collapse=" "),paste0("results/",pheno,".",snp,".",env,".adjCohort.V7.info"))
    write(paste(table(data$disease),collapse=" "),paste0("results/",pheno,".",snp,".",env,".adjCohort.V7.info"),append=T)
    
    if(genoModel==2){
      
      write(paste("this many cases and controls - being analysed:",collapse=" "),paste0("results/",pheno,".",snp,".",env,".adjCohort.V7.info"),append=T)
      write(paste(table(data[ !is.na(data[,"env"]) & !is.na(data$disease) & !is.na(data$sex) & !is.na(data$age) & !is.na(data$genoHEcov) & !is.na(data$id),"disease"]),collapse=" "),paste0("results/",pheno,".",snp,".",env,".adjCohort.V7.info"),append=T)
      
      write(paste("controls have these genotypes:",collapse=" "),paste0("results/",pheno,".",snp,".",env,".adjCohort.V7.info"),append=T)
      write(paste(table(data[ !is.na(data[,"env"]) & !is.na(data$disease) & !is.na(data$sex) & !is.na(data$age) & !is.na(data$genoHEcov) & !is.na(data$id) & data$disease%in%0,"genoPrint"]),collapse=" "),paste0("results/",pheno,".",snp,".",env,".adjCohort.V7.info"),append=T)
      
      write(paste("cases have these genotypes:",collapse=" "),paste0("results/",pheno,".",snp,".",env,".adjCohort.V7.info"),append=T)
      write(paste(table(data[ !is.na(data[,"env"]) & !is.na(data$disease) & !is.na(data$sex) & !is.na(data$age) & !is.na(data$genoHEcov) & !is.na(data$id) & data$disease%in%1,"genoPrint"]),collapse=" "),paste0("results/",pheno,".",snp,".",env,".adjCohort.V7.info"),append=T)
      
      
  } else{
      
      write(paste("this many cases and controls - being analysed:",collapse=" "),paste0("results/",pheno,".",snp,".",env,".adjCohort.V7.info"),append=T)
      write(paste(table(data[ !is.na(data[,"env"]) & !is.na(data$disease) & !is.na(data$sex) & !is.na(data$age) & !is.na(data$genoCov) & !is.na(data$id),"disease"]),collapse=" "),paste0("results/",pheno,".",snp,".",env,".adjCohort.V7.info"),append=T)
      
      write(paste("controls have these genotypes:",collapse=" "),paste0("results/",pheno,".",snp,".",env,".adjCohort.V7.info"),append=T)
      write(paste(table(data[ !is.na(data[,"env"]) & !is.na(data$disease) & !is.na(data$sex) & !is.na(data$age) & !is.na(data$genoCov) & !is.na(data$id) & data$disease%in%0,"genoPrint"]),collapse=" "),paste0("results/",pheno,".",snp,".",env,".adjCohort.V7.info"),append=T)
      
      write(paste("cases have these genotypes:",collapse=" "),paste0("results/",pheno,".",snp,".",env,".adjCohort.V7.info"),append=T)
      write(paste(table(data[ !is.na(data[,"env"]) & !is.na(data$disease) & !is.na(data$sex) & !is.na(data$age) & !is.na(data$genoCov) & !is.na(data$id) & data$disease%in%1,"genoPrint"]),collapse=" "),paste0("results/",pheno,".",snp,".",env,".adjCohort.V7.info"),append=T)
      
  }
    
} else{

    write(paste("this many cases and controls:",collapse=" "),paste0("results/",pheno,".",snp,".adjCohort.V7.info"))
    write(paste(table(data$disease),collapse=" "),paste0("results/",pheno,".",snp,".adjCohort.V7.info"),append=T)
    
    if(genoModel==2){
      
      write(paste("this many cases and controls:",collapse=" "),paste0("results/",pheno,".",snp,".adjCohort.V7.info"))
      write(paste(table(data$disease),collapse=" "),paste0("results/",pheno,".",snp,".adjCohort.V7.info"),append=T)
      
      write(paste("this many cases and controls - being analysed:",collapse=" "),paste0("results/",pheno,".",snp,".adjCohort.V7.info"),append=T)
      write(paste(table(data[ !is.na(data$disease) & !is.na(data$sex) & !is.na(data$age) & !is.na(data$genoHE) & !is.na(data$id),"disease"]),collapse=" "),paste0("results/",pheno,".",snp,".adjCohort.V7.info"),append=T)
      
      write(paste("controls have these genotypes:",collapse=" "),paste0("results/",pheno,".",snp,".adjCohort.V7.info"),append=T)
      write(paste(table(data[ !is.na(data$disease) & !is.na(data$sex) & !is.na(data$age) & !is.na(data$genoHE) & !is.na(data$id) & data$disease%in%0,"genoPrint"]),collapse=" "),paste0("results/",pheno,".",snp,".adjCohort.V7.info"),append=T)
      
      write(paste("cases have these genotypes:",collapse=" "),paste0("results/",pheno,".",snp,".adjCohort.V7.info"),append=T)
      write(paste(table(data[ !is.na(data$disease) & !is.na(data$sex) & !is.na(data$age) & !is.na(data$genoHE) & !is.na(data$id) & data$disease%in%1,"genoPrint"]),collapse=" "),paste0("results/",pheno,".",snp,".adjCohort.V7.info"),append=T)
      
    
  } else{
      
      write(paste("this many cases and controls:",collapse=" "),paste0("results/",pheno,".",snp,".adjCohort.V7.info"))
      write(paste(table(data$disease),collapse=" "),paste0("results/",pheno,".",snp,".adjCohort.V7.info"),append=T)
      
      write(paste("this many cases and controls - being analysed:",collapse=" "),paste0("results/",pheno,".",snp,".adjCohort.V7.info"),append=T)
      write(paste(table(data[ !is.na(data$disease) & !is.na(data$sex) & !is.na(data$age) & !is.na(data$geno) & !is.na(data$id),"disease"]),collapse=" "),paste0("results/",pheno,".",snp,".adjCohort.V7.info"),append=T)
    
      write(paste("controls have these genotypes:",collapse=" "),paste0("results/",pheno,".",snp,".adjCohort.V7.info"),append=T)
      write(paste(table(data[ !is.na(data$disease) & !is.na(data$sex) & !is.na(data$age) & !is.na(data$geno) & !is.na(data$id) & data$disease%in%0,"genoPrint"]),collapse=" "),paste0("results/",pheno,".",snp,".adjCohort.V7.info"),append=T)

      write(paste("cases have these genotypes:",collapse=" "),paste0("results/",pheno,".",snp,".adjCohort.V7.info"),append=T)
      write(paste(table(data[ !is.na(data$disease) & !is.na(data$sex) & !is.na(data$age) & !is.na(data$geno) & !is.na(data$id) & data$disease%in%1,"genoPrint"]),collapse=" "),paste0("results/",pheno,".",snp,".adjCohort.V7.info"),append=T)

    
  }
    

}

print("data looks like this")
print(head(data))


write.table(data,"/home/emil/tmp/runGLMM.data",col=T,row=F,qu=F)

GRM<-read.table("/home/emil/tmp/runGLMM.sXX.txt",as.is=T,h=F)

GRM<-as.matrix(GRM)
rownames(GRM)<-data$id
colnames(GRM)<-data$id

## The Cholesky algorithm will fail if the matrix is not positive definite, 
##c<-chol(GRM)

wald1<-"NULL"
wald2<-"NULL"

if(ifEnv){

    if(all(c("b99","b99_ihit","bbh")%in%colnames(data))){
        print("1")

        if(genoModel==2){

            wald1<-GMMAT::glmm.wald(disease ~ age + sex + genoHEcov + genoHOcov + genoHOcovEnv + env + b99 + b99_ihit + bbh, data=data, id = "id", snps = c("SNPHE"), infile = "/home/emil/tmp/runGLMM.snpHE.tab",infile.ncol.skip = 3,kins = as.matrix(GRM))
            wald2<-GMMAT::glmm.wald(disease ~ age + sex + genoHEcov + genoHOcov + genoHEcovEnv + env + b99 + b99_ihit + bbh, data=data, id = "id", snps = c("SNPHO"), infile = "/home/emil/tmp/runGLMM.snpHO.tab",infile.ncol.skip = 3,kins = as.matrix(GRM))
            
        } else{
        
            wald1<-GMMAT::glmm.wald(disease ~ age + sex + genoCov + env + b99 + b99_ihit + bbh, data=data, id = "id", snps = c("SNP1"), infile = "/home/emil/tmp/runGLMM.snp2.tab",infile.ncol.skip = 3,kins = as.matrix(GRM))

        }
        
    } else if(all(c("b99","b99_ihit")%in%colnames(data))){
        print("2")

        if(genoModel==2){
            
            wald1<-GMMAT::glmm.wald(disease ~ age + sex + genoHEcov + genoHOcov + genoHOcovEnv + env + b99 + b99_ihit, data=data, id = "id", snps = c("SNPHE"), infile = "/home/emil/tmp/runGLMM.snpHE.tab",infile.ncol.skip = 3,kins = as.matrix(GRM))
            wald2<-GMMAT::glmm.wald(disease ~ age + sex + genoHEcov + genoHOcov + genoHEcovEnv + env + b99 + b99_ihit, data=data, id = "id", snps = c("SNPHO"), infile = "/home/emil/tmp/runGLMM.snpHO.tab",infile.ncol.skip = 3,kins = as.matrix(GRM))
            
        } else{
            
            wald1<-GMMAT::glmm.wald(disease ~ age + sex + genoCov + env + b99 + b99_ihit, data=data, id = "id", snps = c("SNP1"), infile = "/home/emil/tmp/runGLMM.snp2.tab",infile.ncol.skip = 3,kins = as.matrix(GRM))
            
        }
        
    } else if(all(c("b99","bbh")%in%colnames(data))){
        print("3")

        if(genoModel==2){
            
            wald1<-GMMAT::glmm.wald(disease ~ age + sex + genoHEcov + genoHOcov + genoHOcovEnv + env + b99 + bbh, data=data, id = "id", snps = c("SNPHE"), infile = "/home/emil/tmp/runGLMM.snpHE.tab",infile.ncol.skip = 3,kins = as.matrix(GRM))
            wald2<-GMMAT::glmm.wald(disease ~ age + sex + genoHEcov + genoHOcov + genoHEcovEnv + env + b99 + bbh, data=data, id = "id", snps = c("SNPHO"), infile = "/home/emil/tmp/runGLMM.snpHO.tab",infile.ncol.skip = 3,kins = as.matrix(GRM))
            
        } else{
                       
            wald1<-GMMAT::glmm.wald(disease ~ age + sex + genoCov + env + b99 + bbh, data=data, id = "id", snps = c("SNP1"), infile = "/home/emil/tmp/runGLMM.snp2.tab",infile.ncol.skip = 3,kins = as.matrix(GRM))
        }
        
    } else if(all(c("b99")%in%colnames(data))){
        print("4")

        if(genoModel==2){
            
            wald1<-GMMAT::glmm.wald(disease ~ age + sex + genoHEcov + genoHOcov + genoHOcovEnv + env + b99, data=data, id = "id", snps = c("SNPHE"), infile = "/home/emil/tmp/runGLMM.snpHE.tab",infile.ncol.skip = 3,kins = as.matrix(GRM))
            wald2<-GMMAT::glmm.wald(disease ~ age + sex + genoHEcov + genoHOcov + genoHEcovEnv + env + b99, data=data, id = "id", snps = c("SNPHO"), infile = "/home/emil/tmp/runGLMM.snpHO.tab",infile.ncol.skip = 3,kins = as.matrix(GRM))
            
        } else{
                        
            wald1<-GMMAT::glmm.wald(disease ~ age + sex + genoCov + env + b99, data=data, id = "id", snps = c("SNP1"), infile = "/home/emil/tmp/runGLMM.snp2.tab",infile.ncol.skip = 3,kins = as.matrix(GRM))

        }
        
    } else if(all(c("b99_ihit")%in%colnames(data))){
        print("5")

                if(genoModel==2){
            
            wald1<-GMMAT::glmm.wald(disease ~ age + sex + genoHEcov + genoHOcov + genoHOcovEnv + env + b99_ihit, data=data, id = "id", snps = c("SNPHE"), infile = "/home/emil/tmp/runGLMM.snpHE.tab",infile.ncol.skip = 3,kins = as.matrix(GRM))
            wald2<-GMMAT::glmm.wald(disease ~ age + sex + genoHEcov + genoHOcov + genoHEcovEnv + env + b99_ihit, data=data, id = "id", snps = c("SNPHO"), infile = "/home/emil/tmp/runGLMM.snpHO.tab",infile.ncol.skip = 3,kins = as.matrix(GRM))
            
        } else{

            wald1<-GMMAT::glmm.wald(disease ~ age + sex + genoCov + env + b99_ihit, data=data, id = "id", snps = c("SNP1"), infile = "/home/emil/tmp/runGLMM.snp2.tab",infile.ncol.skip = 3,kins = as.matrix(GRM))
            
        }
        
    } else if(all(c("bbh")%in%colnames(data))){
        print("6")

        if(genoModel==2){
            
            wald1<-GMMAT::glmm.wald(disease ~ age + sex + genoHEcov + genoHOcov + genoHOcovEnv + env + bbh, data=data, id = "id", snps = c("SNPHE"), infile = "/home/emil/tmp/runGLMM.snpHE.tab",infile.ncol.skip = 3,kins = as.matrix(GRM))
            wald2<-GMMAT::glmm.wald(disease ~ age + sex + genoHEcov + genoHOcov + genoHEcovEnv + env + bbh, data=data, id = "id", snps = c("SNPHO"), infile = "/home/emil/tmp/runGLMM.snpHO.tab",infile.ncol.skip = 3,kins = as.matrix(GRM))
            
        } else{
        
            wald1<-GMMAT::glmm.wald(disease ~ age + sex + genoCov + env + bbh, data=data, id = "id", snps = c("SNP1"), infile = "/home/emil/tmp/runGLMM.snp2.tab",infile.ncol.skip = 3,kins = as.matrix(GRM))
        }
        
    } else{
        print("7")

        if(genoModel==2){
            wald1<-GMMAT::glmm.wald(disease ~ age + sex + genoHEcov + genoHOcov + genoHOcovEnv + env, data=data, id = "id", snps = c("SNPHE"), infile = "/home/emil/tmp/runGLMM.snpHE.tab",infile.ncol.skip = 3,kins = as.matrix(GRM))
            wald2<-GMMAT::glmm.wald(disease ~ age + sex + genoHEcov + genoHOcov + genoHEcovEnv + env, data=data, id = "id", snps = c("SNPHO"), infile = "/home/emil/tmp/runGLMM.snpHO.tab",infile.ncol.skip = 3,kins = as.matrix(GRM))
            
        } else{
            
            wald1<-GMMAT::glmm.wald(disease ~ age + sex + genoCov + env, data=data, id = "id", snps = c("SNP1"), infile = "/home/emil/tmp/runGLMM.snp2.tab",infile.ncol.skip = 3,kins = as.matrix(GRM))
        }
      

    }

    if(genoModel==2){
        write.table(wald1,paste0("results/",pheno,".",snp,".",env,".adjCohort.fullHE.wald.V7.res"),col=T,row=F,qu=F)
        write.table(wald2,paste0("results/",pheno,".",snp,".",env,".adjCohort.fullHO.wald.V7.res"),col=T,row=F,qu=F)                
    } else if(genoModel==1){        
        write.table(wald1,paste0("results/",pheno,".",snp,".",env,".adjCohort.recessive.wald.V7.res"),col=T,row=F,qu=F)        
    } else{
        write.table(wald1,paste0("results/",pheno,".",snp,".",env,".adjCohort.wald.V7.res"),col=T,row=F,qu=F)        
    }
  
} else{
    
    if(all(c("b99","b99_ihit","bbh")%in%colnames(data))){
        print("8")

        if(genoModel==2){
            
            wald1<-GMMAT::glmm.wald(disease ~ age + sex + genoHO + b99 + b99_ihit + bbh, data=data, id = "id", snps = c("SNPHE"), infile = "/home/emil/tmp/runGLMM.snpHE.tab",infile.ncol.skip = 3,kins = as.matrix(GRM))
            wald2<-GMMAT::glmm.wald(disease ~ age + sex + genoHE + b99 + b99_ihit + bbh, data=data, id = "id", snps = c("SNPHO"), infile = "/home/emil/tmp/runGLMM.snpHO.tab",infile.ncol.skip = 3,kins = as.matrix(GRM))
            
        } else{
            
            wald1<-GMMAT::glmm.wald(disease ~ age + sex + b99 + b99_ihit + bbh, data=data, id = "id", snps = c("SNP1"), infile = "/home/emil/tmp/runGLMM.snp2.tab",infile.ncol.skip = 3,kins = as.matrix(GRM))
        }
                
    } else if(all(c("b99","b99_ihit")%in%colnames(data))){
        print("9")

        if(genoModel==2){
            
            wald1<-GMMAT::glmm.wald(disease ~ age + sex + genoHO + b99 + b99_ihit, data=data, id = "id", snps = c("SNPHE"), infile = "/home/emil/tmp/runGLMM.snpHE.tab",infile.ncol.skip = 3,kins = as.matrix(GRM))
            wald2<-GMMAT::glmm.wald(disease ~ age + sex + genoHE + b99 + b99_ihit, data=data, id = "id", snps = c("SNPHO"), infile = "/home/emil/tmp/runGLMM.snpHO.tab",infile.ncol.skip = 3,kins = as.matrix(GRM))
            
        } else{
            
            wald1<-GMMAT::glmm.wald(disease ~ age + sex + b99 + b99_ihit, data=data, id = "id", snps = c("SNP1"), infile = "/home/emil/tmp/runGLMM.snp2.tab",infile.ncol.skip = 3,kins = as.matrix(GRM))
        }

            
    } else if(all(c("b99","bbh")%in%colnames(data))){
        print("10")

        if(genoModel==2){
            
            wald1<-GMMAT::glmm.wald(disease ~ age + sex + genoHO + b99 + bbh, data=data, id = "id", snps = c("SNPHE"), infile = "/home/emil/tmp/runGLMM.snpHE.tab",infile.ncol.skip = 3,kins = as.matrix(GRM))
            wald2<-GMMAT::glmm.wald(disease ~ age + sex + genoHE + b99 + bbh, data=data, id = "id", snps = c("SNPHO"), infile = "/home/emil/tmp/runGLMM.snpHO.tab",infile.ncol.skip = 3,kins = as.matrix(GRM))
            
        } else{

            wald1<-GMMAT::glmm.wald(disease ~ age + sex + b99 + bbh, data=data, id = "id", snps = c("SNP1"), infile = "/home/emil/tmp/runGLMM.snp2.tab",infile.ncol.skip = 3,kins = as.matrix(GRM))
        }
            
    } else if(all(c("b99")%in%colnames(data))){
        print("11")

        if(genoModel==2){
            
            wald1<-GMMAT::glmm.wald(disease ~ age + sex + genoHO + b99, data=data, id = "id", snps = c("SNPHE"), infile = "/home/emil/tmp/runGLMM.snpHE.tab",infile.ncol.skip = 3,kins = as.matrix(GRM))
            wald2<-GMMAT::glmm.wald(disease ~ age + sex + genoHE + b99, data=data, id = "id", snps = c("SNPHO"), infile = "/home/emil/tmp/runGLMM.snpHO.tab",infile.ncol.skip = 3,kins = as.matrix(GRM))
            
        } else{
        
            wald1<-GMMAT::glmm.wald(disease ~ age + sex + b99, data=data, id = "id", snps = c("SNP1"), infile = "/home/emil/tmp/runGLMM.snp2.tab",infile.ncol.skip = 3,kins = as.matrix(GRM))
        }
            
    } else if(all(c("b99_ihit")%in%colnames(data))){
        print("12")
        if(genoModel==2){
            
            wald1<-GMMAT::glmm.wald(disease ~ age + sex + genoHO + b99_ihit, data=data, id = "id", snps = c("SNPHE"), infile = "/home/emil/tmp/runGLMM.snpHE.tab",infile.ncol.skip = 3,kins = as.matrix(GRM))
            wald2<-GMMAT::glmm.wald(disease ~ age + sex + genoHE + b99_ihit, data=data, id = "id", snps = c("SNPHO"), infile = "/home/emil/tmp/runGLMM.snpHO.tab",infile.ncol.skip = 3,kins = as.matrix(GRM))
            
        } else{
            
            wald1<-GMMAT::glmm.wald(disease ~ age + sex + b99_ihit, data=data, id = "id", snps = c("SNP1"), infile = "/home/emil/tmp/runGLMM.snp2.tab",infile.ncol.skip = 3,kins = as.matrix(GRM))
        }
            
    } else if(all(c("bbh")%in%colnames(data))){
        print("13")
        if(genoModel==2){
            
            wald1<-GMMAT::glmm.wald(disease ~ age + sex + genoHO + bbh, data=data, id = "id", snps = c("SNPHE"), infile = "/home/emil/tmp/runGLMM.snpHE.tab",infile.ncol.skip = 3,kins = as.matrix(GRM))
            wald2<-GMMAT::glmm.wald(disease ~ age + sex + genoHE + bbh, data=data, id = "id", snps = c("SNPHO"), infile = "/home/emil/tmp/runGLMM.snpHO.tab",infile.ncol.skip = 3,kins = as.matrix(GRM))
            
        } else{
        
            wald1<-GMMAT::glmm.wald(disease ~ age + sex + bbh, data=data, id = "id", snps = c("SNP1"), infile = "/home/emil/tmp/runGLMM.snp2.tab",infile.ncol.skip = 3,kins = as.matrix(GRM))
        }
            
    } else{
        print("14")
        if(genoModel==2){
            
            wald1<-GMMAT::glmm.wald(disease ~ age + sex + genoHO, data=data, id = "id", snps = c("SNPHE"), infile = "/home/emil/tmp/runGLMM.snpHE.tab",infile.ncol.skip = 3,kins = as.matrix(GRM))
            wald2<-GMMAT::glmm.wald(disease ~ age + sex + genoHE, data=data, id = "id", snps = c("SNPHO"), infile = "/home/emil/tmp/runGLMM.snpHO.tab",infile.ncol.skip = 3,kins = as.matrix(GRM))
            
        } else{

            wald1<-GMMAT::glmm.wald(disease ~ age + sex, data=data, id = "id", snps = c("SNP1"), infile = "/home/emil/tmp/runGLMM.snp2.tab",infile.ncol.skip = 3,kins = as.matrix(GRM))
        }

    }


    if(genoModel==2){
        write.table(wald1,paste0("results/",pheno,".",snp,".adjCohort.fullHE.wald.V7.res"),col=T,row=F,qu=F)
        write.table(wald2,paste0("results/",pheno,".",snp,".adjCohort.fullHO.wald.V7.res"),col=T,row=F,qu=F)                
    } else if(genoModel==1){        
        write.table(wald1,paste0("results/",pheno,".",snp,".adjCohort.recessive.wald.V7.res"),col=T,row=F,qu=F)        
    } else{
        write.table(wald1,paste0("results/",pheno,".",snp,".adjCohort.wald.V7.res"),col=T,row=F,qu=F)        
    }


}
