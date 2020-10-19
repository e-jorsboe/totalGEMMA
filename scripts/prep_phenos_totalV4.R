### based on /home/emil/OMNI5/ricco_scripts/prep_phenos_bimbam_newV2.R

l<-commandArgs(TRUE)

####

trans<-l[1] ## == 2 binary pheno
phe<-l[2] ## phenotype
cov<-l[3]
cohort<-l[4]
out<-l[5]
phenofile<-l[6]
genofile<-l[7]
version<-l[8]
year<-l[9]
env<-l[10]
source("scripts/transfunsV4.R")

#####
## 
## trans<-"1"
## phe<-"bmi"
## cov<-"sex,age,B99_IHIT"
## cohort<-"IHIT"
## out=paste0("multiPheno9_",phe,"_addtrans")
## phenofile="/home/emil/OMNI5/megaChipQC/Greenland_Pheno_23maj2016newBarcode_PAX_WGS_MEGA_inMetabo_inMega.csv"
## genofile="/emc/emil/riccoStuff/gemmaTotal/data/megaChip/megaV2_autosOnlySNPsSwitchedID"
## version=1
## year=2017
## env="NULL"

phe1<-gsub(",","_",phe)

if(grepl(pattern = "\\.mgf$",genofile)){    
    fam<-read.table(gsub("RECESSIVE","ADDITIVE",gsub("mgf$","fam",genofile)),as.is=T)  
} else{
    fam<-read.table(paste0(genofile,".fam"),as.is=T)  
}

fam$orgOrder<-1:nrow(fam)

phenos<-read.csv(phenofile,as.is=T)

## NOT all indidivuals from .fam are pheno file necesarilly!!
int<-intersect(phenos$particid,fam$V1)

phenos<-phenos[phenos$particid %in% fam$V1,]
fam<-fam[ fam$V1%in%phenos$particid,]

## EVEN THOUGH FAM is ordered like this in the first place
fam2<-fam[order(fam$V1),]
phenos<-phenos[order(phenos$particid),]

if(any(grepl("cohort",cov))){

    cohorts<-cbind(as.numeric(phenos$cohort=="IHIT"),as.numeric(phenos$cohort=="B99"),as.numeric(phenos$cohort=="BBH"),as.numeric(phenos$cohort=="B99_IHIT"))
    colnames(cohorts)<-c("IHIT","B99","BBH","B99_IHIT")

    phenos2<-cbind(phenos,cohorts)

    cov<-sub("cohort","B99,BBH,B99_IHIT",cov)
    cov<-unlist(strsplit(cov,","))
}
        
trans<-as.numeric(unlist(strsplit(trans,",")))
phe<-unlist(strsplit(phe,","))

if(env!="NULL"){
    pheCov<-c(phe,cov,env)
} else{
    pheCov<-c(phe,cov)
}

pheCov<-pheCov[pheCov!="" &!is.na(pheCov)]

print(pheCov)

if(!all(c(pheCov)%in%c(names(phenos2)))){
  cat("Both the phenotype and the covs need to be in the phenotype file\n")
  q()
}

if(toupper(cohort)=="ALL"){
    
    ## looks at number rows of covar and sees if there is just one missing
    keep<-!apply(is.na(phenos2[,pheCov,drop=FALSE]),1,any)
    cat("Number of ind with non missing phentypes ",sum(keep),"\n")
    
    if(all(phenos2[keep,"B99"]==0)){
        pheCov<-pheCov[pheCov!="B99"]
        cov<-cov[cov!="B99"]
    }
    
    ## IF no BBH remove from covars
    if(all(phenos2[keep,"BBH"]==0)){
        pheCov<-pheCov[pheCov!="BBH"]
        cov<-cov[cov!="BBH"]
    }
    
    if(all(phenos2[keep,"B99_IHIT"]==0)){
        pheCov<-pheCov[pheCov!="B99_IHIT"]
        cov<-cov[cov!="B99_IHIT"]
    }
    
} else{

    ## we do not adjust for cohort when only analysing one cohort
    
    keep2<-phenos2$cohort2==toupper(cohort)
    keep<-!apply(is.na(phenos2[,pheCov,drop=FALSE]),1,any)
    keep<-keep & keep2
    cat("Number of ind with non missing phentypes ",sum(keep),"\n")
    

    pheCov<-pheCov[pheCov!="B99"]
    pheCov<-pheCov[pheCov!="IHIT"]
    pheCov<-pheCov[pheCov!="BBH"]
    pheCov<-pheCov[pheCov!="B99_IHIT"]
    
    cov<-cov[cov!="B99"]
    cov<-cov[cov!="BBH"]
    cov<-cov[cov!="IHIT"]   
    cov<-cov[cov!="B99_IHIT"]
    
}

if(sum(keep)==0){
    cat("No indis to analyse for\n")
    write("",file=paste0("run",year,"_v",version,"/",phe1,"/",out,".empty"))  
    q()
    
}

if(any(c("marine_mammals_energi","pct_grl_diet","Traditional_diet")%in%cov)){

    keep2<-phenos2[,"energi_ok"]==1
    keep<-keep & keep2

}

##fattyAcidPhenos<-read.table("data/fattyAcidPhenos.list",as.is=T)

##if(any(phe%in%fattyAcidPhenos)){

##  keep2<-phenos2[,"energi_ok"]==1
##  keep<-keep & keep2
##}

phe_all<-list()
phenoIndex<-1

for(q in trans){

    if(q=="0"){
        cat("Untransformed phe is used\n")
        phe_all[[phenoIndex]]<-phenos2[keep,phe[phenoIndex]]        
    } else{
        if(q=="2"){
            cat("transformed phe (sexes pooled) is used\n")
            phe_all[[phenoIndex]]<-randomqnormtransform(phenos2[keep ,phe[phenoIndex]])                                     
        } else{
            if(q=="1"){
                cat("transformed phe (sexes sep) is used\n")
                
                phe_all[[phenoIndex]]<-phenos2[keep,c(phe[phenoIndex],"sex")]
                phe_all[[phenoIndex]][ phe_all[[phenoIndex]]$sex==1, phe[phenoIndex]]<-randomqnormtransform(phenos2[keep & phenos2$sex==1,phe[phenoIndex]]) ### trans within sexes
                phe_all[[phenoIndex]][ phe_all[[phenoIndex]]$sex==2, phe[phenoIndex]]<-randomqnormtransform(phenos2[keep & phenos2$sex==2,phe[phenoIndex]]) ### trans within sexes                
                phe_all[[phenoIndex]]<-phe_all[[phenoIndex]][,1]                
            } else{
                if(q=="3"){
                    cat("cctransformed phe is used\n")
                    ## CCTRANSFORMED ASSUMES THAT THEY ARE ALL 0 and !!!!!!!!!
                    phe_all[[phenoIndex]]<-cctransformV2(phenos2[keep ,phe[phenoIndex]])                  
                } else{
                    cat("Error: trans type not known!\n")
                    q()
                }
            }
        }
    }
    phenoIndex<-phenoIndex+1
}

## fam2 has SAME ordering as phenos
fam2All<-cbind(fam2[keep,1:5],do.call(cbind,phe_all),fam2[keep,"orgOrder"])
## imperateive that .fam file has orignal order from plink files
## as otherwise genotypes are put on wrong
fam2All<-fam2All[ order(fam2All[,ncol(fam2All)]),]
fam2All<-fam2All[,1:(ncol(fam2All)-1)]

## check that B99 and B99_IHIT are not perfectly inversly correlated

print("EMIL:")
print(cov)

if(any(grepl("^B99$",cov)) & any(grepl("B99_IHIT",cov))){

    print("DO I GO HERE - EMIL")
    if(cor(phenos2[keep,"B99"],phenos2[keep,"B99_IHIT"])**2==1){
        cov<-cov[cov!="B99_IHIT"]
    }
}

print("AFTER:")
print(cov)

covars<-cbind(1,phenos2[keep,cov],fam2[keep,"orgOrder"])
covars<-covars[ order(covars[,ncol(covars)]),]
covars<-covars[,1:(ncol(covars)-1)]

if(grepl(pattern = "\\.mgf$",genofile)){

  ## those to cut from .mgf dosage file
  write(paste(fam2[keep,"orgOrder"]+3,collapse = ","),file=paste0("run",year,"_v",version,"/",phe1,"/",out,"_cols.txt"))    
  write.table(fam2All[,1:2],file=paste0("run",year,"_v",version,"/",phe1,"/",out,"_IDs.txt"),row=F,col=F,qu=F)    
  ## pheno, where fam file has to have N phenos - 1 extra columns, for multi phenotype analysis
  write.table(fam2All[,6:ncol(fam2All)],file=paste0("run",year,"_v",version,"/",phe1,"/",out,"_pheno.txt"),row=F,col=F,qu=F)
  ## covariates
  write.table(covars,file=paste0("run",year,"_v",version,"/",phe1,"/",out,".cov"),row=F,col=F,qu=F)    
} else{

  ## those to keep in plink file
  write.table(fam2All[,1:6],file=paste0("run",year,"_v",version,"/",phe1,"/tmp",out,".fam"),row=F,col=F,qu=F)
  ## pheno, where fam file has to have N phenos - 1 extra columns, for multi phenotype analysis
  write.table(fam2All[,1:ncol(fam2All)],file=paste0("run",year,"_v",version,"/",phe1,"/pheno",out,".fam"),row=F,col=F,qu=F)
  ## covariates
  write.table(covars,file=paste0("run",year,"_v",version,"/",phe1,"/",out,".cov"),row=F,col=F,qu=F) 

}

if(env!="NULL"){
    enviroment<-cbind(phenos2[keep,env],fam2[keep,"orgOrder"])
    enviroment<-enviroment[ order(enviroment[,ncol(enviroment)]),]
    enviroment<-enviroment[,1:(ncol(enviroment)-1)]
    write.table(enviroment,file=paste0("run",year,"_v",version,"/",phe1,"/",out,".env"),row=F,col=F,qu=F) 
}
