#!/bin/bash

## do one script for each cohort, like runALL.sh, runIHIT.sh
## run from /emc/emil/riccoStuff/gemmaTotal/
## for more see /emc/emil/riccoStuff/gemmaTotal/README

## this script should have these arguments
COHORT=$1
version=$2
YEAR=$3
## if multiple phenos PHENO1,PHENO2,... then with multiple phenotypes, should also have multiple qTrans 1,3,...
PHENOFILE=$4
## has to be the additive one for dosages, meaning ADDITIVE*.mgf
GENOFILE=$5

## /ricco/porno/emil/imputation/analysisdescription_impu2016_v${version}.txt
DESCRIPTION=$6
CORES=$7
ONLYTRANS=${8:-"0"}
SUFFIX=${9:-NONE}
## list of SNPs to be removed based on QC, rs IDs (optional)
## assumed to be an .Rdata file like this "/net/pontus/pontus/data/anders/shiny/dataGreenland2015c/csvInfo"
## object loaded assumed to be called csvInfo, with passFilter (if passed filter) and snpID (same as GENOFILE) columns
## or can be used for designating which sites (list if IDs) to include for the GSM calculations - if ends with ".gsm"
BADSNPS=${10:-NONE}
## if only want to run with subset of SNPs - has to be list of SNP-IDs that GEMMA can take as parameter
## has to be absolute path to file
SNPS=${11:-NONE}
GSMPLINK=${12:-NONE}

########################

SLP=60

mkdir -p run${YEAR}_v${version}
mkdir -p nohups

if `echo $GENOFILE | grep -q "\\.gz$"`; then
    GENOFILE=`echo $GENOFILE | sed -e 's/\\.gz$//g'`
fi


if `echo $GENOFILE | grep -q "\\.mgf$"`; then
    pigz -d ${GENOFILE}.gz
    pigz -d `echo ${GENOFILE}.gz | sed -e 's/ADDITIVE/RECESSIVE/g'`
fi


## first do additives, because they are the ones where the sXX matrix is calculated
## Run through all phenos in analyses description file and do adduntrans 
while read -a traitinfo; do    
    cnt1=$(top -b -u emil -n1 | grep emil | grep -w "^R$" | wc -l)
    cnt2=$(top -b -u emil -n1 | grep emil | grep plink | wc -l)
    cnt3=$(top -b -u emil -n1 | grep emil | grep gemma | wc -l)
    while [ $(($cnt1 + $cnt2 + $cnt3)) -gt $(($CORES-1)) ]
    do
        sleep ${SLP}
	cnt1=$(top -b -u emil -n1 | grep emil | grep -w "^R$" | wc -l)
	cnt2=$(top -b -u emil -n1 | grep emil | grep plink | wc -l)
	cnt3=$(top -b -u emil -n1 | grep emil | grep gemma | wc -l)
    done

    PHE=${traitinfo[0]}
    PHE=`echo $PHE | sed -e 's/,/ \n/g' | sort | tr '\n' ',' | sed -e 's/ //g' |  sed -e 's/,$//g'`
    PHE1=`echo $PHE | sed -e 's/,/_/g'`

    REC=0
    TRANS=${traitinfo[1]}    
    COV=${traitinfo[2]}

    TRANS_REC="${TRANS}_${REC}"
    
    ## Additive tranformed        
    if `echo $GENOFILE | grep -q "\\.mgf$"`; then    
	OUT=run${YEAR}_v${version}_addtrans_${COHORT}_${PHE1}_dosages
    else	
	OUT=run${YEAR}_v${version}_addtrans_${COHORT}_${PHE1}	
    fi
    
    if [ ${#traitinfo[@]} -gt 3 ]; then
	ENV=${traitinfo[3]}
	OUT="${OUT}_gxe${ENV}"
    else
	ENV="NULL"       
    fi

    if [ "$SUFFIX" != "NONE" ]; then
	OUT="${OUT}_${SUFFIX}"
    fi
    
    if [ "$GSMPLINK" != "NONE" ]; then
	GENOFILE="${GENOFILE}GSMFILE${GSMPLINK}"
    fi
    
    echo $PHE $COV $TRANS_REC $ENV $PHENOFILE $GENOFILE $OUT $BADSNPS $SNPS
    
    ## in order to wait for iniating jobs so can be counted
done < $DESCRIPTION | xargs -n9 -P${CORES} ./runGEMMA_singleRunV4pScore.sh

while read -a traitinfo; do
    cnt1=$(top -b -u emil -n1 | grep emil | grep -w "^R$" | wc -l)
    cnt2=$(top -b -u emil -n1 | grep emil | grep plink | wc -l)
    cnt3=$(top -b -u emil -n1 | grep emil | grep gemma | wc -l)
    
    while [ $(($cnt1 + $cnt2 + $cnt3)) -gt $(($CORES-1)) ]
    do
        sleep ${SLP}
	cnt1=$(top -b -u emil -n1 | grep emil | grep -w "^R$" | wc -l)
	cnt2=$(top -b -u emil -n1 | grep emil | grep plink | wc -l)
	cnt3=$(top -b -u emil -n1 | grep emil | grep gemma | wc -l)	
    done
      
    PHE=${traitinfo[0]}
    PHE=`echo $PHE | sed -e 's/,/ \n/g' | sort | tr '\n' ',' | sed -e 's/ //g' |  sed -e 's/,$//g'`
    PHE1=`echo $PHE | sed -e 's/,/_/g'`
    ## Recessive tranformed
    REC=1
    TRANS=${traitinfo[1]}
    COV=${traitinfo[2]}

    TRANS_REC="${TRANS}_${REC}"
        
    ## .mgf GENOFILES should be ADDITIVE*.mgf and RECESSIVE*.mgf respectively
    GENOFILE=`echo $GENOFILE | sed -e 's/ADDITIVE/RECESSIVE/g'`          
    if `echo $GENOFILE | grep -q "\\.mgf$"`; then 	  
	OUT=run${YEAR}_v${version}_rectrans_${COHORT}_${PHE1}_dosages
    else	  
	OUT=run${YEAR}_v${version}_rectrans_${COHORT}_${PHE1}	  
    fi
    
    if [ ${#traitinfo[@]} -gt 3 ]; then
	ENV=${traitinfo[3]}
	OUT="${OUT}_gxe${ENV}"
    else
	ENV="NULL"       
    fi
        
    if [ "$SUFFIX" != "NONE" ]; then
	OUT="${OUT}_${SUFFIX}"
    fi

    if [ "$GSMPLINK" != "NONE" ]; then
	GENOFILE="${GENOFILE}GSMFILE${GSMPLINK}"
    fi
    
    echo $PHE $COV $TRANS_REC $ENV $PHENOFILE $GENOFILE $OUT $BADSNPS $SNPS
        
done < $DESCRIPTION | xargs -n9 -P${CORES} ./runGEMMA_singleRunV4pScore.sh

if [ $ONLYTRANS != "0" ]; then
    echo "Only transformed analyses will be run!"
    exit 1
fi

while read -a traitinfo; do
    cnt1=$(top -b -u emil -n1 | grep emil | grep -w "^R$" | wc -l)
    cnt2=$(top -b -u emil -n1 | grep emil | grep plink | wc -l)
    cnt3=$(top -b -u emil -n1 | grep emil | grep gemma | wc -l)
    
    while [ $(($cnt1 + $cnt2 + $cnt3))  -gt $(($CORES-1)) ]
    do
        sleep ${SLP}
	cnt1=$(top -b -u emil -n1 | grep emil | grep -w "^R$" | wc -l)
	cnt2=$(top -b -u emil -n1 | grep emil | grep plink | wc -l)
	cnt3=$(top -b -u emil -n1 | grep emil | grep gemma | wc -l)	 
    done

    PHE=${traitinfo[0]}
    ## sorts the phenotypes so cannot have more runs with the same phenotypes (just sorted differently)
    PHE=`echo $PHE | sed -e 's/,/ \n/g' | sort | tr '\n' ',' | sed -e 's/ //g' |  sed -e 's/,$//g'`
    PHE1=`echo $PHE | sed -e 's/,/_/g'`

    REC=0
    ## ok with just 0 if multile phenos, R script will take care of this
    TRANS=`echo ${traitinfo[1]} | sed -e 's/1/0/g' | sed -e 's/2/0/g' | sed -e 's/3/0/g'`     
    COV=${traitinfo[2]}

    TRANS_REC="${TRANS}_${REC}"

    GENOFILE=`echo $GENOFILE | sed -e 's/RECESSIVE/ADDITIVE/g'`
    if `echo $GENOFILE | grep -q "\\.mgf$"`; then     
	OUT=run${YEAR}_v${version}_adduntrans_${COHORT}_${PHE1}_dosages
    else	
	OUT=run${YEAR}_v${version}_adduntrans_${COHORT}_${PHE1}
    fi

    ## if row has more than 3 elements
    if [ ${#traitinfo[@]} -gt 3 ]; then
	ENV=${traitinfo[3]}
	OUT="${OUT}_gxe${ENV}"
    else
	ENV="NULL"       
    fi

    if [ "$SUFFIX" != "NONE" ]; then
	OUT="${OUT}_${SUFFIX}"
    fi
    
    if [ "$GSMPLINK" != "NONE" ]; then
	GENOFILE="${GENOFILE}GSMFILE${GSMPLINK}"
    fi
    
    echo $PHE $COV $TRANS_REC $ENV $PHENOFILE $GENOFILE $OUT $BADSNPS $SNPS
    ## in order to wait for iniating jobs so can be counted
    
done < $DESCRIPTION | xargs -n9 -P${CORES} ./runGEMMA_singleRunV4pScore.sh
## analysis description, like Ida's

while read -a traitinfo; do
    cnt1=$(top -b -u emil -n1 | grep emil | grep -w "^R$" | wc -l)
    cnt2=$(top -b -u emil -n1 | grep emil | grep plink | wc -l)
    cnt3=$(top -b -u emil -n1 | grep emil | grep gemma | wc -l)       
    while [ $(($cnt1 + $cnt2 + $cnt3)) -gt $(($CORES-1)) ]
    do
        sleep ${SLP}
	cnt1=$(top -b -u emil -n1 | grep emil | grep -w "^R$" | wc -l)
	cnt2=$(top -b -u emil -n1 | grep emil | grep plink | wc -l)
	cnt3=$(top -b -u emil -n1 | grep emil | grep gemma | wc -l)
    done
    
    PHE=${traitinfo[0]}
    PHE=`echo $PHE | sed -e 's/,/ \n/g' | sort | tr '\n' ',' | sed -e 's/ //g' |  sed -e 's/,$//g'`
    PHE1=`echo $PHE | sed -e 's/,/_/g'`
    
    ## Recessive untranformed
    REC=1
    ## ok with just 0 if multile phenos, R script will take care of this
    TRANS=`echo ${traitinfo[1]} | sed -e 's/1/0/g' | sed -e 's/2/0/g' | sed -e 's/3/0/g'`
    COV=${traitinfo[2]}

    TRANS_REC="${TRANS}_${REC}"
       
    ## .mgf GENOFILES should be ADDITIVE*.mgf and RECESSIVE*.mgf respectively
    GENOFILE=`echo $GENOFILE | sed -e 's/ADDITIVE/RECESSIVE/g'`    
    if `echo $GENOFILE | grep -q "\\.mgf$"`; then     
	OUT=run${YEAR}_v${version}_recuntrans_${COHORT}_${PHE1}_dosages
    else	
	OUT=run${YEAR}_v${version}_recuntrans_${COHORT}_${PHE1}	
    fi
    
    ## gxe does not work for recessive model
    if [ ${#traitinfo[@]} -gt 3 ]; then
	ENV=${traitinfo[3]}
	OUT="${OUT}_gxe${ENV}"
    else
	ENV="NULL"       
    fi
    
    if [ "$SUFFIX" != "NONE" ]; then
	OUT="${OUT}_${SUFFIX}"
    fi

    if [ "$GSMPLINK" != "NONE" ]; then
	GENOFILE="${GENOFILE}GSMFILE${GSMPLINK}"
    fi
    
    echo $PHE $COV $TRANS_REC $ENV $PHENOFILE $GENOFILE $OUT $BADSNPS $SNPS
    
    ## in order to wait for iniating jobs so can be counted
done < $DESCRIPTION | xargs -n9 -P${CORES} ./runGEMMA_singleRunV4pScore.sh


if `echo $GENOFILE | grep -q "\\.mgf$"`; then
    pigz ${GENOFILE}
    pigz `echo $GENOFILE | sed -e 's/ADDITIVE/RECESSIVE/g'`
fi
