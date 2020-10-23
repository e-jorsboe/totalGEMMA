#!/bin/bash

PHE=$1
COV=$2
TRANS_REC=$3
ENV=$4
PHENOFILE=$5
GENOFILE=$6
OUT=$7
BADSNPS=$8
SNPS=$9
##version=${10}
##YEAR=${11}



if [ "$OUT" == "NONE" ]; then
    ##emil
    echo "I AM HERE!"
    exit 1
fi
   
## set default something!!!
MODEL=`echo $OUT | cut -f3 -d"_"` 
if [ "$ENV" != "NULL" ] && ([ "$MODEL" == "rectrans" ] || [ "$MODEL" == "recuntrans" ]); then
    ## gxe does not work for recessive
    exit 1

    #####################################################
    ## THIS MEANT THAT untrans analyses were not run 
    #####################################################
    
    ##elif [ "$TRANS" == "0" ] && ([ `echo $OUT | cut -f3 -d"_"` == "adduntrans" ] || [ `echo $OUT | cut -f3 -d"_"` == "recuntrans" ]) ; then
    ## if pheno does not have to be transformed (binary pheno) does not need to run untransed model!
    ## AS NOW I GENERATE THE GSM IN THE TRANSED MODEL - WHICH ARE RUN FIRST!!
    ##exit 1
    
fi

GSMGENO="0"

if `echo $GENOFILE | grep -q "GSMFILE" `; then
    GENOFILEtmp=`echo $GENOFILE | sed -e 's/GSMFILE/^/g'`
    GENOFILE1=`echo $GENOFILEtmp | cut -f1 -d"^"`
    GENOFILE2=`echo $GENOFILEtmp | cut -f2 -d"^"`
    GENOFILE=$GENOFILE1
    GSMGENO="1"
fi


## sXX can be calculated regardless model, same for adduntrans (calculated for this), addtrans, recuntrans, rectrans
OUTsXX="`echo $OUT | cut -f1,2 -d"_"`_`echo $OUT | cut -f4- -d"_"`"

TRANS=`echo $TRANS_REC | cut -f1 -d"_"`
REC=`echo $TRANS_REC | cut -f2 -d"_"`

version=`echo $OUT | cut -f2 -d"_" | sed -e 's/v//g'`
YEAR=`echo $OUT | cut -f1 -d"_" | sed -e 's/run//g'`
COHORT=`echo $OUT | cut -f4 -d"_"`

SUFFIX=`echo $OUT | rev | cut -f1 -d"_" | rev`

if [ "$SUFFIX" == "gxe${ENV}" ]; then
    SUFFIX=""
else
    SUFFIX="_${SUFFIX}"
fi

echo "single file"
echo $PHE $COV $TRANS_REC $ENV $PHENOFILE $GENOFILE $OUT $BADSNPS $SNPS $version $YEAR $COHORT
echo "end single file"

PHE1=`echo $PHE | sed -e 's/,/_/g'`
mkdir -p  run${YEAR}_v${version}/$PHE1
> nohups/$OUT

## 0 = Untransed, 2 trans (sexes pooled), 1 trans (sexes sep), 3 made binary to 1 and 2, (0 -> 1  and !=0 -> 2)
Rscript scripts/prep_phenos_totalV4.R $TRANS $PHE $COV $COHORT $OUT $PHENOFILE $GENOFILE $version $YEAR $ENV > nohups/$OUT 2>&1

## if no indis in analysed dataset - it quits
if [ -f run${YEAR}_v${version}/${PHE1}/${OUT}.empty ]; then
    rm run${YEAR}_v${version}/${PHE1}/${OUT}.empty
    echo "no Indis in analyses - quitting!"
    exit 1
fi

if `echo $GENOFILE | grep -q "\\.mgf$"`; then
    ## to create file that STDOUT/STDERR is being written to    
    COLS=`cat run${YEAR}_v${version}/${PHE1}/${OUT}_cols.txt`
    ## cut columns from list of columns to cut

    if [ "$GSMGENO" == "1" ]; then
	cut -f 1,2,3,$COLS -d "," $GENOFILE2 > run${YEAR}_v${version}/${PHE1}/${OUT}gsm.mgf
    fi
    
    cut -f 1,2,3,$COLS -d "," $GENOFILE > run${YEAR}_v${version}/${PHE1}/${OUT}.mgf
    
    ## get correlation/covariance matrix
    if [ "$MODEL" == "addtrans" ]; then
	if `echo $BADSNPS | grep -q "\\.gsm$"`; then

	    if [ "$GSMGENO" == "1" ]; then
		prog/gemma.linux -g run${YEAR}_v${version}/${PHE1}/${OUT}gsm.mgf -p run${YEAR}_v${version}/${PHE1}/${OUT}_pheno.txt -gk 2 -o ${OUTsXX} -maf 0.0 -miss 0.8 -snps $BADSNPS >> nohups/$OUT 2>&1
	    else
		prog/gemma.linux -g run${YEAR}_v${version}/${PHE1}/${OUT}.mgf -p run${YEAR}_v${version}/${PHE1}/${OUT}_pheno.txt -gk 2 -o ${OUTsXX} -maf 0.0 -miss 0.8 -snps $BADSNPS >> nohups/$OUT 2>&1
	    fi
	    
	    K=output/${OUTsXX}.sXX.txt
	    
	else

	    if [ "$GSMGENO" == "1" ]; then
		prog/gemma.linux -g run${YEAR}_v${version}/${PHE1}/${OUT}gsm.mgf -p run${YEAR}_v${version}/${PHE1}/${OUT}_pheno.txt -gk 2 -o ${OUTsXX} -maf 0.05 -miss 0.01 >> nohups/$OUT 2>&1
	    else
	    ## sXX should be NAMED *_dosages because one for all 4 runs	    
		prog/gemma.linux -g run${YEAR}_v${version}/${PHE1}/${OUT}.mgf -p run${YEAR}_v${version}/${PHE1}/${OUT}_pheno.txt -gk 2 -o ${OUTsXX} -maf 0.05 -miss 0.01 >> nohups/$OUT 2>&1
	    fi
		
	    K=output/${OUTsXX}.sXX.txt
	fi	    
    else
	    K=output/${OUTsXX}.sXX.txt		

    fi
    n=`head -n 1 run${YEAR}_v${version}/${PHE1}/${OUT}_pheno.txt | wc -w`    
    ## choose model
    if [ "$REC" == "1" ]; then	
	if [ "$SNPS" != "NONE" ]; then
	    echo "RECModel"
	    prog/gemma.linux -g run${YEAR}_v${version}/${PHE1}/${OUT}.mgf -a `echo $GENOFILE | sed -e 's/mgf$/ann/g' | sed -e 's/RECESSIVE/ADDITIVE/g'` -p run${YEAR}_v${version}/${PHE1}/${OUT}_pheno.txt -k $K -lmm 3 -o ${OUT} -c run${YEAR}_v${version}/${PHE1}/${OUT}.cov -n `seq 1 $n` -maf 0 -miss 1 -r2 1 -snps $SNPS >> nohups/$OUT 2>&1
	else	   	
	    echo "RECModel"
	    prog/gemma.linux -g run${YEAR}_v${version}/${PHE1}/${OUT}.mgf -a `echo $GENOFILE | sed -e 's/mgf$/ann/g' | sed -e 's/RECESSIVE/ADDITIVE/g'` -p run${YEAR}_v${version}/${PHE1}/${OUT}_pheno.txt -k $K -lmm 3 -o ${OUT} -c run${YEAR}_v${version}/${PHE1}/${OUT}.cov -n `seq 1 $n` -maf 0 -miss 1 -r2 1 >> nohups/$OUT 2>&1
	fi
	    
    elif [ "$ENV" != "NULL" ]; then
	if [ "$SNPS" != "NONE" ]; then
	     prog/gemma.linux -g run${YEAR}_v${version}/${PHE1}/${OUT}.mgf -a `echo $GENOFILE | sed -e 's/mgf$/ann/g' | sed -e 's/RECESSIVE/ADDITIVE/g'` -p run${YEAR}_v${version}/${PHE1}/${OUT}_pheno.txt -k $K -lmm 3 -o ${OUT} -c run${YEAR}_v${version}/${PHE1}/${OUT}.cov -n `seq 1 $n` -maf 0 -miss 1 -r2 1 -snps $SNPS -gxe run${YEAR}_v${version}/${PHE1}/${OUT}.env  >> nohups/$OUT 2>&1
	else
	    ## can only do score regression
	    prog/gemma.linux -g run${YEAR}_v${version}/${PHE1}/${OUT}.mgf -a `echo $GENOFILE | sed -e 's/mgf$/ann/g' | sed -e 's/RECESSIVE/ADDITIVE/g'` -p run${YEAR}_v${version}/${PHE1}/${OUT}_pheno.txt -k $K -lmm 3 -o ${OUT} -c run${YEAR}_v${version}/${PHE1}/${OUT}.cov -n `seq 1 $n` -maf 0 -miss 1 -r2 1 -gxe run${YEAR}_v${version}/${PHE1}/${OUT}.env >> nohups/$OUT 2>&1	   
	fi
	
    else	
	if [ "$SNPS" != "NONE" ]; then
	    prog/gemma.linux -g run${YEAR}_v${version}/${PHE1}/${OUT}.mgf -a `echo $GENOFILE | sed -e 's/mgf$/ann/g' | sed -e 's/RECESSIVE/ADDITIVE/g'` -p run${YEAR}_v${version}/${PHE1}/${OUT}_pheno.txt -k $K -lmm 3 -o ${OUT} -c run${YEAR}_v${version}/${PHE1}/${OUT}.cov -n `seq 1 $n` -maf 0 -miss 1 -r2 1 -snps $SNPS >> nohups/$OUT 2>&1	    
	else
	    prog/gemma.linux -g run${YEAR}_v${version}/${PHE1}/${OUT}.mgf -a `echo $GENOFILE | sed -e 's/mgf$/ann/g' | sed -e 's/RECESSIVE/ADDITIVE/g'` -p run${YEAR}_v${version}/${PHE1}/${OUT}_pheno.txt -k $K -lmm 3 -o ${OUT} -c run${YEAR}_v${version}/${PHE1}/${OUT}.cov -n `seq 1 $n` -maf 0 -miss 1 -r2 1 >> nohups/$OUT 2>&1
	fi
    fi
        
    if [ "$SNPS" == "NONE" ]; then
	rm run${YEAR}_v${version}/${PHE1}/${OUT}.mgf
	Rscript scripts/generate_manhattan_pScore.R output/${OUT}.assoc.txt $PHE $OUT $ENV $BADSNPS >> nohups/$OUT 2>&1
	mkdir -p manhattan
	mkdir -p qqplots
	mv ${OUT}_manhattan.png manhattan/
	mv ${OUT}_qqp.png qqplots/
	mv ${OUT}_below5e-08.assoc.txt manhattan/
    fi
    
else    

    if [ "$GSMGENO" == "1" ]; then
	
	plink --bfile $GENOFILE2 --keep run${YEAR}_v${version}/${PHE1}/tmp${OUT}.fam --make-bed --out run${YEAR}_v${version}/${PHE1}/${OUT}gsm >> nohups/$OUT 2>&1
	cp run${YEAR}_v${version}/$PHE1/pheno${OUT}.fam run${YEAR}_v${version}/$PHE1/${OUT}gsm.fam
    fi
	
    plink --bfile $GENOFILE --keep run${YEAR}_v${version}/${PHE1}/tmp${OUT}.fam --make-bed --out run${YEAR}_v${version}/${PHE1}/${OUT} >> nohups/$OUT 2>&1
    mv run${YEAR}_v${version}/$PHE1/pheno${OUT}.fam run${YEAR}_v${version}/$PHE1/${OUT}.fam
    
    if  [ "$MODEL" == "addtrans" ]; then
	if `echo $BADSNPS | grep -q "\\.gsm$"`; then

	    if [ "$GSMGENO" == "1" ]; then
		prog/gemma.linux -bfile run${YEAR}_v${version}/${PHE1}/${OUT}gsm -gk 2 -o ${OUTsXX} -maf 0.0 -miss 0.8 -snps $BADSNPS >> nohups/$OUT 2>&1
		K=output/${OUTsXX}.sXX.txt
	    else
		prog/gemma.linux -bfile run${YEAR}_v${version}/${PHE1}/${OUT} -gk 2 -o ${OUTsXX} -maf 0.0 -miss 0.8 -snps $BADSNPS >> nohups/$OUT 2>&1
		K=output/${OUTsXX}.sXX.txt
	    fi
	else
	    
	    if [ "$GSMGENO" == "1" ]; then
		prog/gemma.linux -bfile run${YEAR}_v${version}/$PHE1/${OUT}gsm -gk 2 -o ${OUTsXX} -maf 0.05 -miss 0.01 >> nohups/$OUT 2>&1
		K=output/${OUTsXX}.sXX.txt
	    else
		prog/gemma.linux -bfile run${YEAR}_v${version}/$PHE1/${OUT} -gk 2 -o ${OUTsXX} -maf 0.05 -miss 0.01 >> nohups/$OUT 2>&1
		K=output/${OUTsXX}.sXX.txt
	    fi
	fi
    else	
	K=output/${OUTsXX}.sXX.txt
    fi
    
    n=`head -n 1 run${YEAR}_v${version}/$PHE1/${OUT}.fam | wc -w`
    ## because first 5 columns of fam file, are not phenos
    n=$((n-5))

    ## choose model
    if [ "$REC" == "1" ]; then
	if [ "$SNPS" != "NONE" ]; then
	    echo "RECModel"
	    prog/gemma_recessive -bfile run${YEAR}_v${version}/$PHE1/${OUT} -k $K -lmm 3 -o $OUT -c run${YEAR}_v${version}/$PHE1/${OUT}.cov -n `seq 1 $n` -maf 0 -miss 1 -r2 1 -snps $SNPS >> nohups/$OUT 2>&1
	else
	    prog/gemma_recessive -bfile run${YEAR}_v${version}/$PHE1/${OUT} -k $K -lmm 3 -o $OUT -c run${YEAR}_v${version}/$PHE1/${OUT}.cov -n `seq 1 $n` -maf 0 -miss 1 -r2 1 >> nohups/$OUT 2>&1
	fi
    elif [ "$ENV" != "NULL" ]; then
	if [ "$SNPS" != "NONE" ]; then
	    prog/gemma.linux -bfile run${YEAR}_v${version}/$PHE1/${OUT} -k $K -lmm 3 -o $OUT -c run${YEAR}_v${version}/$PHE1/${OUT}.cov -n `seq 1 $n` -maf 0 -miss 1 -r2 1 -snps $SNPS -gxe run${YEAR}_v${version}/${PHE1}/${OUT}.env >> nohups/$OUT 2>&1
	else
	    prog/gemma.linux -bfile run${YEAR}_v${version}/$PHE1/${OUT} -k $K -lmm 3 -o $OUT -c run${YEAR}_v${version}/$PHE1/${OUT}.cov -n `seq 1 $n` -maf 0 -miss 1 -r2 1 -gxe run${YEAR}_v${version}/${PHE1}/${OUT}.env >> nohups/$OUT 2>&1
	fi
    else
	if [ "$SNPS" != "NONE" ]; then
	    prog/gemma.linux -bfile run${YEAR}_v${version}/$PHE1/${OUT} -k $K -lmm 3 -o $OUT -c run${YEAR}_v${version}/$PHE1/${OUT}.cov -n `seq 1 $n` -maf 0 -miss 1 -r2 1 -snps $SNPS >> nohups/$OUT 2>&1
	else	    	    	    
	    prog/gemma.linux -bfile run${YEAR}_v${version}/$PHE1/${OUT} -k $K -lmm 3 -o $OUT -c run${YEAR}_v${version}/$PHE1/${OUT}.cov -n `seq 1 $n` -maf 0 -miss 1 -r2 1 >> nohups/$OUT 2>&1
	fi
    fi



    if [ "$SNPS" == "NONE" ]; then
	rm run${YEAR}_v${version}/${PHE1}/${OUT}.bed
	Rscript scripts/generate_manhattan_pScore.R output/${OUT}.assoc.txt $PHE $OUT $ENV $BADSNPS >> nohups/$OUT 2>&1
	mkdir -p manhattan
	mkdir -p qqplots
	mv ${OUT}_manhattan.png manhattan/
	mv ${OUT}_qqp.png qqplots/
	mv ${OUT}_below5e-08.assoc.txt manhattan/
    fi
fi
