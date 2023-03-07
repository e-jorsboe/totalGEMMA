#!/bin/bash
# ----------------QSUB Parameters----------------- #

#PBS -W group_list=emijor
#PBS -A cu_10048

#PBS -m n
#PBS -N part2GemmaArray
#PBS -t 1-856%140
#PBS -l nodes=1:ppn=1
#PBS -l walltime=60:00:00
#PBS -l mem=4gb

#PBS -o logs/
#PBS -e logs/

module load tools
module load R/3.1.2

# ----------------Your Commands------------------- #

sleep $PBS_ARRAYID

cd $PBS_O_WORKDIR

VERSION=3 
PHENOFILE=/home/projects/cu_10048/data/Greenland_Pheno_17nov2015.csv 
FAMFILE=/home/projects/cu_10048/data/datplus_QCed_newIDs_allcohorts_allSNPsJuli2015_geno030_nodups.fam 
ANALYSIS_DESCRIPTION=/home/projects/cu_10048/people/emijor/part2rundescription_impu2016_v3.txt
ABOVE=04
VERSION=3

traitinfo=`sed -n "${PBS_ARRAYID}p" $ANALYSIS_DESCRIPTION`

PHE=`echo $traitinfo | cut -f1 -d" "`
TRANS=`echo $traitinfo | cut -f2 -d" "`
COV=`echo $traitinfo | cut -f3 -d" "`
MODEL=`echo $traitinfo | cut -f4 -d" "`
COHORT=`echo $traitinfo | cut -f5 -d" "`

PBS_ARRAYID2=$((PBS_ARRAYID+856))

echo $PBS_ARRAYID2 $PHE $TRANS $COV $MODEL $COHORT > /home/projects/cu_10048/data/runLogs/run${PBS_ARRAYID2}.txt

OUT=impu2016_v${VERSION}_${PHE}


if [ "$MODEL" -eq "1" ]; then

Rscript /home/projects/cu_10048/people/emijor/prep_phenos_bimbam_IHIT_B99V3.R $PHE $COV 0 $ABOVE $PHENOFILE $FAMFILE ${OUT}_untrans $COHORT > /home/projects/cu_10048/data/nohups/${OUT}_adduntrans_${COHORT}.nohup

COLS=`cat /home/projects/cu_10048/data/generated/${PHE}/newAbove${ABOVE}/${OUT}_untrans_${COHORT}_cols.txt`


elif [ "$MODEL" -eq "2" ]; then

Rscript /home/projects/cu_10048/people/emijor/prep_phenos_bimbam_IHIT_B99V3.R $PHE $COV $TRANS $ABOVE $PHENOFILE $FAMFILE ${OUT}_trans $COHORT > /home/projects/cu_10048/data/nohups/${OUT}_addtrans_${COHORT}.nohup

COLS=`cat /home/projects/cu_10048/data/generated/${PHE}/newAbove${ABOVE}/${OUT}_trans_${COHORT}_cols.txt`


elif [ "$MODEL" -eq "3" ]; then

Rscript /home/projects/cu_10048/people/emijor/prep_phenos_bimbam_IHIT_B99V3.R $PHE $COV 0 $ABOVE $PHENOFILE $FAMFILE ${OUT}_untrans $COHORT > /home/projects/cu_10048/data/nohups/${OUT}_recuntrans_${COHORT}.nohup

COLS=`cat /home/projects/cu_10048/data/generated/${PHE}/newAbove${ABOVE}/${OUT}_untrans_${COHORT}_cols.txt`


elif [ "$MODEL" -eq "4" ]; then

Rscript /home/projects/cu_10048/people/emijor/prep_phenos_bimbam_IHIT_B99V3.R $PHE $COV $TRANS $ABOVE $PHENOFILE $FAMFILE ${OUT}_trans $COHORT > /home/projects/cu_10048/data/nohups/${OUT}_rectrans_${COHORT}.nohup

COLS=`cat /home/projects/cu_10048/data/generated/${PHE}/newAbove${ABOVE}/${OUT}_trans_${COHORT}_cols.txt`

fi



if [ "$MODEL" -eq "1" ]; then

cut -f 1,2,3,$COLS -d "," /home/projects/cu_10048/data/add46xtraR2above04.phased.impute2.mgf > /home/projects/cu_10048/data/generated/${PHE}/newAbove${ABOVE}/${OUT}_adduntrans_${COHORT}.mgf


elif [ "$MODEL" -eq "2" ]; then

cut -f 1,2,3,$COLS -d "," /home/projects/cu_10048/data/add46xtraR2above04.phased.impute2.mgf > /home/projects/cu_10048/data/generated/${PHE}/newAbove${ABOVE}/${OUT}_addtrans_${COHORT}.mgf


elif [ "$MODEL" -eq "3" ]; then

cut -f 1,2,3,$COLS -d "," /home/projects/cu_10048/data/rec46xtraR2above04.phased.impute2.mgf > /home/projects/cu_10048/data/generated/${PHE}/newAbove${ABOVE}/${OUT}_recuntrans_${COHORT}.mgf


elif [ "$MODEL" -eq "4" ]; then

cut -f 1,2,3,$COLS -d "," /home/projects/cu_10048/data/rec46xtraR2above04.phased.impute2.mgf > /home/projects/cu_10048/data/generated/${PHE}/newAbove${ABOVE}/${OUT}_rectrans_${COHORT}.mgf


fi


K=/home/projects/cu_10048/data/output/${OUT}_${COHORT}.sXX.txt

if [ "$MODEL" -eq "1" ]; then

/home/projects/cu_10048/apps/gemma -g /home/projects/cu_10048/data/generated/${PHE}/newAbove${ABOVE}/${OUT}_adduntrans_${COHORT}.mgf -a /home/projects/cu_10048/data/all46xtraR2above04.phased.impute2.ann -p /home/projects/cu_10048/data/generated/${PHE}/newAbove${ABOVE}/${OUT}_untrans_${COHORT}_pheno.txt -k $K -lmm 4 -o ${OUT}_adduntrans_${COHORT} -outdir /home/projects/cu_10048/data/output -c /home/projects/cu_10048/data/generated/${PHE}/newAbove${ABOVE}/${OUT}_untrans_${COHORT}.cov -maf 0 -miss 0.8 >> /home/projects/cu_10048/data/nohups/${OUT}_adduntrans_${COHORT}.nohup


elif [ "$MODEL" -eq "2" ]; then

/home/projects/cu_10048/apps/gemma -g /home/projects/cu_10048/data/generated/${PHE}/newAbove${ABOVE}/${OUT}_addtrans_${COHORT}.mgf -a /home/projects/cu_10048/data/all46xtraR2above04.phased.impute2.ann -p /home/projects/cu_10048/data/generated/${PHE}/newAbove${ABOVE}/${OUT}_trans_${COHORT}_pheno.txt -k $K -lmm 4 -o ${OUT}_addtrans_${COHORT} -outdir /home/projects/cu_10048/data/output -c /home/projects/cu_10048/data/generated/${PHE}/newAbove${ABOVE}/${OUT}_trans_${COHORT}.cov -maf 0 -miss 0.8 >> /home/projects/cu_10048/data/nohups/${OUT}_addtrans_${COHORT}.nohup


elif [ "$MODEL" -eq "3" ]; then

/home/projects/cu_10048/apps/gemma -g /home/projects/cu_10048/data/generated/${PHE}/newAbove${ABOVE}/${OUT}_recuntrans_${COHORT}.mgf -a /home/projects/cu_10048/data/all46xtraR2above04.phased.impute2.ann -p /home/projects/cu_10048/data/generated/${PHE}/newAbove${ABOVE}/${OUT}_untrans_${COHORT}_pheno.txt -k $K -lmm 4 -o ${OUT}_recuntrans_${COHORT} -outdir /home/projects/cu_10048/data/output -c /home/projects/cu_10048/data/generated/${PHE}/newAbove${ABOVE}/${OUT}_untrans_${COHORT}.cov -maf 0 -miss 0.8 >> /home/projects/cu_10048/data/nohups/${OUT}_recuntrans_${COHORT}.nohup


elif [ "$MODEL" -eq "4" ]; then

/home/projects/cu_10048/apps/gemma -g /home/projects/cu_10048/data/generated/${PHE}/newAbove${ABOVE}/${OUT}_rectrans_${COHORT}.mgf -a /home/projects/cu_10048/data/all46xtraR2above04.phased.impute2.ann -p /home/projects/cu_10048/data/generated/${PHE}/newAbove${ABOVE}/${OUT}_trans_${COHORT}_pheno.txt -k $K -lmm 4 -o ${OUT}_rectrans_${COHORT} -outdir /home/projects/cu_10048/data/output -c /home/projects/cu_10048/data/generated/${PHE}/newAbove${ABOVE}/${OUT}_trans_${COHORT}.cov -maf 0 -miss 0.8 >> /home/projects/cu_10048/data/nohups/${OUT}_rectrans_${COHORT}.nohup


fi

if [ "$MODEL" -eq "1" ]; then

rm  /home/projects/cu_10048/data/generated/${PHE}/newAbove${ABOVE}/${OUT}_adduntrans_${COHORT}.mgf

elif [ "$MODEL" -eq "2" ]; then

rm  /home/projects/cu_10048/data/generated/${PHE}/newAbove${ABOVE}/${OUT}_addtrans_${COHORT}.mgf

elif [ "$MODEL" -eq "3" ]; then

rm  /home/projects/cu_10048/data/generated/${PHE}/newAbove${ABOVE}/${OUT}_recuntrans_${COHORT}.mgf

elif [ "$MODEL" -eq "4" ]; then

rm  /home/projects/cu_10048/data/generated/${PHE}/newAbove${ABOVE}/${OUT}_rectrans_${COHORT}.mgf

fi

echo finished $PBS_ARRAYID2 $PHE $TRANS $COV $MODEL $COHORT >> /home/projects/cu_10048/data/runLogs/run${PBS_ARRAYID2}.txt
