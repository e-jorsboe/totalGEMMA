#!/bin/bash

## myregionfile
## chr start end, 
## 1 (chr)
## 1 (start)     
## 10000000 (end)

C=$1

CHR=$(head -2 /home/emil/OMNI5/IMPUTE2/regions/myregionfile${C} | tail -1)
BGN=$(head -3 /home/emil/OMNI5/IMPUTE2/regions/myregionfile${C} | tail -1)
END=$(head -4 /home/emil/OMNI5/IMPUTE2/regions/myregionfile${C} | tail -1)

for i in `seq $BGN 5000000 $END`
do

 
 j=$((i+5000000-1))
 s_j=$((j/1000000))
 s_i=$((i/1000000))
 echo $s_i
 echo $s_j


 mkdir -p /home/emil/OMNI5/megaChipQC/forImpu/imputation/AllMenDiscMega/chr${CHR}

 /home/emil/software/impute_v2.3.2_x86_64_static/impute2 -merge_ref_panels -use_prephased_g -m /home/emil/OMNI5/genetic_map_b37/WGSdata/filledOut2_All_map_chr${CHR}_b37.txt -h /home/emil/GL_WGS/1000gRefPanel/hg19/refPanelIntersectingMarkersV3MenDisc/gbrCEUchbCDXchsRefChr${CHR}overlapMaf001.hap.gz /home/emil/GL_WGS/data/QCpassedAll_hg19picard/byChrMenDisc/phased/All_hg19picard_chr${CHR}_menDiscV2.hap -l /home/emil/GL_WGS/1000gRefPanel/hg19/refPanelIntersectingMarkersV3MenDisc/gbrCEUchbCDXchsRefChr${CHR}overlapMaf001.legend.gz /home/emil/GL_WGS/data/QCpassedAll_hg19picard/byChrMenDisc/phased/All_hg19picard_chr${CHR}_menDiscV2.leg -known_haps_g /home/emil/OMNI5/megaChipQC/forImpu/byChrMenDisc/phased/megaV2_plusCPT1Asnp_indisNotInWGS_chr${CHR}_mac5_MenDisc.haps -strand_g /home/emil/OMNI5/megaChipQC/forImpu/byChrMenDisc/phased/megaV2_plusCPT1Asnp_indisNotInWGS_chr${CHR}_mac5_MenDisc.strand -int ${i} ${j} -Ne 1500 -o /home/emil/OMNI5/megaChipQC/forImpu/imputation/AllMenDiscMega/chr${CHR}/${s_i}_${s_j}_chr${CHR}.AllMenDiscMega.impute2


done
