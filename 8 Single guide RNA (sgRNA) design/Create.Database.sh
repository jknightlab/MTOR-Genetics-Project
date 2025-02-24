#!/bin/bash
#$ -cwd
#$ -V
#$ -N test 
#$ -S /bin/bash 
#$ -q short.qc
#$ -j y
#$ -pe shmem 2 
#$ -b n 


mkdir tmp

/apps/well/java/jdk1.8.0_latest/bin/java -Xmx6000m \
-jar FlashFry-assembly-1.9.3.jar \
 index \
 --tmpLocation ./tmp \
 --database hg38_Cas9NGG_database \
 --reference hg38.fa.gz \
 --enzyme spcas9ngg
