#!/bin/bash
#SBATCH -A {}.prj 
#SBATCH -p short --cpus-per-task 8
#SBATCH -J featurecount            
#SBATCH -o featurecount-%a.o

echo "Started at: "`date`

##### results folder
if [[ ! -e featureCounts ]]; then
mkdir -p featureCounts
fi

GTF=/gtfs/gencode.v31.annotation.gtf.gz

/apps/htseq/subread-1.6.4-Linux-x86_64/bin/featureCounts -T 12 \
-p -s 2 -a $GTF -g gene_name \
--extraAttributes gene_id \
-o featureCounts/featureCounts_Tcell.day0.day4.txt \
Mapping/mapped.hisat2/*bam

echo "Ended at: "`date`

