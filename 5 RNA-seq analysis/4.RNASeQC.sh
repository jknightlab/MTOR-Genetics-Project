#!/bin/bash
#SBATCH -A {}.prj 
#SBATCH -p short --cpus-per-task 1
#SBATCH -J rnaseqc            
#SBATCH -o rnaseqc.o           
#SBATCH -e rnaseqc.e    
#SBATCH --array 1-12:1 

echo "Started at: "`date`

##### results folder
if [[ ! -e RNASeQC ]]; then
mkdir -p RNASeQC
fi

# Get sample ID for this task
SAMPLE_NAME=$(cat sample.key_12.txt | tail -n+${SLURM_ARRAY_TASK_ID} | head -1 | cut -f1 )
echo "SAMPLE_NAME: $SAMPLE_NAME"

GTF=gencode.v31.collapsed.gene.gtf

~/rnaseqc.v2.3.6.linux $GTF \
--mapping-quality=60 \
Mapping/mapped.hisat2/${SAMPLE_NAME}.hisat2.bam \
RNASeQC/

echo "Ended at: "`date`

