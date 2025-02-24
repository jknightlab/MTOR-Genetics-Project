#!/bin/bash
#SBATCH -A jknight.prj 
#SBATCH -p short --cpus-per-task 4      
#SBATCH -J bigwig             
#SBATCH -o bigwig-%a.o          
#SBATCH --array 1-12:1 

module load SAMtools/1.9-foss-2018b
module load deepTools/3.3.1-foss-2018b-Python-3.6.6

# output folder
if [[ ! -e Bigwig.file ]]; then
mkdir Bigwig.file
fi

# Get fastq file names
FASTQ=$(cat sample.key_12.txt | tail -n+${SLURM_ARRAY_TASK_ID} | head -1 | cut -f1 )

###
bamCoverage -b Mapping/mapped.hisat2/${FASTQ}_chr.uniq.mapped.bam \
--normalizeUsing RPKM \
--numberOfProcessors 20 \
--binSize 20 \
--outFileFormat bigwig \
-o Bigwig.file/${FASTQ}_RPKM.bw
