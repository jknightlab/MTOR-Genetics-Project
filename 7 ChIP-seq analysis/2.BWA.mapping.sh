#!/bin/bash
#SBATCH -A jknight.prj 
#SBATCH -p short --cpus-per-task 8      # Partition (queue) - asking for short queue & 4-cpu job
#SBATCH -J bwa             # Job name
#SBATCH -o bwa-%a.o          # Standard error "slurm-%A_%a.out", "%A" is replaced by the job ID and "%a" with the array index.
#SBATCH --array 1-12:1 

echo "------------------------------------------------"
echo `date`: Executing task ${SLURM_ARRAY_TASK_ID} of job ${SLURM_ARRAY_JOB_ID} on `hostname` as user ${USER}
echo "Started at: "`date`
echo "------------------------------------------------"

# load samtools
#module load SAMtools/1.9-foss-2018b
module load SAMtools/1.18-GCC-12.3.0

# output folder
if [[ ! -e Mapping/mapped.bwa-MEM ]]; then
mkdir Mapping/mapped.bwa-MEM
fi

# Get fastq file names
FASTQ=$(cat sample.key.txt | tail -n+${SLURM_ARRAY_TASK_ID} | head -1 | cut -f1 )
FASTQ1=$FASTQ"_1_val_1.fq.gz"
FASTQ1=Mapping/AdapterTrimmed/$FASTQ1
FASTQ2=$FASTQ"_2_val_2.fq.gz"
FASTQ2=Mapping/AdapterTrimmed/$FASTQ2
echo $FASTQ1  $FASTQ2

REF=/well/jknight/users/kwz374/refs/REF.UCSC/hg38.fa
## /apps/well/bwa/0.7.12-gcc4.9.1/bwa index $REF

echo "##### ***Running bwa-MEM on SE SAMPLE: $FASTQ1"
/apps/well/bwa/0.7.12-gcc4.9.1/bwa mem -t 8 \
-M $REF $FASTQ1 | samtools sort > Mapping/mapped.bwa-MEM/${FASTQ}.bam 

#echo "##### ***Running bwa-MEM on PE SAMPLE: $FASTQ1 $FASTQ2"
#/apps/well/bwa/0.7.12-gcc4.9.1/bwa mem -t 12 \
#-M $REF $FASTQ1 $FASTQ2 | samtools sort -@ 12 > Mapping/mapped.bwa-MEM/${FASTQ}.bam 

samtools index -@ 12 Mapping/mapped.bwa-MEM/${FASTQ}.bam

########*****************************************************************#######
echo "###Running MAPQ10 on $FASTQ"
samtools view -@ 12 -q 10 Mapping/mapped.bwa-MEM/${FASTQ}.bam -b > Mapping/mapped.bwa-MEM/${FASTQ}.MAPQ10.bam

echo "###Running picard on $FASTQ"
/apps/well/java/jdk1.8.0_latest/bin/java -Xmx8g -jar \
/apps/well/picard-tools/2.21.1/picard.jar MarkDuplicates \
INPUT= Mapping/mapped.bwa-MEM/${FASTQ}.MAPQ10.bam \
OUTPUT= Mapping/mapped.bwa-MEM/${FASTQ}.MAPQ10.dedup.bam \
METRICS_FILE= Mapping/mapped.bwa-MEM/dedup_bam_metrics_${FASTQ}.txt \
REMOVE_DUPLICATES=true

samtools index -@ 12 Mapping/mapped.bwa-MEM/${FASTQ}.MAPQ10.dedup.bam

echo "Ended at: "`date`
