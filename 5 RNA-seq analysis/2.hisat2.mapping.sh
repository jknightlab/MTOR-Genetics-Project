#!/bin/bash
#SBATCH -A jknight.prj 
#SBATCH -p short --cpus-per-task 4      # Partition (queue) - asking for short queue & 4-cpu job
#SBATCH -J hisat2             # Job name
#SBATCH -o hisat2-%a.o          # Standard error "slurm-%A_%a.out", "%A" is replaced by the job ID and "%a" with the array index.
#SBATCH --array 1-12:1 

echo "------------------------------------------------"
echo `date`: Executing task ${SLURM_ARRAY_TASK_ID} of job ${SLURM_ARRAY_JOB_ID} on `hostname` as user ${USER}
echo "Started at: "`date`
echo "------------------------------------------------"

# output folder
if [[ ! -e Mapping/mapped.hisat2 ]]; then
mkdir Mapping/mapped.hisat2
fi


# Get fastq file names
FASTQ=$(cat sample.key_12.txt | tail -n+${SLURM_ARRAY_TASK_ID} | head -1 | cut -f1 )
FASTQ1="$(echo "${FASTQ}" | sed -e "s|,|_1_val_1.fq.gz,Mapping/AdapterTrimmed/\/|g")"
FASTQ2="$(echo "${FASTQ}" | sed -e "s|,|_2_val_2.fq.gz,Mapping/AdapterTrimmed/\/|g")"
FASTQ1=$FASTQ1"_R1_val_1.fq.gz"
FASTQ2=$FASTQ2"_R2_val_2.fq.gz"
FASTQ1=Mapping/AdapterTrimmed/$FASTQ1
FASTQ2=Mapping/AdapterTrimmed/$FASTQ2
echo $FASTQ1
echo $FASTQ2

# load samtools
ml use -a /apps/eb/skylake/modules/all
ml use -a /apps/eb/2020b/skylake/modules/all
module load SAMtools/1.9-foss-2018b

# The HISAT2 index name.
IDX=/well/jknight/users/kwz374/refs/grch38/genome

echo "Running Hisat2: $FASTQ1 $FASTQ2"

/apps/well/hisat2/2.1.0/hisat2 --threads 12  -x $IDX -1 $FASTQ1 -2 $FASTQ2 | samtools sort -@ 12  > Mapping/mapped.hisat2/${FASTQ}.hisat2.bam 

samtools index -@ 12 Mapping/mapped.hisat2/${FASTQ}.hisat2.bam



##echo "Ended at: "`date`


