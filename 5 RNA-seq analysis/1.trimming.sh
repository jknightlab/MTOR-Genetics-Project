#!/bin/bash
#SBATCH -A {}.prj 
#SBATCH -p short --cpus-per-task 4      
#SBATCH -J trmming            
#SBATCH -o triming-%a.o         
#SBATCH --array 1-12:1 

echo "------------------------------------------------"
echo `date`: Executing task ${SLURM_ARRAY_TASK_ID} of job ${SLURM_ARRAY_JOB_ID} on `hostname` as user ${USER}
echo "Started at: "`date`
echo "------------------------------------------------"


# Set parameters
## source config.sh
if [[ ! -e "Mapping" ]]; then
mkdir Mapping
fi

# Fastq files 
FASTQDIR="~/raw.RNAseq"
echo "Fastq directory: "$FASTQDIR

FASTQ=$(cat sample.key_12.txt | tail -n+${SLURM_ARRAY_TASK_ID} | head -1 | cut -f1 )
echo $FASTQ
FASTQ1=$FASTQDIR/$FASTQ"_R1.fastq.gz"
FASTQ2=$FASTQDIR/$FASTQ"_R2.fastq.gz"
READS=$FASTQ1" "$FASTQ2

ml use -a /apps/eb/skylake/modules/all
ml use -a /apps/eb/2020b/skylake/modules/all
module load Trim_Galore/0.6.2-GCCcore-8.2.0-Java-11

trim_galore --cores 4 --paired --gzip --output_dir Mapping/AdapterTrimmed/ $READS

echo "Ended at: "`date`

