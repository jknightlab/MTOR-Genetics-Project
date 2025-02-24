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

## SE Fastq files 
FASTQDIR="/well/jknight/users/kwz374/MTOR-genetic-project/CD4+.T.cell.5mC-and-5hmC.medIP-Seq/fastq.file"

echo "Fastq directory: "$FASTQDIR

FASTQ=$(cat sample.key.txt | tail -n+${SLURM_ARRAY_TASK_ID} | head -1 | cut -f1 )
echo $FASTQ
FASTQ=$FASTQDIR/$FASTQ".fastq.gz"

module load Trim_Galore/0.6.2-GCCcore-8.2.0-Java-11

trim_galore --cores 16 --gzip --output_dir Mapping/AdapterTrimmed/ $FASTQ

echo "Ended at: "`date`

