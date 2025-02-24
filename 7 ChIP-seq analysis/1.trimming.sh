#!/bin/bash
#SBATCH -A jknight.prj 
#SBATCH -p short --cpus-per-task 4     # Partition (queue) - asking for short queue & 4-cpu job
#SBATCH -J trmming             # Job name
#SBATCH -o triming-%a.o          # Standard error "slurm-%A_%a.out", "%A" is replaced by the job ID and "%a" with the array index.
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
##FASTQDIR="/well/jknight/users/kwz374/MTOR-genetic-project/CD4+.T.cell.5mC-and-5hmC.medIP-Seq/fastq.file"
##FASTQDIR="/well/jknight/users/kwz374/MTOR-genetic-project/CD4+.T.cell.H3K9me3-and-KAP1-ChIP-seq/fastq.file"

#echo "Fastq directory: "$FASTQDIR

#FASTQ=$(cat sample.key.txt | tail -n+${SLURM_ARRAY_TASK_ID} | head -1 | cut -f1 )
#echo $FASTQ
#FASTQ=$FASTQDIR/$FASTQ".fastq.gz"

#module load Trim_Galore/0.6.2-GCCcore-8.2.0-Java-11

#trim_galore --cores 16 --gzip --output_dir Mapping/AdapterTrimmed/ $FASTQ

# PE Fastq files 
FASTQDIR="/well/jknight/users/kwz374/MTOR-genetic-project/CD4+.T.cell.5hmC.medIP-Seq_rep2_PRJEB13837/fastq.file"
echo "Fastq directory: "$FASTQDIR

FASTQ=$(cat sample.key.txt | tail -n+${SLURM_ARRAY_TASK_ID} | head -1 | cut -f1 )
echo $FASTQ
FASTQ1=${FASTQDIR}/${FASTQ}_1.fastq.gz
FASTQ2=${FASTQDIR}/${FASTQ}_2.fastq.gz
READS=$FASTQ1" "$FASTQ2

module load Trim_Galore/0.6.10-GCCcore-12.3.0

trim_galore --cores 16 --paired --gzip --output_dir Mapping/AdapterTrimmed/ $READS

echo "Ended at: "`date`

