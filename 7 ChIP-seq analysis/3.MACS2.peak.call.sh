#!/bin/bash
#SBATCH -A jknight.prj 
#SBATCH -p short --cpus-per-task 4
#SBATCH -J macs5           
#SBATCH -o macs5-%a.o          # Standard error "slurm-%A_%a.out", "%A" is replaced by the job ID and "%a" with the array index.

echo "Started at: "`date`

module load MACS2/2.2.6-foss-2018b-Python-3.6.6

FILE=/well/jknight/users/kwz374/MTOR-genetic-project/CD4+.T.cell.5mC-and-5hmC.medIP-Seq/Mapping/mapped.bwa-MEM/

#--treatment ${FILE}SRR2932482_GSM1936640_Activated_hmC_medIP-Seq_7_12_Homo_sapiens_MeDIP-Seq.MAPQ10.dedup.bam \
#--control ${FILE}SRR2932487_GSM1936645_Activated_input_7_12_Homo_sapiens_ChIP-Seq.MAPQ10.dedup.bam \

#--treatment ${FILE}SRR2932480_GSM1936638_Naive_hmC_medIP-Seq_7_12_Homo_sapiens_MeDIP-Seq.MAPQ10.dedup.bam \
#--control ${FILE}SRR2932484_GSM1936642_Naive_input_7_12_Homo_sapiens_ChIP-Seq.MAPQ10.dedup.bam \

#--treatment ${FILE}SRR2932481_GSM1936639_Activated_mC_medIP-Seq_10_12_Homo_sapiens_MeDIP-Seq.MAPQ10.dedup.bam \
#--control ${FILE}SRR2932488_GSM1936646_Activated_input_10_12_Homo_sapiens_ChIP-Seq.MAPQ10.dedup.bam \

macs2 callpeak \
--treatment ${FILE}SRR2932479_GSM1936637_Naive_mC_medIP-Seq_10_12_Homo_sapiens_MeDIP-Seq.MAPQ10.dedup.bam \
--control ${FILE}SRR2932485_GSM1936643_Naive_input_10_12_Homo_sapiens_ChIP-Seq.MAPQ10.dedup.bam \
--name mc_native.T \
--format BAM \
--keep-dup all \
--gsize hs \
--qvalue 0.01 \
--SPMR \
-B \
--outdir Macs2

## bigwig
REF=/well/jknight/users/kwz374/refs/REF.UCSC/hg38.fa
##
sort -k1,1 -k2,2n Macs2/mc_native.T_treat_pileup.bdg  > Macs2/mc_native.T_treat_pileup_sorted.bdg
##
~/bedGraphToBigWig \
Macs2/mc_native.T_treat_pileup_sorted.bdg \
${REF}.fai Macs2/mc_native.T_treat_pileup_sorted.bw

#### normalised bigwig
echo "##### ***Rmacs2 bdgcmp for fold enrichment"
macs2 bdgcmp \
-t Macs2/mc_native.T_treat_pileup.bdg \
-c Macs2/mc_native.T_control_lambda.bdg \
-o Macs2/mc_native.T_FE.bdg \
-m FE

sort -k1,1 -k2,2n Macs2/mc_native.T_FE.bdg  > Macs2/mc_native.T_FE_sorted.bdg
#
~/bedGraphToBigWig \
Macs2/mc_native.T_FE_sorted.bdg \
${REF}.fai \
Macs2/mc_native.T_FE_sorted_normalized.bw


echo "Ended at: "`date`

