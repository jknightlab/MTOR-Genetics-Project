#!/bin/bash

# Chr    ProbeID GeneticDistance ProbeBp Gene    Orientation PathOfEsd
# Columns are chromosome, probe ID(can be the ID of an exon or a transcript for RNA-seq data),
# genetic distance (can be any arbitary value), physical position, gene ID, gene orientation,
# and PathOfEsd
# (this information will only be used for graphic presentation, please see [SMR plot]). from https://yanglab.westlake.edu.cn/software/smr/#SMR&HEIDIanalysis

INPUT_DIR="esd.fat"
BASE_NAME=$(echo "$INPUT_DIR" | sed 's/^[^.]*\.//')
OUTPUT_DIR="./"        
gffdb="./gencode.v41lift37.basic.annotation.gff3"      
feature_ids="${INPUT_DIR}/feature_ids.txt"
gwas_file="./T2D_cohort2024_subset_gwas_sum_ma"
gene_list="./glist_hg19_refseq.txt"
bfile="./1000genomes_EUR_chr1_dedup"

### 1. make ESD files
echo -e "*** make the file for smr --eqtl-flist"
# Generate esd_file_list.txt using awk
awk -v input_dir="$INPUT_DIR" '
BEGIN { 
    print "Chr\tProbeID\tGeneticDistance\tProbeBp\tGene\tOrientation\tPathOfEsd"; 
}
NR == FNR {
    split($9, ATTR, ";");
    split(ATTR[6], TMPLIST, "[=.]"); 
    GENEID = TMPLIST[2];
    CHROM = $1; 
    sub(/^chr/, "", CHROM);
    PROBEID_DICT[GENEID] = CHROM"\t"GENEID"\t0\t"$5"\t"GENEID"\t"$7;
}
NR != FNR {
    if (PROBEID_DICT[$1]) 
        print PROBEID_DICT[$1]"\t"input_dir"/"$1".esd";
    else 
        next;
}
' "$gffdb" "$feature_ids" | sort -k1,1g -k4,4g > "$INPUT_DIR/esd_file_list.txt"

# Check if esd_file_list.txt is empty
if [[ ! -s "$INPUT_DIR/esd_file_list.txt" ]]; then
    echo "Error: esd_file_list.txt is empty! Please check input files and processing steps." >&2
    exit 1
else
    echo "esd_file_list.txt is successfully generated and contains data."
fi

### 2. Make a BESD file from eQTL summary data in ESD format
echo -e "*** Generate BESD database from eQTL summary ..."
./smr --make-besd --eqtl-flist ${INPUT_DIR}/esd_file_list.txt --out ${OUTPUT_DIR}${BASE_NAME}_eqtl_summary

#module load PLINK/2.00a2.3_x86_64
#plink2 --bfile ~/1kg_EUR/1000genomes_EUR_chr1 --rm-dup force-first --make-bed --out 1000genomes_EUR_chr1_dedup

### 3. run SMR and HEIDI test
./smr --bfile $bfile \
--gwas-summary $gwas_file \
--beqtl-summary ${OUTPUT_DIR}${BASE_NAME}_eqtl_summary \
--out ${OUTPUT_DIR}${BASE_NAME} \
--thread-num 8 


### 4. get the results for plotting
./smr --bfile $bfile \
--gwas-summary $gwas_file \
--beqtl-summary ${OUTPUT_DIR}${BASE_NAME}_eqtl_summary \
--out myplot${BASE_NAME} \
--plot --probe MTOR \
--probe-wind 1000 --gene-list $gene_list

###
### 4. SMR locus and effect plots
module load R/4.3.2-gfbf-2023a 

Rscript - <<EOF
source("./plot/plot_SMR.r") # with fixed line 
plot.resul <- ReadSMRData("./plot/myplot${BASE_NAME}.MTOR.txt")
pdf("SNP_effect.pdf", width = 4, height = 4.5)
SMREffectPlot(data=plot.resul)
dev.off()

pdf("locus_plot.pdf", width = 7, height = 7)
SMRLocusPlot(data=plot.resul,smr_thresh=8.4e-6, heidi_thresh=0.05, plotWindow=1000, max_anno_probe=16)
dev.off()
EOF

# merge the results
awk 'FNR==1 { if (NR==1) { print $0, "eQTL.dataset"; } next; } {print $0, FILENAME}' OFS='\t' *.smr | sed 's/\.smr//g' > merged_output.smr.txt




