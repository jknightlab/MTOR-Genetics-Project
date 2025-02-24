#!/bin/bash

# Define the study and sample group combinations
studies=("BLUEPRINT" "Schmiedel_2018" "Schmiedel_2018" "FUSION" "TwinsUK")
sample_groups=("neutrophil" "CD4_T-cell_anti-CD3-CD28" "CD8_T-cell_anti-CD3-CD28" "adipose_naive" "fat")

# Loop through each study and sample group combination
for i in "${!studies[@]}"; do
    study=${studies[$i]}
    sample_group=${sample_groups[$i]}

    echo "Running analysis for $study and $sample_group..."
    
    # Run the R script with the current study and sample group
    Rscript ~/MTOR.project/smr.analysis/1.smr.file.prep.R "$study" "$sample_group"
done