# MTOR-Genetics-Project
#### _Zhang, P. et al., Context-specific regulatory genetic variation in MTOR dampens Neutrophil-T cell crosstalk in sepsis, modulating disease. (2025)_

<div align="center">
  <img src="output/model.png" alt="Screenshot" width="80%" />
</div>

Sepsis is a heterogeneous clinical syndrome with a high mortality rate and personalised stratification strategies are proposed as essential to successful targeted therapeutics. Here, we characterise genetic variation that modulates MTOR, a critical regulator of metabolism and immune responses in sepsis. The effects are highly context specific, involving a regulatory element that affects MTOR expression in activated T cells with opposite direction of effect in neutrophils. The lead variant significantly interacts with the known sepsis prognostic marker neutrophil-to-lymphocyte ratio, showing activity specific to sepsis endotype and a pleiotropic effect on type 2 diabetes (T2D) risk. Using ex vivo models, we demonstrate that activated T cells promote immunosuppressive sepsis neutrophils through released cytokines, a process dampened by hypoxia and the mTOR inhibitor rapamycin. The G-allele of the variant, associated with decreased risk of T2D, is associated with reduced mTOR signaling in T cells and improved survival in sepsis patients due to pneumonia. We define a novel epigenetic mechanism that fine-tunes MTOR transcription and T cell activity via the variant-containing regulatory element, which exhibits an allelic effect upon vitamin C treatment. Our findings reveal how common genetic variation can interact with disease state/endotype to modulate immune cell-cell communication, providing a patient stratification strategy to inform more effective treatment of sepsis and suggesting putative mechanisms underlying the variable efficacy of vitamin C therapy in sepsis.

## Source of Data:
- Genomic Advances in Sepsis (GAinS) [whole blood gene expression](https://ega-archive.org/datasets/EGAD00001008730) and [genotyping](https://ega-archive.org/datasets/EGAD00001015369) data.
- ATAC-seq data from human primary immune cells including [macrophages](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE172116) | [monocytes](https://zenodo.org/record/8158923) | [neutrophils](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE150018) | [NK and dendritic cells](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE118189) | [CAR T cells](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE168882).
- RNA-seq data in primary CD4+ or CD8+ T cells treated with anti-CD3/CD28 Dynabeads and mTOR inhibitor [Rapamycin](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE129829).
- MeDIP-Seq for [5hmC](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE74850) in CD4+ T cells.
- Histone modification ChIP-seq and CTCF ChIP-seq/ChIA-PET from the [ENCODE](https://www.encodeproject.org/) project.
- Sepsis whole blood [scRNA-seq](https://zenodo.org/records/7924238) data.
- Genotype data from the [UK Biobank](https://www.ukbiobank.ac.uk/) for individuals with confirmed bacterial pneumonia.
- Genotype data from the [1000 Genomes Project](http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/). 
- [The eQTL Catalogue](https://www.ebi.ac.uk/eqtl/) - Expression Quantitative trait loci (eQTL) recomputed from public datasets derived from 75 tissues/cell types and 14 treatments.
- Type 2 diabetes GWAS summary statistics - [dataset1](https://www.diagram-consortium.org/downloads.html) and [dataset2](https://ftp.ncbi.nlm.nih.gov/dbgap/studies/phs001672/analyses/)

## Table of Contents
- [1 Identification of context specific eQTLs](#1-identification-of-context-specific-eqtls)
  - 1.1 Interaction analysis
  - 1.2 Pairwise linkage disequilibrium (R²) for the MTOR locus
  - 1.3 Colocalisation analysis
- [2 Survival analysis](#2-survival-analysis)
- [3 Summary data-based Mendelian randomisation](#3-summary-data-based-mendelian-randomisation)
  - 3.1 Input file preparation
  - 3.2 SMR/HEIDI analysis
- [4 Pairwise fixation index (Fst)](#4-pairwise-fixation-index-fst)
- [5 RNA-seq analysis](#5-rna-seq-analysis)
  - 5.1 RNA-seq analysis pipeline
  - 5.2 PCA & Differential expression analysis
  - 5.3 Cell type deconvolution
  - 5.4 UMAP projection of sepsis whole blood scRNA-seq data 
- [6 ATAC-seq analysis](#6-atac-seq-analysis)
  - 6.1 ATAC-seq analysis pipeline
  - 6.2 PCA & Differential chromatin accessibility analysis
- [7 ChIP-seq analysis](#7-chip-seq-analysis)
  - 7.1 ChIP-seq analysis pipeline
- [8 Single guide RNA (sgRNA) design](#8-single-guide-rna-sgrna-design)

## 1 Identification of context specific eQTLs
#### 1.1 Interaction analysis
eQTL interaction was determined using a linear mixed model as implemented in R package lmerTest: p ~ g + i + g:i + (1|donor), where p is the gene expression corrected for population structure and PEER factors, g is the genotype vector, i is an interaction term for SRS endotypes or log2 transformed NLR, and (1|donor) is a random effect accounting for variability between donors.
```bash
# Interaction between MTOR.SNP*SRS or NLR
Rscript ./1.LMM-model-interaction.R <SRS_data_file> <cell_counts_file>
# Inverted eQTLs between Neutrophils and T cells
Rscript ./2.inverted.SRS-interacting.eQTLs.R
```
#### 1.2 Pairwise linkage disequilibrium (R²) for the MTOR locus
R² was calculated using PLINK with GAinS genotype data (n=1,168) and visualised with LDheatmap. 
```
bash ./3.r2.ld.map.sh
```
#### 1.3 Colocalisation analysis
Colocalisation analysis was performed in a ±200-bp window surrounding the MTOR lead eQTL, using R package coloc with default settings. eQTL summary statistics from healthy tissue and cell types were retrieved from eQTL Caltalog. 
```bash
Rscript 4.1.coloc.R && Rscript 4.2.plot.resul.R
```

## 2 Survival analysis
To assess the association between genetic variants and 28-day mortality, we use Cox proportional-hazards model and logistic regression adjusting for age, sex, and the first seven genotype principal components.
```bash
Rscript ./cox.PH_logistic.regression.R
```

## 3 Summary data-based Mendelian randomisation
We performed the SMR analysis using eQTLs as instrumental variables to identify genes whose expression is associated with T2D risk due to pleiotropy and/or causality. Genes were included in the analysis if they had at least one cis-eQTL (P < 5e⁻⁸) within a 2 Mbp window around GWAS loci, following the default settings of the [SMR tool (v1.3.1)](https://yanglab.westlake.edu.cn/software/smr/#SMR&HEIDIanalysis). The HEIDI (heterogeneity in dependent instruments) test was applied to differentiate functional associations from linkage effects. LD correlation between SNPs was estimated using 1000 Genomes Project data for Europeans.
```bash
# 1 Input file preparation
bash ./1.run.1.sh
"#!/bin/bash
studies=("BLUEPRINT" "Schmiedel_2018" "Schmiedel_2018" "FUSION" "TwinsUK")
sample_groups=("neutrophil" "CD4_T-cell_anti-CD3-CD28" "CD8_T-cell_anti-CD3-CD28" "adipose_naive" "fat")
# Loop through each study and sample group combination
for i in "${!studies[@]}"; do
    study=${studies[$i]}
    sample_group=${sample_groups[$i]}
    echo "Running analysis for $study and $sample_group..."
    # Run the R script with the current study and sample group
    Rscript ./1.smr.file.prep.R "$study" "$sample_group"
done"

# 2 SMR/HEIDI analysis
bash ./2.smr.run.sh
```

## 4 Pairwise fixation index (Fst)

## 5 RNA-seq analysis
#### 5.1 This section describes the RNA-seq analysis pipeline. 
Raw RNA-seq reads were trimmed using Trim Galore (v0.6.2) and aligned to human genome (hg38) using the HISAT2 (v2.1.0). Transcript quantification was performed using featureCounts (v1.6.2) with GENCODE v31 annotations. The bigwig files normalised by RPKM (Reads Per Kilobase per Million mapped reads) were generated using the bamCoverage function of deepTools (version 3.3.1). Differential gene expression analysis was conducted using DESeq2 on raw read counts. 
```
sbatch --wrap="bash 0.get.key.name.sh && bash 1.trimming.sh && bash 2.hisat2.mapping.sh && bash 3.featureCounts.sh && 4.RNASeQC.sh && 5.make.bigwigs.sh"
```
#### 5.2 PCA & Differential expression analysis 
related to Fig.4 & 5
```
parallel ::: "Rscript ./DEseq2_T.cells_fig.4.R" "Rscript ./DEseq2_T.cells_fig.5.R" "Rscript ./DEseq2-CD4-RNAsesq-Rapamycin-PRJNA532911.R"
```
#### 5.3 Cell type deconvolution
Cell-type deconvolution was performed with [CIBERSORTx](https://cibersortx.stanford.edu/) using a reference panel derived from a [sepsis whole blood scRNA-seq dataset](https://pubmed.ncbi.nlm.nih.gov/37095375/). A signature matrix was built by the Create Signature Matrix analysis module with parameters min. expression = 0.25, replicates = 100 and sampling = 0.5. The CIBERSORTx absolute scores of each cell type in bulk samples were then obtained using the mixture file (Bulk RNAseq count matrix normalised by DESeq2), the signature matrix derived from single cell RNA-seq, the single cell reference matrix for S-mode batch correction and with 100 permutations via the Impute Cell Fractions analysis module.
```bash
# association between cell absolute scores and sepsis SRS endotypes
Rscript ./CIBERSORTx_linear.regression_fig.S1.R
```

#### 5.4 UMAP projection of sepsis whole blood scRNA-seq data 
```bash
Rscript ./resul-Seurat.R
```

## 6 ATAC-seq analysis
This section describes the ATAC-seq analysis pipeline.
Sequencing reads for ATAC-seq were aligned to the human genome (hg38) using Bowtie2 (v2.2.5). Data were filtered for quality control using Picard (v2.0.1) and Samtools (v1.9) before peak calling with MACS2 (v2.1.0). Differential peak analysis was performed using DESeq2, considering peaks present in at least 30% of samples. Potential batch effects and/or technical variation were assessed through principal component analysis and incorporated as covariates in the DESeq2 design formula. 

#### 6.1 ATAC-seq analysis pipeline
Detailed scripts are [publicly available](https://pubmed.ncbi.nlm.nih.gov/37388915/).
```bash
# Reads alignment and peak calling 
/ATAC_seq_analysis_slurm_hg38/ATAC_MAIN \
--fastq ./raw.ATACseq/ \
--out ./Results_ATAC \
--min 3 --noidr 8
```
#### 6.2 PCA & Differential chromatin accessibility analysis
```bash
Rscript ./02.ATACseq_DESeq2.R
```

## 7 ChIP-seq analysis
This section describes the analysis for hMeDIP-seq. Raw sequencing reads trimmed using Trim Galore (v0.6.2), and aligned to human genome (hg38) using the BWA-mem alignment algorithm (v0.7.12). The binary alignment and map (BAM) files were filtered to remove reads with a mapping quality score less than 10 and duplicate reads using SAMtools (v1.9) and Picard (v2.21.1). The normalised fold enrichment tracks over the corresponding input controls were generated by using the callpeak function with the --SPMR flag, then passing the bedgraph outputs into the bdgcmp function of MACS2 and the bedGraphToBigWig tool.
```bash
sbatch --wrap="bash 0.get.key.name.sh && bash 1.trimming.sh && bash 2.BWA.mapping.sh && bash 3.MACS2.peak.call.sh"
```

## 8 Single guide RNA (sgRNA) design
We designed and selected top ranked single guide RNA (sgRNA) based on the scoring metrics using [FlashFry](https://bmcbiol.biomedcentral.com/articles/10.1186/s12915-018-0545-0).
```
bash ./Create.Database.sh
bash ./Find.Score.gRNAs.sh
```

## Contact
Ping Zhang (ping.zhang@well.ox.ac.uk)

