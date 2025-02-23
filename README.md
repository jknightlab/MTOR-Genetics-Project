# MTOR-Genetics-Project
#### _P. Zhang et al., Context-specific regulatory genetic variation in MTOR dampens Neutrophil-T cell crosstalk in sepsis, modulating disease. (2025)_

![Screenshot](output/model.png)
### 
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
  - 1.4
- [2 Survival analysis](#2-survival-analysis)
- [3 Summary data-based Mendelian randomisation](#3-summary-data-based-mendelian-randomisation)
  - 3.1 Input file preparation
  - 3.2 SMR/HEIDI analysis
- [4 Pairwise fixation index (Fst)](#4-pairwise-fixation-index-fst)
  - 4.1 ###
  - 4.2 ###
- [5 RNA-seq analysis](#5-rna-seq-analysis)
  - 4.1 RNA-seq analysis pipeline
  - 4.2 PCA & Differential expression analysis
  - 4.3 UMAP projection of sepsis whole blood scRNA-seq data 
  - 4.3 Cell type deconvolution
- [6 ATAC-seq analysis](#6-atac-seq-analysis)
  - 5.1 ATAC-seq analysis pipeline
  - 5.2 PCA & Differential expression analysis
- [7 ChIP-seq analysis](#7-chip-seq-analysis)
  - 6.1 ChIP-seq analysis pipeline
- [8 Single guide RNA (sgRNA) design](#8-single-guide-rna-sgrna-design)

## 1 Identification of context specific eQTLs
#### 1.1 Interaction analysis
#### 1.2 Pairwise linkage disequilibrium (R²) for the MTOR locus
#### 1.3 Colocalisation analysis

## 2 Survival analysis

## 3 Summary data-based Mendelian randomisation
We performed the SMR analysis using eQTLs as instrumental variables to identify genes whose expression is associated with T2D risk due to pleiotropy and/or causality. Genes were included in the analysis if they had at least one cis-eQTL (P < 5e⁻⁸) within a 2 Mbp window around GWAS loci, following the default settings of the SMR tool (v1.3.1). The HEIDI (heterogeneity in dependent instruments) test was applied to differentiate functional associations from linkage effects. LD correlation between SNPs was estimated using 1000 Genomes Project data for Europeans.

## 4 Pairwise fixation index (Fst)

## 5 RNA-seq analysis
This section describes the RNA-seq analysis pipeline.

## 6 ATAC-seq analysis

## 7 ChIP-seq analysis
This section describes the analysis for hMeDIP-seq.

## 8 Single guide RNA (sgRNA) design
We designed and selected top ranked single guide RNA (sgRNA) based on the scoring metrics using FlashFry 
