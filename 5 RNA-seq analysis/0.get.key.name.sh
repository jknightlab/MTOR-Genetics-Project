#!/bin/bash

ls raw.RNAseq/ | awk -F '_R1.fastq.gz|_R2.fastq.gz' '{print $1}' | uniq > sample.key_12.txt
