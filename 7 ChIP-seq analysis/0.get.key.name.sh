#!/bin/bash

ls fastq.file/ | awk -F '.fastq.gz' '{print $1}'  > sample.key.txt

ls fastq.file/ | awk -F '_1.fastq.gz|_2.fastq.gz' '{print $1}' | uniq > sample.key.txt
