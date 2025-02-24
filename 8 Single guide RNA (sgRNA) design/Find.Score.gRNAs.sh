#!/bin/bash

REF=/ref/hg38.fa

bedtools getfasta -bed regions_of_interest.bed \
 -fi $REF -fo pilot.regions.fa

java -Xmx4g \
-jar FlashFry-assembly-1.9.3.jar \
 discover \
 --database hg38_Cas9NGG_database \
 --fasta pilot.regions.fa \
 --output pilot.regions.output.txt

java -Xmx4g -jar FlashFry-assembly-1.9.3.jar \
 score \
 --input pilot.regions.output.txt \
 --output pilot.regions.output.Scored.txt \
 --scoringMetrics doench2014ontarget,doench2016cfd,dangerous,hsu2013,minot,rank \
 --database hg38_Cas9NGG_database
