#!/bin/sh
use .bcftools-1.9
pop=$1
echo "sorting step" 
for file in *_${pop%%}.vcf.gz; do bcftools sort $file -o sorted_${file}; done
for file in *_${pop%%}.vcf.gz; do rm $file; done
echo "concat step"
bcftools concat -o pf7_${pop%%}_qc.vcf.gz sorted_*_${pop%%}.vcf.gz
echo "tidying step"
#for file in sorted_*_${pop%%}.vcf.gz; do rm $file; done
