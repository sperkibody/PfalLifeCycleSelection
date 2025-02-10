#!/bin/sh
#  split
use .vcftools-0.1.14
use .bcftools-1.9
use Tabix
pop=$1
SetDir=$2
BASE1="/seq/plasmodium/Pfal_life_cycle_sel"
BASE="${BASE1%%}/pf7_vcf"

for FILE in $(basename -s .txt ${BASE1%%}/gene_sets_${SetDir%%}/*)
do
	##### VCF file for biallelic SNPs within specified life stage
	vcftools --gzvcf ${BASE%%}/pf7_${pop%%}_qc.vcf.gz \
	--bed ${BASE1%%}/${SetDir%%}_bed/${FILE%%}.bed \
	--remove-indels \
	--mac 1 \
	--remove-filtered-all \
	--recode \
	--recode-INFO-all \
	--out ${BASE%%}/${pop%%}/${FILE%%}.CoreOnly.SNPsOnly

	#### compress and index file
	bgzip ${BASE%%}/${pop%%}/${FILE%%}.CoreOnly.SNPsOnly.recode.vcf
	tabix -p vcf ${BASE%%}/${pop%%}/${FILE%%}.CoreOnly.SNPsOnly.recode.vcf.gz
done
