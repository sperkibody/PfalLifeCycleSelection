#!/bin/sh

#  final_qc.sh
#  remove calls with DP < 5 and het calls
#
use .vcftools-0.1.14
use Tabix
use .gatk-4.2.3.0

pop=$1
BASE="/seq/plasmodium/Pfal_life_cycle_sel"
	
##### FILTER FOR MIN DP 5
vcftools --gzvcf ${BASE%%}/pf7_vcf/pf7_${pop%%}_coi1.vcf.gz --out ${BASE%%}/pf7_vcf/pf7_${pop%%}_qc1 --minDP 5 --recode --recode-INFO-all
	
##### OVERWRITE HET CALLS
gatk VariantFiltration \
-V ${BASE%%}/pf7_vcf/pf7_${pop%%}_qc1.recode.vcf \
-O ${BASE%%}/pf7_vcf/pf7_${pop%%}_qc2.vcf \
--genotype-filter-expression "isHet == 1" \
--genotype-filter-name "isHetFilter"
	
rm ${BASE%%}/pf7_vcf/pf7_${pop%%}_qc1.recode.vcf
	
gatk SelectVariants \
-V ${BASE%%}/pf7_vcf/pf7_${pop%%}_qc2.vcf \
--set-filtered-gt-to-nocall \
-O ${BASE%%}/pf7_vcf/pf7_${pop%%}_qc.vcf
	
rm ${BASE%%}/pf7_vcf/pf7_${pop%%}_qc2.vcf

#### compress and index file
bgzip ${BASE%%}/pf7_vcf/pf7_${pop%%}_qc_coi1.vcf
tabix -p vcf ${BASE%%}/pf7_vcf/pf7_${pop%%}_qc_coi1.vcf.gz
