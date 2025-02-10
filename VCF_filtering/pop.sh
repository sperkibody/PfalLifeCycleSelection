#!/bin/sh
use .bcftools-1.9

#### specify chrom./contig (VCF), study pop, and life stage (gene set)
Chrom=$1
BASE1="/seq/plasmodium/Pfal_life_cycle_sel"
BASE="${BASE1%%}/${InputVCF%%}"
### iterate through all contigs (chromosomes)
#declare -a arr=("cambodia" "drc" "tanzania" "ghana")
declare -a arr=("drc")
for StudyPop in "${arr[@]}"
do
	fn="${BASE1%%}/pf7_vcf/Pf3D7_${Chrom%%}_v3_qc_full.pf7.vcf.gz"
	##### Filter pf7 data for population
	bcftools view -S ${BASE1%%}/pf7_pops/${StudyPop%%}.txt $fn > ${BASE1%%}/pf7_vcf/${Chrom%%}_${StudyPop%%}.vcf.gz
done
