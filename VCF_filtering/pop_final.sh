#!/bin/sh
use .bcftools-1.9

BASE1="/seq/plasmodium/Pfal_life_cycle_sel"

### iterate through all contigs (chromosomes)
declare -a arr=("cambodia" "drc" "tanzania" "ghana")
for StudyPop in "${arr[@]}"
do
	fn="${BASE1%%}/pf7_vcf/pf7_${StudyPop%%}_qc.recode.vcf"
	##### Filter pf7 data for selected samples
	bcftools view -S ${BASE1%%}/pf7_pops/${StudyPop%%}.txt $fn > ${BASE1%%}/pf7_vcf/pf7_${StudyPop%%}_ibd_corrected.vcf
done
