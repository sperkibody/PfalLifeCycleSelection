#!/bin/sh

#  runpi.sh
#  
#
#
use .vcftools-0.1.14
use Tabix
use Anaconda3
pop=$1
SetDir=$2
BASE1="/seq/plasmodium/Pfal_life_cycle_sel"
BASE="${BASE1%%}/pf7_vcf"

for FILE in $(basename -s .txt ${BASE1%%}/gene_sets_${SetDir%%}/*)
do
	python ${BASE1%%}/pnps.py $FILE $SetDir $pop
done
