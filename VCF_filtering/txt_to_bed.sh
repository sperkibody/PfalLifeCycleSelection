#!/bin/sh
#  split
SetDir=$1
BASE1="/seq/plasmodium/Pfal_life_cycle_sel"
mkdir $BASE1/${SetDir%%}_bed
for FILE in $(basename -s .txt ${BASE1%%}/gene_sets_${SetDir%%}/*)
do
	##### VCF file for specified life stage
	grep -F -f $BASE1/gene_sets_${SetDir%%}/${FILE%%}.txt /seq/plasmodium/data/bed/PlasmoDB-31_Pfalciparum3D7_gene_coord_CORE_ONLY.bed > $BASE1/${SetDir%%}_bed/${FILE%%}.bed
done
