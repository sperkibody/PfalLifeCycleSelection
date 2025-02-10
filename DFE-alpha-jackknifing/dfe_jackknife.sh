#!/bin/sh

use Anaconda3
use .gsl-2.6
pop=$1
SetDir=$2
neut=$3
BASE1="/seq/plasmodium/Pfal_life_cycle_sel"
#BASE="${BASE1%%}/pf7_vcf"

# iterate through workflow for all life stages 
for FILE in $(basename -s .txt ${BASE1%%}/gene_sets_${SetDir%%}/*)
do
	echo ${FILE%%}
	while IFS= read -r line || [[ -n "$line" ]]; do
			# pull out gene to exclude in this iteration
			echo $line > dfe_files/jackknife_exclude/${FILE%%}.txt
			echo 'Exclude: ' $line
			# write new sfs
			python ${BASE1%%}/sfs_dfe.py ${FILE} ${SetDir%%} ${pop%%} ${neut%%}
			# write new between-species divergence file (population-agnostic) 
			python ${BASE1%%}/divergence_filter.py ${FILE} ${neut%%}
			# run DFE-alpha executables
			cd ${BASE1%%}/dfe-alpha-release-2.16
			echo 'Running DFE-alpha...'
			# neutral sites
			./est_dfe -c config_expansion/config_sc0_${FILE%%}.txt
			# selected (ns) sites
			./est_dfe -c config_expansion/config_sc1_${FILE%%}.txt
			# estimate adaptive divergence
			./est_alpha_omega -c config_ao_${FILE%%}.txt > final_output/dfe_alpha_results_${FILE%%}_${pop%%}_exclude_${line%%}.txt
			cd ..
	done < ${BASE1%%}/gene_sets_${SetDir%%}/${FILE%%}.txt
done
