#!/bin/sh
use .vcftools-0.1.14
use .bcftools-1.9
use Tabix
use Java-1.8

# use R-4.1
#### specify chrom./contig (VCF), study pop, and life stage (gene set)
pop=$1
SetDir=$2
BASE1="/seq/plasmodium/Pfal_life_cycle_sel"
BASE="${BASE1%%}/pf7_vcf/${pop%%}"

for FILE in $(basename -s .txt ${BASE1%%}/gene_sets_${SetDir%%}/*)
do
	### overwrite old file
	java -Xmx8g -jar /cil/shed/apps/external/snpEff/snpEff-4.1g/snpEff.jar -formatEff -classic -no_shift_hgvs -v Pf3D7v91 ${BASE%%}/${FILE%%}.CoreOnly.SNPsOnly.recode.vcf.gz > ${BASE%%}/${FILE%%}.CoreOnly.SNPsOnly.recode.ann.vcf
	##### NS SNPs only
	cat ${BASE%%}/${FILE%%}.CoreOnly.SNPsOnly.recode.ann.vcf | java -jar /cil/shed/apps/external/snpEff/snpEff-4.1g/SnpSift.jar filter "( EFF[*].EFFECT != 'synonymous_variant' )" > ${BASE%%}/${FILE%%}.CoreOnly.NS_SNPsOnly.recode.vcf
	##### S SNPs only
	cat ${BASE%%}/${FILE%%}.CoreOnly.SNPsOnly.recode.ann.vcf | java -jar /cil/shed/apps/external/snpEff/snpEff-4.1g/SnpSift.jar filter "( EFF[*].EFFECT = 'synonymous_variant' )" > ${BASE%%}/${FILE%%}.CoreOnly.S_SNPsOnly.recode.vcf
	##### S SNPs at FFD sites only 
	# filter for lines in FFD codon change annotation list 
	grep '#CHROM' ${BASE%%}/${FILE%%}.CoreOnly.S_SNPsOnly.recode.vcf > header.vcf
	grep -F -f ${BASE1%%}/ffd_codon_changes.txt ${BASE%%}/${FILE%%}.CoreOnly.S_SNPsOnly.recode.vcf > ${BASE%%}/${FILE%%}.CoreOnly.FFD_SNPsOnly_base.recode.vcf
	cat header.vcf ${BASE%%}/${FILE%%}.CoreOnly.FFD_SNPsOnly_base.recode.vcf > ${BASE%%}/${FILE%%}.CoreOnly.FFD_SNPsOnly.recode.vcf
	rm header.vcf
	rm ${BASE%%}/${FILE%%}.CoreOnly.FFD_SNPsOnly_base.recode.vcf

	bgzip ${BASE%%}/${FILE%%}.CoreOnly.NS_SNPsOnly.recode.vcf
	bgzip ${BASE%%}/${FILE%%}.CoreOnly.S_SNPsOnly.recode.vcf
	bgzip ${BASE%%}/${FILE%%}.CoreOnly.FFD_SNPsOnly.recode.vcf
	tabix -p vcf ${BASE%%}/${FILE%%}.CoreOnly.NS_SNPsOnly.recode.vcf.gz
	tabix -p vcf ${BASE%%}/${FILE%%}.CoreOnly.S_SNPsOnly.recode.vcf.gz
	tabix -p vcf ${BASE%%}/${FILE%%}.CoreOnly.FFD_SNPsOnly.recode.vcf.gz
	rm ${BASE%%}/${FILE%%}.CoreOnly.SNPsOnly.recode.ann.vcf
done

