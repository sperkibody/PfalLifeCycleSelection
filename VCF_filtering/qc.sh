#!/bin/sh
use .bcftools-1.9
# use R-4.1
#### specify chrom./contig (VCF), study pop, and life stage (gene set)
Chrom=$1
BASE1="/seq/plasmodium/Pfal_life_cycle_sel"
BASE="${BASE1%%}/${InputVCF%%}"
### iterate through all contigs (chromosomes)
fn="${BASE1%%}/pf7_vcf/Pf3D7_${Chrom%%}_v3.pf7.vcf.gz"
echo $fn
#### Apply quality filters
bcftools view \
--include 'FILTER="PASS" && N_ALT=1 && CDS==1 && TYPE="snp" && VQSLOD>2.0' \
--output-type z \
--output-file ${BASE1%%}/pf7_vcf/Pf3D7_${Chrom%%}_v3_qc_full.pf7.vcf.gz \
${fn%%}
bcftools index --tbi ${BASE1%%}/pf7_vcf/Pf3D7_${Chrom%%}_v3_qc_full.pf7.vcf.gz
