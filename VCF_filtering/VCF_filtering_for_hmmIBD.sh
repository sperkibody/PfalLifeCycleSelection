use .vcftools-0.1.14
use Tabix

InputVCF=$1
OutputBase=$2

vcftools --gzvcf ${InputVCF} \
	--positions /seq/plasmodium/data/bed/sfs_preferred_Pfal_snps_June2020.txt \
	--remove-filtered-all \
	--maf 0.05 \
	--recode \
	--recode-INFO-all \
	--out ${OutputBase}.for.hmmIBD

perl vcf_to_hmmIBD_input_multiallelic_w_min_depth.pl ${OutputBase}.for.hmmIBD.recode.vcf 5 ${OutputBase}.for.hmmIBD
