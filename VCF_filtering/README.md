## Scripts 
Scripts used to filter Pf7 VCFs. 
- txt_to_bed.sh: converts gene list (e.g. in gene_sets_all.txt) to core gene regions in bed format. This will be used for filtering step in split.sh and can be used on more specific gene sets. 
- pf1.sh: example script to download Pf7 VCFs locally, which are split by chromosome. Needs to be edited to specify chromosome for download (can download up to all 14). 
- bgzip.sh: zips VCF files. 
- concat.sh: script to combine 14 separate Pf7 VCFs into single file for each population dataset.   
- pop.sh: filter for population samples after initial qc (qc.sh). List of population samples can be output by scripts in sample filtering directory.   
- qc.sh: Initial QC variant filtering. Includes filtering of variants with Low VQSLOD, filtering for biallelic SNPs, QC Pass. 
- final_qc.sh: Final QC step for population-level VCFs, including filtering of calls with DP<5 or het calls. 
- VCF_filtering_for_hmmIBD_array.sh: VCF filtering step for input into hmmIBD to identify clonal clusters. 
- pop_final.sh: Additional filtering of samples at the population level after clonal samples have been identified and pruned via hmmIBD and sample filtering script in directory. 
- split.sh: Filter VCF for coding SNPs including variants in genes of interest, specified by bed file. 
- split_snp.sh: Split gene set VCFs into separate VCF files with NS, S, and FFD SNPs only. 
## Data directories
- PlasmoDB-31_Pfalciparum3D7_gene_coord_CORE_ONLY.bed: bed file with genes falling in the core region. This can be subset to include genes of interest listed in txt files. 
- gene_sets_all.txt: Includes all core genes in the MCA SS2 assay. 
- ffd_codon_changes.txt: File listing FFD codon changes (as annotated by SNPSift) for custom script filtering VCF-encoded FFD variants. 
