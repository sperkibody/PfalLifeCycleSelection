## Sample selection: 

- Pf7_fws.txt: Fws estimates per sample from Pf7. 
- *.clusters.0.25.txt: hmmIBD clonal clusters from population samples. Choose single representative sample from clusters in pf7_filtering.
- pf7_pops/*.txt: Lists of sample IDs selected for each population after all filtering steps, including quality filtering, relatedness filtering, and monogenomic sample selection. 

## Gene set selection
- breadth_final.csv: breadth labels from 50% expression and DE approach
- m3_combined.csv: primary gene sets. Supp file 2.
- stage_specific.csv: secondary gene sets. Supp file 3. 

## VCF filtering: 

- PlasmoDB-31_Pfalciparum3D7_gene_coord_CORE_ONLY.bed: bed file with genes falling in the core region. This can be subset to include genes of interest listed in txt files. 
- gene_sets_all.txt: Includes all core genes in the MCA SS2 assay. 
- ffd_codon_changes.txt: File listing FFD codon changes (as annotated by SNPSift) for custom script filtering VCF-encoded FFD variants. 

## Diversity:
- props_adj_breadth.txt: File containing site count estimates from count_ns_syn_sites_per_transcript.pl. Includes preliminary breadth estimates before final DE filter (this column is removed before analysis in statistical analysis and plotting scripts, but acts as a place holder at this stage). File is used for pnps.py statistic generate scripts.
- min_expression_2.5.txt: contains basic set of genes expressed in assay. Useful for specifying genes for use in initial statistic generation files. Alternatively, can use gene_sets_assay/*.txt directory or other desired directory (which can be used to produce output file with genes labeled by gene set).
- divstats/: includes output from runpi.sh and fst.py, before processing into single diversity value file.
- AltTranscriptLengths.csv: Includes transcript lengths for each gene. Used to validate that site counts consider the longest transcript for each gene and adjust counts in 14 special cases where longest transcript was not initially considered.
- helb_seropos.txt: Includes gene IDs for antigens identified by Helb. et al in antibody reactivity screen (https://pubmed.ncbi.nlm.nih.gov/26216993/). Data from supp file. 
- diversity_values.csv: Output of diversity_df_filtering.R for use in analysis and plotting files. Includes diversity statistics for 3,699 genes in each of the four populations. NA values denote cases where statistics or labels could not be calculated (e.g. no variants in population for Fst calculations or gene not detected in assay). Pf7 variant call pass rate, site counts, breadth labels included. Each row corresponds to diversity values in one of the 4 populations for a given gene.
- stage_div_comparisons.csv: meta-analysis results for stage-stage comparisons of primary gene set from stage_analysis.R. 
- fst_results.csv: meta-analysis results for stage-stage comparisons from fst_plots.R.
## Divergence: 
- Pfal_Preichenowi_orthos.csv: Ortholog groups from PlasmoDB comparing <i>P. reichenowi</i> and <i>P. falciparum</i>. May include multiple groupings for single <i>P. falciparum</i> genes.
- Pfal_Ppraefal_orthos.csv: Ortholog groups from PlasmoDB comparing <i>P. praefalciparum</i> and <i>P. falciparum</i>. May include multiple groupings for single <i>P. falciparum</i> genes.
- gene_sets_assays: Detected genes with initial breadth estimates (later overwritten). Includes ~4930 relevant genes detected in assay for the purposes of including those genes in dN/dS analysis (breadth estimates/filenames serve as placeholders and can be replaced as desired to do dN/dS calculations for a different set of labeled genes).
- lifestage_orthologs: directory name that can be used to hold files including ortholog group lists downloaded from PlasmoDB. Also serves as intermediate file directory for ouputs from ortho_match and plasmoDB searches. 
- breadth_divergence_by_gene_praefal.csv: output from run of divergence_deidentified.ipynb including divergence estimates for Pfal/Ppraefal.
- breadth_divergence_by_gene_preich.csv: output from run of divergence_deidentified.ipynb including divergence estimates for Pfal/Preich. 
## DFE-alpha: 
- dfe_files: Directory for tmp files, including sfs_dfe.py output, divergence_filter.py output, and jackknife exclude gene ID for each jackknife run. Includes divergence estimates subset for gene sets of interest (e.g. primary gene sets, m3_final)
- dfe_output_files: includes DFE-alpha output processed as .csv files from .txt files. 
