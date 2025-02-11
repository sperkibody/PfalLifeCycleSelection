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
