## Sample selection: 

- Pf7_fws.txt: Fws estimates per sample from Pf7. 
- *.clusters.0.25.txt: hmmIBD clonal clusters from population samples. Choose single representative sample from clusters in pf7_filtering. 

## VCF filtering: 

- PlasmoDB-31_Pfalciparum3D7_gene_coord_CORE_ONLY.bed: bed file with genes falling in the core region. This can be subset to include genes of interest listed in txt files. 
- gene_sets_all.txt: Includes all core genes in the MCA SS2 assay. 
- ffd_codon_changes.txt: File listing FFD codon changes (as annotated by SNPSift) for custom script filtering VCF-encoded FFD variants. 

## Diversity:

## Divergence: 
- Pfal_Preichenowi_orthos.csv: Ortholog groups from PlasmoDB comparing P. reichenowi and P. falciparum. May include multiple groupings for single <i>P. falciparum</i> genes.
- Pfal_Ppraefal_orthos.csv: Ortholog groups from PlasmoDB comparing P. praefalciparum and P. falciparum. May include multiple groupings for single P. falciparum genes.
## DFE-alpha: 
- dfe_files: Directory for tmp files, including sfs_dfe.py output, divergence_filter.py output, and jackknife exclude gene ID for each jackknife run. Includes divergence estimates subset for gene sets of interest (e.g. primary gene sets, m3_final)
- dfe_output_files: includes DFE-alpha output processed as .csv files from .txt files. 
