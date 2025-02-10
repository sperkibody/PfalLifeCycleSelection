breadth_divergence_by_gene_praefal.csv: output from run of divergence_deidentified.ipynb including divergence estimates for Pfal/Ppraefal.

breadth_divergence_by_gene_preich.csv: output from run of divergence_deidentified.ipynb including divergence estimates for Pfal/Preich. 

gene_sets_assay: directory of genes with detectable expression in the assay (minimal filter). Includes preliminary breadth estimates which are later overwritten by DE filter. Can be used as input to ortho_match.ipynb to include full set of detected assay genes. 

(Fig. S8)
final_stats_breadth.R: Use diversity and divergence estimates with breadth labels to calculate breadth correlations and make breadth plots in supp fig. 8 

(Fig. 2)
dNdSfig2.R: Takes output file from divergence_deidentified.ipynb. Calculates dN/dS from estimates in intermediate .csv file generated. 

divergence_deidentified.ipynb: Takes fasta output from plasmoDB search (after ortho_match.ipynb list generation) and calculates dN and dS for sequences. 

lifestage_orthologs: contains files including ortholog group lists downloaded from PlasmoDB. Also serves as intermediate file directory for ouputs from ortho_match and plasmoDB searches. 

ortho_match.ipynb: Takes input from PlasmoDB (files listing ortholog groups for species pairs). Outputs gene IDs as txt files for input into plasmoDB search strategy, described in file. 
