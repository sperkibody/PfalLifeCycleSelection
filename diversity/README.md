# Diversity analysis 
These files calculate and compare diversity statistics from Pf7 VCFs, after VCF filtering steps. 

count_ns_syn_sites_per_transcript.pl: Generate site counts and proportions for genes. 

min_expression_2.5.txt: contains basic set of genes expressed in assay. Useful for specifying genes for use in initial statistic generation files. 

Pf7_missing_data_calc.ipynb: Calculates missing data proportion in Pf7 for given gene set (input). Combines estimates with site count estimate file from count_ns_syn_sites_per_transcript.pl to make props_adj_breadth.txt. 

props_adj_breadth.txt: File containing site count estimates from count_ns_syn_sites_per_transcript.pl. Includes preliminary breadth estimates before final DE filter (this column is removed before analysis in statistical analysis and plotting scripts, but acts as a place holder at this stage). File is used for pnps.py statistic generate scripts. 

diversity_df_filtering.R: After runpi.sh executes and results are put into data directory, filter and prepare diversity statistics as one dataframe for subsequent analyses. 

(Fig. S2, S13)
pf7_missingness_viz.R: analyze pf7 variant call missingness data across genome. Also looks at proportion of FFD sites across gene set. 

(Fig. 3, S8, Table S1, Table S2)
pnps.py: calculate main diversity statistics from VCF files using scikit-allel. 
stage_analysis.R: Specify gene sets (primary or secondary) at beginning of file. Statistically compare and plot diversity statistics by stage. 
runpi.sh: Shell file for running pnps.py on a specified population. 

Fst (Fig. 4, Table S3)
fst.py: calculate fst using scikit-allel
fst_plots.R: plot and statistically compare Fst results 

Genome-wide length correlations (Fig. S4-5)
genome_wide_correlations.R: examine and calculate genome-wide correlations of diversity/divergence/sequence statisitcs and total coding length.
