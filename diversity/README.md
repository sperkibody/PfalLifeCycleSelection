## Scripts
These files calculate and compare diversity statistics from Pf7 VCFs, after VCF filtering steps. 

- count_ns_syn_sites_per_transcript.pl: Generate site counts and proportions for genes. 
- Pf7_missing_data_calc.ipynb: Calculates missing data proportion in Pf7 for given gene set (input). Combines estimates with site count estimate file from count_ns_syn_sites_per_transcript.pl to make props_adj_breadth.txt. 

(Fig. 3, S9, Table S1, Table S2)
- pnps.py: calculate main diversity statistics from VCF files using scikit-allel.
- runpi.sh: Shell file for running pnps.py on a specified population.
- diversity_df_filtering.R: After runpi.sh executes and results are put into data directory, filter and prepare diversity statistics as one dataframe for subsequent analyses. 
- stage_analysis.R: Specify gene sets (primary or secondary) at beginning of file. Statistically compare and plot diversity statistics by stage. 

(Fig. S1)
- pf7_missingness_viz.R: analyze pf7 variant call missingness data across genome. Also looks at proportion of FFD sites across gene set if desired. 

Fst (Fig. 4, Table S3)
- fst.py: calculate fst using scikit-allel
- fst_plots.R: plot and statistically compare Fst results 

Genome-wide length correlations (Fig. S4-5)
- genome_wide_correlations.R: examine and calculate genome-wide correlations of diversity/divergence/sequence statisitcs and total coding length.

## Data
- props_adj_breadth.txt: File containing site count estimates from count_ns_syn_sites_per_transcript.pl. Includes preliminary breadth estimates before final DE filter (this column is removed before analysis in statistical analysis and plotting scripts, but acts as a place holder at this stage). File is used for pnps.py statistic generate scripts.
- min_expression_2.5.txt: contains basic set of genes expressed in assay. Useful for specifying genes for use in initial statistic generation files. Alternatively, can use gene_sets_assay/*.txt directory or other desired directory (which can be used to produce output file with genes labeled by gene set).
- divstats/: includes output from runpi.sh and fst.py, before processing into single diversity value file.
- AltTranscriptLengths.csv: Includes transcript lengths for each gene. Used to validate that site counts consider the longest transcript for each gene and adjust counts in 14 special cases where longest transcript was not initially considered.
- diversity_values.csv: Output of diversity_df_filtering.R for use in analysis and plotting files. Includes diversity statistics for 3,699 genes in each of the four populations. NA values denote cases where statistics or labels could not be calculated (e.g. no variants in population for Fst calculations or gene not detected in assay). Pf7 variant call pass rate, site counts, breadth labels included. Each row corresponds to diversity values in one of the 4 populations for a given gene.
- stage_div_comparisons.csv: meta-analysis results for stage-stage comparisons of primary gene set from stage_analysis.R. 
- fst_results.csv: meta-analysis results for stage-stage comparisons from fst_plots.R.
