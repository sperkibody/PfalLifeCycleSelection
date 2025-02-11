# DFE-alpha jackknifing
Scripts executing DFE-alpha based on diversity and divergence datasets using a leave-one-out approach. 

- divergence_filter.py: updates divergence estimates with jackknifed gene excluded. Can be modified to change what species is used a reference for divergence data. Also can be modified to use S or FFD sites as neutral. Writes divergence estimates for filtered genes to temporary divergence files in corresponding directory (make sure to create directory for this and specify name/path appropriately).

- sfs_dfe.py: updates dfe estimates with jackknifed gene excluded. Can be modified to change what population data is used as input. Also can be modified to use S or FFD sites as neutral.  

- dfe_jackknife.sh: runs divergence_filter.py, sfs_dfe.py, and dfe-alpha executables via jackknifing approach, excluding one gene in specified gene set at a time. Specific if S or FFD sites are used as neutral. Also specify what population dataset is used. Config files used can also be modified to specify file paths and model (epoch number, generation number, population size at each epoch, etc.) as described in DFE-alpha documentation. 

- dfe_files: Directory for tmp files, including sfs_dfe.py output, divergence_filter.py output, and jackknife exclude gene ID for each jackknife run. Includes divergence estimates subset for gene sets of interest (e.g. primary gene sets, m3_final)

- config_expansion: includes config files used for DFE-alpha analysis. Files are set to run with demographic expansion model. 

- DFE_convert_txt_to_csv.ipynb: Converts output from DFE analysis pipeline into CSV for subsequent graphing and analysis (input and output filenames need to be modified/specified). 

(Fig. 5, S10-12)
- dfe_plots.R: visualizes and summarizes statistics on resulting data. Statistics are compiled into .csv format by DFE_convert_txt_to_csv.ipynb. 
