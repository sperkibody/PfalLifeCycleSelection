library(ggplot2)
library(tidyverse)
library(dplyr)
###############MAIN DIVERSITY STATISTICS##########################
setwd('~')
path="/Users/sap8772/Documents/Neafsey_lab/"
fn_out = paste(path,'diversity_values.csv', sep="")
breadth_labels <- read.csv(paste(path,'de50_combined.csv',sep="")) %>% dplyr::rename(c("GENE"="value", "breadth"="variable"))
props=read.csv(paste(path,'props_adj_breadth.txt',sep = ""), sep='\t')
props$X <- gsub("location=", "", props$COORD)
props$X <- sapply(props$X, function(x) substr(x, 0, nchar(x)-3))
props <- props %>% arrange(GENE, desc(TOTAL_CODING_LENGTH))
#props %>% filter(value <= 0.7 | ag) %>% nrow()
# retain only the longest transcript 
props <- props[!duplicated(dplyr::select(props, GENE)),]
# exclude any gene outside of the core region or with all variants filtered at variant quality level 
all_core <- read.csv(paste(path,'all.bed',sep=""),
                     sep="\t",header=FALSE) 
min_expressed_assay <- read.table(paste(path,'min_expression_2.5.txt',sep=""),sep='\n')
sum(!props$GENE %in% all_core$V6)
sum(!props$GENE %in% min_expressed_assay$V1)
qc_variants_filtered <- c("PF3D7_0115000","PF3D7_0201900") # variants not found in any pop VCFs after variant QC filtering step 
props <- props[props$GENE %in% all_core$V6,]
props <- props[props$GENE %in% min_expressed_assay$V1,]
props <- props[!props$GENE %in% qc_variants_filtered,]

stat_df <- function(stat, pops = c('cambodia', 'drc', 'tanzania', 'ghana'), 
                    args = c('all'),divpath=divpath){
  for(pop in pops){
    setwd(paste(divpath,pop,'/',sep=''))
    #p0 = paste(stage, '_fst.csv', sep='')
    for(i in 1:length(args)){
      print(paste(args[i],'_', stat, '.csv', sep = ''))
      ls_r <- read.csv(paste(args[i],'_', stat, '.csv', sep = ''))
      ls <- data.frame(ls_r)
      if(i==1){
        pidf <- full_join(props,data.frame(ls))
      }
      else{
        ls <- full_join(props,data.frame(ls))
        pidf <- rbind(pidf,ls)}
    }
    pidf$pop <- rep(pop, nrow(pidf))
    if (pop=='cambodia'){ pidf_all <- pidf }
    else{pidf_all <- rbind(pidf_all,pidf)}
  }
  return(pidf_all)
}
### new df 
divpath=paste(path,'divstats_td_fixed/',sep="")
args0=c('all')
# read in all diversity statistic files 
pi_ns_all <- stat_df('pi_ns', args=args0,divpath=divpath) 
pi_s_all <- stat_df('pi_s', args=args0,divpath=divpath) 
pnps_all <- stat_df('pi_ns_div_pi_s', args=args0,divpath=divpath) 
tajimas_d_all <- stat_df('tajimas_d', args=args0,divpath=divpath) 
fst_all <- read.csv(paste(divpath,'all_fst.csv',sep='/'))
# convert negative Fst values (which can be produced by Hudson estimator) to 0 
for(c in 2:7){
  fst_all[,c][fst_all[,c] < 0] <- 0
}
diversity <- full_join(pi_ns_all, pi_s_all)
diversity <- full_join(diversity, pnps_all)
diversity <- full_join(diversity, tajimas_d_all)
diversity <- left_join(diversity, fst_all)
diversity %>% group_by(pop) %>% summarize(n()) %>% View()

# FILTERING STEPS 
# filter out antigens for main analysis 
ags <- read.table(paste(path,'helb_seropos.txt',sep=""))[[1]]
diversity$ag <- diversity$GENE %in% ags
# check for 14 special cases of splice variants and correct statistics based on site count estimates 
coding_length_correction <- read.csv(paste(path,'AltTranscriptLengths.csv',sep=""))
head(coding_length_correction)
cds_missing <- coding_length_correction %>% filter(CDSMissingFromLongest > 0)
div_corr <- diversity[diversity$GENE %in% cds_missing$GENE,] %>% left_join(cds_missing, by="GENE")
cdsratio <- (div_corr$LongestTranscript / div_corr$TotalNonRedundantCDS)
pi_s_corrected <- div_corr$pi_s * cdsratio
pi_ns_corrected <- div_corr$pi_ns * cdsratio
pnps_corrected <- pi_ns_corrected / pi_s_corrected
div_corr$pi_s <- pi_s_corrected
div_corr$pi_ns <- pi_ns_corrected
div_corr$pi_ns_div_pi_s <- pnps_corrected
div_corr <- div_corr[,1:(ncol(div_corr)-ncol(coding_length_correction)+1)]
diversity <- diversity[!diversity$GENE %in% cds_missing$GENE,]
# rejoin with updated genes 
diversity <- rbind(diversity,div_corr)
# make sure duplicates have been removed at GENE level
diversity %>% group_by(pop,GENE) %>% summarize(n=n()) %>% filter(n>1) %>% dplyr::select(GENE)

###############FINISH FILTERING STEPS AND SAVE DF##########################
diversity <- diversity %>% filter(value > 0.7)
diversity <- diversity %>% filter(!ag) 

# add 0s to invariant sites that were filtered from VCFs 
diversity$pi_ns[is.na(diversity$pi_ns)] <- 0
diversity$pi_s[is.na(diversity$pi_s)] <- 0
# update pnps estimates to reflect cases where pi_ns = 0 and pi_s > 0
diversity$pnps <- diversity$pi_ns/diversity$pi_s
# ensure that updated estimate only changed Inf, not numeric values
ggplot(diversity, aes(x=pi_ns_div_pi_s,y=pnps)) + geom_point()
# joining with breadth, filtering out NAs for breadth or pn/ps 
diversity <- left_join(diversity, breadth_labels)
View(diversity)
div_processed <- diversity %>% dplyr::select(-c(variable,ID,X,STOP_CODONS)) %>% 
  dplyr::rename(c("Antigen"="ag", "missing_data_proportion"="value", 
                  "tajimas_d" = "tajimas_d_all", "population"="pop")) %>% 
  dplyr::select(COORD, GENE, NAME, PROP_NS, PROP_SYN, PROP_FFD, TOTAL_CODING_LENGTH, 
                population, pi_ns, pi_s, 
                pnps, tajimas_d, fst.drc.ghana, fst.drc.tanzania, fst.drc.cambodia,
                fst.ghana.tanzania, fst.tanzania.cambodia, fst.ghana.cambodia, 
                missing_data_proportion, breadth) 

write.csv(div_processed, fn_out, row.names = FALSE)

div_processed %>% group_by(population)  %>% summarize(n())
