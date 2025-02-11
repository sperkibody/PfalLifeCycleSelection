## Scripts
pf7_filtering: Choose samples from representative population using Pf7 metadata, Pf7 Fws information (to filter out complex infections). Use clusters from hmmIBD in second part of code to filter samples down to cluster representatives. 
## Data: 
- Pf7_fws.txt: Fws estimates per sample from Pf7. 
- *.clusters.0.25.txt: hmmIBD clonal clusters from population samples. Choose single representative sample from clusters in pf7_filtering.
- pf7_pops/*.txt: Lists of sample IDs selected for each population after all filtering steps, including quality filtering, relatedness filtering, and monogenomic sample selection. 
