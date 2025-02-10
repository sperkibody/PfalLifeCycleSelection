import warnings
warnings.filterwarnings('ignore')
import allel
import pandas as pd
import numpy as np
import sys
import os
os.chdir('/seq/plasmodium/Pfal_life_cycle_sel')
# specify life stage
lifestage = sys.argv[1]
# specify basic gene set, e.g. m3
set_base = sys.argv[2]
# specify population e.g. Dakar
pop = sys.argv[3]
# neutral sites SYN or FFD 
neut_key = sys.argv[4]

print(lifestage)
notfound = []
exclude = []
with open('dfe_files/jackknife_exclude/'+lifestage+'.txt') as f: 
	for line in f:
		exclude.append(line.strip())
#iterate through all genes in bed file to get start and stop positions
props = pd.read_csv('props_adj_breadth.txt', sep = '\t')
# get total counts of S and NS sites from bed files
ns = 0
# neutral sites--either synonymous or FFD synonymous 
neutral = 0
with open(set_base+'_bed/'+lifestage+'.bed') as f:
	for line in f:
		bdat = line.strip().split('\t')
		gene = bdat[5]
		chrom0 = bdat[0]
		start0 = np.min([int(bdat[1]), int(bdat[2])])
		end0 = np.max([int(bdat[1]), int(bdat[2])])
		if (gene in exclude): continue 
		# multiply site counts by missing data rate estimates 
		region0=str(chrom0)+':'+str(start0)+'-'+str(end0)
		entry = props['COORD'].str.contains(region0)
		if not np.sum(entry): 
			print(gene, region0)
			print('Not found')
			continue
		props1 = props[entry].sort_values(by='TOTAL_CODING_LENGTH')
		ns_gene = props1['NS'].values[0]
		missing_data_prop = props1['value'].values[0]
		ns += ns_gene*missing_data_prop
		neutral_gene = props1.loc[:,neut_key].values[0]
		neutral += neutral_gene*missing_data_prop
ns = np.round(ns)
neutral = np.round(neutral) 
print(ns, 'total ns sites and ', neutral, ' total ' + neut_key + ' sites')

def filter_callset(callset, chrom, start, end):
	# chrom mask 
	chrom_mask = callset['variants/CHROM']==chrom
	# mask for position
	pos_mask = np.isin(callset['variants/POS'], np.arange(start, end))
	# combine masks
	gene_mask = ~(np.logical_and(chrom_mask, pos_mask))
	# filter dataset
	return(callset['calldata/GT'][gene_mask])

# create folded sfs files to input into DFE analysis for S and NS sites across entire life stage 
print('Reading NS variants...')
callset_ns = allel.read_vcf('pf7_vcf/' + pop + '/'  + lifestage + '.CoreOnly.NS_SNPsOnly.recode.vcf.gz', fields=['calldata/GT', 'variants/POS', 'variants/CHROM'])
# filter out excluded gene 
callset_ns = filter_callset(callset_ns, chrom0, start0, end0)
gt = allel.GenotypeArray(callset_ns)
ac = gt.count_alleles()
ns_sfs = allel.sfs_folded(ac)
print(len(ns_sfs))
# calculate the number of invariant NS sites + add to variant sites w/ ac 0  
nonseg_count = np.round(ns - np.sum(ns_sfs)) + ns_sfs[0]
ns_sfs[0] = nonseg_count

print('Reading neutral variants...')
if (neut_key=='SYN'): 
	callset_neutral = allel.read_vcf('pf7_vcf/' + pop + '/' + lifestage + '.CoreOnly.S_SNPsOnly.recode.vcf.gz', fields=['calldata/GT', 'variants/POS', 'variants/CHROM'])
else: 
	callset_neutral = allel.read_vcf('pf7_vcf/' + pop + '/' + lifestage + '.CoreOnly.FFD_SNPsOnly.recode.vcf.gz', fields=['calldata/GT', 'variants/POS', 'variants/CHROM'])
# filter out excluded gene 
callset_neutral = filter_callset(callset_neutral, chrom0, start0, end0)
gt = allel.GenotypeArray(callset_neutral)
ac = gt.count_alleles()
neutral_sfs = allel.sfs_folded(ac)
# calculate the number of invariant neutral sites + add to variant sites w/ ac 0 
nonseg_count = np.round(neutral - np.sum(neutral_sfs)) + neutral_sfs[0]
neutral_sfs[0] = nonseg_count

assert len(neutral_sfs)==len(ns_sfs)

# FIX: VCF CODING AS DIPLOID CREATES SFS WHERE K // 2 = 0. PARASITES ARE HAPLOID. 
ind_haploid = np.arange(0, len(ns_sfs), 2)
# note: for est_dfe folded sfs input, second half of vector should be all 0s 
ns_sfs_haploid = np.zeros(len(ns_sfs))
neutral_sfs_haploid = np.zeros(len(ns_sfs))
for i in range(len(ind_haploid)): 
	ns_sfs_haploid[i] = int(ns_sfs[ind_haploid[i]])
	neutral_sfs_haploid[i] = int(neutral_sfs[ind_haploid[i]])

fn = 'dfe_files/' + lifestage + '_dfe.txt' 
dfe_str = '1\n' + str(len(ns_sfs_haploid)-1) + '\n'
for i in ns_sfs_haploid:
	dfe_str += str(int(i)) + ' '
dfe_str += '\n'
for i in neutral_sfs_haploid:
	dfe_str += str(int(i)) + ' '

print(dfe_str)
with open(fn, 'w') as dfe_file:
	dfe_file.write(dfe_str)

