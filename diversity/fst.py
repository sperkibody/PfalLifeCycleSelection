import sys
import allel
import pandas as pd
import numpy as np
import os
os.chdir('/seq/plasmodium/Pfal_life_cycle_sel')

#Calculate pi and theta via scikit-allel
split = sys.argv[1]
lifestage = sys.argv[2]

def get_hudson(pop1, pop2, ls, region0):
	base = 'pf7_vcf/'
	f1 = base + pop1 + '/' + ls + '.CoreOnly.SNPsOnly.recode.vcf.gz'
	f2 = base + pop2 + '/' + ls + '.CoreOnly.SNPsOnly.recode.vcf.gz'
	fields0 = ['variants/POS','variants/CHROM', 'calldata/GT']
	callset1 = allel.read_vcf(f1, region = region0, fields = fields0)
	callset2 = allel.read_vcf(f2, region = region0, fields = fields0)
	if (callset1 is None) | (callset2 is None): return(None)
	vpos1 = callset1["variants/POS"]
	vpos2 = callset2["variants/POS"]
	vc1 = callset1["variants/CHROM"]
	vc2 = callset2["variants/CHROM"]
	# check for variants in both datasets 
	filter_match_1 = np.zeros(len(vpos1))
	filter_match_2 = np.zeros(len(vpos2))
	for i in range(len(vpos1)):
		found = vpos2==vpos1[i]
		if not np.sum(found): next
		else: 
			chromi = vc1[i]
			# loop through to find variant position match where chromosome also matches
			for j,k in zip(np.array(range(len(vpos2)))[found], vc2[found]):
				if chromi==k: 
					filter_match_1[i]=1
					filter_match_2[j]=1
	# filter gt arrays and convert to allele counts 
	gt1 = allel.GenotypeArray(callset1["calldata/GT"])[filter_match_1>0]
	gt2 = allel.GenotypeArray(callset2["calldata/GT"])[filter_match_2>0]
	if (gt1.shape[0]==0) | (gt2.shape[0]==0): return(None)
	ac1 = gt1.count_alleles()
	ac2 = gt2.count_alleles()
	if (ac1 is None) | (ac2 is None): return(None)
	num, den = allel.hudson_fst(ac1, ac2)
	fst = np.sum(num) / np.sum(den)
	return(fst)

#iterate through all genes in bed file to get start and stop positions; calculate pi and theta per gene  
fst_drc_ghana = {}
fst_drc_tanzania = {}
fst_ghana_tanzania = {}
fst_drc_cambodia = {}
fst_ghana_cambodia = {}
fst_tanzania_cambodia = {}
with open(split+'_bed/'+lifestage+'.bed') as f:
	for line in f:  
		bdat = line.strip().split('\t')
		chrom = bdat[0]
		p0 = int(bdat[1])
		p1 = int(bdat[2])
		gene = bdat[5]
		region0=chrom+':'+bdat[1]+'-'+bdat[2]
		print(region0)
		fst_drc_ghana[region0] = get_hudson('drc', 'ghana', lifestage, region0)
		fst_drc_tanzania[region0] = get_hudson('drc', 'tanzania', lifestage, region0)
		fst_ghana_tanzania[region0] = get_hudson('ghana', 'tanzania', lifestage, region0)
		fst_drc_cambodia[region0] = get_hudson('drc', 'cambodia', lifestage, region0)
		fst_ghana_cambodia[region0] = get_hudson('ghana', 'cambodia', lifestage, region0)
		fst_tanzania_cambodia[region0] = get_hudson('tanzania', 'cambodia', lifestage, region0)

fdf = pd.DataFrame({'fst:drc-ghana':pd.Series(fst_drc_ghana),'fst:drc-tanzania':pd.Series(fst_drc_tanzania), 'fst:drc-cambodia':pd.Series(fst_drc_cambodia), 'fst:ghana-tanzania':pd.Series(fst_ghana_tanzania), 'fst:ghana-cambodia':pd.Series(fst_ghana_cambodia), 'fst:tanzania-cambodia':pd.Series(fst_tanzania_cambodia)})
fdf.to_csv('diversity_stats/' + split + '/pf7/' + lifestage + '_fst.csv')