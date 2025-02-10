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

#iterate through all genes in bed file to get start and stop positions
p_s = {}
p_ns = {}
p_ns_p_s = {}
thetas_ns = {}
thetas_s = {}
td_ns = {}
td_s = {}
td_all = {}
notfound = []
props = pd.read_csv('props_adj_breadth.txt', sep = '\t')
with open(set_base+'_bed/'+lifestage+'.bed') as f:
	for line in f:
		bdat = line.strip().split('\t')
		chrom = bdat[0]
		p0 = int(bdat[1])
		p1 = int(bdat[2])
		gene = bdat[5]
		region0=chrom+':'+bdat[1]+'-'+bdat[2]
		print(region0)
		entry = props['COORD'].str.contains(region0)
		if not np.sum(entry): continue
		# pick up relevant transcript, defaulting to 1st (GENE.1) if coordinates are the same for both 
		props1 = props[entry].sort_values(by='TOTAL_CODING_LENGTH')
		ns = props1['NS'].values[0]
		missing_data_prop = props1['value'].values[0]
		ns *= missing_data_prop
		print('ns: ')
		print(ns)
		s = props1['SYN'].values[0]
		s *= missing_data_prop
		callset = allel.read_vcf('pf7_vcf/' + pop + '/'  + lifestage + '.CoreOnly.NS_SNPsOnly.recode.vcf.gz', fields=['variants/POS','calldata/GT', 'calldata/DP'], region = region0)
		if callset is None:
			print('Not found')
			notfound.append(region0)
			pi_ns = 0
		else:
			gt = allel.GenotypeArray(callset['calldata/GT'])
			ac = gt.count_alleles()
			pos = callset['variants/POS']
			#pi_ns = allel.sequence_diversity(pos, ac)
			pi_ns = np.sum(allel.mean_pairwise_difference(ac))/ns
			theta_ns = allel.watterson_theta(pos,ac)
			tajimas_d_ns = allel.tajima_d(ac, pos=pos)
			p_ns[region0] = pi_ns
			thetas_ns[region0] = theta_ns
			td_ns[region0] = tajimas_d_ns
		callset_s = allel.read_vcf('pf7_vcf/' + pop + '/' + lifestage + '.CoreOnly.S_SNPsOnly.recode.vcf.gz', fields=['variants/POS','calldata/GT','calldata/DP'], region = region0)
		if callset_s is None:
			print('Not found')
			notfound.append(region0)
			pi_s = 0
		else:
			gt = allel.GenotypeArray(callset_s['calldata/GT'])
			ac = gt.count_alleles()
			pos = callset_s['variants/POS']
			#pi_s = allel.sequence_diversity(pos, ac)
			pi_s = np.sum(allel.mean_pairwise_difference(ac))/s
			theta_s = allel.watterson_theta(pos,ac)
			tajimas_d_s = allel.tajima_d(ac, pos=pos)
#		if ((pi_s > 0) | (pi_ns > 0)):
			p_s[region0] = pi_s
			p_ns_p_s[region0] = pi_ns/pi_s
			thetas_s[region0] = theta_s
			td_s[region0] = tajimas_d_s

		callset = allel.read_vcf('pf7_vcf/' + pop + '/' + lifestage + '.CoreOnly.SNPsOnly.recode.vcf.gz', fields=['variants/POS','calldata/GT', 'calldata/DP'], region = region0)
		if callset is None:
			print('Not found')
			notfound.append(region0)
		else: 
			gt = allel.GenotypeArray(callset['calldata/GT'])
			pos = callset['variants/POS']
			ac = gt.count_alleles()
			tajimas_d=allel.tajima_d(ac, pos=pos)
			td_all[region0]=tajimas_d
 
psdf = pd.DataFrame.from_dict(p_s, orient="index")
psdf.columns = ['pi_s']
psdf.to_csv('diversity_stats/'+set_base+'/pf7/'+pop+'/'+lifestage+'_pi_s.csv')

pnsdf = pd.DataFrame.from_dict(p_ns, orient="index")
pnsdf.columns = ['pi_ns']
pnsdf.to_csv('diversity_stats/'+set_base+'/pf7/'+pop+'/'+lifestage+'_pi_ns.csv')

pnspsdf = pd.DataFrame.from_dict(p_ns_p_s, orient="index")
pnspsdf.columns = ['pi_ns_div_pi_s']
pnspsdf.to_csv('diversity_stats/'+set_base+'/pf7/'+pop+'/'+lifestage+'_pi_ns_div_pi_s.csv')

tdnsdf = pd.DataFrame.from_dict(td_ns, orient="index")
tdnsdf.columns = ['tajimas_d_ns']
tdnsdf.to_csv('diversity_stats/'+set_base+'/pf7/'+pop+'/'+lifestage+'_tajimas_d_ns.csv')

tdsdf = pd.DataFrame.from_dict(td_all, orient="index")
tdsdf.columns = ['tajimas_d_all']
tdsdf.to_csv('diversity_stats/'+set_base+'/pf7/'+pop+'/'+lifestage+'_tajimas_d.csv')

tnsdf = pd.DataFrame.from_dict(thetas_ns, orient="index")
tnsdf.columns = ['theta_ns']
tnsdf.to_csv('diversity_stats/'+set_base+'/pf7/'+pop+'/'+lifestage+'_theta_ns.csv')

tsdf = pd.DataFrame.from_dict(thetas_s, orient="index")
tsdf.columns = ['theta_s']
tsdf.to_csv('diversity_stats/'+set_base+'/pf7/'+pop+'/'+lifestage+'_theta_s.csv')

nfdf = pd.DataFrame(notfound)
nfdf.to_csv('diversity_stats/'+set_base+'/pf7/'+pop+'/'+lifestage+'_notfound.csv')
