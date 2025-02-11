import warnings
warnings.filterwarnings('ignore')
import pandas as pd 
import sys
import os

os.chdir('/seq/plasmodium/Pfal_life_cycle_sel/dfe_files')
# read in pre-computed divergence file for chosen reference species. Divergence file already includes stage labels.
#div = pd.read_csv('divergence_by_gene_praefal_m3_final.csv')
div = pd.read_csv("divergence_by_gene_preich_m3_final.csv")
# read in life stage and excluded gene for jackknife iteration
lifestage = sys.argv[1].strip()
with open('jackknife_exclude/'+lifestage+'.txt') as f: 
	for line in f:
		excl = line.strip()
# filter divergence data by life stage and exclude gene
totals = div[(div['stage']==lifestage) & (div['Gene ID'] != excl)]
# neutral sites SYN or FFD 
neut_key = sys.argv[2]
if(neut_key=='SYN'): totals = totals.agg({'ns':'sum', 's':'sum', 'ns_sites':'sum', 's_sites':'sum'})
else: totals = totals.agg({'ns':'sum', 'ffd':'sum', 'ns_sites':'sum', 'ffd_sites':'sum'})

# write to divergence files in overwrite mode
fn = 'divergence_files/' + lifestage + '.txt'
with open(fn, 'w') as file:
	file.write('1 ' + str(totals['ns_sites']) + ' ' + str(totals['ns']) + '\n')
	if(neut_key=='SYN'): file.write('0 ' + str(totals['s_sites']) + ' ' + str(totals['s']))
	else: file.write('0 ' + str(totals['ffd_sites']) + ' ' + str(totals['ffd']))
print(fn)
