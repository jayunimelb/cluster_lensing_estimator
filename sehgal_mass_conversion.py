

import numpy as np
import pickle,gzip
 

 # note frequency correction!!
sehgal_90_sz_file = 'data/sehgal_2009_sims/tsz_cib_radio_ksz_extracts/tSZ_extracts_090_M200min1.3_zmin0.25_boxsize100.0am_dx0.5am.pkl.gz'
sehgal_148_sz_file  = 'data/sehgal_2009_sims/tsz_cib_radio_ksz_extracts/tSZ_extracts_148_M200min1.3_zmin0.25_boxsize100.0am_dx0.5am.pkl.gz'

sehgal_90_data = pickle.load(gzip.open(sehgal_90_sz_file))
sehgal_148_data = pickle.load(gzip.open(sehgal_148_sz_file))

keynames = sehgal_90_data.keys()

mass_vals = (np.asarray(keynames)[:,3])*1e14

# convert the mass values to M500 rho_critical from M200 mean



mass_arr = (np.arange(2,6))*1e14

for mm in mass_arr:
	result_dic_90 = {}
	result_dic_148 = {}
	result_dic_file_90 = 'data/sehgal_2009_sims/tsz_cib_radio_ksz_extracts/tSZ_extracts_090_M200min1.3_zmin0.25_boxsize100.0am_dx0.5am_mass_.pkl.gz'
	for i,mval in enumerate(mass_vals):
		factor = (mm/mval)**(5/3.)
		result_dic_90[keynames[i]] = sehgal_90_data[keynames[i]]*factor
		result_dic_148[keynames[i]] = sehgal_148_data[keynames[i]]*factor
	pickle.dump(result_dic_90,gzip.open('','w'))
	pickle.dump(result_dic_148,gzip.open('','w'))




