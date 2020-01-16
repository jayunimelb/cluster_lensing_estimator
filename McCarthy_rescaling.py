import numpy as np
import pickle,gzip
data_file = 'data/tsz_y_cutouts_M200m_v01.pk'

data = pickle.load(open(data_file))

M500 = data['M500c']
ycts = data['y_cutouts']

mass_arr = np.arange(2,6)*1e14

for mm in mass_arr:
	results = {}

	output = np.zeros(ycts.shape)
	
	for i,mval in enumerate(M500):
	
		output [:,:,i]= ycts[:,:,i]*((mm/mval)**(5/3.))
	results['M500c'] = M500		
	results['scaled_y_cutouts'] = output
	pickle.dump(results,open('data/tsz_y_cutouts_M500c_v01_%s_1e14.pk'%(mm/1e14),'w'))
