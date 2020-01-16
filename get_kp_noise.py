import numpy as np
import pickle,gzip
import sys
from pylab import *
inner_radius = float(sys.argv[1])
noises = [1.0,2.0,3.0,5.0]
ideal_std,tsz_std,fltr_std = np.zeros(len(noises)),np.zeros(len(noises)),np.zeros(len(noises))
parent_folder = 'kpa_var'
for i,ns in enumerate(noises): 
	ideal_file = '%s/noise_level_%s_mass_5.0/ideal.pkl.gz'%(parent_folder, ns)
	tsz_file = '%s/noise_level_%s_mass_5.0/tsz.pkl.gz'%(parent_folder, ns)
	filt_file = '%s/noise_level_%s_mass_5.0/tsz_filter.pkl.gz'%(parent_folder, ns)	
	ideal_data = pickle.load(gzip.open(ideal_file))
	tsz_data = pickle.load(gzip.open(tsz_file))
	filt_data =pickle.load(gzip.open(filt_file))
	RADIUS = ideal_data['radius']
	ideal_kappa = ideal_data['KAPPA_QE_ARR'][0]
	tsz_kappa = tsz_data['KAPPA_QE_ARR'][0]
	filt_kappa = filt_data['KAPPA_QE_ARR'][0]
	inds_inner= np.where(RADIUS<=inner_radius)
	ideal_std[i]= np.std(ideal_kappa[inds_inner])
	tsz_std[i] = np.std(tsz_kappa[inds_inner])
	fltr_std[i] = np.std(filt_kappa[inds_inner])

plot(noises,ideal_std,'r')
plot(noises,tsz_std,'r',ls = '--')
plot(noises,fltr_std,'r',ls = '-.')
show()
#from IPython import embed;embed()	

