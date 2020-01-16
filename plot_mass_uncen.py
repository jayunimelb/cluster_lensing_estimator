import numpy as np
import pickle,gzip,glob
import matplotlib.pyplot as plt

parent_dir_sz = 'data_sp/sims/var_ratio/_no_lensing/0.5_tSZ_white_3.0/add_tSZ_1_match_filter_0/arnaud_profile/'
parent_dir_no_sz =  'data_sp/sims/var_ratio/_no_lensing/0.5_no_tSZ_white_3.0/add_tSZ_0_match_filter_0/'
#noise_levels = [1.0,3.0,5.0,7.0,9.0]
mass_levels = [1.0,2.0,3.0,4.0,5.0]
no_sz_uncen = np.zeros(len(mass_levels))
sz_uncen = np.zeros(len(mass_levels))
filt_uncen = np.zeros(len(mass_levels))

def fn_get_avg_kappa(data,inds_inner):
	stkp = data['stacked_kappa_qe']
	#return np.std(stkp[inds_inner])
	data = data['kappa_qe_arr']
	kp_std = np.zeros(len(data))
	for i in range(len(kp_std)):
		kp_std[i] = np.std(data[i][inds_inner])

	return np.mean(kp_std)

clra, cldec = 0., 0.
boxsize,nx,ny,dx,dy = 30, 60 ,60,0.5,0.5

minval, maxval = clra-boxsize/2/60.,  clra+boxsize/2/60.

ra = np.linspace(minval, maxval, nx)
dec = np.linspace(minval, maxval, ny)
RA, DEC = np.meshgrid(ra,dec)
RADIUS = (RA ** 2. + DEC ** 2.)**0.5
RADIUS = RADIUS*60.
inds_inner= np.where(RADIUS<=10.)
y_arr = np.zeros(len(mass_levels))

for i,nn in enumerate(mass_levels):
	
	no_tsz_file =glob.glob('%s/mass_%s/null_test/*'%(parent_dir_no_sz,nn))[0]
	tsz_file =glob.glob('%s/mass_%s/null_test/data*'%(parent_dir_sz ,nn))[0]
	no_sz_data  = pickle.load(gzip.open(no_tsz_file))
	sz_data = pickle.load(gzip.open(tsz_file))
	#from IPython import embed;embed()

	no_sz_sigma = fn_get_avg_kappa(no_sz_data,inds_inner)
	total_sigma = fn_get_avg_kappa(sz_data,inds_inner)
	sz_sigma = np.sqrt(total_sigma**2 - no_sz_sigma**2)
	
	y_arr[i] = (sz_sigma/no_sz_sigma)#**2
	#y_arr[i] = (total_sigma- no_sz_sigma)/no_sz_sigma
plt.plot(np.asarray(mass_levels),y_arr, '*', ls = '--',lw = 3, color = 'b', markersize = 12)

t = ('Experimental noise level 3uK-arcmin')
font = {'weight': 'normal',
'family': 'serif',
'size': 14}
plt.text(1,2, t, fontdict = font)
plt.xlim(0.5,5.5)
plt.ylabel('$\sigma_{sz}/\sigma_{k}$', fontsize  = 24,  fontweight='bold')
plt.xlabel('Mass in ($10^{14} M_{\odot}$)',fontsize  = 14,  fontweight='bold' )
x_arr = np.asarray(mass_levels)
mass_val = np.interp(1.0, y_arr,x_arr)
plt.plot((mass_val, mass_val),(0,1), color = 'k',ls='--', lw = 2.)
plt.plot((0, mass_val),(1,1), color = 'k',ls = '--', lw = 2.)
plt.plot(mass_val, 1, 'ko',markersize = 10)
#plt.axvline(mass_val, ymin = 0., ymax = .25)
#plt.axhline(1, xmin = 0, xmax = mass_val/(5.5))
from IPython import embed;embed()
plt.show()
"""
plt.plot(noise_levels,no_sz_uncen,'*',label = 'ideal case')
plt.plot(noise_levels,sz_uncen,'*', label = 'tSZ in second leg')
plt.plot(noise_levels,filt_uncen,'*', label = 'Template fitting')
legend = plt.legend(loc='upper right', shadow=True, fontsize='small')
legend.get_frame().set_facecolor('#00FFCC')
plt.xlabel("Experimental noise levels (in  uK)")
plt.ylabel('Percentage mass uncertanity')
plt.xlim(0,10)
plt.show()


for i,nn in enumerate(mass_levels)#enumerate(noise_levels):
	no_tsz_file ='%s/0.5_no_tSZ_white_%s/add_tSZ_0_match_filter_0/no_fgs/mass_5.0/results.pkl.gz'%(parent_dir,nn)
	tsz_file =  '%s/0.5_tSZ_white_%s/add_tSZ_1_match_filter_0/arnaud_profile/no_fgs/mass_5.0/results.pkl.gz'%(parent_dir,nn)
	filt_file = '%s/0.5_tSZ_white_%s/add_tSZ_1_match_filter_1_matchedfilterszcorrection/arnaud_profile/no_fgs/mass_5.0/fit_size_5/bigger_box/bounded_params/results.pkl.gz'%(parent_dir,nn)
	no_sz_data  = pickle.load(gzip.open(no_tsz_file))
	sz_data = pickle.load(gzip.open(tsz_file))
	filt_data = pickle.load(gzip.open(filt_file))
	#from IPython import embed;embed()
	no_sz_uncen[i] = no_sz_data['mass_uncen']
	sz_uncen[i] = sz_data['mass_uncen']
	filt_uncen[i] = filt_data['mass_uncen']
"""