import numpy as np
import pickle,gzip,sys
import matplotlib.pyplot as plt
file1 = 'data_sp/sims/des/tsz_fitting/SNR_None_None/0.5_no_tSZ_white_1.0/add_tSZ_0_match_filter_0/Mccarthy_sims/mass_2.0/results.pkl.gz'
file2  = 'data_sp/sims/des/tsz_fitting/SNR_None_None/0.5_tSZ_white_1.0/add_tSZ_1_match_filter_0/Mccarthy_sims/mass_2.0/results.pkl.gz'
#file1 = 'data_sp/sims/des/tsz_fitting/SNR_None_None/0.5_tSZ_white_1.0/add_tSZ_1_match_filter_1_matchedfilterszcorrection/Mccarthy_sims/mass_2.0/fit_size_5/bigger_box/bounded_params/results.pkl.gz'
#file2 = 'data_sp/sims/des/tsz_fitting/SNR_None_None/0.5_tSZ_white_1.0/add_tSZ_1_match_filter_1_matchedfilterszcorrection/Mccarthy_sims/mass_2.0/fit_size_5/bigger_box/bounded_params/amp_guess/results.pkl.gz'
#file3 = 'data_sp/sims/des/tsz_fitting/SNR_None_None/0.5_tSZ_white_1.0/add_tSZ_1_match_filter_1_matchedfilterszcorrection/Mccarthy_sims/mass_2.0/fit_size_5/bigger_box/bounded_params/amp_guess/beam_2_arcmin/results.pkl.gz'
#file3 = 'data_sp/sims/des/tsz_fitting/SNR_None_None/0.5_tSZ_white_1.0/add_tSZ_1_match_filter_1_matchedfilterszcorrection/Mccarthy_sims/mass_2.0/fit_size_5/bigger_box/bounded_params/results.pkl.gz'

file1 = 'data_sp/sims/des/tsz_fitting//SNR_None_None/0.5_tSZ_white_5.0/add_tSZ_1_match_filter_0/arnaud_profile/no_fgs/mass_5.0/results.pkl.gz'
file2 = 'data_sp/sims/des/tsz_fitting//SNR_None_None/0.5_tSZ_white_5.0/add_tSZ_1_match_filter_1_matchedfilterszcorrection/arnaud_profile/no_fgs/mass_5.0/fit_size_5/bigger_box/bounded_params/results.pkl.gz'
data1 = pickle.load(gzip.open(file1))

data2 = pickle.load(gzip.open(file2))
#data3 = pickle.load(gzip.open(file3))
keys1 = data1.keys()
keys2 = data2.keys()
#keys3 = data3.keys()

ax = plt.subplot(121)
for i,key in enumerate(keys1):
	Aarr = data1[key]['Aarr']
	lkhd = data1[key]['Likelihood']
	plt.plot(Aarr/1e14,lkhd,lw = 3, color = 'b')
	if i ==0:
		plt.plot(Aarr/1e14,lkhd,alpha = 0.25, color = 'b', label = 'Combined Likelihood')#label ='Simulations')
	else:
		plt.plot(Aarr/1e14,lkhd,alpha = 0.25, color = 'b')
	if i ==0:
		L1 = lkhd
	else:
		L1 = L1*lkhd
plt.xlim(3,7)
plt.xlabel('Mass$(10^{14}$$ M_{\odot})$', fontsize = 14,fontweight = 'bold');#plt.plot(Aarr/1e14,L1/L1.max(),lw = 3, color = 'k', label = 'Combined likelihood')
plt.ylabel('Normalised Likelihood', fontsize = 14,fontweight = 'bold')
plt.axvline(5.0, color = 'k', lw = 2, ls ='--', label = 'Input mass');plt.legend(loc='upper right', shadow=True, fontsize='small')
plt.title('tSZ induced variance')

ax = plt.subplot(122)
for i,key in enumerate(keys2):
	Aarr = data2[key]['Aarr']
	lkhd = data2[key]['Likelihood']
	plt.plot(Aarr/1e14,lkhd,lw = 3, color = 'b')
	if i ==0:
		plt.plot(Aarr/1e14,lkhd,alpha = 0.25, color = 'b',  label = 'Combined Likelihood')#label ='Simulations')
	else:
		plt.plot(Aarr/1e14,lkhd,alpha = 0.25, color = 'b')
	if i ==0:
		L2 = lkhd
	else:
		L2 = L2*lkhd

plt.xlim(3,7)
plt.xlabel('Mass$(10^{14}$$ M_{\odot})$', fontsize = 14,fontweight = 'bold');#plt.plot(Aarr/1e14,L2/L2.max(),lw = 3, color = 'k', label = 'Combined likelihood')
plt.ylabel('Normalised Likelihood', fontsize = 14,fontweight = 'bold')
plt.axvline(5.0, color = 'k', lw = 2, ls ='--', label = 'Input mass');plt.legend(loc='upper right', shadow=True, fontsize='small')
plt.title('template fitting')
plt.text(60, .025, r'$\mu=100,\ \sigma=15$')
plt.show()
sys.exit()
"""
ax = plt.subplot(313)
for i,key in enumerate(keys3):
	Aarr = data3[key]['Aarr']
	lkhd = data3[key]['Likelihood']
	plt.plot(Aarr/1e14,lkhd,alpha = 0.25, color = 'b')
	if i ==0:
		plt.plot(Aarr/1e14,lkhd,alpha = 0.25, color = 'b', label ='Simulations')
	else:
		plt.plot(Aarr/1e14,lkhd,alpha = 0.25, color = 'b')
	if i ==0:
		L3 = lkhd
	else:
		L3 = L3*lkhd

plt.xlim(0,4)
plt.xlabel('Mass$(10^{14}$$ M_{\odot})$', fontsize = 14,fontweight = 'bold');plt.plot(Aarr/1e14,L3/L3.max(),lw = 3, color = 'k', label = 'Combined likelihood with fitting')
plt.ylabel('Normalised Likelihood', fontsize = 14,fontweight = 'bold')
plt.axvline(2.0, color = 'k', lw = 3, ls ='--', label = 'Input mass');plt.legend(loc='upper right', shadow=True, fontsize='small')

from IPython import embed;embed()
sys.exit()

plt.show()
"""