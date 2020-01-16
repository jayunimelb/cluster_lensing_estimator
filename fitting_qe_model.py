"""



"""

import numpy as np
import modules, sys, glob, os, pickle,gzip
scl_cmb = modules.scl_cmb

sims = scl_cmb.simulations()

# get the covariance matrix using JK on data 
def fn_select_clusters(kappa_qe_arr, CLUS_IDENTIFIER, rich1, rich2, z1 = None, z2 = None):

	passed_inds = []
	for kcnt, cluskey in enumerate( CLUS_IDENTIFIER ):
		ra, dec, z_val, rich, weight = cluskey

		passed = 0
		if rich >= rich1 and rich<rich2:
			passed = 1

		if z1<>None:
			if z_val>=z1: 
				passed = 1
			else: 
				passed = 0

		if z2>None:
			if z_val<z1: 
				passed = 1
			else: 
				passed = 0

		if passed: passed_inds.append(kcnt)

	if len(passed_inds) == 0: return None, None

	passed_inds = np.asarray( passed_inds )

	return kappa_qe_arr[passed_inds], CLUS_IDENTIFIER[passed_inds]


ipfolder, noofsims = sys.argv[1], int(sys.argv[2])

try:
	minrich, maxrich = sys.argv[3], sys.argv[4]
except:
	minrich, maxrich = 20, 101

fnamearr = glob.glob('%s/st*'%(ipfolder))
cov_file = '%s/kappa_COV_%s_JK.pkl.gz' %(ipfolder,noofsims)#, fname.split('/')[-1])

# make sure all the data files have same number of clusters
kappa_qe_arr = None
if not os.path.exists(cov_file):
	print '\n\n getting kappa_COV using a JK approch\n'
	for fname in fnamearr:
		kappadic = pickle.load(gzip.open(fname,'rb'))
		if not 'kappa_qe_arr' in kappadic:
			continue
		kappa_qe_arr = np.asarray( kappadic['kappa_qe_arr'] )
		CLUS_IDENTIFIER = np.asarray( kappadic['CLUS_IDENTIFIER'] )
		break
	if sims.is_seq(kappa_qe_arr):
		totalclus = len(kappa_qe_arr)
		kappa_qe_arr, CLUS_IDENTIFIER = fn_select_clusters(kappa_qe_arr, CLUS_IDENTIFIER, minrich, maxrich)

		kappa_COV = sims.fn_get_COV_from_JK(kappa_qe_arr,CLUS_IDENTIFIER,noofsims,[ (minrich, maxrich) ] )
		pickle.dump(kappa_COV, gzip.open(cov_file,'wb'), protocol = 2)
	else:
		print 'No kappa_qe_arr for JK cov'


## Generting models

param_dict = kappadic['param_dict']
boxsize = param_dict['boxsize']
dx,dy = 0.5,0.5

nx,ny = int(boxsize/dx), int(boxsize/dy)
mapparams = nx,ny,dx,dy
clra,cldec = 0,0
boxsize = 100
minval, maxval = clra-boxsize/2/60.,  clra+boxsize/2/60.
ra = np.linspace(minval, maxval, nx)
dec = np.linspace(minval, maxval, ny)
RA, DEC = np.meshgrid(ra,dec)
Bl  = sims.fn_beam_stuff(mapparams, use_beam = 1., exp_beam =1., return_beam = 1)
if param_dict['add_TF'] == 0.:
	TWODTF = sims.fn_get_HPF(mapparams, ideal = 1)[0]
minA, maxA, delA = 0, 8.0, 0.05
Aarr = np.arange(minA,maxA+delA,delA) * 1e14
opfolder = '%s/model'
kappa_model_file = '%s/stacked_kappa_model_minrich%s_maxrich%s.pkl.gz' %(opfolder, minrich, maxrich)%(ipfolder)
if not os.path.exists(opfolder): os.system('mkdir %s' %(opfolder))
if os.path.exists(kappa_model_file):
	kappa_model_dic = pickle.load(gzip.open(kappa_model_file, 'rb'))
else:
	kappa_model_dic = {}

for acnt, AA in enumerate(Aarr):
	keyname = AA/1e14
	if keyname not in kappa_model_dic: 
		rad_profile_info = [binsize, minbin, maxbin]
		KAPPA_MODELS_rad_profiles = sims.fn_get_STACKED_KAPPA_MODEL_FN_MASS(AA, alpha_fit, RA, DEC, CLUS_IDENTIFIER,param_dict, mapparams, Bl = sims.Bl, TWODTF = TWODTF, kappa_model_fixed = 0,rad_profile_info = rad_profile_info, use_TF_beam = 1, perform_offsets_correction = 0)
		kappa_model_dic[keyname] = KAPPA_MODELS_rad_profiles

pickle.dump(kappa_model_dic, gzip.open(kappa_model_file, 'wb'), protocol = 2)




from IPython import embed;embed()
sys.exit()

