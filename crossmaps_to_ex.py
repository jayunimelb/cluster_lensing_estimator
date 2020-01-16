
import numpy as np, os, sys

produce_maps = int(sys.argv[1])

'''
noiseval = 6.0 #sptpol 150ghz maps
rsarr = [739494, 570842, 381439]
comps = ['nocross_map', '2am1.2ambeam', '5am1.2ambeam']
paramsfiles = ['params_planck_r_0.0_2015_cosmo_lensed_LSS_nocrossmap.txt', 'params_planck_r_0.0_2015_cosmo_lensed_LSS_crossmap_2am1.2ambeam.txt', 'params_planck_r_0.0_2015_cosmo_lensed_LSS_crossmap_5am1.2ambeam.txt']
tszarr = [0, 0, 0]

#comps = ['nocross_map_gradcut_1500_with_tSZ', 'nocross_map_gradcut_1000_with_tSZ']
#paramsfiles = ['params_planck_r_0.0_2015_cosmo_lensed_LSS_nocrossmap_gradcut_1500.txt', 'params_planck_r_0.0_2015_cosmo_lensed_LSS_nocrossmap_gradcut_1000.txt']
#tszarr = [1,1]

#comps = ['nocross_map', 'nocross_map_gradcut_1500', 'nocross_map_gradcut_1000', '2am1.2ambeam', '5am1.2ambeam']
#paramsfiles = ['params_planck_r_0.0_2015_cosmo_lensed_LSS_nocrossmap.txt', 'params_planck_r_0.0_2015_cosmo_lensed_LSS_nocrossmap_gradcut_1500.txt', 'params_planck_r_0.0_2015_cosmo_lensed_LSS_nocrossmap_gradcut_1000.txt', 'params_planck_r_0.0_2015_cosmo_lensed_LSS_crossmap_2am1.2ambeam.txt', 'params_planck_r_0.0_2015_cosmo_lensed_LSS_crossmap_5am1.2ambeam.txt']

#comps = ['nocross_map', 'nocross_map_gradcut_1500_with_tSZ', 'nocross_map_gradcut_1000_with_tSZ', '2am1.2ambeam', '5am1.2ambeam']
#paramsfiles = ['params_planck_r_0.0_2015_cosmo_lensed_LSS_nocrossmap.txt', 'params_planck_r_0.0_2015_cosmo_lensed_LSS_nocrossmap_gradcut_1500.txt', 'params_planck_r_0.0_2015_cosmo_lensed_LSS_nocrossmap_gradcut_1000.txt', 'params_planck_r_0.0_2015_cosmo_lensed_LSS_crossmap_2am1.2ambeam.txt', 'params_planck_r_0.0_2015_cosmo_lensed_LSS_crossmap_5am1.2ambeam.txt']

#comps = ['nocross_map_gradcut_1500_with_tSZ_with_offsets']
#paramsfiles = ['params_planck_r_0.0_2015_cosmo_lensed_LSS_nocrossmap_gradcut_1500_with_ofssets.txt']
#tszarr = [1]
"""
comps = ['nocross_map', '2am1.2ambeam', '5am1.2ambeam', '5am1.2ambeam_low_noise']
#paramsfiles = ['params_planck_r_0.0_2015_cosmo_lensed_LSS_nocrossmap.txt', 'params_planck_r_0.0_2015_cosmo_lensed_LSS_crossmap_5am1.2ambeam.txt', 'params_planck_r_0.0_2015_cosmo_lensed_LSS_crossmap_5am1.2ambeam_low_noise.txt']
tszarr = [0]
"""

comps = ['nocross_map', '2am1.2ambeam', '5am1.2ambeam','5am1.2ambeam_high_noise']
paramsfiles = ['params_planck_r_0.0_2015_cosmo_lensed_LSS_nocrossmap.txt', 'params_planck_r_0.0_2015_cosmo_lensed_LSS_crossmap_2am1.2ambeam.txt', 'params_planck_r_0.0_2015_cosmo_lensed_LSS_crossmap_5am1.2ambeam.txt','params_planck_r_0.0_2015_cosmo_lensed_LSS_crossmap_5am1.2ambeam_high_noise.txt']
tszarr = [0,0, 0, 0]

labels_beam = ['1.2^{\\prime}, 1.2^{\\prime}', '2^{\\prime},1.2^{\\prime}', '5^{\\prime},1.2^{\\prime}', '5^{\\prime},1.2^{\\prime}']
labels_noise = ['6, 6', '24, 6', '55, 6', '100, 6']

"""
comps = ['nocross_map_gradcut_1500_with_tSZ', 'nocross_map_gradcut_1000_with_tSZ']
paramsfiles = ['params_planck_r_0.0_2015_cosmo_lensed_LSS_nocrossmap_gradcut_1500.txt', 'params_planck_r_0.0_2015_cosmo_lensed_LSS_nocrossmap_gradcut_1000.txt']
tszarr = [1,1]
"""

comps = ['2am1.2ambeam_with_tSZ_cleaning']
paramsfiles = ['params_planck_r_0.0_2015_cosmo_lensed_LSS_crossmap_2am1.2ambeam_with_tSZ_cleaning.txt']
tszarr = [1]

#comps = ['nocross_map', '2am1.2ambeam', '5am1.2ambeam','5am1.2ambeam_high_noise', 'nocross_map_gradcut_1500_with_tSZ', '2am1.2ambeam_with_tSZ_cleaning']
#paramsfiles = ['params_planck_r_0.0_2015_cosmo_lensed_LSS_nocrossmap.txt', 'params_planck_r_0.0_2015_cosmo_lensed_LSS_crossmap_2am1.2ambeam.txt', 'params_planck_r_0.0_2015_cosmo_lensed_LSS_crossmap_5am1.2ambeam.txt','params_planck_r_0.0_2015_cosmo_lensed_LSS_crossmap_5am1.2ambeam_high_noise.txt','params_planck_r_0.0_2015_cosmo_lensed_LSS_nocrossmap_gradcut_1500.txt', 'params_planck_r_0.0_2015_cosmo_lensed_LSS_crossmap_2am1.2ambeam_with_tSZ_cleaning.txt']
#tszarr = [0,0, 0, 0, 1, 1]
labels_beam = ['1.2^{\\prime}, 1.2^{\\prime}', '2^{\\prime},1.2^{\\prime}', '5^{\\prime},1.2^{\\prime}', '5^{\\prime},1.2^{\\prime}', '1.2^{\\prime}', '2^{\\prime},1.2^{\\prime}']
labels_noise = ['6, 6', '24, 6', '55, 6', '100, 6', '6, 6', '24, 6']

comps = ['nocross_map', '2am1.2ambeam', '5am1.2ambeam', '2am1.2ambeam_with_tSZ_cleaning', 'nocross_map_gradcut_1500_with_tSZ']
paramsfiles = ['params_planck_r_0.0_2015_cosmo_lensed_LSS_nocrossmap.txt', 'params_planck_r_0.0_2015_cosmo_lensed_LSS_crossmap_2am1.2ambeam.txt', 'params_planck_r_0.0_2015_cosmo_lensed_LSS_crossmap_5am1.2ambeam.txt','params_planck_r_0.0_2015_cosmo_lensed_LSS_crossmap_2am1.2ambeam_with_tSZ_cleaning.txt', 'params_planck_r_0.0_2015_cosmo_lensed_LSS_nocrossmap_gradcut_1500.txt']
tszarr = [0,0, 0, 1]
labels_beam = ['1.2^{\\prime}, 1.2^{\\prime}', '2^{\\prime},1.2^{\\prime}', '5^{\\prime},1.2^{\\prime}', '2^{\\prime},1.2^{\\prime}', '1.2^{\\prime}']
labels_noise = ['6, 6', '24, 6', '55, 6', '24, 6', '6, 6']

comps = ['sptpol_map_like_sims/sptpollike_with_tSZ_cleaning']
paramsfiles = ['sptpol_map_like_sims/params_planck_r_0.0_2015_cosmo_lensed_LSS_crossmap_sptpol_with_FG_tSZ_cleaning.txt']
labels_beam = ['2^{\\prime},1.2^{\\prime}']

comps = ['sptpol_map_like_sims/sptpollike_with_tSZ_cleaning_v2']
paramsfiles = ['sptpol_map_like_sims/params_planck_r_0.0_2015_cosmo_lensed_LSS_crossmap_sptpol_with_FG_tSZ_cleaning.txt']
labels_beam = ['2^{\\prime},1.2^{\\prime}']

#comps = ['sptpol_map_like_sims/sptpollike_with_tSZ_cleaning_v2', 'sptpol_map_like_sims/sptpollike_with_tSZ_cleaning_v3']
#paramsfiles = ['sptpol_map_like_sims/params_planck_r_0.0_2015_cosmo_lensed_LSS_crossmap_sptpol_with_FG_tSZ_cleaning.txt', 'sptpol_map_like_sims/params_planck_r_0.0_2015_cosmo_lensed_LSS_crossmap_sptpol_FG_sims.txt']
#labels_beam = ['sptpol', 'sptpol']

tszarr = [1, 1]
noiseval = 'sptpol'
labels_noise = ['sptpol', 'sptpol']
noofsims = 2
totalclus = 1900

which_exp = 'sptpol_des'

comps = ['sptpol_map_like_sims/sptpollike_with_tSZ_cleaning_v4_rich_30_40']
paramsfiles = ['sptpol_map_like_sims/params_planck_r_0.0_2015_cosmo_lensed_LSS_crossmap_sptpol_FG_sims.txt']

comps = ['sptpol_map_like_sims/sptpollike_with_tSZ_cleaning_v5_rich_30_40']
paramsfiles = ['sptpol_map_like_sims/params_planck_r_0.0_2015_cosmo_lensed_LSS_crossmap_sptpol_FG_sims.txt']

comps = ['sptpol_map_like_sims/sptpollike_with_tSZ_cleaning_v6_rich_30_40']
paramsfiles = ['sptpol_map_like_sims/params_planck_r_0.0_2015_cosmo_lensed_LSS_crossmap_sptpol_FG_sims.txt']

comps = ['sptpol_map_like_sims/sptpollike_with_tSZ_cleaning_v7_rich_30_40']
paramsfiles = ['sptpol_map_like_sims/params_planck_r_0.0_2015_cosmo_lensed_LSS_crossmap_sptpol_FG_sims.txt']

#comps = ['sptpol_map_like_sims/sptpollike_with_tSZ_cleaning_v4_rich_30_40_no_mask']
#paramsfiles = ['sptpol_map_like_sims/params_planck_r_0.0_2015_cosmo_lensed_LSS_crossmap_sptpol_FG_sims.txt']

noiseval = 'sptpol';labels_beam = ['2.0, 1.2'];labels_noise = ['sptpol']
tszarr = [1]
noofsims = 2
totalclus = 500
which_exp = 'sptpol_des_rich'

"""
comps = ['2am1.2ambeam', 'sptpol_map_like_sims/sptpollike_with_tSZ_cleaning_v2']
paramsfiles = ['params_planck_r_0.0_2015_cosmo_lensed_LSS_nocrossmap.txt', 'sptpol_map_like_sims/params_planck_r_0.0_2015_cosmo_lensed_LSS_crossmap_sptpol_with_FG_tSZ_cleaning.txt']
tszarr = [0,1]
totalclus = 1900
labels_beam = ['2^{\\prime},1.2^{\\prime}', '2^{\\prime},1.2^{\\prime}',]
labels_noise = ['24, 6', 'sptpol']
which_exp = 'sptpol_des'
"""


#### for SPT-SZ EBAX

comps = ['offsets_v2/withtsz_plus_mask']
paramsfiles = ['sptpol_map_like_sims/params_planck_r_0.0_2015_cosmo_lensed_LSS_crossmap_sptsz_FG_sims.txt']
tszarr = [1]

"""
comps = ['offsets_v2/notsz_plus_mask']
paramsfiles = ['sptpol_map_like_sims/params_planck_r_0.0_2015_cosmo_lensed_LSS_crossmap_sptsz_FG_sims_notSZ.txt']
tszarr = [0]

comps = ['offsets/notsz_nomask']
paramsfiles = ['sptpol_map_like_sims/params_planck_r_0.0_2015_cosmo_lensed_LSS_crossmap_sptsz_FG_sims_notSZ_nomask.txt']
tszarr = [0]
"""
#noofsims = 2

noiseval = 18.;labels_beam = ['2.0, 1.2'];labels_noise = ['sptpol']
totalclus = 3500
which_exp = 'sptsz_des'
import glob
#rsarr = glob.glob('data/sims/sptsz_des/crossmaps/offsets/withtsz_plus_mask/SNR_None_None/0.5_tSZ_white_18.0/stacked*')
#rsarr = glob.glob('data/sims/sptsz_des/crossmaps/%s/SNR_None_None/0.5_*/stacked*' %(comps[0]))
#rsarr = map(lambda x: int( x.split('/')[-1].split('_')[-1].replace('.pkl.gz','') ), rsarr)
rsarr = [92633, 145598, 397327, 351012, 911917, 46576, 21648, 739494, 345804, 805436,  966089, 570842, 208142, 737186, 94932, 391048, 76613, 360918, 675753, 524681, 508428, 915807, 244716, 730139, 103136]

#rsarr = [76613, 94932, 103136, 208142, 244716, 360918, 391048, 508428, 524681, 570842, 675753, 730139, 737186, 915807]
noofsims = len(rsarr)
#### for SPT-SZ EBAX
'''


##########################################################################################
##########################################################################################
#20180111 - SPTpol x DES RM year 3 SNR check

noiseval = 'sptpol'# 150ghz maps
rsarr = [739494, 570842, 381439, 253783, 18472, 76613, 94932, 103136, 208142, 244716, 360918, 391048, 508428, 524681,
5708, 6753, 139, 7186, 91537, 3649, 7988, 204503, 471022, 5361, 562017]#
#rsarr=[570842, 675753, 730139, 737186, 915807]

#comps = ['sptpol_map_like_sims/year3_lgt_5_sims']
#paramsfiles = ['sptpol_map_like_sims/params_planck_r_0.0_2015_cosmo_lensed_LSS_crossmap_sptpol_FG_sims_year3_lgt_5_sims.txt']
#tszarr = [1]

comps = ['sptpol_map_like_sims/DES_RM_like_sample/year3_lgt_5']
paramsfiles = ['sptpol_map_like_sims/DES_RM_sims_testing/params_planck_r_0.0_2015_cosmo_lensed_LSS_crossmap_sptpol_FG_notSZ_sims_year3_lgt_5_sims.txt']
tszarr = [0]

comps = ['sptpol_map_like_sims/DES_RM_like_sample/year3_lgt_5_with_FG']
paramsfiles = ['sptpol_map_like_sims/DES_RM_sims_testing/params_planck_r_0.0_2015_cosmo_lensed_LSS_crossmap_sptpol_FG_notSZ_sims_year3_lgt_5_sims_with_FG.txt']
tszarr = [0]


"""
#20180122_1845 for Sanjay's SNR
comps = ['sptpol_map_like_sims/DES_RM_like_sample/year3_lgt_5_Gaussbeam']
paramsfiles = ['sptpol_map_like_sims/DES_RM_sims_testing/params_planck_r_0.0_2015_cosmo_lensed_LSS_crossmap_sptpol_FG_notSZ_sims_year3_lgt_5_sims_Gaussbeam.txt']
tszarr = [0]
"""
'''
comps = ['sptpol_map_like_sims/DES_RM_like_sample/normal']
paramsfiles = ['sptpol_map_like_sims/DES_RM_sims_testing/params_planck_r_0.0_2015_cosmo_lensed_LSS_crossmap_sptpol_FG_notSZ_sims.txt']
tszarr = [0]
'''

"""
#20180122_1705
comps = ['sptpol_map_like_sims/DES_RM_like_sample/normal_Gaussbeam']
paramsfiles = ['sptpol_map_like_sims/DES_RM_sims_testing/params_planck_r_0.0_2015_cosmo_lensed_LSS_crossmap_sptpol_FG_notSZ_sims_Gaussbeam.txt']
tszarr = [0]
"""
"""
comps = ['sptpol_map_like_sims/DES_RM_like_sample/normal_nomask']
paramsfiles = ['sptpol_map_like_sims/DES_RM_sims_testing/params_planck_r_0.0_2015_cosmo_lensed_LSS_crossmap_sptpol_FG_notSZ_sims_nomask.txt']
tszarr = [0]

comps = ['sptpol_map_like_sims/DES_RM_like_sample/normal_nomask_Gaussbeam']
paramsfiles = ['sptpol_map_like_sims/DES_RM_sims_testing/params_planck_r_0.0_2015_cosmo_lensed_LSS_crossmap_sptpol_FG_notSZ_sims_nomask_Gaussbeam.txt']
tszarr = [0]
"""

#comps = ['sptpol_map_like_sims/DES_RM_like_sample/normal_nomask_samenoise_samebeam']
#paramsfiles = ['sptpol_map_like_sims/DES_RM_sims_testing/params_planck_r_0.0_2015_cosmo_lensed_LSS_crossmap_sptpol_FG_notSZ_sims_nomask_samenoise_samebeam.txt']
#tszarr = [0]

'''
comps = ['sptpol_map_like_sims/DES_RM_like_sample/1.2am_beam']
paramsfiles = ['sptpol_map_like_sims/DES_RM_sims_testing/params_planck_r_0.0_2015_cosmo_lensed_LSS_crossmap_sptpol_FG_notSZ_sims_1.2am_beam.txt']
tszarr = [0]

comps = ['sptpol_map_like_sims/DES_RM_like_sample/no_filtering']
paramsfiles = ['sptpol_map_like_sims/DES_RM_sims_testing/params_planck_r_0.0_2015_cosmo_lensed_LSS_crossmap_sptpol_FG_notSZ_sims_no_filtering.txt']
tszarr = [0]

comps = ['sptpol_map_like_sims/DES_RM_like_sample/no_filtering_1.2am_beam']
paramsfiles = ['sptpol_map_like_sims/DES_RM_sims_testing/params_planck_r_0.0_2015_cosmo_lensed_LSS_crossmap_sptpol_FG_notSZ_sims_no_filtering_1.2am_beam.txt']
tszarr = [0]

comps = ['sptpol_map_like_sims/DES_RM_like_sample/no_filtering_1.2am_beam_no_crossmaps']
paramsfiles = ['sptpol_map_like_sims/DES_RM_sims_testing/params_planck_r_0.0_2015_cosmo_lensed_LSS_NOcrossmap_sptpol_FG_notSZ_sims_no_filtering_1.2am_beam.txt']
tszarr = [0]

comps = ['sptpol_map_like_sims/DES_RM_like_sample/no_filtering_1.2am_beam_no_crossmaps_withmask']
paramsfiles = ['sptpol_map_like_sims/DES_RM_sims_testing/params_planck_r_0.0_2015_cosmo_lensed_LSS_NOcrossmap_sptpol_FG_notSZ_sims_no_filtering_1.2am_beam_withmask.txt']
tszarr = [0]
'''

totalclus = 3500
which_exp = 'sptpol_des'
### noofsims = 5

comps = ['sptpol_map_like_sims/SPTpol_SZ_sample']
paramsfiles = ['sptpol_map_like_sims/DES_RM_sims_testing/params_planck_r_0.0_2015_cosmo_lensed_LSS_crossmap_sptpol_FG_notSZ_sims.txt']
totalclus = 350
which_exp = 'sptpol'
tszarr = [0]


##########################################################################################
##########################################################################################

minSNR, maxSNR = None, None
if produce_maps:
	###for n in range(noofsims):
	start, end = int(sys.argv[2]), int( sys.argv[3] )
	for n in range(start, end):
		#### rsval = np.random.randint(1e6)
		rsval = rsarr[n]
		for i in range(len( comps) ):
			#cmd = 'python s1_sims_QE_tester.py %s 0.5 %s %s %s sptpol_des params/crossmaps/%s crossmaps/%s 0 0' %(totalclus, rsval, noiseval, tszarr[i], paramsfiles[i], comps[i])
			cmd = 'python s1_2_sims_QE_tester_crossmaps.py %s 0.5 %s %s %s %s params/crossmaps/%s crossmaps/%s 0 0 %s %s' %(totalclus, rsval, noiseval, tszarr[i], which_exp, paramsfiles[i], comps[i], minSNR, maxSNR)
			if i<len(comps)-1 or len(comps) == 1:
				cmd = 'nohup %s &' %(cmd)
			print cmd#;sys.exit()
			os.system(cmd)

		'''
		cmd = 'nohup python s1_sims_QE_tester.py %s 0.5  %s %s 0 sptpol_des params/crossmaps/ crossmaps/%s 0 0 &' %(totalclus, rsval, noiseval, comps[1])
		os.system(cmd)

		cmd = 'python s1_sims_QE_tester.py %s 0.5  %s %s 0 sptpol_des params/crossmaps/ crossmaps/%s 0 0' %(totalclus, rsval, noiseval, comps[2])
		os.system(cmd)
		'''
	sys.exit()


##########################################################################################
##########################################################################################
plot_maps = 1
import scipy.optimize as optimize
def fitting_func(p, p0, X, Y, MAP, lbounds = None, ubounds = None, fixed = None, return_fit = 0):
	'''
	if hasattr(fixed, '__len__'):
	p[fixed] = p0[fixed]

	if hasattr(lbounds, '__len__'):
	linds = abs(p)<abs(lbounds)
	p[linds] = lbounds[linds]

	if hasattr(ubounds, '__len__'):
	uinds = abs(p)>abs(ubounds)
	p[uinds] = ubounds[uinds]
	'''

	def fn_gaussian(p, xp, yp):
		#rota = np.radians(p[6])
		#xp = xp * np.cos(rota) - yp * np.sin(rota)
		#yp = xp * np.sin(rota) + yp * np.cos(rota)
		g = p[0]+p[1]*np.exp( -(((p[2]-xp)/p[4])**2+ ((p[3]-yp)/p[5])**2)/2.)

		return g

	if not return_fit:
		return np.ravel(fn_gaussian(p, X, Y) - MAP)
	else:
		return fn_gaussian(p, X, Y)

if plot_maps:
	from pylab import *
	import pickle, gzip, glob#, modules
	from matplotlib import rc;rc('text', usetex=True);rc('font', weight='bold');matplotlib.rcParams['text.latex.preamble'] = [r'\boldmath']
	from matplotlib.font_manager import fontManager, FontProperties

	sbpl = 1
	labfsval = 7
	for n in range(noofsims):
		for i in range(len(comps)):

			#folderpath = 'data/sims/sptpol_des/crossmaps/v1/%s/0.5*white_%s' %(comps[i], noiseval)
			#folderpath = 'data/sims/sptpol_des/crossmaps/%s/0.5*white_%s' %(comps[i], noiseval)
			folderpath = 'data/sims/%s/crossmaps/%s/SNR_%s_%s/' %(which_exp,comps[i], minSNR, maxSNR)
			if not os.path.exists(folderpath):
				folderpath = 'data/sims/%s/crossmaps/%s/' %(which_exp,comps[i])
			folderpath = '%s/0.5*' %(folderpath)
			#print os.path.exists(folderpath)#;sys.exit()
			### files = glob.glob('%s/*stacked*%s.pkl.gz' %(folderpath, rsarr[n]))
			files = glob.glob('%s/*stacked*.pkl.gz' %(folderpath))[n:n+1]

			sn = 10
			for fcnt, fname in enumerate(files):
				print fname
				#if fcnt <> n: continue
				#if fcnt == n: continue
				kappadic = pickle.load(gzip.open(fname,'rb'))
				kappa_qe = kappadic['stacked_kappa_qe']
				dx = dy = kappadic['reso_arcmin']
				add_noise = kappadic['add_noise']
				#expnoiselevel = [kappadic['expnoiselevel'][0]]
				boxsize = kappadic['boxsize']
				nx, ny = kappa_qe.shape

				clra, cldec = 0., 0.

				minval, maxval = clra-boxsize/2/60.,  clra+boxsize/2/60.
				#minval, maxval = 0,  clra+boxsize/60.
				ra = dec = np.linspace(minval, maxval, nx)

				RA, DEC = np.meshgrid(ra,dec)
				"""
				RADEC = [RA, DEC]
				radprofile = sims.fn_radial_profile(kappa_qe, RADEC, bin_size = 3., minbin = 0., maxbin = 30.)
				#radprofile = sims.fn_radial_profile(kappa_qe, RADEC)
				"""
	
				"""
				smooth_beam_arcmins = 1.
				fwhm = np.radians(smooth_beam_arcmins/60.)
				mapparams = [nx, ny, dx, dy]
				lx, ly = sims.get_lxly(mapparams)
				EL = (lx ** 2. + ly ** 2.)**0.5	
				LPF = sims.gauss_beam(fwhm, EL)[0]
				kappa_qe = np.fft.ifft2( np.fft.fft2(kappa_qe) * LPF).real
				"""

				e1, e2 = int(nx/2 - sn/dx/2), int(nx/2 + sn/dx/2)
				kappa_qe = kappa_qe[e1:e2, e1:e2]
				RA = RA[e1:e2, e1: e2] * 60.
				DEC = DEC[e1:e2, e1: e2] * 60.

				#fit Gaussian to the map to get the centroid
				height, amp = 0,np.max(kappa_qe)
				x_cen, y_cen = 0.,0.
				wx, wy = 1., 1. #arcmins
				rot = 0. #degrees
				p0 = [height, amp, x_cen, y_cen, wx, wy, rot]
				p1, pcov, infodict, errmsg, success = optimize.leastsq(fitting_func, p0[:], args=(p0, RA, DEC, kappa_qe), full_output=1)
				xcen, ycen = p1[2], p1[3]
				centroid_arr.append((xcen, ycen))

				if (1):
					MAP_FIT = fitting_func(p1,p1,RA,DEC,kappa_qe,return_fit = 1)

					minval, maxval = None, None#-0.05, 0.2#0.18
					interpval = 'bicubic'#'None'
					ax = subplot(noofsims,len(comps),sbpl); 
					css = imshow(kappa_qe, origin='lower', extent = [np.min(RA), np.max(RA), np.min(DEC), np.max(DEC)], cmap = cm.jet, interpolation = interpval, vmin = minval, vmax = maxval);
					plot(xcen, ycen, 'k+', ms = 8.); text(xcen*2., ycen*2., r'\textbf{(%.3f, %.3f)}' %(xcen, ycen), fontsize = 10, color = 'k')
					#minval, maxval = css.get_clim()
					#cbar = colorbar()
					#cbar = colorbar(ticks = [minval, maxval], shrink = 0.7, orientation='horizontal')#, format='%0.2f')#;plot(0.,0.,'k+', ms = 10.)#;xlim(0, new_nx); ylim(0, new_ny);
					cbar = colorbar(shrink = 0.7, orientation='horizontal')#, format='%0.2f')#;plot(0.,0.,'k+', ms = 10.)#;xlim(0, new_nx); ylim(0, new_ny);
					#cbar.ax.set_yticklabels(['$%.2f$' %(minval), '$%.2f$' %(maxval)])

					cbar.ax.tick_params(labelsize=labfsval+2, length = 0.1)#, pad = -5.)#, labelleft = 1)
					#textlab = '%s' %(comps[i]).replace('_', '\_').replace('am1.2', 'am x 1.2')

					textlab = '$Beams = %s$\n$\\Delta T = %s\ \\mu K^{\\prime}$' %(labels_beam[i], labels_noise[i])
					if n == 0:
						title(textlab, fontsize = 8)
					#axis('off')
					setp(ax.get_xticklabels(), visible=False);setp(ax.get_yticklabels(), visible=False);  
					ax.xaxis.set_ticks_position('none');ax.yaxis.set_ticks_position('none');

					grid(1, ls = 'solid', lw = .5)
					if i == 0:
						ylabel('$Realisation\ %s$' %(n+1), fontsize = 10)
					"""
					title(fname.split('/')[-1])
					subplot(122)#, xscale = 'log'); #plot(radprofile[:,0], radprofile[:,1], 'ro-');
					errorbar(radprofile[:,0],radprofile[:,1], yerr = [radprofile[:,2],radprofile[:,2]], color='r', marker = 'o')#, lw = 2., ls = 'solid')
					#xlim(0, 10.)
					show()
					"""
				sbpl += 1
	show()
	#savefig('reports/20171109_crossmap_results/sims_expectations.png', bbox_inches = 'tight', pad_inches = 0.1)
	#sys.exit()

	#plot centroid array
	cinds = np.linspace(10, 255, len(centroid_arr))
	colorarr = [cm.jet(int(c)) for c in cinds]
	shift = []
	for cntr, (x,y) in enumerate(centroid_arr):
		shift.append( (x**2. + y**2.)**0.5 )
		plot(x, y, '+', color = colorarr[cntr])
		print x, y
		grid(1, ls = 'solid', lw = .5)
	xlim(-5.,5.); ylim(-5.,5.)
	xlabel(r'\textbf{arcmin}', fontsize = 14)
	ylabel(r'\textbf{arcmin}', fontsize = 14)
	show()

