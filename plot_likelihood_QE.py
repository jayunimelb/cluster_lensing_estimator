'''
Input: input folder with simulation of data with stacked kappa, KAPPA_COV matrix or stacked kappa from which cov is calculated


Output : likelihood plot

'''


import numpy as np
import pickle,gzip
import numpy as np, pickle, gzip, sys, glob, time, os, scipy as sc
sys.path.append('modules')
import modules
import sky_local
#import modules.radialProfile as radialProfile
from pylab import *

def fn_calc_likelihood(MAP, MODEL, C, modelmass, simcnt, onedfit = 0):

	#imshow(C);colorbar();show();quit()
	Cinv = sc.linalg.pinv2(C)
	C = np.mat(C)

	sign, logdetval = np.linalg.slogdet(C)
	logdetval = logdetval * sign

	#t1 = -.5 * npixels * np.log( 2 * np.pi )
	#t2 = -0.5 * logdetval

	if onedfit:
		#raprfmap = sims.fn_radial_profile(MAP, RADEC)
		raprfmap = sims.fn_radial_profile(MAP, RADEC, bin_size=binsize, minbin=0.0, maxbin=maxbin)
		d = raprfmap[:,1].flatten() ### - np.mean(raprfmap[:,1].flatten())

		#raprfmodel = sims.fn_radial_profile(MODEL, RADEC)
		raprfmodel = sims.fn_radial_profile(MODEL, RADEC, bin_size=binsize, minbin=0.0, maxbin=maxbin)
		dp = raprfmodel[:,1].flatten() ### - np.mean(raprfmodel[:,1].flatten())

	else:

		d= MAP.flatten()
		dp = MODEL.flatten()


	'''
	subplot(221);imshow(MAP);colorbar();title(np.mean(MAP))
	subplot(222);imshow(MODEL);colorbar();title(np.mean(MODEL))
	subplot(223)
	plot(raprfmap[:,0], d, 'ro')
	colorinds = np.linspace(0,255,len(Marr)); colorarr = np.asarray( [cm.jet( int(c) ) for c in colorinds] )
	ccurrcolind = np.where(Marr == MM)[0][0];colorval = colorarr[ccurrcolind]
	plot(raprfmap[:,0], dp, color = colorval)#;title('SIM = %s; Model mass = %s' %(simcnt, modelmass))
	show()
	'''

	d = d-dp	
	t3 = -0.5 * np.asarray( np.dot(d.T, np.dot( Cinv, d ))).squeeze()

	logLval = t3# + t2

	return logLval
sims = modules.scl_cmb.simulations()
kov_data = pickle.load(gzip.open('/Users/sanjaykumarp/transfer/spt/stacked_kappa_sims_for_cov_file_single_cluster.pkl.gz_500sims'))

nx,ny = 200,200; dx,dy = 0.5,0.5 #!!! Hard coded
mapparams = [nx,ny,dx,dy]
boxsize = nx * dx
clra, cldec = 0., 0.
minval, maxval = clra-boxsize/2/60.,  clra+boxsize/2/60.
ra = dec = np.linspace(minval, maxval, nx)
RA, DEC = np.meshgrid(ra,dec)
RADEC = [RA, DEC]
binsize = 1.0;maxbin = 10.
RADPROFILES = np.asarray( map(lambda x: sims.fn_radial_profile(x, RADEC, bin_size=binsize, minbin=0.0, maxbin=maxbin), kov_data) )
RADPRF = RADPROFILES[:,:,1]
totbins = np.shape(RADPRF)[1]
noofsims = kov_data.shape[0]
kappa_COV = sims.calcCov(RADPRF, noofsims, npixels = totbins)
gauss = 1
if gauss:
	kappa_COV = pickle.load(gzip.open('gauss_files/SNR_None_None/0.5_no_tSZ_sptpol/kappa_COV_100sims_10.0am_1.0dx_10cluters.pkl.gz'))
	#kappa_COV = kappa_COV/10.
#paramfile = 'srini_files/SNR_None_None/0.5_no_tSZ_sptpol/params_planck_r_0.0_2015_cosmo_lensed_LSS_crossmap_sptpol_FG_notSZ_sims_year3_lgt_5_sims_used.txt'
paramfile = 'gauss_files/SNR_None_None/0.5_no_tSZ_sptpol/params_planck_r_0.0_2015_cosmo_lensed_LSS_crossmap_sptpol_FG_notSZ_sims_year3_lgt_5_sims_Gaussbeam_used.txt'
params = np.recfromtxt(paramfile,usecols=[0],delimiter = '=')
paramvals = np.recfromtxt(paramfile,usecols=[1],delimiter = '=')
param_dict = {}
cc = 3.0
for p,pval in zip(params,paramvals):
        tmp = pval.strip()
        try:
                float(tmp)
                if tmp.find('.')>-1:
                        param_dict[p.strip()] = float(tmp)
                else:
                        param_dict[p.strip()] = int(tmp)
        except:
                if tmp == 'None':
                        param_dict[p.strip()] = None
                elif tmp[0] == '[':
                        #param_dict[p.strip()] = np.asarray( map(tmp[1:-1].split(',')) )
                        dummy = tmp[1:-1].split(',')[0]
                        try:
                                param_dict[p.strip()] = np.asarray( map(float, tmp[1:-1].split(',')) )
                        except:                         
                                param_dict[p.strip()] = tmp[1:-1].split(',')
                else:
                        param_dict[p.strip()] = tmp.replace('\'','')

sims._ini_params(param_dict)

#fnamearr = glob.glob('srini_files/SNR_None_None/0.5_no_tSZ_sptpol/stac*pkl.gz')
fnamearr = glob.glob('gauss_files/SNR_None_None/0.5_no_tSZ_sptpol/stac*pkl.gz')

TWODTF = sims.fn_get_EBAX_2016_anal(mapparams, l1=sims.inidic['l1'], l2=sims.inidic['l2'], l3=sims.inidic['l3'])[0]
A_arr = np.arange(0,5,0.2)*1e14
onedfit = 1
det_sig = np.zeros(len(fnamearr))
catalog = '/Users/sanjaykumarp/transfer/spt/map_150.pkl.gz_cutouts_150ghz_no_downsampling_year3_lgt5_redmapper_clusters_TQUforMLE'
cluskeys = sorted(pickle.load(gzip.open(catalog))['cutouts'].keys())
#imshow(kappa_COV/35.);colorbar();show();quit()
for fcnt,fname in enumerate(fnamearr):
	logLarr = np.zeros(len(A_arr))
	sim_data = pickle.load(gzip.open(fname))
	for i,AA in enumerate(A_arr):
		if i <2:
			continue
		if 0==1:#os.path.exists('stacked_model_kappa_A_%s.pkl.gz'%(AA)):
			stacked_kappa = pickle.load(gzip.open('stacked_model_kappa_A_%s.pkl.gz'%(AA)))
		else:
			kappa_model = np.zeros((len(sim_data['CLUS_IDENTIFIER']), sim_data['stacked_kappa_qe'].shape[0],sim_data['stacked_kappa_qe'].shape[0]))
			for j,clus_iden in enumerate(sim_data['CLUS_IDENTIFIER']):
				zz = clus_iden[2]
				rich = clus_iden[3]
				MM = sims.fn_rich_mass_M17(rich,zz, A = AA)
				h = param_dict['h']
				cc = sims.c_Duffy(MM,zz,h)
				param_dict = sim_data['param_dict']
				sims.fn_lensing_ini(mapparams,param_dict, RA, DEC, [clra], [cldec], [MM], [cc], [zz], param_dict['z_lss'], param_dict['mass_def'], param_dict['rho_def'], truncate_kappa_at_radius = param_dict['truncate_kappa_at_radius'])
				kappa_true = sims.KAPPA
				lx, ly = sims.get_lxly(mapparams)
				L = (lx**2. + ly**2.)**0.5
				sims.exp_beam = param_dict['exp_beam']
				sims.fn_beam_stuff(mapparams)
				print sims.exp_beam;sys.exit()
				ngaus = int( np.sqrt(8. * np.log(2.)) / np.radians(sims.exp_beam/60.) )
				try:
					ngaus = int( np.sqrt(8. * np.log(2.)) / np.radians(sims.exp_beam/60.) )
				except:
					ngaus = int( np.sqrt(8. * np.log(2.)) / np.radians(sims.inidic['exp_beam']/60.) )
				above_beam_scale = np.where(L>=ngaus)
				kappa_true_fft = np.fft.fft2(kappa_true)
				kappa_true_fft[above_beam_scale] = 0.					
				kappa_true = (kappa_true_fft * TWODTF)
				kappa_true = np.fft.ifft2(kappa_true_fft).real
				kappa_model[j] = kappa_true
			stacked_kappa = np.mean(kappa_model,axis =0)
			pickle.dump(stacked_kappa,gzip.open('stacked_model_kappa_A_%s.pkl.gz'%(AA),'w'))
		logL = fn_calc_likelihood(sim_data['stacked_kappa_qe'], stacked_kappa, kappa_COV/350., AA, fcnt, onedfit = onedfit)
		logLarr[i] = logL
	max_mass = A_arr[np.where(logLarr == max(logLarr))[0]]
	''' #!!! to show radial profiles

	kappa_max_mass = pickle.load(gzip.open('stacked_model_kappa_%s.pkl.gz'%(float('%.1g' % max_mass)),'rb'))
	RADPROFILES_model = np.asarray( map(lambda x: sims.fn_radial_profile(x, RADEC, bin_size=binsize, minbin=0.0, maxbin=maxbin), [kappa_max_mass]) )
	RADPRF_model = RADPROFILES_model[:,:,1]
	RADPROFILES_data = np.asarray( map(lambda x: sims.fn_radial_profile(x, RADEC, bin_size=binsize, minbin=0.0, maxbin=maxbin), [sim_data['stacked_kappa_qe']]) )
	RADPRF_data = RADPROFILES_data[:,:,1]
	'''
	det_sig[fcnt] = np.sqrt(max(logLarr)-logLarr[0])
	print det_sig[fcnt]
	tmp = logLarr - max(logLarr)
	L = np.exp(tmp); L/=max(L)
	plot(A_arr,L)

show()
