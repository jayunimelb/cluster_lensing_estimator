'''
Input: input folder with simulation of data with stacked kappa, KAPPA_COV matrix or stacked kappa from which cov is calculated


Output : likelihood plot

'''

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

	d = d-dp	
	t3 = -0.5 * np.asarray( np.dot(d.T, np.dot( Cinv, d ))).squeeze()

	logLval = t3# + t2

	return logLval

import numpy as np
import pickle,gzip
import numpy as np, pickle, gzip, sys, glob, time, os, scipy as sc
sys.path.append('modules')
import modules
import sky_local
from pylab import *
sims = modules.scl_cmb.simulations()
ip_folder = sys.argv[1]

fnamearr = glob.glob('%s/stacked*_03500_clusters*'%(ip_folder))
paramfile = glob.glob('%s/*used.txt'%(ip_folder))[0]
kappa_file = glob.glob('%s/kappa_COV*100cluters.pkl.gz'%(ip_folder))[0]
kappa_COV = pickle.load(gzip.open(kappa_file))
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
nx,ny = 200,200; dx,dy = 0.5,0.5 #!!! Hard coded
mapparams = [nx,ny,dx,dy]
boxsize = nx * dx
clra, cldec = 0., 0.
minval, maxval = clra-boxsize/2/60.,  clra+boxsize/2/60.
ra = dec = np.linspace(minval, maxval, nx)
RA, DEC = np.meshgrid(ra,dec)
RADEC = [RA, DEC]
binsize = 1.0;maxbin = 10.
A_arr = np.arange(0,6.0,0.1)*1e14
onedfit = 1
det_sig = np.zeros(len(fnamearr))
nx,ny = 200,200; dx,dy = 0.5,0.5 #!!! Hard coded
mapparams = [nx,ny,dx,dy]
boxsize = nx * dx
clra, cldec = 0., 0.
minval, maxval = clra-boxsize/2/60.,  clra+boxsize/2/60.
ra = dec = np.linspace(minval, maxval, nx)
RA, DEC = np.meshgrid(ra,dec)
RADEC = [RA, DEC]
binsize = 1.0
maxbin = 10.0
TWODTF = sims.fn_get_EBAX_2016_anal(mapparams, l1=sims.inidic['l1'], l2=sims.inidic['l2'], l3=sims.inidic['l3'])[0]
kappa_file = glob.glob('%s/kappa_COV*100cluters.pkl.gz'%(ip_folder))[0]
kappa_COV = pickle.load(gzip.open(kappa_file))

for fcnt,fname in enumerate(fnamearr):
	logLarr = np.zeros(len(A_arr))
	sim_data = pickle.load(gzip.open(fname))
	for i,AA in enumerate(A_arr):
		if os.path.exists('%s/test_stacked_model_kappa_A_%s.pkl.gz'%(ip_folder,AA)):
			stacked_kappa = pickle.load(gzip.open('%s/test_stacked_model_kappa_A_%s.pkl.gz'%(ip_folder,AA)))

		else:
			kappa_model = np.zeros((len(sim_data['CLUS_IDENTIFIER']), sim_data['stacked_kappa_qe'].shape[0],sim_data['stacked_kappa_qe'].shape[0]))
			for j,clus_iden in enumerate(sim_data['CLUS_IDENTIFIER']):
				zz = clus_iden[2]
				rich = clus_iden[3]
				MM = sims.fn_rich_mass_M17(rich,zz, A = AA)
				h = param_dict['h']
				cc = sims.c_Duffy(MM,zz,h)
				param_dict['rho_def'] = 'mean'
				param_dict = sim_data['param_dict']	
				sims.exp_beam = param_dict['exp_beam']
				sims.inidic = param_dict
				sims.fn_lensing_ini(mapparams,param_dict, RA, DEC, [clra], [cldec], [MM], [cc], [zz], param_dict['z_lss'], param_dict['mass_def'], param_dict['rho_def'], truncate_kappa_at_radius = param_dict['truncate_kappa_at_radius'])
				kappa_true = sims.KAPPA
				lx, ly = sims.get_lxly(mapparams)
				L = (lx**2. + ly**2.)**0.5
				ngaus = int( np.sqrt(8. * np.log(2.)) / np.radians(sims.exp_beam/60.) )
				above_beam_scale = np.where(L>=ngaus)
				kappa_true_fft = np.fft.fft2(kappa_true)
				kappa_true_fft[above_beam_scale] = 0.					
				TWODTF = sims.fn_get_EBAX_2016_anal(mapparams, l1=sims.inidic['l1'], l2=sims.inidic['l2'], l3=sims.inidic['l3'])[0]	
				kappa_true = np.fft.ifft2((kappa_true_fft * TWODTF)).real
				kappa_model[j] = kappa_true

			
			stacked_kappa = np.mean(kappa_model,axis =0)
						
			pickle.dump(stacked_kappa,gzip.open('%s/test_stacked_model_kappa_A_%s.pkl.gz'%(ip_folder,AA),'w'))
		
		logL = fn_calc_likelihood(sim_data['stacked_kappa_qe'], stacked_kappa, kappa_COV/35., AA, fcnt, onedfit = onedfit)
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
	
	tmp = logLarr - max(logLarr)
	L = np.exp(tmp); L/=max(L)
	plot(A_arr,L)
print np.mean(det_sig)
#axvline(x = 2.35*1e14)
show()
