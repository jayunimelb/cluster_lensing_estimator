'''
Input: input folder with simulation of data with stacked kappa, KAPPA_COV matrix or stacked kappa from which cov is calculated


Output : likelihood plot

'''
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
		passed_inds = np.asarray( passed_inds )
	return kappa_qe_arr[passed_inds], CLUS_IDENTIFIER[passed_inds]

def fn_ang_dist(ip_x1, ip_y1, ip_x2, ip_y2):
	x1 = np.radians(ip_x1)
	y1 = np.radians(ip_y1)
	x2 = np.radians(ip_x2)
	y2 = np.radians(ip_y2)
	ang_dist_rad = np.arccos(np.sin(y1)*np.sin(y2)+(np.cos(y1)*np.cos(y2)*np.cos(x1 - x2)))
	ang_dist_deg = np.degrees(ang_dist_rad)
	return ang_dist_deg * 60.

def get_kappa_model(param_dict,clus_iden, norm=2.35*1e14, alpha_fit = 1.12):
	zz,rich = clus_iden[2],clus_iden[3]
	param_dict = sim_data['param_dict']
	h = param_dict['h']
	MM = sims.fn_rich_mass_M17(rich,zz, A= norm)
	cc = sims.c_Duffy(MM,zz,h)
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
	return kappa_true


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

psrcfile = 'data/500sqdeg_bleem/quick_mm_point_source_file_150.0GHz_6.00000mJy.txt'
psrcdata = np.loadtxt(psrcfile,usecols=[1,2,5],skiprows = 10)
psrcra,psrcdec,psrcflux = np.asarray(zip(*psrcdata))
params = np.recfromtxt(paramfile,usecols=[0],delimiter = '=')
paramvals = np.recfromtxt(paramfile,usecols=[1],delimiter = '=')
param_dict = {}

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


use_data = 1

if use_data:
	fnamearr = glob.glob('%s/data_*'%(ip_folder))

	
#if new_des_year3_full:
#	fnamearr = glob.glob('%s/*full/data*'%(ip_folder)) 
#if new_des_year3_vl:
	#fnamearr = glob.glob('%s/*vl/data*'%(ip_folder))
sims._ini_params(param_dict)
A_min, A_max = float(sys.argv[2]), float(sys.argv[3])
A_arr = np.arange(A_min,A_max,0.1)*1e14
onedfit = 1
ptsrc_cuts = 0
rich_cuts = 0
rich_min = 24
min_sep_ptsrc = 17 
#mean_field = pickle.load(gzip.open('%s/MEANFIELD_500sims_10.0am_1.0dx_100cluters.pkl.gz'%(ip_folder)))# doesn't work now (have to change for sims)


if use_data:
	if not param_dict['perform_mean_sub']:
		mean_field = pickle.load(gzip.open('/data19/sri/QE/2016_11/sims/sptpol_des/crossmaps/sptpol_map_like_sims/data/year3_lgt_5/SNR_20.0_1000.0/0.5_no_tSZ_sptpol/data_stacked_kappa_21829_clusters.pkl.gz_20180123_212904_redmapper_randoms'))
#!!! Check param dict and sims.inidic make sure they are the same!!!
det_sig = np.zeros(len(fnamearr))
for fcnt,fname in enumerate(fnamearr):
	logLarr = np.zeros(len(A_arr))
	sim_data = pickle.load(gzip.open(fname))
	sim_data['param_dict'] = param_dict 
	nx,ny = sim_data['stacked_kappa_qe'].shape
	nx,ny = nx/2, ny/2
	dx,dy = sim_data['reso_arcmin'],sim_data['reso_arcmin']
	clra, cldec = 0., 0.
	boxsize = nx * dx
	minval, maxval = clra-boxsize/2/60.,  clra+boxsize/2/60.
	ra = dec = np.linspace(minval, maxval, nx)
	RA, DEC = np.meshgrid(ra,dec)
	RADEC = [RA, DEC]
	binsize = 1.0;maxbin = 10.
	mapparams = [nx,ny,dx,dy]
	TWODTF = sims.fn_get_EBAX_2016_anal(mapparams, l1=sims.inidic['l1'], l2=sims.inidic['l2'], l3=sims.inidic['l3'])[0] ## make sure you pass the correct dictionary file
	for i,AA in enumerate(A_arr):
		# include mkdir command here to create kappa model
		if os.path.exists('%s/model_kappa_norm/test_stacked_model_kappa_A_%s.pkl.gz'%(ip_folder,AA)):
			stacked_kappa = pickle.load(gzip.open('%s/model_kappa_norm/test_stacked_model_kappa_A_%s.pkl.gz'%(ip_folder,AA)))
		else:
			kappa_model = np.zeros((len(sim_data['CLUS_IDENTIFIER']), nx,ny))
			for j,clus_iden in enumerate(sim_data['CLUS_IDENTIFIER']):
				kappa_model[j] = get_kappa_model(sim_data['param_dict'], clus_iden, norm = AA)
			stacked_kappa = np.mean(kappa_model,axis =0)
			pickle.dump(stacked_kappa,gzip.open('%s/model_kappa_norm/test_stacked_model_kappa_A_%s.pkl.gz'%(ip_folder,AA),'w'))
		#if 1==1: ## just to dump the model and not to calculate the log likelihood
		#continue 
		## though it says sim_data if use_data = 1 then it is DATA!!
 		if use_data  and  ptsrc_cuts:
			new_kappa_qe_arr = []
			kappa_qe_arr = sim_data['kappa_qe_arr']
			CLUS_IDENTIFIER = sim_data['CLUS_IDENTIFIER']
			for keycnt, keyname in enumerate( CLUS_IDENTIFIER):
				ra, dec, z, conf, weight = keyname
				currra, currdec = keyname[0:2]
				ang_dist_arcmins = fn_ang_dist(psrcra, psrcdec, np.tile(currra,len(psrcdec)), np.tile(currdec,len(psrcdec)) )
				if min(ang_dist_arcmins) <= min_sep_ptsrc:
					continue
				new_kappa_qe_arr.append(kappa_qe_arr[keycnt])
			total_clus = len(new_kappa_qe_arr)
			new_kappa_qe_arr = np.asarray(new_kappa_qe_arr)
			logL = fn_calc_likelihood(np.mean(new_kappa_qe_arr,axis =0)-mean_field['stacked_kappa_qe'], stacked_kappa, kappa_COV/35., AA, fcnt, onedfit = onedfit)
		elif use_data and rich_cuts:
			new_kappa_qe_arr = []
			kappa_qe_arr = sim_data['kappa_qe_arr']
			CLUS_IDENTIFIER = sim_data['CLUS_IDENTIFIER']
			for keycnt, keyname in enumerate( CLUS_IDENTIFIER):
				ra, dec, z, conf, weight = keyname
				if z < rich_min:
					continue
				new_kappa_qe_arr.append(kappa_qe_arr[keycnt])
			total_clus = len(new_kappa_qe_arr)
			new_kappa_qe_arr = np.asarray(new_kappa_qe_arr)
			logL = fn_calc_likelihood(np.mean(new_kappa_qe_arr,axis =0)-mean_field['stacked_kappa_qe'], stacked_kappa, kappa_COV/35., AA, fcnt, onedfit = onedfit)
		else:
			#from IPython import embed;embed()
			total_clus = len(sim_data['CLUS_IDENTIFIER'])
			
			#imshow(sim_data['stacked_kappa_qe'][90:110,90:110]);colorbar();show()#;quit()
			imshow(stacked_kappa);colorbar();title('%s'%(AA));show()#;quit()
			logL = fn_calc_likelihood(sim_data['stacked_kappa_qe'], stacked_kappa, kappa_COV/15., AA, fcnt, onedfit = onedfit)
			total_clus = sim_data['totalclus']
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
	output = {}
	print fname
	output[fname] = np.zeros((2,len(L)))
	output[fname][0,:] = A_arr
	output[fname][1,:] = L
	pickle.dump(output,gzip.open('%s/lkhd_output_data_%s_total_clus.pkl.gz'%(ip_folder,total_clus),'w'))
 	#plot(A_arr,L)
exit()
show()
