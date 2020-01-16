"""

get kappa for simulations

"""


import numpy as np
import modules,pickle,gzip,time,os

scl_cmb = modules.scl_cmb
sims = scl_cmb.simulations()
sims_data_file = 'data_sp/sptsz/min_variance/simulated_lensed_cmb_cutouts_500_clusters_150GHz_Arnaud_tSZ_fg_noise.pkl.gz'
sims_data_file_90 = 'data_sp/sptsz/min_variance/lensed_cmb_cutouts_500_clusters_90GHz_Arnaud_tSZ_fg_noise.pkl.gz'
sims_data_file_220 = 'data_sp/sptsz/min_variance/lensed_cmb_cutouts_500_clusters_220GHz_Arnaud_tSZ_fg_noise.pkl.gz'
lensed_cmb_data = pickle.load(gzip.open(sims_data_file))
lensed_cmb_data_90 = pickle.load(gzip.open(sims_data_file_90))
lensed_cmb_data_220 = pickle.load(gzip.open(sims_data_file_220))

cutouts = lensed_cmb_data['cutouts']

paramfile = 'plnk_sz_grad_maps.txt'
params = np.recfromtxt(paramfile,usecols=[0],delimiter = '=')
paramvals = np.recfromtxt(paramfile,usecols=[1],delimiter = '=')
param_dict = {}
import time
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


def fn_apodize(MAP, mapparams, mask = 'square', just_return_mask = 1):

	import numpy as np, scipy.ndimage as ndimage, scipy.signal as signal

	nx, ny, dx, dy = mapparams

	start =  time.time()
	pix = dx
	radius = (nx * pix)/10. #(nx * pix)/3. #changed on 20171206
	npix_cos = int(radius/pix)
	ker=np.hanning(npix_cos)
	ker2d=np.asarray( np.sqrt(np.outer(ker,ker)) )

	##from IPython import embed; embed()
	MASKf=np.zeros((nx,ny))
	minval, maxval = -(nx*pix)/2,  (nx*pix)/2
	x = y = np.linspace(minval, maxval, nx)
	X, Y = np.meshgrid(x,y)
	xc, yc = 0., 0.
	if mask == 'circle':
		if padding_zeros>1:
			radius = (nx * dx/2) - (nx/3) ###(nx/50) #changed on 20171206
		else:
			radius = (nx * dx/2) - (nx/50) ###(nx/50) #changed on 20171206
		inds=np.where((X-xc)**2. + (Y-yc)**2. <= radius**2.) #all in arcmins
	elif mask == 'square':
		radius = (nx * dx/2) - 1.
		inds=np.where((abs(X)<=radius) & (abs(Y)<=radius)) #all in arcmins

	MASKf[inds]=1.

	apodMASKf=ndimage.convolve(MASKf, ker2d)#, mode='wrap')
	apodMASKf/=apodMASKf.max()

	###imshow(apodMASKf);colorbar();show()

	return apodMASKf
## overwriting params dict file here
param_dict['fit_FG'] =  'all_Gaussian'
sims.inidic = param_dict
padding_zeros = 1
boxsize = param_dict['boxsize']
dx,dy = 0.5,0.5
nx,ny = int(boxsize/dx), int(boxsize/dy)
simmapparams = nx,ny,dx,dy
Dlfile_len = param_dict['Dlfile_len']
Dlfile_unlen = param_dict['Dlfile_unlen']
Dls_len = np.loadtxt(Dlfile_len,usecols=[0,1])
Dls_unlen = np.loadtxt(Dlfile_unlen,usecols=[0,1])
timedatestr='%s' %(time.strftime('%Y%m%d_%H%M%S', time.localtime()))
log_file = 'tmp/logs_%s.txt' %(timedatestr)
if os.path.exists(log_file):
	os.system('rm %s' %(log_file))
sims._ini_log_file(log_file)
## below parameters are hard coded, based on simulations
sims.just_SZfreemap = 0
sims.exp_beam = 1
sims.noise_present = 1
sims.expnoiselevel = [18.]
sims.sehgal_fg_added = 0
sims.matched_filter = 0
lgmca_sims = 0
sim_lensed_cmb_cutouts = pickle.load(gzip.open('data_sp/sptsz/min_variance/lensed_cmb_cutouts_500_clusters.pkl.gz'))
sims.cmb_noise_randomseedval = 1
ilc_method = 1
sims.ilc_method = ilc_method
cilc_var = np.zeros(len(cutouts.keys()))
for i,key in enumerate(cutouts.keys()):
	apodMASK = fn_apodize(cutouts[key][0], simmapparams, mask = 'circle', just_return_mask = 1)
	cutouts[key] = cutouts[key]#/1e6
	#OBSMAP = cutouts[key]*apodMASK
	# simulate lgmca maps
	if lgmca_sims:
		# convolve with 5 arcmin beam and noise of 45uK'
		# no tSZ or foregrounds 
		mapparams = simmapparams
		Bl= sims.fn_beam_stuff( mapparams, use_beam = 1, nu = 150, return_beam = 1, exp_beam = 5, beam_error = 0)[0]
		sims.Bl_gradient_map = Bl
		sim_lensed_cmb_cutouts[i,0] = np.fft.ifft2( np.fft.fft2( sim_lensed_cmb_cutouts[i,0] ) * Bl ).real
		CMB_noise = sims.fn_tf_noise(np.asarray([sim_lensed_cmb_cutouts[i]]), mapparams, 1,  add_TF = param_dict['add_TF'], nu =90,add_noise = 1,  cluscnt = i,expnoiselevel =[45.]) 
		CMB_noise = CMB_noise/1e6
		plnk_map = CMB_noise[0]*apodMASK
		# should be over written
		sims.inidic['cross_maps'] = 1
		sims.inidic['add_TF_grad'] = 1

		sims.expnoiselevel = [18.]
		sims.use_lgmca = 1
	if ilc_method:
	
		ilc_wghts = np.asarray([ 0.17010394,  0.82393362,  0.00596244])
		cilc_wghts = np.asarray([-1.33470857,  2.22896331,  0.10574526])
		data_150 = cutouts[key] 
		data_90 =  lensed_cmb_data_90[i]
		data_220 = lensed_cmb_data_220[i]
		final_map =  ilc_wghts[0]*data_90 + ilc_wghts[1]*data_150 + ilc_wghts[2]*data_220
		#final_map = final_map/1e6
		if i !=499:
			cilc_var[i] = np.std(final_map)
			continue
		else:
			cilc_var[i] = np.std(final_map)
			from IPython import embed;embed()
		final_map = final_map/1e6
		OBSMAP = final_map*apodMASK

 	KAPPA_QE = sims.fn_get_kappa_QE(OBSMAP, simmapparams, Dls_len, Dls_unlen, tszfree = param_dict['perform_tSZ_removal'], OBSMAP2 = plnk_map, use_data =0)
	kappa_real = (KAPPA_QE.real - np.mean(KAPPA_QE.real))
	if i ==0:
		stacked_kappa_qe = kappa_real
	else:
		stacked_kappa_qe = stacked_kappa_qe + kappa_real	
	stacked_kappa_qe = stacked_kappa_qe/500.
if not lgmca_sims:
	pickle.dump(stacked_kappa_qe,gzip.open('data_sp/sptsz/min_variance/stacked_kappa_150Ghz.pkl.gz','w'))
if lgmca_sims and ilc_method:
	pickle.dump(stacked_kappa_qe,gzip.open('data_sp/sptsz/min_variance/stacked_kappa_150Ghz_vs_plnck_map_cilc.pkl.gz','w'))

