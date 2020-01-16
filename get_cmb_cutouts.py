"""
load lensed_cmb cutotus 

add tSZ + add TF and noise +  

"""



import numpy as np
import pickle,gzip
import modules,time,os,sys

scl_cmb = modules.scl_cmb
sims = scl_cmb.simulations()
paramfile = 'plnk_sz_grad_maps.txt'
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



lensed_cmb_file = 'data_sp/sptsz/min_variance/lensed_cmb_cutouts_500_clusters.pkl.gz'
lensed_cmb_cutouts = pickle.load(gzip.open(lensed_cmb_file,'r'))
lensed_cmb_90,lensed_cmb_150, lensed_cmb_220= np.zeros(lensed_cmb_cutouts.shape), np.zeros(lensed_cmb_cutouts.shape),np.zeros(lensed_cmb_cutouts.shape)
totalclus = len(lensed_cmb_cutouts)
## adding tSZ (note it has CIB, radio and other foregrounds)
tSZ_file = 'data_sp/sptsz/min_variance/tSZ_CIB_radio_sehgal_sims.pkl.gz'
tSZ= pickle.load(gzip.open(tSZ_file,'r'))
tSZ_90,tSZ_150 = tSZ['90Ghz'], tSZ['150Ghz']

boxsize = param_dict['boxsize']
dx,dy = 0.5,0.5
nx, ny = int(boxsize/dx), int(boxsize/dy)
mapparams = [nx, ny, dx, dy]
add_noise = 1
timedatestr='%s' %(time.strftime('%Y%m%d_%H%M%S', time.localtime()))
log_file = 'tmp/logs_%s_%s.txt' %(timedatestr,2)
if os.path.exists(log_file):
	os.system('rm %s' %(log_file))
sims._ini_log_file(log_file)
# setting it arbitararilty for now check later whether it will bias the results
sims.cmb_noise_randomseedval = 200
sims.inidic = param_dict
Arnaud_profile  =1
noise_90, noise_150, noise_220 = [40.], [18.], [220.] # in uK'
sims.tqulen = 1
for i in range(totalclus):


	## for 90 and 150 Ghz cutouts adding tSZ from Sehgal sims and it has other foregrounds  
	
	
	if not Arnaud_profile:
		lensed_cmb_90[i,0] = lensed_cmb_cutouts[i,0] + tSZ_90[i]
		lensed_cmb_150[i,0]= lensed_cmb_cutouts[i,0] + tSZ_150[i]

	## for Arnaud's profile just add tSZ
	else:
		nx_smaller = ny_smaller = 60
		s, e = int( (nx-nx_smaller)/2 ), int( (nx+nx_smaller)/2 )
		simmapparams_smaller = [nx_smaller, ny_smaller, dx, dy]
		tSZ_150= sims.fn_get_Arnaud_ymap_from_Nikhel_Gupta(simmapparams_smaller, [4e14], [0.7], cosmo_params_dict = param_dict, dx_finer = 0.1, mass_def = param_dict['mass_def'], rho_def = param_dict['rho_def'],random_seed = i)[0]
		tSZ_90 = sims.fn_get_Arnaud_ymap_from_Nikhel_Gupta(simmapparams_smaller, [4e14], [0.7], cosmo_params_dict = param_dict, dx_finer = 0.1, mass_def = param_dict['mass_def'], rho_def = param_dict['rho_def'],random_seed = i, nu = 90e9)[0]
		beta_model_90, beta_model_150 = np.zeros( (nx, ny) ), np.zeros( (nx, ny) )
		beta_model_90[s:e, s:e],beta_model_150[s:e, s:e] = tSZ_90,tSZ_150
		## add foregrounds for arnaud profile
		fg_90 = sims.fn_get_foreground_power('all_Gaussian', mapparams, 1, perform_lensing = 0, nu1=90, nu2=90, ign_DG_RG = 0)[0] * 1e6
		from IPython import embed;embed()
		fg_150 = sims.fn_get_foreground_power('all_Gaussian', mapparams, 1, perform_lensing = 0, nu1=150, nu2=150, ign_DG_RG = 0)[0]*1e6
		
		lensed_cmb_90[i,0] = lensed_cmb_cutouts[i,0]  + fg_90[0] + beta_model_90
		lensed_cmb_150[i,0] = lensed_cmb_cutouts[i,0]  + fg_150[0] + beta_model_150
	# as 220 Ghz channel doesn't have tSZ, just adding george et al foregrounds
	lensed_cmb_220[i,0] = lensed_cmb_cutouts[i,0] + sims.fn_get_foreground_power('all_Gaussian', mapparams, 1, perform_lensing = 0, nu1=220, nu2=220, ign_DG_RG = 0)[0]
	## now beam is convolved --- should be added later but for now as beam doesn't vary between frequencies so ignored
	
	Bl= sims.fn_beam_stuff( mapparams, use_beam = 1, nu = 150, return_beam = 1, exp_beam = 1, beam_error = 0)
	lensed_cmb_90[i,0]  = np.fft.ifft2( np.fft.fft2( lensed_cmb_90[i,0]  ) * Bl ).real
	lensed_cmb_150[i,0]  = np.fft.ifft2( np.fft.fft2( lensed_cmb_150[i,0]  ) * Bl ).real
	lensed_cmb_220[i,0]  = np.fft.ifft2( np.fft.fft2( lensed_cmb_220[i,0]  ) * Bl ).real
	lensed_cmb_90[i] = sims.fn_tf_noise(np.asarray([lensed_cmb_90[i]]), mapparams, 1,  add_TF = param_dict['add_TF'], nu =90,add_noise = add_noise,  cluscnt = i,expnoiselevel =noise_90) 

	lensed_cmb_150[i] = sims.fn_tf_noise(np.asarray([lensed_cmb_150[i]]), mapparams, 1, add_TF = param_dict['add_TF'],nu = 150, add_noise = add_noise,  cluscnt = i,expnoiselevel =noise_150) 
	lensed_cmb_220[i] = sims.fn_tf_noise(np.asarray([lensed_cmb_220[i]]), mapparams, 1, add_TF = param_dict['add_TF'],nu = 220, add_noise = add_noise,  cluscnt = i, expnoiselevel =noise_220) 

if not Arnaud_profile:
	pickle.dump(lensed_cmb_90,gzip.open('data_sp/sptsz/min_variance/lensed_cmb_cutouts_%s_clusters_90GHz_fg_noise.pkl.gz'%(totalclus),'w'))
	pickle.dump(lensed_cmb_150,gzip.open('data_sp/sptsz/min_variance/lensed_cmb_cutouts_%s_clusters_150GHz_fg_noise.pkl.gz'%(totalclus),'w'))
	pickle.dump(lensed_cmb_220,gzip.open('data_sp/sptsz/min_variance/lensed_cmb_cutouts_%s_clusters_220GHz_fg_noise.pkl.gz'%(totalclus),'w'))
else:
	pickle.dump(lensed_cmb_90,gzip.open('data_sp/sptsz/min_variance/lensed_cmb_cutouts_%s_clusters_90GHz_Arnaud_tSZ_fg_noise.pkl.gz'%(totalclus),'w'))
	pickle.dump(lensed_cmb_150,gzip.open('data_sp/sptsz/min_variance/lensed_cmb_cutouts_%s_clusters_150GHz_Arnaud_tSZ_fg_noise.pkl.gz'%(totalclus),'w'))
	pickle.dump(lensed_cmb_220,gzip.open('data_sp/sptsz/min_variance/lensed_cmb_cutouts_%s_clusters_220GHz_Arnaud_tSZ_fg_noise.pkl.gz'%(totalclus),'w'))
