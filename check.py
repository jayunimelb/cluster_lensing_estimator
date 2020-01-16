

import numpy as np 
import pickle,gzip,sys
import modules
import scipy as sc

scl_cmb = modules.scl_cmb

sims = scl_cmb.simulations()
def twoD_Gaussian( (x, y), amplitude,xo, yo, sigma_x, sigma_y, offset):
	"""
	convert fwhm to sigma
	"""

	a = 1./(2*sigma_x**2) 
	c = 1./(2*sigma_y**2)
	
	g = offset + amplitude*np.exp( - (a*((x-xo)**2)  + c*((y-yo)**2)))
		
	return g.ravel()


nx, ny = 200,200
dx = dy = 0.5
boxsize = nx * dx
mapparams = [nx, ny, dx, dy]
clra, cldec = 0., 0.
minval, maxval = clra-boxsize/2/60.,  clra+boxsize/2/60.
ra = dec = np.linspace(minval, maxval, nx)
RA, DEC = np.meshgrid(ra,dec)
RADEC = [RA, DEC]
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
sims._ini_params(param_dict)


def twoD_Gaussian(self, (x, y), amplitude,xo, yo, sigma_x, sigma_y, offset):
	a = 1./(2*sigma_x**2) 
	c = 1./(2*sigma_y**2)
	g = offset + amplitude*np.exp( - (a*((x-xo)**2)  + c*((y-yo)**2)))
	return g.ravel()


def fn_get_model(mm,param_dict,mapparams):
	nx,ny,dx,dy = mapparams
	nx_smaller = ny_smaller = 60
	s, e = int( (nx-nx_smaller)/2 ), int( (nx+nx_smaller)/2 )
	simmapparams_smaller = [nx_smaller, ny_smaller, dx, dy]
	beta_model_smaller = sims.fn_get_Arnaud_ymap_from_Nikhel_Gupta(simmapparams_smaller, [mm], [0.7], cosmo_params_dict = param_dict, dx_finer = 0.1, mass_def = param_dict['mass_def'], rho_def = param_dict['rho_def'], sigma_logY_scatter =None)[0]
	beta_model = np.zeros( (ny, nx) )
	beta_model[s:e, s:e] = beta_model_smaller
	beta_model_smaller_90 = sims.fn_get_Arnaud_ymap_from_Nikhel_Gupta(simmapparams_smaller, [mm], [0.7], cosmo_params_dict = param_dict, dx_finer = 0.1, mass_def = param_dict['mass_def'], rho_def = param_dict['rho_def'],nu =90e9, sigma_logY_scatter =None)[0]
	beta_model_90 =  np.zeros( (ny, nx) )
	
	beta_model_90[s:e, s:e] = beta_model_smaller_90
	beta_model_90 = np.fft.ifft2( np.fft.fft2(beta_model_90 ) * Bl ).real
	beta_model = np.fft.ifft2( np.fft.fft2(beta_model) * Bl ).real
	sub_model = beta_model_90- beta_model
		
	model_1d = sims.fn_radial_profile(sub_model,  RADEC, bin_size=1, minbin=0.0, maxbin=10)[:,1]
	return sub_model
simmapparams = 200,200,0.5,0.5


Bl = sims.fn_beam_stuff(simmapparams,use_beam = 1,exp_beam = 1.0,return_beam = 1)[0]
tsz_arr = []
data_mass= 1.20*1e14

mapparams = 200, 200, 0.5,0.5
data_model = fn_get_model(data_mass, param_dict, mapparams)

from IPython import embed;embed()


for i in range(1000):
	nn = 3. *np.sqrt(2)
	noise = sims.fn_get_noise(mapparams, expnoiselevel =[nn])
	noise_data = data_model + noise[0]
	tsz_arr.append(noise_data)
	fwhm_x,fwhm_y = 2,2 
	sigma_x = fwhm_x / np.sqrt(8. * np.log(2.))
	sigma_y = fwhm_y / np.sqrt(8. * np.log(2.))
	fit_size,epsilon,offset, = 5,0.00001,0
	param_bounds = ([-400,fit_size -1,fit_size -1, sigma_x-epsilon,sigma_y-epsilon, -10],[100,fit_size +1,fit_size+1, sigma_x +epsilon+2,sigma_y +epsilon+2, 10]) 
	initial_guess = amp, x0,y0, sigma_x,sigma_y,offset
tsz_arr  = np.asarray(tsz_arr)
# calculate cov matrix

rf= []
no_of_sims = 1000
for i in range(no_of_sims):
	RADPROFILES = sims.fn_radial_profile(tsz_arr[i], RADEC, bin_size=1, minbin=0.0, maxbin=10)[:,1]
	rf.append(RADPROFILES)	
xx = np.asarray(rf)
xx = xx.T
cov_mat = np.cov(xx,xx)[10:, 0:10]
FULL_COV = cov_mat *(no_of_sims-1)
cov_mat  = FULL_COV
minA, maxA, delA = 0, 8.0, 0.1		
mass_arr = np.arange(minA,maxA+delA,delA) * 1e14
det_arr = []

model = {}

lx, ly = sims.get_lxly(mapparams)
l2d = np.sqrt(lx**2. + ly**2.)
below_l_g_cut = np.where(l2d < 1000)
above_l_g_cut = np.where(l2d >14000)
#data_fft = np.fft.fft2(data_model)
#data_fft[below_l_g_cut] = 0
#data_fft[above_l_g_cut] = 0
#data_model = np.fft.ifft2(data_fft).real


for i,arr in enumerate(tsz_arr[400:800]):

	tsz_fft = np.fft.fft2(tsz_arr[i])
	tsz_fft[below_l_g_cut] = 0
	tsz_fft[above_l_g_cut] = 0
	data_2d = np.fft.ifft2(tsz_fft).real
 	data = sims.fn_radial_profile(data_2d, RADEC, bin_size=1, minbin=0.0, maxbin=10)[:,1]
	d = data.flatten()
	kappa_COV = np.mat(cov_mat)
	Cinv = sc.linalg.pinv2(kappa_COV)
	correction_factor = (float(no_of_sims) - len(d) - 1.)/no_of_sims
	Cinv = Cinv * correction_factor
	logLarr =[]
	for bb, mm in enumerate(mass_arr):
		
		try:
			model_1d = model[mm] 
		except:
			nx_smaller = ny_smaller = 60
			s, e = int( (nx-nx_smaller)/2 ), int( (nx+nx_smaller)/2 )
			simmapparams_smaller = [nx_smaller, ny_smaller, dx, dy]
			beta_model_smaller = sims.fn_get_Arnaud_ymap_from_Nikhel_Gupta(simmapparams_smaller, [mm], [0.7], cosmo_params_dict = param_dict, dx_finer = 0.1, mass_def = param_dict['mass_def'], rho_def = param_dict['rho_def'], random_seed = bb, sigma_logY_scatter =None)[0]
			beta_model = np.zeros( (ny, nx) )
			beta_model[s:e, s:e] = beta_model_smaller
			beta_model_smaller_90 = sims.fn_get_Arnaud_ymap_from_Nikhel_Gupta(simmapparams_smaller, [mm], [0.7], cosmo_params_dict = param_dict, dx_finer = 0.1, mass_def = param_dict['mass_def'], rho_def = param_dict['rho_def'],nu =90e9,  random_seed = bb,sigma_logY_scatter =None)[0]
			beta_model_90 =  np.zeros( (ny, nx) )
		
			beta_model_90[s:e, s:e] = beta_model_smaller_90
			beta_model_90 = np.fft.ifft2( np.fft.fft2(beta_model_90 ) * Bl ).real
			beta_model = np.fft.ifft2( np.fft.fft2(beta_model) * Bl ).real
			sub_model = beta_model_90- beta_model
			sub_fft = np.fft.fft2(sub_model)
			sub_fft[below_l_g_cut] = 0
			sub_fft[above_l_g_cut] = 0
			sub_model = np.fft.ifft2(sub_fft).real
			model_1d = sims.fn_radial_profile(sub_model,  RADEC, bin_size=1, minbin=0.0, maxbin=10)[:,1]
			model[mm] = model_1d
		#from IPython import embed;embed()
		diff_arr = d - model_1d
		logLval = -0.5 * np.asarray( np.dot(diff_arr.T, np.dot( Cinv, diff_arr ))).squeeze()
		logLarr.append( logLval )
		
	

	det_val = sims.fn_get_det_significance(mass_arr,logLarr)
	det_arr.append(det_val)
	tmp = logLarr - max(logLarr); L = np.exp(tmp); L = L/max(L)
	
res = {}
res['det_arr'] = det_arr
pickle.dump(res,gzip.open('sz_det_lkhd_%s_filtered_%s_sims.pkl.gz'%(data_mass,no_of_sims),'w'))


