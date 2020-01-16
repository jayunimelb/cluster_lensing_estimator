


import numpy as np
import modules
scl_cmb = modules.scl_cmb
sims = scl_cmb.simulations()
sims.T_cmb = 2.73

# read params file
paramfile = 'params_150GHz_sims_no_mask.txt'
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
param_dict['T_cmb'] = 2.73
sims._ini_params(param_dict)




nx,ny,dx,dy = 200,200,0.5,0.5
simmapparams= mapparams = nx,ny,dx,dy 
#!! adds_tSZ
M_200_cl, z_L_cl  = [3e14], [0.7]
nx_smaller = ny_smaller = 60
##print (nx-nx_smaller)/2
s, e = int( (nx-nx_smaller)/2 ), int( (nx+nx_smaller)/2 )
simmapparams_smaller = [nx_smaller, ny_smaller, dx, dy]
beta_model_smaller = sims.fn_get_Arnaud_ymap_from_Nikhel_Gupta(simmapparams_smaller, M_200_cl, z_L_cl, cosmo_params_dict = param_dict, dx_finer = 0.1, mass_def = param_dict['mass_def'], rho_def = param_dict['rho_def'])[0]
beta_model = np.zeros( (ny, nx) )
beta_model[s:e, s:e] = beta_model_smaller

## add noise to it for now white noise
expnoiselevel = [2.]
NOISE = sims.fn_get_noise(mapparams, expnoiselevel)
final_map = NOISE[0] + beta_model
#dxdy = dx * sims.arcmins2radians * dy  * sims.arcmins2radians	
#beta_model_lmat =  dxdy*np.fft.fft2(beta_model)
#beta_model_lmat =beta_model_lmat/np.max(beta_model_lmat)
#from IPython import embed;embed()
final_map = final_map/1e6
# calculate the noise power
noise_level = expnoiselevel[0]
DeltaT = noise_level*1e-6*(1./60.)*(np.pi/180.) 
lx, ly = sims.get_lxly(mapparams)
l2d = np.sqrt(lx**2. + ly**2.)
N_lmat = (DeltaT**2.) + np.zeros(l2d.shape)
Bl_for_filter = sims.fn_beam_stuff(mapparams, exp_beam = 2., return_beam = 1)
Bl_for_filter = Bl_for_filter**2
dxdy = dx * sims.arcmins2radians * dy  * sims.arcmins2radians
tsz_lmat = dxdy*np.fft.fft2(final_map)
weight_inv_var_lmat = 1.# /((N_lmat))
T_inv_var_lmat = weight_inv_var_lmat*tsz_lmat
T_inv_var_mat = (1/dxdy)*np.fft.ifft2(Bl_for_filter*T_inv_var_lmat)
from IPython import embed;embed()