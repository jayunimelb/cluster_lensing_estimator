

import numpy as np
import sys, pickle, gzip, os, time, glob, numpy as np, argparse
import modules
sys.path.append('modules')

import scl_cmb
sims = modules.scl_cmb.simulations()

paramfile = 'plnk_sz_grad_maps.txt'
params = np.recfromtxt(paramfile,usecols=[0],delimiter = '=')
paramvals = np.recfromtxt(paramfile,usecols=[1],delimiter = '=')
param_dict = {}
masss_val = 2*1e14
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

nx_smaller = ny_smaller = 60
nx, ny,dx,dy = 200, 200,0.5,0.5
s, e = int( (nx-nx_smaller)/2 ), int( (nx+nx_smaller)/2 )

for bb in range(1000):
	simmapparams_smaller = [nx_smaller, ny_smaller, dx, dy]
	beta_model_smaller = sims.fn_get_Arnaud_ymap_from_Nikhel_Gupta(simmapparams_smaller, [masss_val], [0.7], cosmo_params_dict = param_dict, dx_finer = 0.1, mass_def = param_dict['mass_def'], rho_def = param_dict['rho_def'], random_seed = bb, sigma_logY_scatter =None)[0]
	beta_model = np.zeros( (ny, nx) )
	beta_model[s:e, s:e] = beta_model_smaller
	beta_model_smaller_90 = sims.fn_get_Arnaud_ymap_from_Nikhel_Gupta(simmapparams_smaller, [masss_val], [0.7], cosmo_params_dict = param_dict, dx_finer = 0.1, mass_def = param_dict['mass_def'], rho_def = param_dict['rho_def'],nu =90e9,  random_seed = bb,sigma_logY_scatter =None)[0]
	beta_model_90 =  np.zeros( (ny, nx) )
	beta_model_90[s:e, s:e] = beta_model_smaller_90
	noise_level = 2*np.sqrt(2)
	noise_sim = 

from IPython import embed;embed()