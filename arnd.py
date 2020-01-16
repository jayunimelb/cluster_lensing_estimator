"""
Checking arnaud profile fitting

"""
import numpy as np
import modules
scl_cmb = modules.scl_cmb
nx_smaller = ny_smaller = 60
dx = dy = 0.5
simmapparams_smaller = [nx_smaller, ny_smaller, dx, dy]
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
sims._ini_params(param_dict)
bb = 0
masss_val = 2e14
beta_model_smaller = sims.fn_get_Arnaud_ymap_from_Nikhel_Gupta(simmapparams_smaller, [masss_val], [0.7], cosmo_params_dict = param_dict, dx_finer = 0.5, mass_def = param_dict['mass_def'], rho_def = param_dict['rho_def'], random_seed = bb)[0]