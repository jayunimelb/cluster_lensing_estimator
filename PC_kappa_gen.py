"""
input: takes in mass(in 1e14 solar mass), redshift of a cluster, and params file

output: lensing convergence  profile in output.pkl.gz

ex: python PC_kappa_gen.py 2 0.7
"""

import modules.scl_cmb as scl_cmb,sys
import numpy as np, pickle,gzip
clus_mass = float(sys.argv[1])*1e14
clus_redshift = float(sys.argv[2])

sims = scl_cmb.simulations()

# read the params file

paramfile = 'params/params.txt'
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

nx,ny,dx,dy = 200,200,0.5,0.5
simmapparams = [nx,ny,dx,dy]
boxsize = 100
clra, cldec = 0., 0.
minval, maxval = clra-boxsize/2/60.,  clra+boxsize/2/60.
ra = np.linspace(minval, maxval, nx)
minval, maxval =cldec-boxsize/2/60.,  cldec+boxsize/2/60.
dec = np.linspace(minval, maxval, ny)
RA, DEC = np.meshgrid(ra,dec)
clus_ra = [clra]
clus_dec = [cldec]
M_200_cl = [clus_mass]
z_L_cl = [clus_redshift]
c_200_cl = [sims.c_Duffy(M,z, param_dict['h']) for (M,z) in zip(M_200_cl, z_L_cl)]

sims.fn_lensing_ini(simmapparams,param_dict, RA, DEC, clus_ra, clus_dec, M_200_cl, c_200_cl, z_L_cl, param_dict['z_lss'], param_dict['mass_def'], param_dict['rho_def'], truncate_kappa_at_radius = param_dict['truncate_kappa_at_radius'])
from IPython import embed;embed()
pickle.dump(sims.KAPPA,gzip.open('output.pkl.gz','w'))
