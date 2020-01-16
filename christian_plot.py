import numpy as np
import modules,sys
scl_cmb = modules.scl_cmb

sims = scl_cmb.simulations()

mapparams = 100, 100, 0.25,0.25

nx,ny,dx,dy = mapparams
boxsize = nx*dx
clra,cldec = 0,0
minval, maxval = clra-boxsize/2/60.,  clra+boxsize/2/60.

ra = np.linspace(minval, maxval, nx)
dec = np.linspace(minval, maxval, ny)
RA, DEC = np.meshgrid(ra,dec)
clus_ra = [clra]
clus_dec = [cldec]
mass,redshift = float(sys.argv[1]), float(sys.argv[2]) 
M_200_cl = [mass*1e14]
z_L_cl = [redshift]

paramfile  = 'plnk_sz_grad_maps.txt'
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
c_200_cl = [sims.c_Duffy(M,z, param_dict['h']) for (M,z) in zip(M_200_cl, z_L_cl)]
from pylab import *
mass = np.arange(1,5.5,0.5)
mass = mass *1e14	
RADEC = [RA,DEC]
import matplotlib.pyplot as plt
jet= plt.get_cmap('jet')
colors = iter(jet(np.linspace(0,1,10)))
redshift =np.arange(0.3, 1.1,0.1)
"""
for i,mm in enumerate(mass):
	M_200_cl = [mm]
	c_200_cl = [sims.c_Duffy(M,z, param_dict['h']) for (M,z) in zip(M_200_cl, z_L_cl)]
	sims.fn_lensing_ini(mapparams,param_dict, RA, DEC, clus_ra, clus_dec, M_200_cl, c_200_cl, z_L_cl, param_dict['z_lss'], param_dict['mass_def'], param_dict['rho_def'], truncate_kappa_at_radius = param_dict['truncate_kappa_at_radius'])
	kp = sims.KAPPA
	rpf = sims.fn_radial_profile(kp, RADEC, bin_size = 0.5)
	plt.plot(rpf[:,0],rpf[:,1]/rpf[:,1].max(),label = 'mass:%s'%(mm), color = next(colors))
"""
M_200_cl = [3.5 *1e14]
for i,zz in enumerate(redshift):
	z_L_cl = [zz]
	c_200_cl = [sims.c_Duffy(M,z, param_dict['h']) for (M,z) in zip(M_200_cl, z_L_cl)]
	sims.fn_lensing_ini(mapparams,param_dict, RA, DEC, clus_ra, clus_dec, M_200_cl, c_200_cl, z_L_cl, param_dict['z_lss'], param_dict['mass_def'], param_dict['rho_def'], truncate_kappa_at_radius = param_dict['truncate_kappa_at_radius'])
	kp = sims.KAPPA
	rpf = sims.fn_radial_profile(kp, RADEC, bin_size = 0.5)
	plt.plot(rpf[:,0],rpf[:,1]/rpf[:,1].max(),label = 'redshift:%s'%(zz), color = next(colors))

legend = plt.legend(loc='upper right', shadow=True, fontsize='small')
t = ('Fixed mass of $3.5 * 10^{14} M_{\odot}$')
plt.text(4,0.8, t)

from IPython import embed;embed()