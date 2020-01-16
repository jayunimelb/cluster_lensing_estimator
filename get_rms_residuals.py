import numpy as np
import pickle,gzip
import modules
import scipy.optimize as opt
scl_cmb = modules.scl_cmb

sims = scl_cmb.simulations()

mapparams = 200,200,0.5,0.5
nx,ny,dx,dy =mapparams
cen_x,cen_y = 100,100

def fn_gaussian((x,y),amplitude):
	x0, y0 = (10-1)/2., (10-1)/2.
	data =  amplitude*np.exp( - (a*((x-x0)**2)  + c*((y-y0)**2)))
	return data.ravel()

def fn_fit_gaussian(map1):
	x,y = np.meshgrid(np.arange(10), np.arange(10))

	initial_guess = map1[100,100]
	param_bounds = ([-500],[0])
	map_fit = map1[95:105,95:105]
	popt, pcov = opt.curve_fit(fn_gaussian, (x, y), map.flatten())
	data = fn_gaussian((x,y),*popt)

	return data


# convolve it by experimental beam
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


Bl_gradient_map = sims.fn_beam_stuff(mapparams, use_beam = param_dict['use_beam'], nu = 150, return_beam = 1, exp_beam = 1.7)


Sehgal_file = 'data/sehgal_2009_sims/tsz_cib_radio_ksz_extracts/tSZ_extracts_090_M200min1.3_zmin0.25_boxsize100.0am_dx0.5am.pkl.gz'
sehgal_dic = pickle.load(gzip.open(Sehgal_file))
keynames = sehgal_dic.keys()
M200m_sehgal_sims = np.asarray(keynames)[:,3]
redshift_sehgal_sims = np.asarray(keynames)[:,2]
from IPython import embed;embed()
no_of_clusters = 100
Mass = 2*1e14
zmin = 0.25
percent_mass_tolerance = 0.05
sehgal_residual = []
for  nn in range(no_of_clusters):
	mm =Mass/1e14
	mtoldel = mm * percent_mass_tolerance
	closestinds = np.where( ( M200m_sehgal_sims>=mm-mtoldel) & ( M200m_sehgal_sims<=mm+mtoldel) & (redshift_sehgal_sims>=zmin))[0]
	selkey = keynames[closestinds[np.random.randint(len(closestinds))]]
	sehgal_sz = sehgal_dic[selkey]
	sehgal_sz_smoothed = np.fft.ifft2(np.fft.fft2(sehgal_sz)*Bl_gradient_map).real
	fitted_sehgal_map =  fn_fit_gaussian(sehgal_sz_smoothed[0])
	from IPython import embed;embed()
	sehgal_residual.append(sehgal_sz_smoothed - fitted_sehgal_map)


from IPython import embed;embed()

Daisuke_file = 'daisuke_sims'



residual_cutout = []



