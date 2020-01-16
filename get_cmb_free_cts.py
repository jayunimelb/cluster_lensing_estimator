import numpy as np
import pickle,gzip
import modules,sys

scl_cmb = modules.scl_cmb
sims = scl_cmb.simulations()
clusterstuff = scl_cmb.cluster_stuff()
mapparams = 200,200,0.5,0.5
nx,ny,dx,dy = mapparams
import scipy.optimize as opt
def twoD_Gaussian((x, y), amplitude, xo, yo, sigma_x, sigma_y, theta, offset):
	xo = float(xo)
	yo = float(yo)    	
	a = (np.cos(theta)**2)/(2*sigma_x**2) + (np.sin(theta)**2)/(2*sigma_y**2)
	b = -(np.sin(2*theta))/(4*sigma_x**2) + (np.sin(2*theta))/(4*sigma_y**2)
	c = (np.sin(theta)**2)/(2*sigma_x**2) + (np.cos(theta)**2)/(2*sigma_y**2)
	g = offset + amplitude*np.exp( - (a*((x-xo)**2) + 2*b*(x-xo)*(y-yo) + c*((y-yo)**2)))
	return g.ravel()


beam = pickle.load(gzip.open('data_sp/sptsz/data/spt/sz_beam.pkl.gz'))
keys = beam.keys()

b150 = beam[('2009', '150')]
b90 = beam[('2010', '90')]
Bl_90 = sims.el1d_to_EL2D(b90,mapparams)
Bl_150 = sims.el1d_to_EL2D(b150,mapparams)

Bl_tSZ_for_150= Bl_90/Bl_150

Bl_tSZ_for_150[Bl_tSZ_for_150 == np.inf] = 0.
Bl_tSZ_for_150[np.where(np.isnan(Bl_tSZ_for_150))] = 0.
boxsize = 100
clra,cldec = 0,0
minval, maxval = clra-boxsize/2/60.,  clra+boxsize/2/60.
ra = np.linspace(minval, maxval, nx)
dec = np.linspace(minval, maxval, ny)
RA, DEC = np.meshgrid(ra,dec)
RADIUS = np.sqrt(RA**2 + DEC**2.)
RADIUS *= 60.
ysz = pickle.load(gzip.open('sptsz_ymap.pkl.gz_cutoutssptsz_clusters_dx0.5am_cutoutsize100.0am'))
keys = ysz.keys()
from pylab import *
for i,key in enumerate(keys):
	tsz_map= ysz[key][0]
	map_fit = tsz_map[96:104,96:104]
	amp = np.mean(tsz_map[RADIUS<=2.])
	x0, y0 = 4,4
	x,y = np.arange(8), np.arange(8)
	x, y = np.meshgrid(x, y)
	sigma_x,sigma_y = 1,1
	theta = 0
	offset = 0
	initial_guess = amp, x0, y0, sigma_x, sigma_y, theta, offset
	try:
		popt, pcov = opt.curve_fit(twoD_Gaussian, (x, y), map_fit.flatten(), p0=initial_guess)
	except:
		p0 =initial_guess
		popt = p0
	data_fitted = twoD_Gaussian((x, y), *popt)
	fit = data_fitted.reshape(8,8)
	full_fit = np.zeros((nx,ny))
	full_fit[96:104,96:104] = fit

	if i >10 and i <15:
		subplot(121);imshow(map_fit);colorbar()
		subplot(122);imshow(fit);colorbar();show()

		#imshow(-2.5*temp[90:110,90:110]);colorbar;show()


sys.exit()





"""
data_150 = pickle.load(gzip.open('data_sp/sptsz/data/spt/sz_cutouts.pkl.gz'))
cuts_150 = data_150['cutouts']
data_90 = pickle.load(gzip.open('data_sp/sptsz/data/spt/sz_cutouts_90GHz.pkl.gz'))
cuts_90 = data_90['cutouts']
from pylab import *
keys = cuts_90.keys()
for i,key in enumerate(keys):
	MAP_150 = cuts_150[key][0]
	MAP_150 = np.fft.ifft2( np.fft.fft2(MAP_150) * Bl_tSZ_for_150 ).real
	MAP_90 = cuts_90[key][0]
	if i >10 and i <15:
		xx = MAP_90 - MAP_150
		subplot(131);imshow(MAP_90[90:110,90:110]);axvline(x=10);axhline(y =10);colorbar()
		subplot(132);imshow(MAP_150[90:110,90:110]);axvline(x=10);axhline(y =10);colorbar()
		subplot(133);imshow(xx[90:110,90:110]);axvline(x=10);axhline(y =10);colorbar()
		show()
"""