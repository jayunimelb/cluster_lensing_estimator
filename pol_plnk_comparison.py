
import numpy as np
import pickle,gzip, healpy as H
import sys
import numpy as np, healpy as hp, pickle,gzip
import modules,sys
from pylab import *
import scipy.ndimage as ndimage

def rotate_tf(tf_map_center, ra0, dec0, ra, dec, tf_rot_angle = None):
	'''
	Based on code by Srini Raghunathan
	tf_map_center is transfer funciton evaluated at the center of the map
	ra0, dec0 specify the map center
	ra, dec are the coordinates at which we want to evaluate the TF
	'''
	if tf_rot_angle == None:
		tf_rot_angle = get_tf_rot_angle(ra,dec,ra0,dec0) #radians
	tf_rotated = ndimage.interpolation.rotate(tf_map_center, tf_rot_angle*(180./np.pi), reshape = False, mode = 'nearest')

	return tf_rotated

def get_tf_rot_angle(ra_orig, dec_orig, ra0_orig, dec0_orig):
	"""
	Author: Srini Raghunathan
	Modified: Eric Baxter

	this is only valid for ZEA (proj 5) projection
	http://arxiv.org/pdf/1111.7245v1.pdf - section 5.4 for reference
	(although note that formulae therein seem to have an error)
	this is roughly the angle required to convert the fourier space filtering to real space filtering
	Should match output of get_proj5_angle.pro in the SPT analysis repo
	"""

	ra_val = np.radians(ra_orig)
	dec_val = np.radians(dec_orig)

	#do the 90deg. subtraction for declinations
	dec_val = np.radians(90.) - dec_val
	dec0 = np.radians(90.) - np.radians(dec0_orig)
	ra0 = np.radians(ra0_orig)

	#nr - numerator; dr - denominator; _1 - first term; _2 - second term, etc. in the equation.
	gamma_dr = 1 + ( np.cos(dec0) * np.cos(dec_val) ) + ( np.sin(dec0) * np.sin(dec_val) * np.cos(ra0 - ra_val) )
	gamma = 0.5 / gamma_dr

	A_1 = gamma * np.sin(dec0) * np.sin(ra_val - ra0) * ( ( np.sin(dec0)*np.cos(dec_val) ) - ( np.sin(dec_val)*np.cos(dec0)*np.cos(ra_val - ra0) ) )
	A_2 = np.cos(dec0) * np.sin(ra_val-ra0)
	A = A_1 + A_2

	B_1 = gamma * np.sin(dec0) * (np.sin(ra_val-ra0))**2. * np.sin(dec_val)
	B_2 = np.cos(ra_val-ra0)

	B = B_1 + B_2

	alpha = np.arctan(-1.*A/B)
	return alpha #in radians

def fn_pick_cutout(x,y,MAPS,cutout_size,map_resol,MAPSHAPE, perform_checks = 0):

	nx = ny = cutout_size/map_resol
	x1,x2 = int(x-nx/2),int(x+nx/2)
	y1,y2 = int(y-ny/2),int(y+ny/2)

	CUTOUT = MAPS[y1:y2,x1:x2]

	return CUTOUT


sky = modules.sky_local
sims = modules.scl_cmb.simulations()
d1 = pickle.load(gzip.open('data/map_150.pkl.gz_cutouts_150ghz_no_downsampling_y3_v6.4.21_redmapper_clusters_vl'))
cts = d1['cutouts']
keys = cts.keys()
smica_beam, sptpol_beam = 5., 1.2


simmapparams = 200,200,0.5,0.5
plnk_map_file = 'data_sp/sptsz/data/LGMCA/WPR2_CMB_muK.fits'
lgmca_map = hp.read_map(plnk_map_file,verbose =0)

lgmca_map_hr = hp.pixelfunc.ud_grade(lgmca_map,nside_out = 2048 *4)
Bl_sptpol = sims.fn_beam_stuff(simmapparams, use_beam = 2, return_beam = 1)
Bl_sptpol = sims.fn_beam_stuff(simmapparams, use_beam = 1, exp_beam = sptpol_beam, return_beam = 1)[0]
Bl_planck = sims.fn_beam_stuff(simmapparams, use_beam = 1, exp_beam = smica_beam, return_beam = 1)[0]
Bl_ratio = Bl_planck/Bl_sptpol
dx = 0.5
xs_cutout = 200
TWODTF = sims.fn_get_EBAX_2016_anal(simmapparams, l1=500, l2=400, l3=15000)[0]
TF = sims.fn_get_HPF(simmapparams, minel=8, maxel=2048, ideal=0, all_iso=1)[0]
i = 0
ra0,dec0 = 180.,0.
nside = 2048
temp = np.zeros((xs_cutout,xs_cutout))
while i  == 0:

	from IPython import embed;embed()
	pol_map = cts[keys[i]][0]
	
	ra, dec = keys[i][0],keys[i][1]
	sptpol_gnom = H.gnomview(lgmca_map, rot = [ra, dec,180.], xsize = xs_cutout, coord = ['G','C'], reso = dx, return_projected_map = 1)
	close()
	tf_rot_angle = get_tf_rot_angle(ra,dec,ra0,dec0)
	tf_rot_angle *= -1
	sptpol_gnom = rotate_tf(sptpol_gnom, ra0, dec0, ra, dec, tf_rot_angle = tf_rot_angle)
	boxsize = xs_cutout 
	map_cen = [0.0, -57.5]
	ra_dec_cen = np.asarray([ra,dec])
	for k in range(boxsize):
		for j in range(boxsize):
			ra,dec = sky.pix2Ang([j,k], ra_dec_cen, dx, temp.shape, proj=0, wrap=True)
			l,b = H.rotator.euler(ra,dec,1)
			PP=H.pixelfunc.ang2pix(nside*4,np.radians(90.-b),np.radians(l))
			temp[j,k] = lgmca_map_hr[PP]
	temp = temp/1e6
	sptpol_gnom = sptpol_gnom/1e6
	map1= np.fft.ifft2( np.fft.fft2(pol_map) * Bl_ratio * TWODTF*TF).real
	map2 = np.fft.ifft2( np.fft.fft2(sptpol_gnom) * TWODTF* TF).real
	map3 = np.fft.ifft2( np.fft.fft2(temp) * TWODTF* TF).real

	subplot(132);imshow(map1*1e6,vmin =-100,vmax = 100);colorbar();title('')#;show()
	subplot(131);imshow(map2*1e6,vmin =-100,vmax = 100);colorbar()#;show()
	subplot(133);imshow(map3*1e6,vmin =-100,vmax = 100);colorbar();title('my version');show()
	i =0
