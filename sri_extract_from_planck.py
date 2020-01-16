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
########################################################################
########################################################################
########################################################################

import matplotlib; matplotlib.use('Agg')
import healpy as H, numpy as numpy, os, sys, glob, pickle, gzip, argparse
import modules.scl_cmb as scl_cmb
import modules.sky_local as sky_local
sims = scl_cmb.simulations()
from pylab import *
import scipy.ndimage as ndimage

########################################################################
### read  planck map
print '\n\t\t read  planck map\t'
fname = 'data/LGMCA/WPR2_CMB_muK.fits'
fname = 'data_sp/sptsz/data/LGMCA/WPR2_CMB_muK.fits'
planck_healpix_MAP = H.read_map(fname, verbose = 0)
nside = 2048
########################################################################


parser = argparse.ArgumentParser(description='')
parser.add_argument('-which_rand', dest='which_rand', action='store', help='which_rand', type=str, default=None)
parser.add_argument('-redmapper', dest='redmapper', action='store', help='redmapper', type=int, default=1)
parser.add_argument('-sptsz', dest='sptsz', action='store', help='sptsz', type=int, default=0)
parser.add_argument('-cutout_size', dest='cutout_size', action='store', help='cutout_size', type=float, default=100.)
parser.add_argument('-dx', dest='dx', action='store', help='dx', type=float, default=0.5)

args = parser.parse_args()

args_keys = args.__dict__
for kargs in args_keys:
	param_value = args_keys[kargs]
	if isinstance(param_value, str):
		cmd = '%s = "%s"' %(kargs, param_value)
	else:
		cmd = '%s = %s' %(kargs, param_value)
	exec(cmd)

'''
print '\n\t\t extract planck map on sptpol patch\t'
boxsize, dx = 2500.,0.5#2deg; 0.5arcmin pixelsize
xs = boxsize/dx
ra0, dec0 = 0.0, -57.5
#ra0, dec0 = 0.0, -59.033333
map_centre = np.asarray( [ra0, dec0] )
planck_MAP_EX = H.gnomview(planck_healpix_MAP, rot = [ra0,dec0,180.], xsize = xs, coord = ['G', 'C'], reso = dx, return_projected_map = 1)
close()
mapshape = planck_MAP_EX.shape
'''
########################################################################
testing = 1
proc_sptpol_map = 1
if proc_sptpol_map and testing:
	fname = 'data/CR_maps_20170910/final_maps/20180124_downsample_x2_no_masking_v2/RADEC.pkl.gz'
	mapname = 'data/CR_maps_20170910/final_maps/20180124_downsample_x2_no_masking_v3/map_150_dsfactor2.pkl.gz'
	maskname = 'data/sanjay_maps_201705xx/total_mask.pkl.gz'
	fname = 'data_sp/sptsz/data/spt/catalog.pkl.gz'


	RADEC = pickle.load(gzip.open(fname, 'rb'))
	SPTpolmap = pickle.load(gzip.open(mapname, 'rb'))[0]
	mask = pickle.load(gzip.open(maskname, 'rb'))
	mask = sims.downsample_map(mask, 2)
	SPTpolmap *= mask

	RA, DEC = RADEC
	nopixels = H.nside2npix(nside)
	sptpol_healpix_MAP=np.zeros(nopixels)
	HIT=np.zeros(nopixels)
	for r in range(np.shape(RA)[0]):
		#l, b = H.rotator.euler(RA[r,:],DEC[r,:],1) #radec to lb
		#PP=H.pixelfunc.ang2pix(nside,np.radians(90.-b),np.radians(l))
		PP=H.pixelfunc.ang2pix(nside,np.radians(90.-DEC[r,:]),np.radians(RA[r,:]))
		sptpol_healpix_MAP[PP]+=SPTpolmap[r,:]
		HIT[PP]+=1
	sptpol_healpix_MAP[HIT>0.]/=HIT[HIT>0.]
	sptpol_healpix_MAP *= 1e6
	'''
	### H.mollview(MAP);show()
	ra0_sptpol, dec0_sptpol = 0.0, -57.5
	SPTpol_MAP_EX = H.gnomview(sptpol_healpix_MAP, rot = [ra0_sptpol, dec0_sptpol,180.], xsize = xs, coord = ['C'], reso = dx, return_projected_map = 1)
	close()
	'''
	## subplot(121);imshow(planck_MAP_EX);colorbar();
	## subplot(122);imshow(SPTpol_MAP_EX);colorbar();show()
########################################################################

########################################################################
print '\n\t\t get TF and beams for comparison with SPTpol map\t'
xs_cutout = cutout_size/dx
simmapparams = [int(xs_cutout), int(xs_cutout), dx, dx]
params_file = '/home/sri/analysis/2016_11/QE_SPTpol/params/crossmaps/sptpol_map_like_sims/params_planck_r_0.0_2015_cosmo_lensed_LSS_crossmap_sptpol_FG.txt'
params_file = 'params_planck_r_0.0_2015_cosmo_lensed_LSS_crossmap_sptpol_planck_SZmaps.txt'
param_dict = sims.fn_get_param_dict_from_paramfile(params_file)
TWODTF = sims.fn_get_EBAX_2016_anal(simmapparams, l1=param_dict['l1'], l2=param_dict['l2'], l3=param_dict['l3'])[0]
#TF = sims.fn_get_EBAX_2016_anal(simmapparams, l1=300., l2=param_dict['l2'], l3=param_dict['l3'])[0]
lmin,lmax = 250, 2000
TF = sims.fn_get_HPF(simmapparams, minel=lmin, maxel=lmax, ideal=0, all_iso=1)[0]

smica_beam, sptpol_beam = 5., 1.2
Bl_planck = sims.fn_beam_stuff(simmapparams, use_beam = 1, exp_beam = smica_beam, return_beam = 1)[0]
Bl_sptpol = sims.fn_beam_stuff(simmapparams, use_beam = 1, exp_beam = sptpol_beam, return_beam = 1)[0]
Bl_ratio = Bl_planck/Bl_sptpol

lmin,lmax = 350, 5000
TF = sims.fn_get_HPF(simmapparams, minel=lmin, maxel=lmax, ideal=0, all_iso=1)[0]
########################################################################

print '\n\t\t read RM cluster catalogue\t'
if redmapper:
	cutouts_file = 'data/CR_maps_20170910/final_maps/20180124_downsample_x2_no_masking_v4/map_150.pkl.gz_cutouts_150ghz_no_downsampling_y3_v6.4.21_lgt5_redmapper_clusters_vl'
elif sptsz:
	 cutouts_file = 'data_sp/sptsz/data/spt/sz_cutouts.pkl.gz'
	
else:
	cutouts_file = 'data/CR_maps_20170910/final_maps/20180124_downsample_x2_no_masking_v4/map_150.pkl.gz_cutouts_150ghz_no_downsamplingsptpol_clusters'

if which_rand<>None:
	cutouts_file = '%s_randoms_%s' %(cutouts_file.replace('y3_v6.4.21_lgt5','year3_lgt5').replace('_clusters_vl',''), which_rand)


cutouts = pickle.load(gzip.open(cutouts_file, 'rb'))['cutouts']
cutout_keys = cutouts.keys()
raRMclus, decRMclus = np.asarray(cutout_keys)[:,0], np.asarray(cutout_keys)[:,1]


ra0,dec0 = 180.,0. #healpix map centre
LGMCA_DIC = {}
tol_frac_zero_inds = 0.01
calib_fac_150 = 0.9097
minval, maxval = -100, 100.


if testing: 
	total = 0
	minrich = 80.

switch_backend('Agg')
for kcnt, (r,d, key) in enumerate(zip(raRMclus, decRMclus, cutout_keys)):
	
	print kcnt, key

	SPTPOL_EX = cutouts[key][0] * 1e6 * calib_fac_150

	if testing:
		if key[3]<minrich: continue

		SPTPOL_EX_flatten = SPTPOL_EX.flatten()
		zero_inds = np.where( (abs(SPTPOL_EX_flatten) == 0.) )[0]
		frac_zero_inds = len(zero_inds) * 1. / len(SPTPOL_EX_flatten.flatten())
		if frac_zero_inds > tol_frac_zero_inds: continue #more than xx% zeros - ignore
	r = r.astype(np.float)
	d = d.astype(np.float)
#	from IPython import embed;embed()
	tf_rot_angle = get_tf_rot_angle(r,d,ra0,dec0)

	tf_rot_angle *= -1
	planck_gnom = H.gnomview(planck_healpix_MAP, rot = [r, d,180.], xsize = xs_cutout, coord = ['G','C'], reso = dx, return_projected_map = 1)
	close()
	#planck_gnom = rotate_tf(planck_gnom, ra0, dec0, r, d, tf_rot_angle = tf_rot_angle)
	if testing:
		MAP1 = np.fft.ifft2( np.fft.fft2(planck_gnom) * TF * TWODTF).real

		sptpol_gnom = H.gnomview(sptpol_healpix_MAP, rot = [r, d,180.], xsize = xs_cutout, coord = ['C'], reso = dx, return_projected_map = 1)
		close()
		sptpol_gnom = rotate_tf(sptpol_gnom, ra0, dec0, r, d, tf_rot_angle = tf_rot_angle)
		MAP2 = np.fft.ifft2( np.fft.fft2(sptpol_gnom) * Bl_ratio * TF).real

		MAP3 = np.fft.ifft2( np.fft.fft2(SPTPOL_EX) * Bl_ratio * TF).real

	LGMCA_DIC[key] = [planck_gnom] #[planck_gnom, sptpol_gnom, MAP1, MAP2]

	if testing:
		subplot(221);
		imshow(MAP1, vmin = minval, vmax = maxval, origin = 'lower');colorbar()
		contour(MAP2);title(key, fontsize = 10);
		subplot(222); imshow(MAP3, vmin = minval, vmax = maxval, origin = 'lower');colorbar();#title(tfa)
		subplot(223); imshow(MAP1, vmin = minval, vmax = maxval, origin = 'lower');colorbar();#title(tfa)
		subplot(224); imshow(MAP2, vmin = minval, vmax = maxval, origin = 'lower');colorbar()
		show()

		total += 1

		if total > 10: break

if redmapper:
	opf = 'data/CR_maps_20170910/final_maps/20180124_downsample_x2_no_masking_v4/LGMCA.pkl.gz_cutouts_y3_v6.4.21_lgt5_redmapper_clusters_vl'
elif sptsz:
	opf = 'data_sp/sptsz/data/LGMCA/sir_lgmca_cutout.pkl.gz'
else:
	opf = 'data/CR_maps_20170910/final_maps/20180124_downsample_x2_no_masking_v4/LGMCA.pkl.gz_cutoutssptpol_clusters'

if which_rand<>None:
	opf = '%s_randoms_%s' %(opf, which_rand)

print opf

pickle.dump(LGMCA_DIC, gzip.open(opf, 'wb'), protocol = 2)
sys.exit()


"""
everything else below is archived
"""

























for kcnt, (r,d, key) in enumerate(zip(raRMclus, decRMclus, cutout_keys)):

	if key[3]<100: continue

	SPTPOL_EX = cutouts[key][0] * 1e6 * calib_fac_150
	SPTPOL_EX_flatten = SPTPOL_EX.flatten()
	zero_inds = np.where( (abs(SPTPOL_EX_flatten) == 0.) )[0]
	frac_zero_inds = len(zero_inds) * 1. / len(SPTPOL_EX_flatten.flatten())
	if frac_zero_inds > tol_frac_zero_inds: continue #more than xx% zeros - ignore

	pixels = sky_local.ang2Pix([r,d], map_centre, dx, mapshape, proj = 5)
	y, x = pixels

	planck_gnom = fn_pick_cutout(x,y,planck_MAP_EX,cutout_size,dx, mapshape) #removed checks on 20171204
	#planck_gnom = rotate_tf(planck_gnom, ra0, dec0, r, d)
	planck_gnom = np.fft.ifft2( np.fft.fft2(planck_gnom) * TF * TWODTF).real
	MAP1 = np.fft.ifft2( np.fft.fft2(planck_gnom) * TF * TWODTF).real

	tf_rot_angle = get_tf_rot_angle(r,d,ra0,dec0)
	sptpol_gnom = fn_pick_cutout(x,y,SPTpol_MAP_EX,cutout_size,dx, mapshape) #removed checks on 20171204
	## sptpol_gnom = rotate_tf(sptpol_gnom, ra0, dec0, r, d, tf_rot_angle = tf_rot_angle)
	MAP2 = np.fft.ifft2( np.fft.fft2(sptpol_gnom) * Bl_ratio * TF).real
	MAP3 = np.fft.ifft2( np.fft.fft2(SPTPOL_EX) * Bl_ratio * TF).real

	subplot(221);
	imshow(MAP1, vmin = minval, vmax = maxval, origin = 'lower');colorbar()
	contour(MAP2);title(key, fontsize = 10);
	subplot(222); imshow(MAP3, vmin = minval, vmax = maxval, origin = 'lower');colorbar();#title(tfa)
	subplot(223); imshow(MAP1, vmin = minval, vmax = maxval, origin = 'lower');colorbar();#title(tfa)
	subplot(224); imshow(MAP2, vmin = minval, vmax = maxval, origin = 'lower');colorbar()
	show()
	sys.exit()

	subplot(121);imshow(MAP1, vmin = minval, vmax = maxval, interpolation = 'bicubic');colorbar()
	subplot(122);imshow(MAP2, vmin = minval, vmax = maxval, interpolation = 'bicubic');colorbar()
	title(key)
	show();sys.exit()

	if kcnt > 25: break
sys.exit()


for kcnt, (r,d, key) in enumerate(zip(raRMclus, decRMclus, cutout_keys)):	
	'''

	if 1==1:#planck_map:
		planck_gnom = fn_pick_cutout(x,y,planck_MAP_EX,cutout_size,dx, mapshape) #removed checks on 20171204
		planck_gnom = rotate_tf(planck_gnom, ra0, dec0, r, d)

		minval, maxval = -100, 100.
		MAP1 = np.fft.ifft2( np.fft.fft2(planck_gnom) * TF * TWODTF).real
		MAP2 = np.fft.ifft2( np.fft.fft2(SPTPOL_EX) * Bl_ratio * TF).real
		#subplot(221);imshow(MAP1, vmin = minval, vmax = maxval);colorbar()
		#subplot(222);imshow(MAP2, vmin = minval, vmax = maxval);colorbar()


	else:

		print 'extracting from sptpol'
		sptpol_gnom = fn_pick_cutout(x,y,SPTpol_MAP_EX,cutout_size,dx, mapshape) #removed checks on 20171204
		### sptpol_EX = rotate_tf(sptpol_gnom, ra0, dec0, r, d)

		minval, maxval = -100, 100.
		MAP1 = np.fft.ifft2( np.fft.fft2(sptpol_gnom) * TF).real
		MAP2 = np.fft.ifft2( np.fft.fft2(SPTPOL_EX) * TF).real
		#subplot(221);imshow(MAP1, vmin = minval, vmax = maxval);colorbar()
		#subplot(222);imshow(MAP2, vmin = minval, vmax = maxval);colorbar()
	'''

	'''
	imshow(MAP1, vmin = minval, vmax = maxval, origin = 'lower');colorbar()
	contour(MAP2);title(key);show()
	#subplot(221);imshow(MAP1, vmin = minval, vmax = maxval);colorbar()
	#subplot(222);imshow(MAP2, vmin = minval, vmax = maxval);colorbar();show()
	'''


	if kcnt > 25: break

sys.exit()












##sys.exit()


#cutouts_file = 'data/sanjay_maps_201705xx/final_maps/20171027_downsample_x2_tSZfree_plus_150/map_150.pkl.gz_cutouts_150ghz_no_downsamplingredmapper_clusters'
#cutouts_file = 'data/sanjay_maps_201705xx/final_maps/20171027_downsample_x2_tSZfree_plus_150_v1/map_150.pkl.gz_cutouts_150ghz_no_downsamplingredmapper_clusters'
#cutouts_file = 'data/sanjay_maps_201705xx/final_maps/20171027_downsample_x2_tSZfree_plus_150_v1_no_masking/map_150.pkl.gz_cutouts_150ghz_no_downsamplingredmapper_clusters'
#cutouts_file = 'data/sanjay_maps_201705xx/final_maps/20171027_downsample_x2_tSZfree_plus_150_v1_no_masking/map_150.pkl.gz_cutouts_150ghz_no_downsamplingsptpol_clusters'
#cutouts_file = 'data/sanjay_maps_201705xx/final_maps/20171027_downsample_x2_tSZfree_plus_150_v2_no_masking/map_150.pkl.gz_cutouts_150ghz_no_downsamplingsptpol_clusters'

cutouts_file = 'data/CR_maps_20170910/final_maps/20180124_downsample_x2_no_masking_v3/map_150.pkl.gz_cutouts_150ghz_no_downsampling_y3_v6.4.21_lgt5_redmapper_clusters_vl'
cutouts = pickle.load(gzip.open(cutouts_file, 'rb'))['cutouts']
cutout_keys = cutouts.keys()
raRMclus, decRMclus = np.asarray(cutout_keys)[:,0], np.asarray(cutout_keys)[:,1]

'''
H.mollview(MAP, coord = ['G','C'])
H.projplot(raRMclus, decRMclus, coord = ['C'], lonlat = 1, color='r', marker = '.', ls = 'None')
show();sys.exit()
'''

### convert sptpol to helapix coords
proc_sptpol_map = 1
if proc_sptpol_map:
	fname = 'data/CR_maps_20170910/final_maps/20180124_downsample_x2_no_masking_v2/RADEC.pkl.gz'
	mapname = 'data/CR_maps_20170910/final_maps/20180124_downsample_x2_no_masking_v3/map_150_dsfactor2.pkl.gz'
	maskname = 'data/sanjay_maps_201705xx/total_mask.pkl.gz'
	RADEC = pickle.load(gzip.open(fname, 'rb'))
	SPTpolmap = pickle.load(gzip.open(mapname, 'rb'))[0]
	mask = pickle.load(gzip.open(maskname, 'rb'))
	mask = sims.downsample_map(mask, 2)
	SPTpolmap *= mask

	RA, DEC = RADEC
	nopixels = H.nside2npix(nside)
	MAP=np.zeros(nopixels)
	HIT=np.zeros(nopixels)
	for r in range(np.shape(RA)[0]):
		#l, b = H.rotator.euler(RA[r,:],DEC[r,:],1) #radec to lb
		#PP=H.pixelfunc.ang2pix(nside,np.radians(90.-b),np.radians(l))
		PP=H.pixelfunc.ang2pix(nside,np.radians(90.-DEC[r,:]),np.radians(RA[r,:]))
		MAP[PP]+=SPTpolmap[r,:]
		HIT[PP]+=1
	MAP[HIT>0.]/=HIT[HIT>0.]
	MAP *= 1e6
	H.mollview(MAP);show()
	ra0, dec0 = 0.0, -57.5
	SPTpol_MAP_EX = H.gnomview(MAP, rot = [ra0,dec0,180.], xsize = xs, coord = ['C'], reso = dx, return_projected_map = 1)
	#SPTpol_MAP_EX = H.cartview(MAP, rot = [ra0,dec0,180.], xsize = xs, coord = ['C'], return_projected_map = 1)
	close()
	imshow(SPTpol_MAP_EX);colorbar();show()
###sys.exit()



boxsize, dx = 100., 0.5#2deg; 0.5arcmin pixelsize
xs = boxsize/dx

simmapparams = [int(xs), int(xs), dx, dx]
params_file = '/Users/sraghunathan/Research/SPTPol/analysis/2016_11/QE_SPTpol/params/crossmaps/sptpol_map_like_sims/params_planck_r_0.0_2015_cosmo_lensed_LSS_crossmap_sptpol_FG.txt'
param_dict = sims.fn_get_param_dict_from_paramfile(params_file)
TWODTF = sims.fn_get_EBAX_2016_anal(simmapparams, l1=param_dict['l1'], l2=param_dict['l2'], l3=param_dict['l3'])[0]
#TF = sims.fn_get_EBAX_2016_anal(simmapparams, l1=300., l2=param_dict['l2'], l3=param_dict['l3'])[0]
lmin,lmax = 250, 2000
TF = sims.fn_get_HPF(simmapparams, minel=lmin, maxel=lmax, ideal=0, all_iso=1)[0]

smica_beam, sptpol_beam = 5., 1.2
Bl_planck = sims.fn_beam_stuff(simmapparams, use_beam = 1, exp_beam = smica_beam, return_beam = 1)[0]
Bl_sptpol = sims.fn_beam_stuff(simmapparams, use_beam = 1, exp_beam = sptpol_beam, return_beam = 1)[0]
Bl_ratio = Bl_planck/Bl_sptpol

planck_map = 1
LGMCA_DIC = {}
tol_frac_zero_inds = 0.01
calib_fac_150 = 0.9097
for kcnt, (r,d, key) in enumerate(zip(raRMclus, decRMclus, cutout_keys)):

	SPTPOL_EX = cutouts[key][0] * 1e6 * calib_fac_150
	SPTPOL_EX_flatten = SPTPOL_EX.flatten()
	zero_inds = np.where( (abs(SPTPOL_EX_flatten) == 0.) )[0]
	frac_zero_inds = len(zero_inds) * 1. / len(SPTPOL_EX_flatten.flatten())
	if frac_zero_inds > tol_frac_zero_inds: continue #more than xx% zeros - ignore

	#if kcnt >10: continue
	print r, d, key

	if planck_map:
		### l, b = H.rotator.euler(r,d,1) #radec to lb
		EX = H.gnomview(MAP, rot = [r,d,180.], xsize = xs, coord = ['G', 'C'], reso = dx, return_projected_map = 1)
		close()



		minval, maxval = -100, 100.
		MAP1 = np.fft.ifft2( np.fft.fft2(EX) * TF * TWODTF).real
		MAP2 = np.fft.ifft2( np.fft.fft2(SPTPOL_EX) * Bl_ratio * TF).real
		subplot(221);imshow(MAP1, vmin = minval, vmax = maxval);colorbar()
		subplot(222);imshow(MAP2, vmin = minval, vmax = maxval);colorbar()

	else:
		#l, b = H.rotator.euler(r,d,1) #radec to lb
		#EX = H.gnomview(MAP, rot = [l, b], xsize = xs, coord = ['G'], reso = dx, return_projected_map = 1)

		EX = H.gnomview(MAP, rot = [r, d, 180.], xsize = xs, coord = ['C'], reso = dx, return_projected_map = 1)
		close()
		minval, maxval = -100, 100.
		MAP1 = np.fft.ifft2( np.fft.fft2(EX) * TF).real
		MAP2 = np.fft.ifft2( np.fft.fft2(SPTPOL_EX) * TF).real
		subplot(221);imshow(MAP1, vmin = minval, vmax = maxval);colorbar()
		subplot(222);imshow(MAP2, vmin = minval, vmax = maxval);colorbar()

	cls_sptpol = sims.fn_plot_pow_spec(simmapparams, [MAP1])[0][0]
	cls_cross = sims.fn_plot_pow_spec(simmapparams, [MAP1, MAP2], cross_power = 1)[0]
	subplot(223);
	loglog(cls_cross[:,0], cls_cross[:,1], color = 'k');
	loglog(cls_sptpol[:,0], cls_sptpol[:,1], color = 'r');ylim(1e-5,1e-1);xlim(lmin,lmax)
	title(key, fontsize = 8)
	show()

	if kcnt>10:
		break

sys.exit()




'''
simmapparams = [int(xs), int(xs), dx, dx]
smica_beam = 5.
Bl_planck = sims.fn_beam_stuff(simmapparams, use_beam = 1, exp_beam = smica_beam, return_beam = 1)[0]
Bl_sptpol = sims.fn_beam_stuff(simmapparams, use_beam = 2, return_beam = 1)
Bl_dummy = np.copy(Bl_planck) *0. + 1.
params_file = '/Users/sraghunathan/Research/SPTPol/analysis/2016_11/QE_SPTpol/params/crossmaps/sptpol_map_like_sims/params_planck_r_0.0_2015_cosmo_lensed_LSS_crossmap_sptpol_FG.txt'
param_dict = sims.fn_get_param_dict_from_paramfile(params_file)
TWODTF = sims.fn_get_EBAX_2016_anal(simmapparams, l1=param_dict['l1'], l2=param_dict['l2'], l3=param_dict['l3'])[0]
#TWODTF = sims.fn_get_HPF(simmapparams, minel = 500., maxel = 2048)[0]


def fn_beam_deconv(MAP, Bl):

    bad = np.where(Bl==0.)
    #Bl[bad] = 1.
    Blinv = 1./Bl
    Blinv[bad] = 0.

    highly_filtered = np.where((Blinv > 1.0e8) | (Blinv < 0) )
    print highly_filtered;sys.exit()
    Blinv[highly_filtered] = 0.0

    return np.fft.ifft2( np.fft.fft2(MAP) * Blinv).real
'''

LGMCA_DIC = {}
for kcnt, (r,d, key) in enumerate(zip(raRMclus, decRMclus, cutout_keys)):
	#if kcnt >10: continue
	print r, d, key

	l, b = H.rotator.euler(r,d,1) #radec to lb
	EX = H.gnomview(MAP, rot = [l,b], xsize = xs, reso = dx, return_projected_map = 1)
	close()
	

	#EX = H.gnomview(MAP, rot = [r,d], coord = ['G','C'], xsize = boxsize, reso = dx, return_projected_map = 1)
	clf();

	SPTPOL_EX = cutouts[key][0] * 1e6
	#SPTPOL_EX = cutouts.values()[np.random.randint(len(cutouts))][0] * 1e6
	#SPTPOL_EX = sims.downsample_map(SPTPOL_EX)
	
	#make SPTpol and Planck similar
	calib_factor = 0.9#/Bl_sptpol
	#SPTPOL_EX = fn_beam_deconv(SPTPOL_EX, Bl_sptpol)
	#EX = fn_beam_deconv(EX, Bl_planck)

	#new_beam_for_sptpol = Bl_sptpol #Bl_planck/Bl_sptpol
	#new_beam_for_sptpol[new_beam_for_sptpol == np.inf] = 0.

	new_beam_for_planck = Bl_dummy #Bl_sptpol#/Bl_planck
	new_beam_for_planck[new_beam_for_planck == np.inf] = 0.

	#TF
	EX = np.fft.ifft2( np.fft.fft2(EX) * TWODTF * new_beam_for_planck).real
	SPTPOL_EX = np.fft.ifft2( np.fft.fft2(SPTPOL_EX) * TWODTF).real/calib_factor
	#SPTPOL_EX = np.fft.ifft2( np.fft.fft2(SPTPOL_EX) * TWODTF * new_beam_for_sptpol).real/calib_factor
	#SPTPOL_EX = np.fft.ifft2( np.fft.fft2(SPTPOL_EX) * TWODTF * new_beam_for_sptpol).real/calib_factor
	
	minval, maxval = -100, 100.
	subplot(221);imshow(EX, vmin = minval, vmax = maxval);colorbar()
	subplot(222);imshow(SPTPOL_EX, vmin = minval, vmax = maxval);colorbar()

	cls_planck = sims.fn_plot_pow_spec(simmapparams, [EX])[0][0]
	cls_sptpol = sims.fn_plot_pow_spec(simmapparams, [SPTPOL_EX])[0][0]
	cls_cross = sims.fn_plot_pow_spec(simmapparams, [EX, SPTPOL_EX], cross_power = 1)[0]

	ax = subplot(223, xscale = 'log', yscale = 'log')
	errorbar(cls_planck[:,0]-10., cls_planck[:,1], yerr= cls_planck[:,2], color ='r', label = 'Planck')
	errorbar(cls_sptpol[:,0], cls_sptpol[:,1], yerr= cls_sptpol[:,2], color ='g', label = 'SPTpol')
	errorbar(cls_cross[:,0]+10., cls_cross[:,1], yerr = cls_cross[:,2], color ='k', label = 'Cross');
	legend(loc=1)
	show()#;sys.exit()


	LGMCA_DIC[key] = EX
	#clf();imshow(EX, vmin = minval, vmax = maxval);colorbar(); #H.mollview(MAP, min = minval, max = maxval);H.projplot(r,d, coord = 'C', lonlat = 1, color='r', marker = 'o', ls = 'None');show()
	#show()

	#EX1 = H.gnomview(MAP, rot = [l,b], xsize = xs, reso = dx, return_projected_map = 1)
	#EX2 = H.gnomview(MAP, rot = [r,d], coord = ['G','C'], xsize = boxsize, reso = dx, return_projected_map = 1)
	#subplot(121);imshow(EX1);colorbar(); subplot(122);imshow(EX2);colorbar();show()

op_file = '%s_planckgradientcutouts' %(cutouts_file)
print op_file
pickle.dump(LGMCA_DIC, gzip.open(op_file, 'wb'), protocol = 2)
sys.exit()
