import numpy as np, healpy as hp, pickle,gzip
import modules,sys
from pylab import *
sky = modules.sky_local


sims = modules.scl_cmb.simulations()

plnk_map_file = 'data_sp/sptsz/data/LGMCA/WPR2_CMB_muK.fits'
lgmca_map = hp.read_map(plnk_map_file,verbose =0)
nside = 2048
lgmca_map_hi_reso = hp.pixelfunc.ud_grade(lgmca_map, nside_out = 2048*4)
new_map_reso = 0.5
nx,ny,dx,dy = 10,10,2,2
mapparams = nx,ny,dx,dy
temp_vals = np.zeros((nx,ny))
x,y = np.arange(0,20,2),np.arange(0,20,2)
Y, X = np.meshgrid(y,x)
boxsize =200
dx = 2 # arcmin resolution
temp = np.zeros((boxsize,boxsize))
try:
	ra_dec_cen = np.asarray([float(sys.argv[1]),float(sys.argv[2])])
	print "passed co-ordinates"
except:
	ra_dec_cen = np.asarray([0,0])
simmapparams = boxsize,boxsize,new_map_reso,new_map_reso
smica_beam, sptpol_beam = 5., 1.2
Bl_planck = sims.fn_beam_stuff(simmapparams, use_beam = 1, exp_beam = smica_beam, return_beam = 1)[0]
Bl_sptpol = sims.fn_beam_stuff(simmapparams, use_beam = 1, exp_beam = sptpol_beam, return_beam = 1)[0]
Bl_ratio = Bl_planck/Bl_sptpol
i =1
while i >0:
	sz_data = pickle.load(gzip.open('data_sp/sptsz/data/spt/sz_cutouts.pkl.gz'));sz_cts = sz_data['cutouts'];keys = sz_cts.keys()
	from IPython import embed; embed()
	
	ra_dec_cen = np.asarray([key[0],key[1]])
	for i in range(boxsize):
		for j in range(boxsize):
			ra,dec = sky.pix2Ang([j,i], ra_dec_cen, new_map_reso, temp.shape, proj=0, wrap=True)
			l,b = hp.rotator.euler(ra,dec,1)
			PP=hp.pixelfunc.ang2pix(nside*4,np.radians(90.-b),np.radians(l))
			temp[i,j] = lgmca_map_hi_reso[PP]

	simmapparams = boxsize,boxsize,new_map_reso,new_map_reso
	import modules;sims = modules.scl_cmb.simulations()
	TWODTF = sims.fn_get_EBAX_2016_anal(simmapparams, l1=500, l2=400, l3=15000)[0]
	TF = sims.fn_get_HPF(simmapparams, minel=8, maxel=2048, ideal=0, all_iso=1)[0]
	
	lgmca_gnom  = hp.gnomview(lgmca_map, flip = 'geo',rot = [ra_dec_cen[0], ra_dec_cen[1],-90], xsize = boxsize, coord = ['G','C'], reso = new_map_reso, return_projected_map = 1)
	sz_map = np.fft.ifft2( np.fft.fft2(sz_cts[key][0]) * Bl_ratio * TF).real
	plnk_map = np.fft.ifft2( np.fft.fft2(temp) * TWODTF*TF).real
	gnom_map = np.fft.ifft2( np.fft.fft2(lgmca_gnom) * TF*TWODTF).real
	subplot(131);imshow(sz_map*1e6,vmin =-100, vmax =100);colorbar();title('sz_map')
	subplot(132);imshow(plnk_map,vmin=-100, vmax =100);colorbar();title('interpolation map')
	subplot(133);imshow(gnom_map,vmin =-100, vmax =100);colorbar();title('gnom view')
	i =1


"""
from IPython import embed;embed()
lgmca_gnom  = hp.gnomview(lgmca_map, flip = 'geo',rot = [ra_dec_cen[0], ra_dec_cen[1],-90], xsize = boxsize, coord = ['G','C'], reso = dx, return_projected_map = 1)
close()
subplot(121);imshow(temp);colorbar()
subplot(122);imshow(lgmca_gnom);colorbar();show()
#imshow(temp - lgmca_gnom);colorbar();show()
"""