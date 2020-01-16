def fn_apodize(MAP, mapparams, apodmask = 'circle', edge_apod = 2, min_kernel_size = 10):

	"""
	MAP = input map to be apodized
	mapparams = [nx, ny, dx, dy] (nx, ny is normally MAP.shape)
	mask = circlular / square
	edge_apod = how many arcmins to apod from the edge. You can play with this number
	"""

	import numpy as np, scipy.ndimage as ndimage

	nx, ny, dx, dy = mapparams

        radius = (nx * dx)/20
	if radius < min_kernel_size: 
		radius = min_kernel_size

        ker=np.hanning(int(radius/dx))
        ker2d=np.asarray( np.sqrt(np.outer(ker,ker)) )

        MASKf=np.zeros((nx,ny))

	#the below grid is in arcmins - note the dx factor
	minval, maxval = -(nx*dx)/2,  (nx*dx)/2
	x = y = np.linspace(minval, maxval, nx)
	X, Y = np.meshgrid(x,y)
	xc, yc = 0., 0.

	if apodmask == 'circle':
		radius = (nx * dx/2) - edge_apod
		inds=np.where((X-xc)**2. + (Y-yc)**2. <= radius**2.) #all in arcmins

	elif apodmask == 'square':
		radius = (nx * dx/2) - edge_apod
	        inds=np.where((abs(X)<=radius) & (abs(Y)<=radius)) #all in arcmins

        MASKf[inds]=1.

	apodMASKf=ndimage.convolve(MASKf, ker2d)#, mode='wrap')
	apodMASKf/=apodMASKf.max()

	#imshow(apodMASKf);colorbar();show();quit()

	return apodMASKf * MAP

def fn_psource_mask(MAP, mapparams, coords, disc_rad = None, min_kernel_size = 10):

	"""
	MAP - input map
	mapparams = [nx, ny, dx, dy] (nx, ny is normally MAP.shape)
	coords - (x,y) pixel coords
	disc_rad - x arcmins to apodize around each point source
		roughly disc_rad/dx will give you the right arcmins. 
	"""

	import numpy as np, scipy.ndimage as ndimage


	x, y = np.arange(np.shape(MAP)[1]), np.arange(np.shape(MAP)[0])
	X, Y = np.meshgrid(x,y)

	#create a mask full of ones and make holes at the location of point sources
	MASK = np.ones(MAP.shape)
	for (i,j) in coords:
		inds=np.where((X-i)**2. + (Y-j)**2. <= disc_rad**2.)
		MASK[inds]=0.

	#imshow(MASK);colorbar();show();quit()

	#create a hanning kernel
	nx, ny, dx, dy = mapparams
        radius = (nx * dx)/20
	if radius < min_kernel_size: 
		radius = min_kernel_size

	ker=np.hanning(int(radius)+1)
	ker2d=np.sqrt(np.outer(ker,ker))

	#imshow(ker2d);colorbar();show();quit()

	#convolve
	MASK=ndimage.convolve(MASK, ker2d)
	MASK/=np.max(MASK)

	#imshow(MASK);colorbar();show();quit()

	return MAP * MASK

###############################################################
###############################################################

import numpy as np
from pylab import *

#map specs
nx, ny, dx, dy = 400, 400, 0.5, 0.5 #200 arcmin map
mapparams = [nx, ny, dx, dy]

# dummy map
MAP = np.random.randn(nx,ny)

#edge apodization
MAP = fn_apodize(MAP, mapparams)

#subplot(121);imshow(MAP);colorbar()
#subplot(122);imshow(MAP_apod);colorbar();show()

#point source apodization
#first create some dummy coordinates
totalsources = 5
coords = [(np.random.randint(nx),np.random.randint(ny)) for i in range(totalsources)]

disc_rad = 6. #6 * dx = 3 arcmin radius in this case
disc_rad_arcmin = disc_rad / dx 
MAP_psrc_masked = fn_psource_mask(MAP, mapparams, coords, disc_rad = disc_rad_arcmin)

subplot(121);imshow(MAP);colorbar()
subplot(122);imshow(MAP_psrc_masked);colorbar();show()




