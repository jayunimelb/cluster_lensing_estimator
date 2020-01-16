"""
### Author: Srinivasan Raghunathan
### CHECK - to be corrected
### Date start: 26.02.2016
### Last udpate: 27.03.2016

"""

class scipy_useful_functions():

	def __init__(self):
		pass

	def _centered(self, arr, newsize):
	    # Return the center newsize portion of the array.
	    newsize = asarray(newsize)
	    currsize = array(arr.shape)
	    startind = (currsize - newsize) // 2
	    endind = startind + newsize
	    myslice = [slice(startind[k], endind[k]) for k in range(len(endind))]
	    return arr[tuple(myslice)]

	def _next_regular(self,target):
	    """
	    Find the next regular number greater than or equal to target.
	    Regular numbers are composites of the prime factors 2, 3, and 5.
	    Also known as 5-smooth numbers or Hamming numbers, these are the optimal
	    size for inputs to FFTPACK.
	    Target must be a positive integer.
	    """
	    if target <= 6:
		return target

	    # Quickly check if it's already a power of 2
	    if not (target & (target-1)):
		return target

	    match = float('inf')  # Anything found will be smaller
	    p5 = 1
	    while p5 < target:
		p35 = p5
		while p35 < target:
		    # Ceiling integer division, avoiding conversion to float
		    # (quotient = ceil(target / p35))
		    quotient = -(-target // p35)

		    # Quickly find next power of 2 >= quotient
		    try:
		        p2 = 2**((quotient - 1).bit_length())
		    except AttributeError:
		        # Fallback for Python <2.7
		        p2 = 2**(len(bin(quotient - 1)) - 2)

		    N = p2 * p35
		    if N == target:
		        return N
		    elif N < match:
		        match = N
		    p35 *= 3
		    if p35 == target:
		        return p35
		if p35 < match:
		    match = p35
		p5 *= 5
		if p5 == target:
		    return p5
	    if p5 < match:
		match = p5
	    return match


	def make_apod(self, MAP, mapparams, fwhm_wght=10., wdth_bord=30., fwhm_apod=15., maxthresh=None, avgthresh=None):
		""" construct an apodization mask, taking this map as a set of weights. the process is
		   (1) smooth the map, with a full-width-at-half-maximum (fwhm) given by fwhm_weight (in arcmin).
		   (2) threshold the smoothed weights, as a percentage of the smoothed maximum (maxthresh) and/or of the smoothed average (avgthresh).
		   (3) remove all pixels within a distance of wdth_bord (in arcmin) of any pixels which have been thresholded to zero.
		   (4) apply a gaussian apodization, with fwhm given by fwhm_apod. 
		"""
		import scipy.ndimage

		nx, ny, dx, dy = mapparams
		reso_arcmin = (dx)# * 180.*60./np.pi)

		smoothwt  = scipy.ndimage.gaussian_filter( MAP, fwhm_wght / reso_arcmin / 2*np.sqrt(2*np.log(2)) )

		threshwt  = np.ones( smoothwt.shape )
		if maxthresh != None:
		    threshwt[ np.where(smoothwt / smoothwt.flatten().max() < maxthresh) ] = 0.0
		if avgthresh != None:
		    threshwt[ np.where(smoothwt / smoothwt.flatten()[np.where(smoothwt.flatten() != 0)].avg() < avgthresh) ] = 0.0

		npix_bord = 2.*int(wdth_bord/reso_arcmin)
		
		xs, ys    = np.meshgrid( np.linspace(-1., 1., npix_bord), np.linspace(-1., 1., npix_bord) )
		kern_bord = np.ones( (npix_bord, npix_bord) )
		kern_bord[ np.where( (xs**2 + ys**2) >= 1. ) ] = 0.0

		imshow(threshwt);colorbar();show()
		quit()


		bordwt = scipy.ndimage.minimum_filter( threshwt, footprint=kern_bord )

		return rmap( self.nx, self.dx, ny=self.ny, dy=self.dy,
			     map=scipy.ndimage.gaussian_filter( bordwt, fwhm_apod / reso_arcmin / 2*np.sqrt(2*np.log(2)) ) )



	def pad_ft(self, a, npad=2):
		""" pad a 2D Fourier transform (produced by np.fft2) with zeros, useful when performing convolutions to ensure that the result is not aliased.
		   * a    = 2D complex Fourier transform.
		   * npad = fractional size to pad. npad=2 will double the size of the Fourier transform in each dimension.
		"""
		if npad==1: return a
		nx,ny = a.shape
		p = np.zeros([nx*npad,ny*npad], dtype=a.dtype)
		p[0:nx,0:ny] = np.fft.fftshift(a)
		p = np.roll(np.roll(p,-nx/2,axis=0),-ny/2,axis=1)
		return p

	def unpad_ft(self, a, npad=2):
		""" un-pad an array in Fourier-space, removing the additional zeros added by 'pad_ft' """
		if npad==1: return a
		nx_pad,ny_pad = a.shape
		nx = int(nx_pad/npad); ny=int(ny_pad/npad)
		return np.roll(np.roll(
		    (np.roll(np.roll(a,nx/2,axis=0),ny/2,axis=1)[0:nx,0:ny]),
		    nx/2,axis=0),ny/2,axis=1)

	def convolve_padded(self, f, g, npad=2):
		""" convolve two 2D complex Fourier transforms 'f' and 'g', using padding by a factor of npad to avoid aliasing.
		   returns r(L) = \int{d^2\vec{l}} f(l) g(L-l).
		"""
		return (self.unpad_ft(np.fft.fft2(
			np.fft.ifft2(self.pad_ft(f,npad=npad)) *
			np.fft.ifft2(self.pad_ft(g,npad=npad))),
			     npad=npad)*npad**2)



class simulations():

	def __init__(self):

		import numpy as np
		self.Tcmb = 2.73 #K
		
		self.degrees2radians = np.pi / 180.
		self.arcmins2radians = self.degrees2radians / 60.
		#self.exp_beam = 1.2 #arcmins
		self.tqulen = 3

		self.inidic = {}
		self.covdic = {}
		self.npixels = None
		self.noofsims = 30000
		

	def _ini_log_file(self, log_file):

		self.log_file = log_file


	def _ini_params(self, addtoinidic):

		if len(self.inidic.keys()) == 0:
			self.inidic = addtoinidic
		else:
			for keys in sorted( addtoinidic.keys() ):
				self.inidic[keys] = addtoinidic[keys]

	def is_seq(self,o):
		return hasattr(o, '__len__')

	def is_seq_of_seq(self,o):
		if not self.is_seq(o):
			return False
		for s in o:
			if not self.is_seq(s):
				return False
			return True

	def get_lxly(self,mapparams):

		import numpy as np

		nx, ny, dx, dy = mapparams
		dx *= self.arcmins2radians
		dy *= self.arcmins2radians

		#print 'hi', dx, dy

		lx, ly = np.meshgrid( np.fft.fftfreq( nx, dx ), np.fft.fftfreq( ny, dy ) )# * np.pi
		lx *= 2* np.pi
		ly *= 2* np.pi

		return lx, ly

	def gauss_beam(self, fwhm, EL, pol=1): #from Healpy
		"""
		Imported from healpy gauss_beam()
		fhwm = experiment fwhm in radians
		EL = 2d \ell matrix
		pol = 1 #includes pol.factor
		"""

		import numpy as np

		sigma = fwhm / np.sqrt(8. * np.log(2.))
		sigma2 = sigma ** 2
		if not self.is_seq(EL):			
			ell = np.arange(EL + 1)
		else:
			ell = EL

		g = np.exp(-.5 * EL**2. * sigma2)

		#beam_lmat = np.exp(-(beam_fwhm**2.)*(l2d**2.)/(16.*np.log(2.)))

		#from pylab import *
		#imshow(g);colorbar();show();quit()

		if not pol: # temperature-only beam
			return g
		else: # polarization beam
			# polarization factors [1, 2 sigma^2, 2 sigma^2, sigma^2]
			pol_factor = np.exp([0., 2*sigma2, 2*sigma2])#, sigma2]) #ignoring T
			G = np.asarray([g * scale for scale in pol_factor])
			return G

        def fn_get_HPF(self, mapparams, ideal = 0, minel = None, maxel = 9000, beamel_cutoff = 0): #Based on Eric's code

                import numpy as np

                nx, ny, dx, dy = mapparams
                dx *= self.arcmins2radians
                dy *= self.arcmins2radians

                lx, ly = self.get_lxly(mapparams)
                L = np.sqrt(lx**2. + ly**2.)
                transfer_lmat =1.+ np.zeros(L.shape)

                if ideal:
                        transfer_lmat = np.asarray( [transfer_lmat for ii in range(self.tqulen)] )
                        return transfer_lmat

		if minel<> None:
		        iso_hipass = np.where(L < minel)
		        transfer_lmat[iso_hipass] = 0.

		if maxel<> None:
		        iso_lopass = np.where(L > maxel)
		        transfer_lmat[iso_lopass] = 0.

		if beamel_cutoff:
		        beam_lowpass = np.where(L> self.elbeam)
		        transfer_lmat[beam_lowpass] = 0.0

		#subplot(111);pcolor(lx,ly,transfer_lmat);colorbar();show();quit()
		
	   	#scan_hipass = np.where( (np.abs(lx) < minel) )
                #transfer_lmat[scan_hipass] = 0.

                transfer_lmat = np.asarray( [transfer_lmat for ii in range(self.tqulen)] )

                return transfer_lmat


	def fn_get_HPF_old(self, mapparams, minel = 300, maxel = 11500): #Based on Eric's code

		import numpy as np

		nx, ny, dx, dy = mapparams
		dx *= self.arcmins2radians
		dy *= self.arcmins2radians

		lx, ly = self.get_lxly(mapparams)
		L = np.sqrt(lx**2. + ly**2.)
		transfer_lmat =1.+ np.zeros(L.shape)
		iso_hipass = np.where(L < minel)
		lmax_nyquist = 2.*np.pi*1./dx/2.*(nx-1.)/nx
		nyquist_lopass = np.where(L > lmax_nyquist)
		transfer_lmat[iso_hipass] = 0.
		transfer_lmat[nyquist_lopass] = 0.0

		#scan_hipass = np.where(np.abs(lx) < minel)
		#scan_lopass = np.where(np.abs(lx) > maxel)
		#transfer_lmat[scan_hipass] = 0.
		#transfer_lmat[scan_lopass] = 0.

		return transfer_lmat

	def Cls2CLS(self,Cls,mapparams): #based on E. Baxter's code

		import numpy as np

		nx, ny, dx, dy = mapparams
		dx *= self.arcmins2radians
		dy *= self.arcmins2radians

		# find how many different Cls were passed
		noofrows = len(Cls)
		if self.is_seq_of_seq(Cls):
			noofcols = np.shape(Cls)[1]
			Cls_tot = noofcols - 1 #first is els, then Cls
		else:
			els = np.arange(2,noofrows)
			Cls = [els,Cls]
			Cls_tot = 1

		lx, ly = self.get_lxly(mapparams)
		L = np.sqrt(lx**2. + ly**2.)


		# processing Cls now
		CLS = np.zeros( (Cls_tot,L.shape[0],L.shape[1]) )
		for clcnt in range(Cls_tot):
			CLS[clcnt,:,:] = np.interp(L.flatten(), Cls[:,0], Cls[:,clcnt+1], right = 0.).reshape(L.shape)

		return CLS

	def Dls2map(self, Dls, mapparams, nosims = 10, passing_Dls = 1, CMB_outputscale = 1, cls_lensed_lss = 1, no_sims_required = 0, random_seed = None): #Handles both T and P.
		"""
		1. Random Gaussian CMB simulations
		2. Supports both T and P
		3. Beam convolution, transfer matrix conv., adding exp. noise level - all available
		
		parameters:
		1. Dls = [ells, DlTT] (or) [ells, DlTT, DlEE, DlBB]

		2. passing_Dls = 1
		if passing_Dls == 0:
			Cls = Dls
		else:
			Cls will be converted to Dls using
			Dl = ( l(l+1) * Cl ) / (2 * pi) - from CAMB

		3. mapparams = [nx, ny, dx, dy] - dimensions and ang. resolution in arcmins

		4. nosims = number of simulations to be performed

		5. exp_beam = beam_FWHM in arcmins

		6. exp_tmat = experiment transfer matrix

		7. exp_noise_level = noise level of experiment in uk-arcmin
		"""

		import numpy as np, time

                if not random_seed == None:
                        np.random.seed(random_seed)

		###########################################################################
		#first check if only T or P is supplied as well; also obtain els
		noofrows = len(Dls)
		if self.is_seq_of_seq(Dls):
			els = np.asarray(Dls[:,0])
			noofcols = np.shape(Dls)[1]
			Dls_tot = noofcols - 1 #first is els, then Dls
			if Dls_tot > 1:
				only_temp = 0
		else:
			els = np.arange(2,noofrows)
			Dls = [els,Dls]
		###########################################################################
		#Dls to Cls
		if passing_Dls:
			if CMB_outputscale == 1:
				Cls = ( self.Tcmb**2. * Dls * 2 * np.pi ) / ( els[:,None] * (els[:,None] + 1) )
			else:
				Cls = ( Dls * 2 * np.pi ) / ( els[:,None] * (els[:,None] + 1) )
			Cls[:,0] = els
		else:
			Cls = Dls

		#exp_var = np.sum( ( (2 * Cls[:,0] + 1) / (4 * np.pi) ) * Cls[:,1] ) #Eq. 11.37 Dodelson
		#print exp_var;quit()


		###########################################################################
		#Cls2map
		nx, ny, dx, dy = mapparams
		dx *= self.arcmins2radians
		dy *= self.arcmins2radians

		#scalefac = np.sqrt((nx * ny) * (dx * dy))
		scalefac = np.sqrt(1./ (dx * dy))

		CLS = self.Cls2CLS(Cls,mapparams)
		if cls_lensed_lss:
			self.CLS = CLS #dump it for later usage
		else:
			self.CLS_unlen = CLS

		if no_sims_required:
			return None


		CLS = np.sqrt(CLS)

		lx, ly = self.get_lxly(mapparams)
		angle = 2.*np.arctan2(lx, -ly)

		########################################################################
		########################################################################
		# processing Cls now
		if Dls_tot>1:
			TFFT, EFFT, BFFT = CLS
		else:
			TFFT = CLS

		#T = np.fft.ifft2( TFFT * scalefac ).real

		if Dls_tot>1:
			FORQ = ( np.cos(angle) * EFFT - np.sin(angle) * BFFT )
			FORU = ( np.sin(angle) * EFFT + np.cos(angle) * BFFT )

			FORQ[np.isnan(FORQ)] = 1e-10
			FORU[np.isnan(FORU)] = 1e-10

			#Q = np.fft.ifft2( FORQ ).real * scalefac
			#U = np.fft.ifft2( FORU ).real * scalefac
			#CAMBMAP = np.asarray([T,Q,U])
			CAMBMAP = np.asarray([TFFT,FORQ,FORU]) * scalefac

			self.tqulen = 3

		else:
			#CAMBMAP = np.asarray([T])
			CAMBMAP = np.asarray([TFFT]) * scalefac
			self.tqulen = 1

		#CAMBMAP = np.asarray([TFFT]) * scalefac
		#self.tqulen = 1

		#CAMBMAP_FFT = np.fft.fft2(CAMBMAP) #go to fourier space for easy convolution
		CAMBMAP_FFT = np.copy(CAMBMAP)

		logline = '\t\tDls converted to maps'
		logfile = open(self.log_file,'a')
		logfile.writelines('%s\n' %(logline))
		logfile.close()
		print logline

		start = time.time()
		logline = '\t\tstarting sims'
		logfile = open(self.log_file,'a')
		logfile.writelines('%s\n' %(logline))
		logfile.close()
		print logline

		#TWODHPF = self.fn_get_HPF(mapparams)

		#GNOISE = []
		SIMMAPS = np.zeros( (nosims, self.tqulen, nx, ny) )
		for simcnt in range(nosims):
			if simcnt % 1000 == 0 and simcnt > 0:
				logline = '\t\t\t\tsimno: %s' %(simcnt)
				logfile = open(self.log_file,'a')
				logfile.writelines('%s\n' %(logline))
				logfile.close()
				print logline
			#TMP = np.fft.ifft2( CAMBMAP_FFT * TWODHPF * np.fft.fft2( np.random.randn(nx,ny) ) ).real
			TMP = np.fft.ifft2( CAMBMAP_FFT * np.fft.fft2( np.random.randn(nx,ny) ) ).real

			#from pylab import *
			#imshow(TMP[0,:,:]);colorbar();show();quit()

			#for bbb in range(self.tqulen):
			#	TMP[bbb,:,:] -= np.mean( TMP[bbb,:,:].ravel() )

			SIMMAPS[simcnt,:,:,:] = TMP

		#GNOISE = np.asarray(GNOISE)
		#SIMMAPS = np.asarray(GNOISE)
		end = time.time()
		#print SIMMAPS.shape

		logline = '\t\tGaussian sims done. time take = %s seconds' %(end-start)
		logfile = open(self.log_file,'a')
		logfile.writelines('%s\n' %(logline))
		logfile.close()
		print logline

		return SIMMAPS


	def fn_perform_lensing_nfw(self, SIMMAPS, mapparams, param_dict, RA, DEC, clus_ra, clus_dec, M_200_list, c_200_list, z_L_list, z_lss, mass_def, rho_def, return_kappa=0, show_plots = 0, smooth_kappa_resol = None):

		#import EBX_modules.nfw_kappa_funcs as EBX_nfw
		import nfw_kappa_funcs as EBX_nfw
		import numpy as np, time
		from scipy import interpolate as intrp

		#MAP RESOL, EL stuffs
		nx, ny, dx, dy = mapparams
		dx *= self.arcmins2radians
		dy *= self.arcmins2radians
		lx, ly = self.get_lxly(mapparams)
		L = np.sqrt(lx**2. + ly**2.)

		#get the kappa map from EBX
		KAPPA = EBX_nfw.get_NFW_kappa_fullmap(param_dict, RA, DEC, clus_ra, clus_dec, M_200_list, c_200_list, z_L_list, z_lss, mass_def, rho_def)

		#smoothing kappa
		if smooth_kappa_resol <> None:
			dx_ori = mapparams[2]
			degrade_fac = smooth_kappa_resol / dx_ori

			print KAPPA.shape
			KAPPA = self.downsample_map(KAPPA, degrade_fac)
			print KAPPA.shape
			quit()		

		#from pylab import *
		#css=imshow(KAPPA);css.set_clim(0.,.7);colorbar();show();quit()

		#make sure kappa does't have nan, inf, etc.
		smallnumber = min(abs(KAPPA.ravel()))
		KAPPA[np.isnan(KAPPA)] = smallnumber
		KAPPA[np.isinf(KAPPA)] = smallnumber

		PHI_FFT = -2. * dx * dy * np.fft.fft2(KAPPA)/(L**2)
		#make sure kappa does't have nan, inf, etc.
		PHI_FFT[np.isnan(PHI_FFT)] = smallnumber
		PHI_FFT[np.isinf(PHI_FFT)] = smallnumber

		DEF_X    = np.fft.ifft2(-1j * PHI_FFT * lx) / ( dx * dy )
		DEF_Y    = np.fft.ifft2(-1j * PHI_FFT * ly) / ( dx * dy )

		logline = '\t\tkappa, phi obtained'
		logfile = open(self.log_file,'a')
		logfile.writelines('%s\n' %(logline))
		logfile.close()
		print logline

		if show_plots:
			from pylab import *
			subplot(221);imshow(KAPPA.real);colorbar();subplot(222);imshow(PHI_FFT.real);colorbar();
			subplot(223);imshow(DEF_X.real);colorbar();subplot(224);imshow(DEF_Y.real);colorbar();show();quit()

		#Angular info on sky
		theta_x_list = np.array(range(1,nx+1)) * dx - dx * 0.5 * (nx - 1.)
		if nx<>ny:
			theta_y_list = np.array(range(1,ny+1))*dy-dy*0.5*(ny - 1.)
			theta_x,theta_y = np.meshgrid(theta_x_list, theta_y_list)
		else:
			theta_x = np.tile(theta_x_list,(nx,1))
			theta_y = np.transpose(theta_x)

		to_evaluate_unlensed_theta_x = (theta_x + DEF_X).flatten()
		to_evaluate_unlensed_theta_y = (theta_y + DEF_Y).flatten()

		SIMMAPS_L = np.copy(SIMMAPS)
		nosims = SIMMAPS_L.shape[0]
		tqulen = self.tqulen

		logline = '\t\tlensing to be performed now - interpolation actually'
		logfile = open(self.log_file,'a')
		logfile.writelines('%s\n' %(logline))
		logfile.close()
		print logline

		start = time.time()
		for cnt in range(nosims):
			if cnt % 1000 == 0 and cnt>0:
				logline = '\t\t\t\tsimno: %s' %(cnt)
				logfile = open(self.log_file,'a')
				logfile.writelines('%s\n' %(logline))
				logfile.close()
				print logline

			for tqucnt in range(tqulen):
				SIMMAPS_L[cnt, tqucnt] = intrp.RectBivariateSpline( theta_y[:,0], theta_x[0,:], SIMMAPS[cnt,tqucnt], kx=5, ky=5).ev(to_evaluate_unlensed_theta_y, to_evaluate_unlensed_theta_x).reshape([ny,nx])
		end = time.time()

		logline = '\t\tlensing complete. time take = %s seconds' %(end-start)
		logfile = open(self.log_file,'a')
		logfile.writelines('%s\n' %(logline))
		logfile.close()
		print logline

		if not return_kappa:
			return SIMMAPS_L
		else:
			return SIMMAPS_L, KAPPA

	def fn_get_phi_thetas(self, mapparams, KAPPA):

		import numpy as np

		nx, ny, dx, dy = mapparams
		lx, ly = self.get_lxly(mapparams)
		L = np.sqrt(lx**2. + ly**2.)

		PHI_FFT = -2. * dx * dy * np.fft.fft2(KAPPA)/(L**2)
		#make sure kappa does't have nan, inf, etc.
		smallnumber = min(abs(KAPPA.ravel()))
		PHI_FFT[np.isnan(PHI_FFT)] = smallnumber

		DEF_X    = np.fft.ifft2(-1j * PHI_FFT * lx) / ( dx * dy )
		DEF_Y    = np.fft.ifft2(-1j * PHI_FFT * ly) / ( dx * dy )

		#Angular info on sky
		theta_x_list = np.array(range(1,nx+1)) * dx - dx * 0.5 * (nx - 1.)
		if nx<>ny:
			theta_y_list = np.array(range(1,ny+1))*dy-dy*0.5*(ny - 1.)
			theta_x,theta_y = np.meshgrid(theta_x_list, theta_y_list)
		else:
			theta_x = np.tile(theta_x_list,(nx,1))
			theta_y = np.transpose(theta_x)

		to_evaluate_unlensed_theta_x = (theta_x + DEF_X).flatten()
		to_evaluate_unlensed_theta_y = (theta_y + DEF_Y).flatten()

		return PHI_FFT, DEF_X, DEF_Y, theta_x, theta_y, to_evaluate_unlensed_theta_x, to_evaluate_unlensed_theta_y

	def downsample_map(self, data, N=8): #from N.Whitehorn
		from numpy import average, split
		width = data.shape[0]
		height= data.shape[1]
		return average(split(average(split(data, width // N, axis=1), axis=-1), height // N, axis=1), axis=-1)

	def fn_perform_lensing_nfw_finer(self, SIMMAPS, mapparams, param_dict, RA, DEC, clus_ra, clus_dec, M_200_list, c_200_list, z_L_list, z_lss, mass_def, rho_def,return_kappa=0, truncate_kappa_at_radius = .5, min_theta_cutoff_arcmins = 3., only_return_lensing_stuffs = 0, ipfac = 5):

		"""
		truncate_kappa_at_radius [r_200] = truncate KAPPA at this radius - see http://arxiv.org/abs/1401.1216
		"""

		#import EBX_modules.nfw_kappa_funcs as EBX_nfw
		import nfw_kappa_funcs as EBX_nfw
		import numpy as np, time
		from scipy import interpolate as intrp

		## make a copy of original RA, DEC
		RA_ori = np.copy(RA)
		DEC_ori = np.copy(DEC)

		#MAP RESOL, EL stuffs
		nx, ny, dx, dy = mapparams
		nx = nx * ipfac
		ny = ny * ipfac
		dx = dx / ipfac
		dy = dy / ipfac

		mapparams_finer = [nx, ny, dx, dy]

		dx *= self.arcmins2radians
		dy *= self.arcmins2radians
		lx, ly = self.get_lxly(mapparams_finer)
		L = np.sqrt(lx**2. + ly**2.)

		#get the kappa map from EBX
		minra, maxra = min(RA.ravel()),max(RA.ravel())
		mindec, maxdec = min(DEC.ravel()),max(DEC.ravel())

		ra = np.linspace(minra, maxra, nx)
		dec = np.linspace(mindec, maxdec, ny)

		RA, DEC = np.meshgrid(ra,dec)

		KAPPA = EBX_nfw.get_NFW_kappa_fullmap(param_dict, RA, DEC, clus_ra, clus_dec, M_200_list, c_200_list, z_L_list, z_lss, mass_def, rho_def)

		#make sure kappa does't have nan, inf, etc.
		smallnumber = min(abs(KAPPA.ravel()))
		KAPPA[np.isnan(KAPPA)] = smallnumber

		#from pylab import *
		#imshow(KAPPA);colorbar();show();quit()

		logline = '\t\tkappa, phi obtained'
		logfile = open(self.log_file,'a')
		logfile.writelines('%s\n' %(logline))
		logfile.close()
		print logline

		#truncate kappa at certain radius
		#http://arxiv.org/abs/1401.1216 using 1 * r_200 for now
		if truncate_kappa_at_radius <> None: ##tochange - for list of M, z
			import modules.scl_cosmo as scl_cosmo
	
			cosmo = scl_cosmo.scl_cosmo_stuffs()
			cosmo.fn_update_cosmo_params(param_dict)
			r_200 = cosmo.fn_r200(M_200_list[0],z_L_list[0]) #metres
			distances = cosmo.fn_cosmo_dist(z_L_list[0]) #metres
			d_A = distances[3] #angular diameter distance

			theta_cutoff_arcmins_actual = np.degrees(r_200/d_A) * 60. #arcmins

			#print M_200_list[0], z_L_list[0], theta_cutoff_arcmins, int(theta_cutoff_arcmins * truncate_kappa_at_radius) + 1 , min_theta_cutoff_arcmins

			theta_cutoff_arcmins_actual = int(theta_cutoff_arcmins_actual * truncate_kappa_at_radius) + 1 
			if theta_cutoff_arcmins_actual< min_theta_cutoff_arcmins:
				theta_cutoff_arcmins = min_theta_cutoff_arcmins
			else:
				theta_cutoff_arcmins = theta_cutoff_arcmins_actual

			#apodMASK = self.fn_kappa_cutoff_mask(KAPPA, RA * 60., DEC * 60., np.degrees(dx/ipfac) * 60., theta_cutoff_arcmins)
			apodMASK = self.fn_kappa_cutoff_mask(KAPPA, RA * 60., DEC * 60., np.degrees(dx) * 60., theta_cutoff_arcmins)

			KAPPA = KAPPA * apodMASK

			logline = '\t\tApplying apodisation to KAPPA at max(%s * r_200/d_A = %s, %s) = %s arcmins' %(truncate_kappa_at_radius, theta_cutoff_arcmins_actual, min_theta_cutoff_arcmins, theta_cutoff_arcmins)
			logfile = open(self.log_file,'a')
			logfile.writelines('%s\n' %(logline))
			logfile.close()
			print logline
			
			"""
			from pylab import *
			subplot(131);pcolor(RA * 60., DEC * 60., KAPPA);colorbar();xlim(-10.,10.); ylim(-10., 10.)
			KAPPA_2 = KAPPA * apodMASK
			subplot(132);pcolor(RA * 60., DEC * 60., KAPPA_2);colorbar();xlim(-10.,10.); ylim(-10., 10.)
			subplot(133);pcolor(RA * 60., DEC * 60., KAPPA - KAPPA_2);colorbar();xlim(-10.,10.); ylim(-10., 10.)
			show();quit()
			"""
			
		PHI_FFT, DEF_X, DEF_Y, theta_x, theta_y, to_evaluate_unlensed_theta_x, to_evaluate_unlensed_theta_y = self.fn_get_phi_thetas(mapparams_finer,KAPPA)

		KAPPA_DEG = intrp.RectBivariateSpline( DEC[:,0], RA[0,:], KAPPA ).ev(DEC_ori.ravel(), RA_ori.ravel()).reshape(RA_ori.shape)
		PHI_FFT_deg, DEF_X_deg, DEF_Y_deg, theta_x_deg, theta_y_deg, to_evaluate_unlensed_theta_x_deg, to_evaluate_unlensed_theta_y_deg = self.fn_get_phi_thetas(mapparams,KAPPA_DEG)

		"""
		if only_return_lensing_stuffs:
			return KAPPA, PHI_FFT, DEF_X, DEF_Y

		show_plots = 1
		if show_plots:
			from pylab import *
			#subplot(221);contourf(RA*60., DEC*60., KAPPA.real);colorbar();subplot(222);contourf(RA*60., DEC*60.,DEF_X.real);colorbar();subplot(223);contourf(RA*60., DEC*60.,DEF_Y.real);colorbar();show();quit()#subplot(224);contourf(RA*60., DEC*60., DEF_X.real + DEF_Y.real);show();quit()
			subplot(221);contourf(RA*60., DEC*60., KAPPA.real);colorbar();subplot(222);contourf(RA*60., DEC*60.,DEF_X.real);colorbar();subplot(223);contourf(RA*60., DEC*60.,DEF_Y.real);colorbar();show();quit()#subplot(224);contourf(RA*60., DEC*60., DEF_X.real + DEF_Y.real);show();quit()
		"""

		SIMMAPS_L = np.copy(SIMMAPS)
		noofsims = SIMMAPS_L.shape[0]
		tqulen = self.tqulen

		logline = '\t\tlensing to be performed now - interpolation actually'
		logfile = open(self.log_file,'a')
		logfile.writelines('%s\n' %(logline))
		logfile.close()
		print logline

		start = time.time()
		for cnt in range(noofsims):
			if cnt % 100 == 0 and cnt>0:
				logline = '\t\t\t\tsimno: %s' %(cnt)
				logfile = open(self.log_file,'a')
				logfile.writelines('%s\n' %(logline))
				logfile.close()
				print logline

			for tqucnt in range(tqulen):

				'''
				UNLENSED_TMP = intrp.RectBivariateSpline( theta_y_deg[:,0], theta_x_deg[0,:], SIMMAPS[cnt,tqucnt] ).ev(theta_y.flatten(), theta_x.flatten()).reshape(RA.shape)
				TMP = intrp.RectBivariateSpline( theta_y[:,0], theta_x[0,:], UNLENSED_TMP ).ev(to_evaluate_unlensed_theta_y, to_evaluate_unlensed_theta_x).reshape(RA.shape)

				'''
				#sparser grid
				#TMP_SPARSE = intrp.RectBivariateSpline( theta_y_deg[:,0], theta_x_deg[0,:], SIMMAPS[cnt,tqucnt] ).ev(to_evaluate_unlensed_theta_y_deg, to_evaluate_unlensed_theta_x_deg).reshape(RA_ori.shape)

				TMP = intrp.RectBivariateSpline( theta_y_deg[:,0], theta_x_deg[0,:], SIMMAPS[cnt,tqucnt] ).ev(to_evaluate_unlensed_theta_y, to_evaluate_unlensed_theta_x).reshape(RA.shape)

				#degrade them to normal resol.
                                TMP_DEGRADED = intrp.RectBivariateSpline( to_evaluate_unlensed_theta_y.reshape(RA.shape)[:,0], to_evaluate_unlensed_theta_x.reshape(RA.shape)[0,:], TMP ).ev(to_evaluate_unlensed_theta_y_deg, to_evaluate_unlensed_theta_x_deg).reshape(RA_ori.shape)

				###TMP_DEGRADED = self.downsample_map(TMP,ipfac)

				'''
				from pylab import *
				subplot(221);css=imshow(TMP_SPARSE);colorbar();
				subplot(222);css=imshow(TMP_DEGRADED);colorbar();

				start = 20
				ex1, ex2, ey1, ey2 = start, len(TMP_DEGRADED)-start,start,len(TMP_DEGRADED)-start

				subplot(223);css=imshow((SIMMAPS[cnt,tqucnt] - TMP_SPARSE)[ex1:ex2, ey1:ey2]);colorbar()
				subplot(224);css=imshow((SIMMAPS[cnt,tqucnt] - TMP_DEGRADED)[ex1:ex2, ey1:ey2]);colorbar();show();quit()

				'''
				'''
				from pylab import *
				subplot(231);imshow(TMP);colorbar()
				subplot(232);imshow(TMP_DEGRADED);colorbar()
				subplot(233);imshow(TMP_SPARSE);colorbar()
				#subplot(234);imshow(UNLENSED_TMP - TMP);colorbar()
				s, e = 10, 90
				subplot(234);css = imshow(TMP_SPARSE[s:e,s:e] - TMP_DEGRADED[s:e,s:e]);colorbar()
				subplot(235);css = imshow(SIMMAPS[cnt,tqucnt][s:e,s:e] - TMP_SPARSE[s:e,s:e]);colorbar()
				subplot(236);css = imshow(SIMMAPS[cnt,tqucnt][s:e,s:e] - TMP_DEGRADED[s:e,s:e]);colorbar()
				show();quit()
				'''

				SIMMAPS_L[cnt, tqucnt] = TMP_DEGRADED

		end = time.time()
		
		logline = '\t\tlensing complete. time take = %s seconds' %(end-start)
		logfile = open(self.log_file,'a')
		logfile.writelines('%s\n' %(logline))
		logfile.close()
		print logline

		if not return_kappa:
			return SIMMAPS_L
		else:
			return SIMMAPS_L, KAPPA

	def fn_gaussian_1d(self,xarr,yarr):

		import scipy.optimize as optimize
		import scipy.special

		xarr=np.asarray(xarr)
		yarr=np.asarray(yarr)

		fitfunc = lambda p, x: p[1]*(np.exp(-(x-p[2])**2/(2*p[3]**2)))
		errfunc = lambda p, x, y: fitfunc(p, x) - y #minimization

		p0 = [0.,np.max(yarr),np.mean(xarr),np.std(xarr)] # Initial guess for the parameters
		p1, success = optimize.leastsq(errfunc, p0[:], args=(xarr, yarr))

		fit=fitfunc(p1,xarr)

		return fit


	def fn_get_exp_noise(self, mapparams, noiselevel): #noiselevel in uK-arcmin

		import numpy as np

		nx, ny, dx, dy = mapparams
		dx *= self.arcmins2radians
		dy *= self.arcmins2radians
		npixels = nx * ny

		noiselevel = np.atleast_1d(noiselevel)
		totaliter = len(noiselevel)

		NOISE = np.zeros( (totaliter, nx,ny) )

		for cnt in range(totaliter):
			DeltaT = noiselevel[cnt]*1e-6*(1./60.)*(np.pi/180.) #K-radian
			#NOISE[cnt,:,:] = np.random.standard_normal([nx,ny]) * DeltaT/dx
			NOISE[cnt,:,:] = np.random.normal(loc=0.0, scale=DeltaT/dx, size=[nx,ny])

			'''
			dummynoise = NOISE[cnt,:,:].ravel()
			hists,binedges=np.histogram(dummynoise,bins = 50)
			gaussfit=self.fn_gaussian_1d(binedges[0:-1],hists)
			clf();hist(dummynoise,bins = 50,color='r')
			plot(binedges[:-1],gaussfit,'k')
			show();quit()
			'''

		return NOISE
		

	def fn_beam_tf(self, SIMMAPS, mapparams, beamfwhmarcmins, nosims = 10, reqdbox = None, use_beam=1, add_TF=0, exp_noise_level=None):

		########################################################################
		########################################################################
		import numpy as np, scipy.ndimage as ndimage
		import time

		scfns = scipy_useful_functions()

		#MAP RESOL, EL stuffs
		nx, ny, dx, dy = mapparams
		dx *= self.arcmins2radians
		dy *= self.arcmins2radians
		lx, ly = self.get_lxly(mapparams)
		L = np.sqrt(lx**2. + ly**2.)

		if use_beam:
			self.exp_beam = beamfwhmarcmins
			exp_beam = self.exp_beam #arcmins
			exp_beam *= self.arcmins2radians

			self.elbeam = np.sqrt(8. * np.log(2.)) / exp_beam

			EL = np.sqrt(lx**2. + ly**2.)
			if self.tqulen>1:
				Bl = self.gauss_beam(exp_beam,EL) #healpy based #in fourier space (or) el space
			else:
				Bl = self.gauss_beam(exp_beam,EL, pol=0) #healpy based #in fourier space (or) el space

		logline = '\t\tperform beam (+TF) convolution to sims'
		logfile = open(self.log_file,'a')
		logfile.writelines('%s\n' %(logline))
		logfile.close()
		print logline

		start = time.time()
		
		self.npixels = nx * ny

		if not reqdbox == None:
			ex1, ex2, ey1, ey2 = reqdbox
			nx, ny = ex2-ex1, ey2-ey1
			self.npixels = nx * ny

		tqulen = self.tqulen

		#SIMS --> #([nosims], [3 -> T,Q,U], [nx, ny])
		SIMS = np.zeros( (nosims, tqulen, nx, ny) )
		for cnt in range(nosims):

			if cnt % 1000 == 0 and cnt>0:
				#logline = '\t\t\t\tsimno: %s for %s' %(cnt, typestr[ii])
				logline = '\t\t\t\tsimno: %s' %(cnt)
				logfile = open(self.log_file,'a')
				logfile.writelines('%s\n' %(logline))
				logfile.close()
				print logline

			MAP = SIMMAPS[cnt,:]

			#MAP = scfns.make_apod(MAP[0], mapparams, maxthresh= 10)#, fwhm_wght=10., wdth_bord=10., fwhm_apod=10., maxthresh=None, avgthresh=None)
			#MAP -= np.mean(MAP)
			#MAP = self.fn_apodize(MAP,mapparams)

			map_in_fft = 0
			#add beam
			if use_beam:

				'''
				#TAKE FFT of CAMB MAP for convolution
				#print 'TAKE FFT of CAMB MAP for convolution'
				MAP_FFT = np.fft.fft2( MAP )
				MAP_FFT =  MAP_FFT * Bl #still in fourier space
				DUMMY1 = np.fft.ifft2(MAP_FFT)[0].real
				subplot(131);imshow(DUMMY1);colorbar()
	
				'''

                                MAP_FFT = np.fft.fft2( MAP )

				MAP_FFT = MAP_FFT * Bl
				#MAP_FFT = self.fn_convolve_using_QL(MAP_FFT, Bl)

				'''
				DUMMY2 = np.fft.ifft2(MAP_FFT).real
				subplot(132);css=imshow(DUMMY2)#;css.set_clim(min(MAP.ravel()), max(MAP.ravel()))
				subplot(133);css=imshow(DUMMY1 - DUMMY2)#;css.set_clim(min(MAP.ravel()), max(MAP.ravel()))
				colorbar();show();quit()
				'''
				
			if map_in_fft:
				MAP = np.fft.ifft2(MAP_FFT).real #fourier to real space after beam, TF
				
			#clf();subplot(131);imshow(np.copy(MAP[0].real));colorbar()
			if exp_noise_level <> None:
				NOISE = self.fn_get_exp_noise(mapparams, exp_noise_level)[0]
				MAP =  MAP + NOISE

				#subplot(132);imshow(NOISE.real);colorbar();subplot(133);imshow(MAP[0] + NOISE);colorbar();show();quit()

				'''
				MAP -= np.mean(MAP.ravel())
				clf()
				for ccc in range(2):
					if ccc == 0:
						dummy = MAP[0].ravel()
					else:
						dummy = NOISE[0].ravel()
					hists,binedges=np.histogram(dummy,bins = 50)
					gaussfit=self.fn_gaussian_1d(binedges[0:-1],hists)
					hist(dummy,bins = 50, histtype='step')
					plot(binedges[:-1],gaussfit,'k')
				show();quit()
				'''


			#subplot(132);imshow(MAP[0]);colorbar()
			if add_TF:
				TWODTF = self.fn_get_HPF(mapparams, ideal = 1)#, minel=None, maxel=10000)#, minel = 300)#, maxel = 3000)#, minel = 555)#, ideal = 1) #still analytic one
                                MAP_FFT = np.fft.fft2( MAP )

				MAP_FFT = MAP_FFT * TWODTF
				#MAP_FFT = self.fn_convolve_using_QL(MAP_FFT, TWODTF)

                                MAP = np.fft.ifft2( MAP_FFT ).real

			'''
			noisespec1d = self.fn_get_powerspectrum(NOISE)
			pspec1d = self.fn_get_powerspectrum(MAP[0])

			ax =  subplot(111, yscale='log')
			plot(noisespec1d,'r');plot(pspec1d,'k');show();quit()
			'''

			#subplot(133);imshow(MAP[0]);colorbar();show();quit()

			if not reqdbox == None:
				MAP = MAP[:, ex1:ex2, ey1:ey2]

			SIMS[cnt,:,:,:] = MAP

		end = time.time()

		logline = '\t\tmap cutouts extracted. time taken = %s seconds. next is COV calc.' %(end-start)
		logfile = open(self.log_file,'a')
		logfile.writelines('%s\n' %(logline))
		logfile.close()
		print logline

		return SIMS

	def fn_apodize(self,MAP,mapparams):

		import numpy as np, scipy.ndimage as ndimage, scipy.signal as signal

		nx, ny, dx, dy = mapparams

		pix = dx
	        radius = (nx * pix)/10.
	        npix_cos = int(radius/pix)
                ker=np.hanning(npix_cos)
	        ker2d=np.asarray( np.sqrt(np.outer(ker,ker)) )

		#npix_bord = 2.*int(wdth_bord/reso_arcmin)
		#xs, ys    = np.meshgrid( np.linspace(-1., 1., npix_bord), np.linspace(-1., 1., npix_bord) )

	        MASKf=np.zeros((nx,ny))
		minval, maxval = -(nx*pix)/2,  (nx*pix)/2
		x = y = np.linspace(minval, maxval, nx)
		X, Y = np.meshgrid(x,y)
		xc, yc = 0., 0.
		radius = (nx * dx/2) - 1.
		#inds=np.where((X-xc)**2. + (Y-yc)**2. <= radius**2.) #all in arcmins
	        inds=np.where((abs(X)<=radius) & (abs(Y)<=radius)) #all in arcmins
	        MASKf[inds]=1.

		apodMASKf=ndimage.convolve(MASKf, ker2d)#, mode='wrap')
		apodMASKf/=apodMASKf.max()

                #ker=signal.hanning(nx, sym=True)
	        #ker2d=np.asarray( np.sqrt(np.outer(ker,ker)) )
		#apodMASKf = np.copy(ker2d)
		#apodMASKf/=apodMASKf.max()

		apodMASKf = np.asarray( [apodMASKf for iii in range(len(MAP)) ] )

		#css=imshow(apodMASKf[0] * MAP[0]);css.set_clim(min(MAP[0].ravel()),max(MAP[0].ravel()));colorbar();show();quit()

		return apodMASKf * MAP


	def fn_convolve_using_QL(self,f,g): #both in fourier space

		import numpy as np

		scfns = scipy_useful_functions()

		if len(g)<>len(f):
			g = np.asarray([g for ii in range(len(f))] )

		array_reshaped = 0
		if len(f.shape)==2:
			array_reshaped = 1
			f = np.asarray([f] )
			g = np.asarray([g] )

		totiter = len(f)
		DUMMY_F = np.zeros(f.shape)

		for mm in range(totiter):
			#pad before convlution
			f_padded = scfns.pad_ft(f[mm])
			g_padded = scfns.pad_ft(g[mm])

			conv =  f_padded * g_padded #convolution

			conv_unpad = scfns.unpad_ft(conv) #unpad now

			DUMMY_F[mm] = np.copy(conv_unpad)

		if array_reshaped:
			return DUMMY_F[0]
		else:
			return DUMMY_F


	def fn_get_powerspectrum(self, MAP):
		
		import numpy as np
		import modules.radialProfile as radialProfile

		MAPFFT = np.fft.fftshift( np.fft.fft2(MAP) )
		MAPPSD = np.abs( MAPFFT )**2
		pspec1d = radialProfile.azimuthalAverage(MAPPSD)

		return pspec1d

	def fn_get_kappa_QE(self, OBSMAP, mapparams, Dls_len, Dls_unlen, noise_level):

		import numpy as np
		import modules.qe_funcs as qe
		###Parameters of QE### from EBX
		#l cuts
		nlens = 60000 #babbloo changing to 27535
		nbeam = 60000
		#what is this? -related to pixel size????
		ngaus = 27535
		nlensmin = 2
		#Filter for gradient
		#7777777777777777777777777777777777777777777777777777777
		ngrad = 2000#99999999999.#1500
		#what is this?
		#nclfilt = 5400


		self.Dls2map(Dls_len, mapparams, cls_lensed_lss = 1, no_sims_required = 1)
		self.Dls2map(Dls_unlen, mapparams, cls_lensed_lss = 0, no_sims_required = 1)

		ell_len = Dls_len[:,0]

		C_l_len_lmat = self.CLS[0]
		C_l_unl_lmat = self.CLS_unlen[0]

		lx, ly = self.get_lxly(mapparams)
		l2d = np.sqrt(lx**2. + ly**2.)
		transfer_lmat = self.fn_get_HPF(mapparams, ideal = 1)[0]
		#transfer_lmat = self.fn_get_HPF(mapparams, ideal = 0, minel=None, maxel=10000)[0]

		nx, ny, dx, dy = mapparams
		dxdy = dx * self.arcmins2radians * dy  * self.arcmins2radians

		beam_fwhm = self.exp_beam*(1./60.)*np.pi/180.     
		#beam_lmat = np.fft.fftshift( np.exp(-(beam_fwhm**2.)*(l2d**2.)/(16.*np.log(2.))) )
		beam_lmat = np.exp(-(beam_fwhm**2.)*(l2d**2.)/(16.*np.log(2.)))

		DeltaT = noise_level*1e-6*(1./60.)*(np.pi/180.) #K-radian
		N_lmat = (DeltaT**2.) + np.zeros(l2d.shape)

		beam_lmat_transfer_lmat = self.fn_convolve_using_QL(beam_lmat, transfer_lmat)
		N_beam_lmat = N_lmat/((beam_lmat_transfer_lmat)**2.)
		#N_beam_lmat = N_lmat/(beam_lmat**2.)

		"""
		subplot(221);imshow(beam_lmat);colorbar()
		subplot(222);imshow(transfer_lmat);colorbar()
		subplot(223);imshow(N_beam_lmat);colorbar()
		show();quit()
		"""

		#Gradient weighting and inverse variance weighting
		print '\t\tGenerating gradient and inv var weights...'
		#Gradient weights (W^TT)
		l_G = ngrad
		weight_gradient_lmat = C_l_unl_lmat/(C_l_len_lmat + N_beam_lmat)
		above_lg = np.where(l2d > l_G)
		weight_gradient_lmat[above_lg] = 0.0

		#Inverse variance weights (W^T)
		weight_inv_var_lmat = l2d /(C_l_len_lmat + N_beam_lmat)

		'''
		subplot(321);imshow(np.fft.ifftshift(np.fft.ifft2(C_l_unl_lmat).real)*1e12);colorbar()
		subplot(322);imshow(np.fft.ifftshift(np.fft.ifft2(C_l_len_lmat).real)*1e12);colorbar()
		subplot(323);imshow(np.fft.ifftshift(np.fft.ifft2(weight_gradient_lmat).real));colorbar()
		subplot(324);imshow(np.fft.ifftshift(np.fft.ifft2(C_l_unl_lmat * weight_gradient_lmat).real)*1e12);colorbar();
		subplot(325);imshow(np.fft.ifftshift(np.fft.ifft2(weight_inv_var_lmat).real));colorbar()
		subplot(326);imshow(np.fft.ifftshift(np.fft.ifft2(C_l_unl_lmat * weight_inv_var_lmat).real));colorbar();

		show();quit()
		'''

		'''
		np.random.seed(100)
		scalefac = np.sqrt(1./ (dx*self.arcmins2radians * dy*self.arcmins2radians))
		RANDOMREAL = np.fft.fft2( np.random.randn(nx,ny) ) * scalefac
		LENMAP = np.fft.ifft2(C_l_len_lmat**.5 * RANDOMREAL ).real
		UNLENMAP = np.fft.ifft2(C_l_unl_lmat**.5 * RANDOMREAL ).real 
		subplot(321);imshow(LENMAP);colorbar()
		subplot(322);imshow(UNLENMAP - LENMAP);colorbar()
		subplot(323);imshow(LENMAP - np.fft.ifft2( np.fft.fft2(LENMAP) * weight_gradient_lmat).real);colorbar()
		subplot(324);imshow(UNLENMAP - np.fft.ifft2( np.fft.fft2(UNLENMAP) * weight_gradient_lmat).real);colorbar()
		subplot(325);imshow(LENMAP - np.fft.ifft2( np.fft.fft2(LENMAP) * weight_gradient_lmat * weight_inv_var_lmat).real);colorbar()
		subplot(326);imshow(UNLENMAP - np.fft.ifft2( np.fft.fft2(UNLENMAP) * weight_gradient_lmat * weight_inv_var_lmat).real);colorbar()
		show();quit()
		'''

		'''
		subplot(131);imshow(OBSMAP);colorbar()
		subplot(132);imshow( np.fft.ifft2( np.fft.fft2(OBSMAP) * weight_gradient_lmat ).real );colorbar()
		subplot(133);imshow( np.fft.ifft2( np.fft.fft2(OBSMAP) * weight_inv_var_lmat ).real );colorbar()
		show();quit()
		'''


		lmax_inv_var = nlens+ngrad# - ensures that we're not using ells beyond range of input C_ls
		above_inv_var_lmax = np.where(l2d > lmax_inv_var)
		weight_inv_var_lmat[above_inv_var_lmax] = 0.0

		#normalization of quadratic estimator
		print '\t\tGenerating QE normalization...'
		norm_lmat = qe.get_norm(lx, ly, C_l_unl_lmat, weight_gradient_lmat, weight_inv_var_lmat, dxdy)
		#norm_lmat = norm_lmat*0. + 1.
		#imshow(norm_lmat.real);colorbar()
		#show();quit()

		#kappa_qe = qe.get_kappa_babbloo(lx, ly, weight_gradient_lmat, weight_inv_var_lmat, OBSMAP, norm_lmat, dxdy, ngaus)
		kappa_qe = qe.get_kappa(lx, ly, weight_gradient_lmat, weight_inv_var_lmat, OBSMAP, norm_lmat, dxdy, ngaus)

		#subplot(121);imshow(kappa_qe.real);colorbar()
		#subplot(122);imshow(kappa_qe_babbloo.real);colorbar();show()
		#show();quit()

		return kappa_qe


from pylab import *

