"""
### Author: Srinivasan Raghunathan
### CHECK - to be corrected
### Date start: 26.02.2016

# should be reorganised with several classes
class simulations
class lensing
class mapgrid
class likelihood
"""

import numpy as np, scipy as sc, os, sys, time
from scipy.interpolate import interp1d
from scipy import interpolate as intrp

class simulations():

	def __init__(self):

		self.Tcmb = 2.73 #K
		self.no_times_sigma_blown = 0 
		self.quiet = 0		
		self.degrees2radians = np.pi / 180.
		self.arcmins2radians = self.degrees2radians / 60.
		self.cmb_150ghz = []
		self.cmb_90ghz = []
		self.cmb_150ghz_after_fit = []
		self.fits = []
		self.amp_sz = []
		self.Bl = None
		#self.exp_beam = 1. #arcmins
		self.exp_beam_90ghz = 1.7 #arcmins
		#self.exp_beam_90ghz = 1.0 #arcmins
		self.cmbrandomseedval = None
		self.tqulen = 3
		self.cmb_noise_randomseedval = None
		self.amp_initial= []
		self.inidic = {}
		self.sz_map = []
		#to store cov infos.
		self.covdic = {}
		self.cmb_90 = None
		self.locdic = {'T':0,'Q':1,'U':2,'E':3,'B':4}

		self.npixels = None
		self.small_number = 1e-5
		self.covfiles = None

		self.map_centre = [0.0, -59.033333]
		self.el_lowpass = 11500		

		#TF
		self.tf_file = 'data/500sqdeg_henning/center_tf_for_sr.npy'
		self.TWODTF = None

		self.perform_EB = 0

		self.C_FG = None
		self.NOISE = None
		self.C_NOISE = None
		self.whichmap = None
#		self.expnoiselevel = np.asarray([5.5,7.7,7.7]) #uK-arcmin
#		self.expnoiselevel = np.asarray([1e-3,1e-3,1e-3]) #uK-arcmin

		"""
		Reichardt et al. 2012: http://arxiv.org/pdf/1111.0932v2.pdf
		Clustered and Poisson DFSG terms
		"""

		self.poisson_power = 7.54 # +/-0.38 uK2
		self.clustered_power = 6.25 # +/-0.52 uK2
		self.clustered_comp_gamma = 0.8
		self.all_clusters_at_same_redshift = 0

		#for tSZ
		self.add_tSZ_cov = 0 #setting this to 1 will add COV_tSZ to the COV matrix during likelihood computation

		self.C_for_regular = None

		#others
		self.inidic['perform_extra_mul'] = 0

		#tSZ stuff
		self.beta_model_90ghz = None
		self.beta_model_150ghz = None
		self.C_tSZ = None
		self.shaw_theta_int = 2.5 #arcmins
		self.shaw_theta_c = 0.75 #arcmins
		self.shaw_beta = 1.
		self.sigma_x = []
		self.sigma_y = []
		self.offset_x = []
		self.offset_y = []
		self.theta = []
		self.amplitude = []
		self.amp_offset = []
		self.stacked_fitted_cmb = np.zeros((200,200))# hard_coded

	def _ini_extra_facs(self, extra_division_factor_for_evals_det):

		self.extra_division_factor_for_evals_det = extra_division_factor_for_evals_det
		

	def _ini_single_params(self,param_name, param_val):
		cmd = '%s = %s' %(param_name, param_val)
		exec(cmd)

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

	def fn_get_param_dict_from_paramfile(self, paramfile):
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

		return param_dict

	def fn_get_lxly_az_angle(self,lx,ly):

		return 2*np.arctan2(lx, -ly)
		#return 2.*np.arctan2(ly, lx)


	def get_lxly(self,mapparams):

		nx, ny, dx, dy = mapparams
		dx *= self.arcmins2radians
		dy *= self.arcmins2radians

		lx, ly = np.meshgrid( np.fft.fftfreq( nx, dx ), np.fft.fftfreq( ny, dy ) )
		lx *= 2* np.pi
		ly *= 2* np.pi

		return lx, ly

	def fn_convert_crit_mean(self, mass, z, param_dict, crit_to_mean = 1):

		from astropy.cosmology import FlatLambdaCDM

		if 'ombh2' not in param_dict:
			param_dict['ombh2'] = 0.022
		omb0 = param_dict['ombh2']/param_dict['h']**2.

		cosmo = FlatLambdaCDM(H0 = param_dict['h']*100., Om0 = param_dict['omega_m'], Ob0 = omb0)
		rho_crit_z = cosmo.critical_density(z).value
		rho_mean_z = cosmo.Om(z)*cosmo.critical_density(z).value
		factor = rho_mean_z/rho_crit_z

		if not crit_to_mean: factor = 1./factor

		return mass * factor

	def fn_random_sampler(self, x, y, howmanysamples = 100000, burn_in = 5000):
		import scipy.integrate as integrate
		import scipy.interpolate as interpolate

		#plot(x, y);show()
		norm = integrate.simps(y, x) #area under curve for norm
		y = y/norm #normalise dn/dM here

		cdf = np.asarray([integrate.simps(y[:i+1], x[:i+1]) for i in xrange(len(x))])
		#print len(cdf)
		#plot(cdf);show()
		cdf_inv = interpolate.interp1d(cdf, x)

		random_sample = cdf_inv(np.random.rand(howmanysamples))
		#plot(random_sample);show()#;sys.exit()

		return random_sample[burn_in:]	

	def fn_get_width_from_sampling(self, x, likelihood_curve):#, sigma_value = [1.]):
		randoms = self.fn_random_sampler(x, likelihood_curve)
		mean_mass = x[np.argmax(likelihood_curve)]
		low_err = mean_mass - np.percentile(randoms, 16.)
		high_err = np.percentile(randoms, 84.) - mean_mass
		return mean_mass, low_err, high_err


	def fn_radial_profile(self, Z, XY, bin_size = 1., minbin = 0., maxbin = 10., to_arcmins = 1):

		Z = np.asarray(Z)
		if XY == None:
			Y, X = np.indices(image.shape)
		else:
			X, Y = XY

		#RADIUS = np.hypot(X,Y) * 60.
		RADIUS = ( (X**2. + Y**2.) ** 0.5 )
		if to_arcmins: RADIUS *= 60.
		#imshow(RADIUS);colorbar();show();quit()

		#ind = np.argsort(RADIUS.flat)
		#print ind;quit()

		binarr=np.arange(minbin,maxbin,bin_size)
		radprf=np.zeros((len(binarr),3))

		hit_count=[]
		#imshow(RADIUS);colorbar()

		for b,bin in enumerate(binarr):
		        ind=np.where((RADIUS>=bin) & (RADIUS<bin+bin_size))
		        radprf[b,0]=(bin+bin_size/2.)
			hits = len(np.where(abs(Z[ind])>0.)[0])

			if hits>0:
				radprf[b,1]=np.sum(Z[ind])/hits
				radprf[b,2]=np.std(Z[ind])
		        hit_count.append(hits)


		hit_count=np.asarray(hit_count)
		std_mean=np.sum(radprf[:,2]*hit_count)/np.sum(hit_count)
		errval=std_mean/(hit_count)**0.5
		radprf[:,2]=errval

		#errorbar(radprf[:,0], radprf[:,1], yerr = [radprf[:,2], radprf[:,2]], color= 'g', marker = 'o')
		#show();quit()
		
		return radprf

	#def fn_offsets_correction(self, mapparams, MODEL, param_dict, redshift = None, fcen = 0.78, CLUS_IDENTIFIER = None):
	def fn_offsets_correction(self, mapparams, MODEL, param_dict, redshift = None, richval = None, totalclus = None, fcen_default_des = 0.78, CLUS_IDENTIFIER = None):

		"""
		fcen = 0.78 #page 14 of https://arxiv.org/pdf/1601.00621.pdf (Rykoff 2016)
		"""
		#20180109: add offsets correction here
		MODEL_FFT = np.fft.fft2( MODEL )

		try:
			fcen = self.inidic['fcen']
		except:
			fcen = fcen_default_des

		try:
			lncmis = self.inidic['lncmis']
			R_lambda = (richval/100.)**0.2 / self.inidic['h']
			offset_scatter = np.exp(lncmis) * R_lambda
			Roff_dist = np.random.rayleigh(offset_scatter, totalclus) ## * theta
			y, x = np.histogram(Roff_dist, bins = 100)
			D_a_sigma_c = self.fn_random_sampler(x[:-1], y, howmanysamples = 1, burn_in = 0)[0]
			'''
			D_a_sigma_c = Roff_dist[np.random.randint(len(Roff_dist))]
			'''
		except:
			D_a_sigma_c = 0.42/param_dict['h'] #0.42h-1 Mpc from 1707.09369

		fmis = 1 - fcen
		if redshift==None:
			CLUS_REDSHIFTS = CLUS_IDENTIFIER[:,2]
			redshift = np.median( CLUS_REDSHIFTS )

		from astropy.cosmology import FlatLambdaCDM
		ombh2 = 0.022
		omb0 = ombh2/param_dict['h']**2.
		cosmo = FlatLambdaCDM(H0 = param_dict['h']*100., Om0 = param_dict['omega_m'], Ob0 = omb0)

		D_a = (cosmo.comoving_distance(redshift)/(1.+redshift) ).value

		lx, ly = self.get_lxly(mapparams)
		L = np.sqrt( lx ** 2. + ly ** 2. )
		sigma_c = D_a_sigma_c/D_a
		##from IPython import embed; embed()
		HOW_MUCH_SMEARING = ( fcen + fmis * np.exp( -0.5 * sigma_c**2. * L**2. ) )
		MODEL_MISCENTERED_FFT = HOW_MUCH_SMEARING * MODEL_FFT
		MODEL_MISCENTERED = np.fft.ifft2( MODEL_MISCENTERED_FFT ).real

		"""
		subplot(131);imshow(MODEL, vmin = -0.02, vmax = 0.2);colorbar()
		subplot(132);imshow(MODEL_MISCENTERED, vmin = -0.02, vmax = 0.2);colorbar()
		subplot(133);imshow(MODEL - MODEL_MISCENTERED);colorbar()
		show()
		"""

		MODEL = np.copy( MODEL_MISCENTERED )

		##imshow(MODEL, interpolation = 'bicubic', cmap = cmap_planck);axhline(30.); axvline(30.);colorbar();show();sys.exit()

		return MODEL

	def fn_get_STACKED_KAPPA_MODEL_FN_MASS(self, AA, alpha_fit, RA, DEC, CLUS_IDENTIFIER, param_dict, mapparams, Bl = None, TWODTF = None, clra = 0., cldec = 0., profile_name = 'NFW', kappa_model_fixed = 0, null_test = 0, lambda_pivot = 30., z_pivot = 0.5, beta_fit = 0.18, theta_max = -1., weights = None, use_TF_beam = 0, rad_profile_info = None, perform_offsets_correction = 1):

		zvals = np.asarray( map(lambda x: x[2], CLUS_IDENTIFIER) )
		richvals = np.asarray( map(lambda x: x[3], CLUS_IDENTIFIER) )
		if len( np.unique(zvals) ) == 1 and len( np.unique(richvals) ) == 1: 
			kappa_model_fixed = 1
		
		
		STACKED_KAPPA_MODEL_ARR = []
		#now loop over all clusters to get a stacked kappa model
		for kcnt, cluskey in enumerate(CLUS_IDENTIFIER):

			if not null_test:
				try:
					ra, dec, z_val, rich, weight, weights_norm = cluskey
				except:
					ra, dec, z_val, rich, weight = cluskey

			else:
				ra, dec = cluskey
				z_val = 0.01#55 #some median value

			#Use M-lambda relation to get M here
			#model: M = A * (lambda / lambda_pivot)**alpha * ( (1+z) / (1+z_pivot) )**beta
			if not kappa_model_fixed:
				
				MM = AA * (rich/lambda_pivot)**alpha_fit * ( (1+z_val) / (1+z_pivot) )**beta_fit
			else:
				MM = AA

			cc = self.c_Duffy(MM, z_val, param_dict['h'], profile_name = profile_name)

			#print MM, AA, z_val, cc

			'''
			#print kcnt, cluskey		
			logline = '\tkcnt = %s, cluskey = %s, dervied concentration = %.3f' %(kcnt, cluskey, cc)
			logfile = open(log_file,'a');logfile.writelines('%s\n' %(logline));logfile.close()
			print logline
			'''



			#from IPython import embed;embed()
			nx,ny,dx,dy = mapparams
			boxsize = param_dict['boxsize']
			#from IPython import embed;embed()
			
			nx= ny = int(boxsize/dx)
			
			minval, maxval = clra-boxsize/2/60.,  clra+boxsize/2/60.
			ra = np.linspace(minval, maxval, nx)

			minval, maxval = cldec-boxsize/2/60.,  cldec+boxsize/2/60.
			dec = np.linspace(minval, maxval, ny)

			RA, DEC = np.meshgrid(ra,dec)
			mapparams = nx,ny,dx,dy
			simmapparams = 200,200,0.5,0.5
			
			
			self.fn_lensing_ini(mapparams,param_dict, RA, DEC, [clra], [cldec], [MM], [cc], [z_val], param_dict['z_lss'], param_dict['mass_def'], param_dict['rho_def'], truncate_kappa_at_radius = param_dict['truncate_kappa_at_radius'], profile_name = profile_name, theta_max = theta_max)
			
			kappa_model = self.KAPPA
			

			if use_TF_beam:

				try:
					above_beam_scale = np.where(L>=ngaus)
				except:
					lx, ly = self.get_lxly(mapparams)
					L = (lx**2. + ly**2.)**0.5
					ngaus = int( np.sqrt(8. * np.log(2.)) / np.radians(self.exp_beam/60.) )
					above_beam_scale = np.where(L>=ngaus)
				'''
				beam_tf_lmat = np.copy(Bl*TWODTF)
				to_deconvolve_lmat = np.copy(beam_tf_lmat)
				bad = np.where(to_deconvolve_lmat==0.)
				to_deconvolve_lmat[bad]= 1.
				deconvolve_lmat = 1./to_deconvolve_lmat
				deconvolve_lmat[bad] = 0.
				highly_filtered = np.where( (deconvolve_lmat==0) | (deconvolve_lmat>=1e10) | (deconvolve_lmat<0.) )
				'''
				
				beam_tf_decov = 1./np.copy(Bl * TWODTF)
				highly_filtered= np.where( (beam_tf_decov==np.inf) | (beam_tf_decov>=1e10) | (beam_tf_decov<0.) )
				kappa_model_fft = np.fft.fft2(kappa_model)
				kappa_model_fft[above_beam_scale] = 0.
				kappa_model_fft[highly_filtered] = 0.
				kappa_model = np.fft.ifft2( kappa_model_fft ).real
				kappa_model = kappa_model[nx/2 - 30: nx/2+30,ny/2-30:ny/2 + 30]

			#from IPython import embed;embed()
			#from IPython import embed;embed()
			#sys.exit()
	
			if perform_offsets_correction:
				kappa_model_ori = np.copy(kappa_model)
				###from IPython import embed; embed()
				kappa_model = self.fn_offsets_correction(mapparams, kappa_model, self.inidic, redshift = z_val, richval = rich, totalclus = self.totalclus_for_analysis)
				if (0):
					if AA==0:continue
					##from IPython import embed; embed()
					subplot(121);imshow(kappa_model);colorbar()
					kappa_model = self.fn_offsets_correction(mapparams, kappa_model, self.inidic, redshift = z_val, richval = rich)
					subplot(122);imshow(kappa_model);colorbar();show();sys.exit()


			STACKED_KAPPA_MODEL_ARR.append( kappa_model )

			if kappa_model_fixed == 1 and (len( np.unique(zvals) ) == 1 and len( np.unique(richvals) ) == 1):
				return kappa_model

		STACKED_KAPPA_MODEL_ARR = np.asarray( STACKED_KAPPA_MODEL_ARR )
		if self.is_seq(rad_profile_info): #20180611: then return radial profile vectors of all cluster models
			binsize, minbin, maxbin = rad_profile_info
			RADEC = [RA,DEC]
			raprf_models = np.asarray( map(lambda x: self.fn_radial_profile(x, RADEC, bin_size=binsize, minbin=0.0, maxbin=maxbin)[:,1], STACKED_KAPPA_MODEL_ARR ) )

			return raprf_models


		#average of all clusters
		STACKED_KAPPA_MODEL = np.mean( STACKED_KAPPA_MODEL_ARR, axis = 0)
		if self.is_seq(weights):
			STACKED_KAPPA_MODEL_INV_VAR_WEIGHTED_SUM = np.sum( np.asarray( map(lambda x, y: x * y, STACKED_KAPPA_MODEL_ARR, weights) ), axis = 0)
			WEIGHTS_SUM = np.sum(weights)
			STACKED_KAPPA_MODEL_INV_VAR_WEIGHTED_MEAN = STACKED_KAPPA_MODEL_INV_VAR_WEIGHTED_SUM / np.sum(weights)

		#return STACKED_KAPPA_MODEL
		STACKED_KAPPA_MODEL_SUM = np.sum( STACKED_KAPPA_MODEL_ARR, axis = 0)
		totalsummed = len(STACKED_KAPPA_MODEL_ARR)

		if self.is_seq(weights):		
			return STACKED_KAPPA_MODEL, STACKED_KAPPA_MODEL_SUM, totalsummed, STACKED_KAPPA_MODEL_INV_VAR_WEIGHTED_MEAN, STACKED_KAPPA_MODEL_INV_VAR_WEIGHTED_SUM, WEIGHTS_SUM
		else:
			return STACKED_KAPPA_MODEL, STACKED_KAPPA_MODEL_SUM, totalsummed


	def gauss_beam(self, fwhm, EL, pol=1): #from Healpy
		"""
		Imported from healpy gauss_beam()
		fhwm = experiment fwhm in radians
		EL = 2d \ell matrix
		pol = 1 #includes pol.factor
		"""

		sigma = fwhm / np.sqrt(8. * np.log(2.))
		sigma2 = sigma ** 2
		if not self.is_seq(EL):			
			ell = np.arange(EL + 1)
		else:
			ell = EL

		g = np.exp(-.5 * EL**2. * sigma2)

		if not pol: # temperature-only beam
			return g
		else: # polarization beam
			# polarization factors [1, 2 sigma^2, 2 sigma^2, sigma^2]
			#pol_factor = np.exp([0., 2*sigma2, 2*sigma2, 2*sigma2, 2*sigma2])#, sigma2]) #ignoring T
			pol_factor_tmp = np.tile(2*sigma2,self.tqulen)
			pol_factor_tmp[0] = 0.
			pol_factor = np.exp(pol_factor_tmp)#, sigma2]) #ignoring T
			G = np.asarray([g * scale for scale in pol_factor])
			return G

        def rotate_tf(self, for_rot, ra, dec, ra0 = None, dec0 = None, in_elspace = 1, return_angle = 0, inv_rot_angle = 0):
                '''
                map_for_rot is transfer funciton / map evaluated at the center of the map
                ra0, dec0 specify the map center
                ra, dec are the coordinates at which we want to evaluate the TF
                '''
                from scipy import ndimage          

                if not in_elspace:
                        for_rot = np.fft.fft2(for_rot)

		else:
			for_rot = np.fft.fftshift(for_rot)

                tf_rot_angle = self.get_tf_rot_angle(ra,dec,ra0,dec0) #radians
		#if inv_rot_angle:
		#	tf_rot_angle = 360.- tf_rot_angle

                rotated_real = ndimage.interpolation.rotate(for_rot.real, np.degrees(tf_rot_angle), reshape = False, mode = 'reflect')
                rotated_imag = ndimage.interpolation.rotate(for_rot.imag, np.degrees(tf_rot_angle), reshape = False, mode = 'reflect')

		rotated = rotated_real + 1j * rotated_imag

                if not in_elspace:
                        rotated = np.fft.ifft2( rotated )
		else:
			rotated = np.fft.fftshift(rotated)

		if not return_angle:
	                return rotated.real
		else:
	                return rotated.real, np.degrees(tf_rot_angle)

        def get_tf_rot_angle(self, ra, dec, ra0 = None, dec0 = None, proj = 5):
                """
                Author: Srini Raghunathan
                Modified: Eric Baxter
                
                this is only valid for ZEA (proj 5) projection
                http://arxiv.org/pdf/1111.7245v1.pdf - section 5.4 for reference
                (although note that formulae therein seem to have an error)
                this is roughly the angle required to convert the fourier space filtering to real space filtering
                Should match output of get_proj5_angle.pro in the SPT analysis repo
                """

                if ra0 == None:
                        ra0, dec0 = np.radians(self.map_centre)
		else:
                        ra0, dec0 = np.radians(ra0), np.radians(dec0)
                
                ra = np.radians(ra)
                dec = np.radians(dec)
                
                #do the 90deg. subtraction for declinations
                dec = np.radians(90.) - dec
                dec0 = np.radians(90.) - dec0

		#SPT idl module has ra - ra0 instead - so changing that
		tmp, tmp0 = ra, ra0
		ra0, ra = tmp, tmp0
		
                
	        #nr - numerator; dr - denominator; _1 - first term; _2 - second term, etc. in the equation.
	        gamma_dr = 1 + ( np.cos(dec0) * np.cos(dec) ) + ( np.sin(dec0) * np.sin(dec) * np.cos(ra0 - ra) )
	        gamma = 0.5 / gamma_dr
	        
	        A_1 = gamma * np.sin(dec0) * np.sin(ra - ra0) * ( ( np.sin(dec0)*np.cos(dec) ) - ( np.sin(dec)*np.cos(dec0)*np.cos(ra - ra0) ) )
	        A_2 = np.cos(dec0) * np.sin(ra-ra0)
	        A = A_1 + A_2

	        B_1 = gamma * np.sin(dec0) * (np.sin(ra-ra0))**2. * np.sin(dec)
	        B_2 = np.cos(ra-ra0)
	        B = B_1 + B_2

	        #alpha = np.arctan(-1.*A/B)
		alpha = np.arctan2(A,B)

                return alpha #in radians                


	def fn_get_HPF(self, mapparams, minel = 500, maxel = 1000, ideal = 0, all_iso = 0, just_lx = 1): #Based on Eric's code


		if minel == None:
			minel = 500
		if maxel == None:
			maxel = 1000

		nx, ny, dx, dy = mapparams
		dx *= self.arcmins2radians
		dy *= self.arcmins2radians

		lx, ly = self.get_lxly(mapparams)
		L = np.sqrt(lx**2. + ly**2.)
		transfer_lmat =1.+ np.zeros(L.shape)
		
		if ideal:

			transfer_lmat = np.asarray( [transfer_lmat for ii in range(self.tqulen)] )
			return transfer_lmat

		if all_iso:
			iso_hipass = np.where(L < minel)
			iso_lowpass = np.where(L > maxel)
			transfer_lmat[iso_hipass] = 0.
			transfer_lmat[iso_lowpass] = 0.

		elif just_lx:
			scan_hipass = np.where( (np.abs(lx) < minel) )
			scan_lopass = np.where( (np.abs(lx) > maxel) )
			transfer_lmat[scan_hipass] = 0.
			transfer_lmat[scan_lopass] = 0.
		else:
			iso_hipass = np.where(L < minel)
			scan_hipass = np.where( (np.abs(lx) < minel) )
			scan_lopass = np.where( (np.abs(lx) > maxel) )
			transfer_lmat[iso_hipass] = 0.
			transfer_lmat[scan_hipass] = 0.
			transfer_lmat[scan_lopass] = 0.

		lmax_nyquist = 2.*np.pi*1./dx/2.*(nx-1.)/nx
		nyquist_lopass = np.where(L > lmax_nyquist)
		transfer_lmat[nyquist_lopass] = 0.0

		#transfer_lmat = np.fft.fftshift(transfer_lmat)
		#imshow(np.fft.fftshift(transfer_lmat), extent = [np.min(lx), np.max(lx), np.min(ly), np.max(ly)] );colorbar();show();quit()

		transfer_lmat = np.asarray( [transfer_lmat for ii in range(self.tqulen)] )

		return transfer_lmat

	def fn_get_EBAX_2016_anal(self, mapparams, l1 = 500, l2 = 400, l3 = 15000):


		nx, ny, dx, dy = mapparams
		dx *= self.arcmins2radians
		dy *= self.arcmins2radians

		lx, ly = self.get_lxly(mapparams)
		L = np.sqrt(lx**2. + ly**2.)
		transfer_lmat =1.+ np.zeros(L.shape)

		ideal_F =1.+ np.zeros(L.shape)
		if l1<>None:
			F1 = np.exp( -(l1/L)**6. )
		else:
			F1 = np.copy(ideal_F)
		if l2<>None:
			F2 = np.exp( -(l2/lx)**6. )
		else:
			F2 = np.copy(ideal_F)

		if l3<>None:
			F3 = np.exp( -(lx/l3)**6. )
		else:
			F3 = np.copy(ideal_F)


		transfer_lmat = F1 * F2 * F3

		transfer_lmat = np.asarray( [transfer_lmat for ii in range(self.tqulen)] )

		return transfer_lmat

	def fn_get_SPTpol_TWODTF(self, mapparams, TWODTF_FULL = None, tf_file = None, dx_full = None):

		"""
		This function extracts the 2dTF for the SPTpol centre field
		1. 2dTF for SPTpol field provided by Tom Crawford
		2. it is the average of TF over different fields
		3. ideal one wants to rotate this based on ra, dec position
		4. 
		.
		.
		.
		n. look into testing/tf_stuff.py
		"""
		nx,ny,dx,dy = mapparams
		if not self.is_seq(TWODTF_FULL):
			if tf_file == None:
				tf_file = self.tf_file
		
			#read the TF
			TWODTF_FULL = np.fft.fftshift( np.load(tf_file) )

		ny_full, nx_full = TWODTF_FULL.shape
		if dx_full == None:
			dx_full = dy_full = dx
		else:
			dy_full = dx_full
		mapparams_full = [nx_full, ny_full, dy_full, dx_full]
		lx_full, ly_full = self.get_lxly(mapparams_full)
		lx, ly = self.get_lxly(mapparams)

		lx_full = np.fft.fftshift( lx_full )
		ly_full = np.fft.fftshift( ly_full )
		lx = np.fft.fftshift( lx )
		ly = np.fft.fftshift( ly )

		#this interpolation method is a lot better than extraction using ifft
		TWODTF = np.fft.fftshift( intrp.RectBivariateSpline( ly_full[:,0],lx_full[0,:],TWODTF_FULL,kx = 5,ky = 5).ev(ly.flatten(),lx.flatten()).reshape([ny,nx]) )

		'''
		#pcolor(lx, ly, TWODTF);colorbar();show();quit()
		imshow(TWODTF,extent=[np.min(lx),np.max(lx),np.min(ly),np.max(ly)]);colorbar();
		xlim(-15e3,15e3);ylim(-15e3,15e3);show();quit()
		'''

		#finally np.tile this guy to match self.tqulen
		TWODTF = np.asarray( [TWODTF for iii in range(self.tqulen) ] )

		return TWODTF

		#ifft to get the real TF and pick the centre
		TWODTF_FULL_REAL = np.fft.ifftshift(np.fft.ifft2(TWODTF_FULL))

		#extracting the centre TF
		#cen_y, cen_x = np.asarray(tf.shape)/2
		import sky_local
		pixels = sky_local.ang2Pix(self.map_centre, self.map_centre, 0.5, TWODTF.shape, proj = 5., round=True, bin_center_zero=True, return_validity=False, use_c_code=False)
		cen_y, cen_x = pixels
		cen_y, cen_x = cen_y[0], cen_x[0]

		#central region will be
		nx,ny,dx,dy = mapparams
		y1, y2 = cen_y - ny/2, cen_y +nx/2
		x1, x2 = cen_x - ny/2, cen_x +nx/2

		#EXTRACTING NOW
		TWODTF_REAL = np.fft.ifftshift(TWODTF_FULL_REAL[y1:y2, x1:x2])
		
		#fft back 
		TWODTF = np.fft.fft2(TWODTF_REAL)

		#plot and make sure it looks okay
		lx,ly = self.get_lxly(mapparams)
		#print np.min(lx),np.max(lx),np.min(ly),np.max(ly)
		imshow(np.fft.fftshift(TWODTF.real),extent = [np.min(lx),np.max(lx),np.min(ly),np.max(ly)]);colorbar();show()#;quit()
		
		#finally np.tile this guy to match self.tqulen
		TWODTF = np.asarray( [TWODTF for iii in range(self.tqulen) ] )
		print TWODTF.shape;quit()

		return TWODTF


	def Cls2CLS(self,Cls,mapparams): #based on E. Baxter's code

		nx, ny, dx, dy = mapparams
		dx *= self.arcmins2radians
		dy *= self.arcmins2radians

		# find how many different Cls were passed
		noofrows = len(Cls)
		if self.is_seq_of_seq(Cls):
			noofcols = np.shape(Cls)[1]
			Cls_tot = noofcols - 1 #first is els, then Cls
		else:
			els = np.arange(noofrows)
			Cls = np.array( [els,Cls] ).T
			Cls_tot = 1

		lx, ly = self.get_lxly(mapparams)
		L = np.sqrt(lx**2. + ly**2.)

		#loglog(Cls[:,0], (Cls[:,0] * (Cls[:,0]+1)) * Cls[:,1])
		#loglog(Cls[:,0], (Cls[:,0] * (Cls[:,0]+1)) * abs(Cls[:,2]))
		#loglog(Cls[:,0], (Cls[:,0] * (Cls[:,0]+1)) * abs(Cls[:,3]))
		#show();quit()

		# processing Cls now
		CLS = np.zeros( (Cls_tot,L.shape[0],L.shape[1]) )
		for clcnt in range(Cls_tot):
			#CLS[clcnt,:,:] = np.interp(L.flatten(), Cls[:,0], abs(Cls[:,clcnt+1]), right = 0.).reshape(L.shape)
			CLS[clcnt,:,:] = np.interp(L.flatten(), Cls[:,0], abs(Cls[:,clcnt+1]), right = 0.).reshape(L.shape) 
			#CHECK - changed from abs(Cls) to Cls - should be okay and removed right=0.
		return CLS

	def Dls2map(self, Dls, mapparams, passing_Dls = 1, CMB_outputscale = 1, return_Clmat = 0): #Handles both T and P.

		#takes Dls as input and converts into CAMBMAPFFT

		###########################################################################
		#first check if only T or P is supplied as well; also obtain els
		noofrows = len(Dls)
		self.mapparams = mapparams
		self.Dls = Dls

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
		self.els = els
		'''
		clf()
		ax = subplot(111, yscale='log', xscale='log')
		plot(Dls[:,0], Dls[:,1]  * 1e12, 'k')
		plot(Dls[:,0], Dls[:,2]  * 1e12, 'r')
		plot(Dls[:,0], Dls[:,3]  * 1e12, 'g')

		show();quit()
		'''

		if passing_Dls:
			if CMB_outputscale == 1:
				Cls = ( self.Tcmb**2. * Dls * 2 * np.pi ) / ( els[:,None] * (els[:,None] + 1) )
			else:
				Cls = ( Dls * 2 * np.pi ) / ( els[:,None] * (els[:,None] + 1) )
			Cls[:,0] = els
		else:
			Cls = Dls

		self.Cls = Cls
		###########################################################################
		#Cls2map
		nx, ny, dx, dy = mapparams
		dx *= self.arcmins2radians
		dy *= self.arcmins2radians

		#self.scalefac = np.sqrt((nx * ny) / (dx * dy))
		self.scalefac = np.sqrt(1./ (dx * dy))

		CLS = self.Cls2CLS(Cls,mapparams)

		if return_Clmat:
			return CLS

		#### Wiener filter by Wayne Hu
		#self.WIENER_FILTER_WHU = self.fn_weighting_using_Wiener_filter_as_Wayne_Hu(mapparams, CLS)

		self.CLS = CLS #just store for future use
		CLS = np.sqrt(CLS)
                for ccc in range(len(CLS)):
                        CLS[ccc][np.isnan(CLS[ccc])] = 0.


		lx, ly = self.get_lxly(mapparams)
		angle = self.fn_get_lxly_az_angle(lx,ly)
		########################################################################
		########################################################################
		# processing Cls now
		if Dls_tot>1:
			#TFFT, EFFT, BFFT = CLS
			TFFT, EFFT, BFFT, TEFFT = CLS
		else:
			TFFT = CLS

		#T = np.fft.ifft2( TFFT * self.scalefac ).real

		if Dls_tot>1:
			FORQ = ( np.cos(angle) * EFFT - np.sin(angle) * BFFT )
			FORU = ( np.sin(angle) * EFFT + np.cos(angle) * BFFT )

			"""removed on 20161201
			FORQ[np.isnan(FORQ)] = 1e-10
			FORU[np.isnan(FORU)] = 1e-10
			"""

			#CAMBMAP_FFT = np.asarray([TFFT,FORQ,FORU]) * self.scalefac
			#self.cambmapvars = ['T','Q','U']
			#CAMBMAP_FFT = np.asarray([TFFT,FORQ,FORU,EFFT, BFFT]) * self.scalefac
			#self.cambmapvars = ['T','Q','U','E','B']
			CAMBMAP_FFT = np.asarray([TFFT, EFFT, BFFT, TEFFT]) * self.scalefac
			self.cambmapvars = ['T','E','B','TE']

		else:
			CAMBMAP_FFT = np.asarray(TFFT) * self.scalefac
			self.cambmapvars = ['T']


		self.CAMBMAP_FFT = CAMBMAP_FFT
		#CAMBMAP_FFT = np.asarray([TFFT])

		###########################################################################
		
		if  len(CAMBMAP_FFT)>1:
			self.tqulen = len(CAMBMAP_FFT) - 1
		else:
			self.tqulen = len(CAMBMAP_FFT)

		logline = '\t\tDls converted to maps'
		logfile = open(self.log_file,'a')
		logfile.writelines('%s\n' %(logline))
		logfile.close()
		#print logline

		return


	def Dls2map_v1(self, Dls, mapparams, passing_Dls = 1, CMB_outputscale = 1, return_Clmat = 0): #Handles both T and P.

		#takes Dls as input and converts into CAMBMAPFFT

		###########################################################################
		#first check if only T or P is supplied as well; also obtain els
		noofrows = len(Dls)
		self.mapparams = mapparams
		self.Dls = Dls

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
		self.els = els
		'''
		clf()
		ax = subplot(111, yscale='log', xscale='log')
		plot(Dls[:,0], Dls[:,1]  * 1e12, 'k')
		plot(Dls[:,0], Dls[:,2]  * 1e12, 'r')
		plot(Dls[:,0], Dls[:,3]  * 1e12, 'g')

		show();quit()
		'''

		if passing_Dls:
			if CMB_outputscale == 1:
				Cls = ( self.Tcmb**2. * Dls * 2 * np.pi ) / ( els[:,None] * (els[:,None] + 1) )
			else:
				Cls = ( Dls * 2 * np.pi ) / ( els[:,None] * (els[:,None] + 1) )
			Cls[:,0] = els
		else:
			Cls = Dls

		###########################################################################
		#Cls2map
		nx, ny, dx, dy = mapparams
		dx *= self.arcmins2radians
		dy *= self.arcmins2radians

		#self.scalefac = np.sqrt((nx * ny) / (dx * dy))
		self.scalefac = np.sqrt(1./ (dx * dy))

		CLS = self.Cls2CLS(Cls,mapparams)

		if return_Clmat:
			return CLS

		if CLS.shape[0]>3:
			CLS = CLS[0:3]

		#### Wiener filter by Wayne Hu
		#self.WIENER_FILTER_WHU = self.fn_weighting_using_Wiener_filter_as_Wayne_Hu(mapparams, CLS)

		self.CLS = CLS #just store for future use

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

		#T = np.fft.ifft2( TFFT * self.scalefac ).real

		if Dls_tot>1:
			FORQ = ( np.cos(angle) * EFFT - np.sin(angle) * BFFT )
			FORU = ( np.sin(angle) * EFFT + np.cos(angle) * BFFT )

			FORQ[np.isnan(FORQ)] = 1e-10
			FORU[np.isnan(FORU)] = 1e-10

			CAMBMAP_FFT = np.asarray([TFFT,FORQ,FORU]) * self.scalefac
			self.cambmapvars = ['T','Q','U']
			#CAMBMAP_FFT = np.asarray([TFFT,FORQ,FORU,EFFT, BFFT]) * self.scalefac
			#self.cambmapvars = ['T','Q','U','E','B']
			#CAMBMAP_FFT = np.asarray([TFFT, EFFT, BFFT]) * self.scalefac
			#self.cambmapvars = ['T','E','B']

		else:
			CAMBMAP_FFT = np.asarray(TFFT) * self.scalefac
			self.cambmapvars = ['T']


		self.CAMBMAP_FFT = CAMBMAP_FFT
		#CAMBMAP_FFT = np.asarray([TFFT])

		###########################################################################
		

		self.tqulen = len(CAMBMAP_FFT)

		logline = '\t\tDls converted to maps'
		logfile = open(self.log_file,'a')
		logfile.writelines('%s\n' %(logline))
		logfile.close()
		print logline

		return

		noofsims = 10
		remove_largescale_modes = 0
		for simcnt in range(noofsims):
			if simcnt % 5000 == 0 and simcnt > 0:
				logline = '\t\t\t\tsimno: %s' %(simcnt)
				logfile = open(self.log_file,'a')
				logfile.writelines('%s\n' %(logline))
				logfile.close()
				print logline

			DUMMY = np.fft.fft2( np.random.randn(nx,ny) )

			if remove_largescale_modes: #removes large scale modes initially - they will not be lensed

				if zero_padded_conv:
					TMP_FFT_1 = CAMBMAP_FFT * DUMMY				
					TMP_FFT = self.fn_convolve_using_QL(TMP_FFT_1[0], TWODTF)
					TMP = np.fft.ifft2( TMP_FFT ).real
				else:
					TMP = np.fft.ifft2( CAMBMAP_FFT * TWODTF * DUMMY ).real
					#imshow(TMP[0,0]);colorbar();show();quit()
				
			else: #large scale modes are present and will be lensed
				TMP = np.fft.ifft2( CAMBMAP_FFT * DUMMY ).real
				#clf();imshow(TMP[0]);colorbar();show();quit()

			#CHECK this
			#TMP = TMP - np.mean(TMP)
			for tqucnt in range(self.tqulen):
				TMP[tqucnt] = TMP[tqucnt] - np.mean(TMP[tqucnt])

			'''
			subplot(131);imshow(TMP[0] * 1e6);colorbar()#;show();quit()
			subplot(132);imshow(TMP[1] * 1e6);colorbar()
			subplot(133);imshow(TMP[2] * 1e6);colorbar();show();quit()#looks okay here
			'''

			#print np.mean(TMP[0,0].ravel())

			#print 'Map variance is', np.var(TMP[0].ravel())#;quit()

			clf()
			#making sure E and B look good
			for cc in range(self.tqulen):
				ax = subplot(2,3,cc+1)
				css=imshow(TMP[cc] * 1e6);colorbar()

			Q, U = np.fft.fft2(TMP[1]),np.fft.fft2(TMP[2])
			E = np.fft.ifft2( np.cos(angle) * Q + np.sin(angle) * U ).real
			B = np.fft.ifft2( -np.sin(angle) * Q + np.cos(angle) * U ).real
			ax = subplot(2,3,5);css=imshow(E * 1e6);colorbar()
			ax = subplot(2,3,6);css=imshow(B * 1e6);colorbar()
			show();quit()

			###########################################################################
			"""
			#clf();ax=subplot(111,yscale='log')#, xscale='log')
			#plot(els,els * (els + 1 ) * Cls[:,1]); plot(els,els * (els + 1 ) * Cls[:,2]);plot(els,els * (els + 1 ) * Cls[:,3])
			### get PHI
			elsclphi = np.loadtxt('data/output_spt_r_0.07_scalCls.dat',usecols = [0,4])
			els, clphi = np.asarray( zip(*elsclphi) )
			clphi /= els**4.
			elsclphi[:,1] = clphi

			CLPHI = self.Cls2CLS(elsclphi,mapparams)
			CLPHI = np.sqrt(CLPHI)
			CAMBPHI_FFT = CLPHI  * self.scalefac

			TWODTF = self.fn_get_HPF(mapparams, minel = 400, maxel = 25500)[0]
			#PHI = np.fft.ifft2( CAMBPHI_FFT[0] * DUMMY * TWODTF[0]).real
			PHI_FFT = CAMBPHI_FFT[0] * DUMMY# * TWODTF[0]

			#now lens E with PHI
			E = TMP[1] * 1e6
			LENS_B = self.calc_lensing_b_first_order(mapparams, E, PHI_FFT, None )
			#subplot(121);imshow(E);colorbar()
			#subplot(122);imshow(PHI);colorbar();show();quit()
			"""
			###########################################################################

			#debug stuff - see if the CMB maps look okay here
			show_cmb_plots = 0
			if show_cmb_plots:
				TMP *= 1e6
				print min(TMP[0,0].ravel()), max(TMP[0,0].ravel())
				if self.tqulen == 1: #only T map
					subplot(121);css=imshow(TMP[0,0]);colorbar()
					subplot(122);hist(TMP[0,0].ravel(),bins=nx,color='g',histtype='step');show();quit()
				else:
					for ppp in range(self.tqulen): #T,Q,U,E,B maps 
						ax = subplot(2,3,ppp+1)
						css=imshow(TMP[ppp]);colorbar()
						if ppp>2:
							clf();css=imshow(TMP[ppp]);colorbar()
							X, Y = np.meshgrid(np.arange(nx),np.arange(ny))
							#UN = VN = np.arctan2(TMP[1],TMP[2]) * 0.5
							quiver(X, Y, UN, VN,headwidth=0,headlength=0)
							show();quit()

						ax = subplot(2,3,ppp+1+self.tqulen)
						hist(TMP[ppp].ravel(),bins=nx,color='g',histtype='step')
					show();quit()

			SIMMAPS[simcnt,:,:,:] = TMP

		end = time.time()

		logline = '\t\tGaussian sims done. time take = %s seconds' %(end-start)
		logfile = open(self.log_file,'a')
		logfile.writelines('%s\n' %(logline))
		logfile.close()
		print logline

		self.SIMMAPS = SIMMAPS


		return SIMMAPS

	def downsample_map(self, data, N=2): #from N.Whitehorn
		from numpy import average, split
		''' original from N.WHitehorn
		width = data.shape[0]
		height= data.shape[1]
		return average(split(average(split(data, width // N, axis=1), axis=-1), height // N, axis=1), axis=-1)
		'''
		height, width = data.shape
		return average(split(average(split(data, width // N, axis=1), axis=-1), height // N, axis=1), axis=-1)


	def fn_lensing_ini(self,mapparams,param_dict, RA, DEC, clus_ra, clus_dec, M_200_list, c_200_list, z_L_list, z_lss, mass_def, rho_def, truncate_kappa_at_radius = .5, min_theta_cutoff_arcmins = 3., only_return_lensing_stuffs = 0, profile_name = 'NFW', theta_max = -1.):

		"""
		truncate_kappa_at_radius [r_200] = truncate KAPPA at this radius - see http://arxiv.org/abs/1401.1216
		"""

		#import EBX_modules.nfw_kappa_funcs as EBX_nfw
		import nfw_kappa_funcs as EBX_nfw

		#MAP RESOL, EL stuffs
		nx, ny, dx, dy = mapparams
		dx *= self.arcmins2radians
		dy *= self.arcmins2radians
		lx, ly = self.get_lxly(mapparams)
		L = np.sqrt(lx**2. + ly**2.)

		L[L==0] = self.small_number

		############################################################################################################
		############################################################################################################
		"""
		#added on 20161231 - Einasto profile. Will add more profiles later
		"""

		#get the kappa map from EBX
		### KAPPA = EBX_nfw.get_NFW_kappa_fullmap(param_dict, RA, DEC, clus_ra, clus_dec, M_200_list, c_200_list, z_L_list, z_lss, mass_def, rho_def)

		KAPPA, NFW_DIC = EBX_nfw.get_NFW_kappa_fullmap(param_dict, RA, DEC, clus_ra, clus_dec, M_200_list, c_200_list, z_L_list, z_lss, mass_def, rho_def, ret_nfw_dic = 1, profile = profile_name, theta_max = theta_max)
		#get the kappa map from EBX
		#css=imshow(KAPPA);css.set_clim(0.,1.);colorbar();show();quit()

		#css=imshow(KAPPA);colorbar();show();quit()
		#make sure kappa does't have nan, inf, etc.
		KAPPA[np.isnan(KAPPA)] = self.small_number

		'''
		logline = '\t\tkappa, phi obtained for %s profile' %(profile_name)
		logfile = open(self.log_file,'a')
		logfile.writelines('%s\n' %(logline))
		logfile.close()
		print logline
		'''

		#truncate kappa at certain radius
		#http://arxiv.org/abs/1401.1216 using 1 * r_200 for now
		#truncate_kappa_at_radius = 1.
		if truncate_kappa_at_radius <> None: ##tochange - for list of M, z
			import modules.scl_cosmo as scl_cosmo

			cosmo = scl_cosmo.scl_cosmo_stuffs()
			cosmo.fn_update_cosmo_params(param_dict)
			r_200 = cosmo.fn_r200(M_200_list[0],z_L_list[0]) #metres
			distances = cosmo.fn_cosmo_dist(z_L_list[0]) #metres
			d_A = distances[3] #angular diameter distance
			theta_cutoff_arcmins_actual = np.degrees(r_200/d_A) * 60. #arcmins

			print r_200, distances[1], np.degrees(r_200/d_A) * 60.;quit()

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
			
		KAPPA[np.isnan(KAPPA)] = self.small_number
		KAPPA[KAPPA == 0.] = self.small_number
		self.KAPPA = KAPPA
		
		PHI_FFT = -2. * dx * dy * np.fft.fft2(KAPPA)/(L**2)
		#make sure kappa does't have nan, inf, etc.
		PHI_FFT[np.isnan(PHI_FFT)] = self.small_number

		DEF_X    = np.fft.ifft2(-1j * PHI_FFT * lx) / ( dx * dy )
		DEF_Y    = np.fft.ifft2(-1j * PHI_FFT * ly) / ( dx * dy )


		if only_return_lensing_stuffs:
			return KAPPA, PHI_FFT, DEF_X, DEF_Y

		#Angular info on sky
		theta_x_list = np.array(range(1,nx+1)) * dx - dx * 0.5 * (nx - 1.)
		if nx<>ny:
			theta_y_list = np.array(range(1,ny+1))*dy-dy*0.5*(ny - 1.)
			theta_x,theta_y = np.meshgrid(theta_x_list, theta_y_list)
		else:
			theta_x = np.tile(theta_x_list,(nx,1))
			theta_y = np.transpose(theta_x)

		'''
		show_plots = 0
		if show_plots:
			from pylab import *
			subplot(221);contourf(np.degrees(theta_x)*60., np.degrees(theta_y)*60., KAPPA.real);colorbar();subplot(222);contourf(np.degrees(theta_x)*60., np.degrees(theta_y)*60.,DEF_X.real);colorbar();subplot(223);contourf(np.degrees(theta_x)*60., np.degrees(theta_y)*60.,DEF_Y.real);colorbar();show();quit()
		'''

		to_evaluate_unlensed_theta_x = (theta_x + DEF_X).flatten().real
		to_evaluate_unlensed_theta_y = (theta_y + DEF_Y).flatten().real

		self.theta_x = theta_x
		self.theta_y = theta_y
		self.to_evaluate_unlensed_theta_x = to_evaluate_unlensed_theta_x
		self.to_evaluate_unlensed_theta_y = to_evaluate_unlensed_theta_y

		return

	def fn_beam_stuff(self, mapparams, use_beam = None, nu = 150, return_beam = 0, exp_beam = None, beam_error = 0):

		#############################################
		#############################################
		####beam stuff
                nx, ny, dx, dy = mapparams
                dx *= self.arcmins2radians
                dy *= self.arcmins2radians
                lx, ly = self.get_lxly(mapparams)
		EL = np.sqrt(lx**2. + ly**2.)

		if use_beam == None:
			use_beam = self.inidic['use_beam']

		Bl = None
		if use_beam == 1: #ideal Gaussian beam

			if exp_beam == None:
				if nu == 150:
					exp_beam = self.exp_beam
				elif nu == 90:
					exp_beam = self.exp_beam_90ghz

			exp_beam *= self.arcmins2radians
			if self.tqulen>1:
				Bl = self.gauss_beam(exp_beam,EL) #healpy based #in fourier space (or) el space
			else:
				Bl = self.gauss_beam(exp_beam,EL, pol=0) #healpy based #in fourier space (or) el space

		elif use_beam == 2: #SPTpol el space beam
			
			"""
			https://pole.uchicago.edu/sptpol/index.php/500dEETE_Beams#5.2F16.2F16
			/data57/jhenning/500dEETE_products/sptpol_beam_20160515.pkl
			"""
			if nu == 150:
				Bl = self.fn_read_SPTpol_el_space_beam(mapparams, beam_error = beam_error)
			elif nu == 90:
				Bl = self.fn_read_SPTpol_el_space_beam(mapparams, nu = 90, beam_error = beam_error)
				#print 'SPTpol beam for 90GHz is not implemented yet. Aborting now ..'
				#sys.exit(0)

		elif use_beam == 3: #SPTpol venus beam
			Bl = self.fn_read_SPTpol_venus_beam(mapparams)

		elif use_beam == 4: #ACTPol venus beam
			Bl = self.fn_read_ACTPol_el_space_beam(mapparams)

		'''
		Bl_gauss = self.gauss_beam(np.radians(self.exp_beam/60.),EL)[0]

		subplot(121);css=imshow(np.fft.fftshift( Bl ),extent = [np.min(lx),np.max(lx), np.min(ly), np.max(ly)]);css.set_clim(0.,1.);colorbar();title('SPTpol')#;grid(1,ls='solid')
		#subplot(121);css=contourf(lx,ly,Bl);css.set_clim(0.,1.);colorbar();title('SPTpol')#;grid(1,ls='solid')
		#contour(lx,ly, Bl_gauss )
		xlim(-16000,16000); ylim(-16000,16000)
		#nx, ny, dx, dy = mapparams
		#subplot(223);imshow(np.fft.ifftshift( np.fft.ifft2(Bl).real ), extent = [-nx * dx, nx * dx, -ny*dy, ny * dy]);colorbar();title('SPTpol ifft');grid(1,ls='solid')

		subplot(122);css=imshow(np.fft.fftshift( Bl_gauss ),extent = [np.min(lx),np.max(lx), np.min(ly), np.max(ly)]);css.set_clim(0.,1.);colorbar();title('Gaussian')#;grid(1,ls='solid')
		xlim(-16000,16000); ylim(-16000,16000)
		#nx, ny, dx, dy = mapparams
		#subplot(224);imshow(np.fft.ifftshift( np.fft.ifft2(Bl_gauss).real ), extent = [-nx * dx, nx * dx, -ny*dy, ny * dy]);grid(1,ls='solid')
		#colorbar();title('Gauss ifft')

		show();quit()
		print  self.inidic['use_beam'];quit()
		'''

		try:
			logline = '\t\tInitialize / get beam now. Requested beam type = %s; nu = %s' %(use_beam, nu)
			logfile = open(self.log_file,'a')
			logfile.writelines('%s\n' %(logline));logfile.close()
			print logline
		except AttributeError as errvals:
			pass

		if return_beam:
			return Bl	
		else:
			self.Bl = Bl

		#############################################
		#############################################

	def fn_get_foreground_power(self, which_fg, mapparams, noofsims, els = None, reqdbox = None,random_seed = None, return_cov = 0, perform_lensing = 0, nu1 = 150, nu2 = None, also_P = 0, pol_frac_per_cent = 0.16, just_Cls = 0, ign_DG_RG = 1, assert_fg_labels = 1, special = None):

		nx, ny, dx, dy = mapparams
		dx *= self.arcmins2radians
		dy *= self.arcmins2radians
		npixels = nx * ny
		self.scalefac = np.sqrt(1./ (dx * dy))

		which_fg = which_fg.replace('_Gaussian','')


		#assert which_fg == 'DG-Cl' or which_fg == 'DG-Po' or which_fg == 'RG' or which_fg == 'all' or which_fg == 'kSZ' or which_fg == 'tSZ'
		#possible_fgs = ['DG-Cl', 'DG-Po', 'RG', 'all', 'tSZ', 'kSZ']
		possible_fgs = ['all', 'tSZ', 'kSZ','DG-Cl','DG-Po','RG','tSZ-CIB','Total','CMB']
		#possible_fgs = ['all','DG-Cl','DG-Po']
		if assert_fg_labels:
			assert which_fg in possible_fgs
		
		try:
			self.george_fg_dic
		except:
			self.george_fg_dic = readsav('data/foreground/george_plot_bestfit_line.sav')

		if nu2 == None: nu2 = nu1

		if nu1 == 90: nu1 = 95
		if nu2 == 90: nu2 = 95
		
		fg_nus = np.asarray( [(95,95), (95,150), (95,220), (150,150), (150,220), (220,220)] )
		
		fg_labels = self.george_fg_dic['ml_dl_labels']
		fg_nu_ind = np.where((fg_nus[:,0] == nu1) & (fg_nus[:,1] == nu2))[0][0]

		fg_els = self.george_fg_dic['ml_l']
		if which_fg <> 'all':
			fg_lab_ind = np.where(fg_labels == which_fg)[0]
			fg_Dls = self.george_fg_dic['ml_dls'][fg_nu_ind][fg_lab_ind][0]
		else:
			fg_Dls_all = self.george_fg_dic['ml_dls'][fg_nu_ind]
			#clf();ax=subplot(111,xscale='log',yscale='log');plot(fg_els, self.george_fg_dic['ml_dls'][fg_nu_ind][1], label = 'CMB')
			for fgcnt, fg in enumerate(possible_fgs):

				########################
				#20170413 - to allow Gaussian realisations of DG, RG when needed. Normally they are set to zero with ign_DG_RG
				#if fg == 'all' or fg == 'DG-Cl' or fg == 'DG-Po' or fg == 'RG' or fg == 'tSZ-CIB':
				#	continue
				if fg == 'all'  or fg == 'tSZ-CIB': continue
				fg_lab_ind = np.where(fg_labels == fg)[0]
				#plot(fg_els, self.george_fg_dic['ml_dls'][fg_nu_ind][fg_lab_ind][0], label = fg)

				if ign_DG_RG and (fg == 'DG-Cl' or fg == 'DG-Po' or fg == 'RG'): continue
				if fg == 'Total' or fg == 'CMB': continue
				########################
				try:
				
					fg_Dls += self.george_fg_dic['ml_dls'][fg_nu_ind][fg_lab_ind][0]
				except:
					fg_Dls = self.george_fg_dic['ml_dls'][fg_nu_ind][fg_lab_ind][0]

				#print fg, fg_lab_ind
				#print self.george_fg_dic['ml_dls'][fg_nu_ind][fg_lab_ind][0]

			#plot(fg_els, fg_Dls, label='all');legend(loc = 4);xlim(1,12e3);ylim(1e-4,1e4);show();quit()
		
		
		#print which_fg;quit()
		fg_Cls = ( fg_Dls * 2 * np.pi ) / ( fg_els * (fg_els + 1) ) * 1e-12
		

		if just_Cls:
			return [fg_els,fg_Cls]

		#alpha = -1.2;amp = np.interp(0.,fg_els, fg_Cls);plaw = amp * fg_els ** alpha;loglog(fg_els, fg_Cls, 'r.');loglog(fg_els, plaw, 'k');show();quit()
		#loglog(fg_els, fg_Cls_pol, 'r')

		#Cls2map
		if also_P:
			fg_Cls_pol = fg_Cls * pol_frac_per_cent/100.
			fg_Cls = np.asarray( [fg_els, fg_Cls, fg_Cls_pol, fg_Cls_pol] ).T
		else:
			fg_Cls = np.asarray( [fg_els, fg_Cls] ).T
		
		self.fg_Cls = fg_Cls
		fg_CLS = self.Cls2CLS(fg_Cls,mapparams)
		FOREGROUND_MAP_FFT = np.sqrt(fg_CLS) * self.scalefac
		
		if random_seed == None:
			random_seed = self.inidic['random_seed_val_for_cov']

		'''
		if noofsims>1:
			np.random.seed(random_seed) #for covariance calculation
			print 'random seed step in fn_get_foreground_power... aborting';sys.exit()
		'''

		if perform_lensing:
			#pull initialized variables here
			theta_x, theta_y = self.theta_x, self.theta_y
			to_evaluate_unlensed_theta_x, to_evaluate_unlensed_theta_y = self.to_evaluate_unlensed_theta_x, self.to_evaluate_unlensed_theta_y		

		if self.is_seq(reqdbox):
			ex1,ex2,ey1,ey2 = reqdbox
			SIMMAPS_FG = np.zeros( (noofsims, len(FOREGROUND_MAP_FFT), ex2-ex1, ey2-ey1) )
		else:			
			SIMMAPS_FG = np.zeros( (noofsims, len(FOREGROUND_MAP_FFT), ny, nx) )

		#percentage of FG - residual
		try:
			residual_fg = self.run_s2_dict['residual_fg']
		except:
			residual_fg = 1.

		if which_fg == 'kSZ': residual_fg = 1. #foreground cleaning not possible for kSZ

		for cnt in range(noofsims):

                        if cnt % 5000 == 0 and cnt>0:
                                logline = '\t\t\tGetting COV for FG (which_pop = %s) simno: %s' %(which_fg, cnt)
                                logfile = open(self.log_file,'a')
                                logfile.writelines('%s\n' %(logline))
                                logfile.close()
                                #print logline
                        
                       	#print np.random.randn(nx)[0], nu1
			np.random.seed(self.cmbrandomseedval)
			FG_DUMMY = np.fft.fft2( np.random.randn(ny,nx) )

			TMP = np.fft.ifft2( FOREGROUND_MAP_FFT * FG_DUMMY ).real
			TMP *= residual_fg

			for tqucnt in range(len(TMP)):
				TMP[tqucnt] -= np.mean(TMP[tqucnt])

			if perform_lensing: #then lens the foregrounds as well
				for tqucnt in range(len(TMP)):
					TMP[tqucnt] = intrp.RectBivariateSpline( theta_y[:,0], theta_x[0,:], TMP[tqucnt], kx = 5, ky = 5).ev(to_evaluate_unlensed_theta_y, to_evaluate_unlensed_theta_x).reshape([ny,nx])

			#subplot(111);imshow(TMP[0]);title('FG');colorbar();show()#;quit()
			#self.FG_radio_tmp = TMP[0]

			if noofsims == 1: SIMMAPS_FG[cnt,:] = TMP; continue

			if nu1 == 150:
				TMP = np.fft.ifft2( self.Bl[0] * np.fft.fft2(TMP) ).real #assuming T and P beams are the same
			elif nu1 == 90:
				TMP = np.fft.ifft2( self.Bl_90_GHz[0] * np.fft.fft2(TMP) ).real #assuming T and P beams are the same

			if self.inidic['degrade_resol_fac']>1:
				print "r";quit()
				TMP_DS = np.zeros( (len(TMP), self.simmapparams[0], self.simmapparams[1]) )
				for tqucnt in range(len(TMP)):
					TMP_DS[tqucnt] = self.downsample_map(TMP[tqucnt],self.inidic['degrade_resol_fac'])
				TMP = TMP_DS

			TMP = np.fft.ifft2( self.TWODTF[0] * np.fft.fft2(TMP) ).real
			if self.is_seq(reqdbox): #if generating sims to get COV then use Beam / TF as well
				ex1,ex2,ey1,ey2 = reqdbox
				SIMMAPS_FG[cnt,:] = TMP[:,ex1:ex2,ey1:ey2]
			else:

				try:
					SIMMAPS_FG[cnt] = TMP
				except:
					SIMMAPS_FG = np.zeros( (noofsims, len(FOREGROUND_MAP_FFT), nx/self.inidic['degrade_resol_fac'], ny/self.inidic['degrade_resol_fac']) ) #just T for now
					SIMMAPS_FG[cnt] = TMP

			''' #20170215: changed to include beam here even for sims (for tSZ subtraction stuff)
			if self.is_seq(reqdbox): #if generating sims to get COV then use Beam / TF as well
				ex1,ex2,ey1,ey2 = reqdbox
				#subplot(121);imshow(TMP[0]);colorbar()

				TMP = np.fft.ifft2( self.Bl[0] * np.fft.fft2(TMP) ).real #assuming T and P beams are the same

				if self.inidic['degrade_resol_fac']>1:
					TMP_DS = np.zeros( (len(TMP), self.simmapparams[0], self.simmapparams[1]) )
					for tqucnt in range(len(TMP)):
						TMP_DS[tqucnt] = self.downsample_map(TMP[tqucnt],self.inidic['degrade_resol_fac'])
					TMP = TMP_DS

				TMP = np.fft.ifft2( self.TWODTF[0] * np.fft.fft2(TMP) ).real

				#subplot(122);imshow(TMP[0]);colorbar();show();quit()
				SIMMAPS_FG[cnt,:] = TMP[:,ex1:ex2,ey1:ey2]
			else:
				SIMMAPS_FG[cnt,:] = TMP
			'''

		if not return_cov: #return maps
			return SIMMAPS_FG
		else:
			COV_FG = {}
			FG_MAPFORCOV = SIMMAPS_FG[:,0] * 1e6
			npixels = FG_MAPFORCOV.shape[1] * FG_MAPFORCOV.shape[2]
			COV_FG['T'] = self.calcCov(FG_MAPFORCOV, noofsims, npixels)
			#imshow(COV_FG['T']);colorbar();show();quit()
			return COV_FG

		##################################################################
		##################################################################


	def fn_get_foreground_DFSG_power(self, mapparams, noofsims, els = None, reqdbox = None,random_seed = None, return_cov = 0, perform_lensing = 0):

		nx, ny, dx, dy = mapparams
		dx *= self.arcmins2radians
		dy *= self.arcmins2radians
		npixels = nx * ny
		self.scalefac = np.sqrt(1./ (dx * dy))

		if els == None:
			try:
				els = self.els
			except:
				els = np.arange(2,3e5)

		normel = 3000.
		for forecnt in range(2):
			if forecnt == 0:
				Dls = np.zeros(len(els)) + (self.poisson_power * 1e-12 * (els/normel)**2.)
				#Cls = Dls / (els)**2.
			else:
				D0 = self.clustered_power * 1e-12
				clus_el_cutoff = 1500 #Baxter et al.
				Dls = D0 * (els/normel)**self.clustered_comp_gamma
				#Dls[els<=clus_el_cutoff] = D0
				#Cls = Dls# / (els)**self.clustered_comp_gamma

			Cls = ( Dls * 2 * np.pi ) / ( els * (els + 1) )

			#Cls2map
			CLS = self.Cls2CLS(Cls,mapparams)
			CLS = np.sqrt(CLS) * self.scalefac

			#semilogy(els,Dls*1e12);vlines(3000,min(Dls*1e12),max(Dls*1e12));xlim(1e3,1e4)

			try:
				FOREGROUND_MAP_FFT += CLS
			except:
				FOREGROUND_MAP_FFT = CLS

		#show();quit()

		if random_seed == None:
			random_seed = self.inidic['random_seed_val_for_cov']

		if noofsims>1:
			np.random.seed(random_seed)

		if perform_lensing:
			#pull initialized variables here
			theta_x, theta_y = self.theta_x, self.theta_y
			to_evaluate_unlensed_theta_x, to_evaluate_unlensed_theta_y = self.to_evaluate_unlensed_theta_x, self.to_evaluate_unlensed_theta_y		

		#TWODTF = self.fn_get_HPF(mapparams, minel = self.inidic['min_el'], maxel = self.inidic['max_el'])[0]
		if self.is_seq(reqdbox):
			ex1,ex2,ey1,ey2 = reqdbox
			SIMMAPS_FG = np.zeros( (noofsims, 1, ex2-ex1, ey2-ey1) ) #just T for now
		else:			
			SIMMAPS_FG = np.zeros( (noofsims, 1, ny, nx) ) #just T for now

		for cnt in range(noofsims):

                        if cnt % 5000 == 0 and cnt>0:
                                logline = '\t\t\t\tGetting COV for DSFG. simno: %s' %(cnt)
                                logfile = open(self.log_file,'a')
                                logfile.writelines('%s\n' %(logline))
                                logfile.close()
                                print logline


			DUMMY = np.fft.fft2( np.random.randn(nx,ny) )
			#TMP = np.fft.ifft2( FOREGROUND_MAP_FFT * DUMMY * TWODTF ).real
			TMP = np.fft.ifft2( FOREGROUND_MAP_FFT * DUMMY ).real
			TMP = TMP - np.mean(TMP)

			#clf();subplot(111);imshow(TMP[0]);title('FG');colorbar();show();quit()

			if perform_lensing: #then lens the foregrounds as well
				TMP[0] = intrp.RectBivariateSpline( theta_y[:,0], theta_x[0,:], TMP[0], kx = 5, ky = 5).ev(to_evaluate_unlensed_theta_y, to_evaluate_unlensed_theta_x).reshape([ny,nx])

			#subplot(122);imshow(TMP[0]);title('FG lensed');colorbar();show();quit()

			if self.is_seq(reqdbox): #if generating sims to get COV then use Beam / TF as well
				ex1,ex2,ey1,ey2 = reqdbox
				#subplot(121);imshow(TMP[0]);colorbar()
				TMP = np.fft.ifft2( self.TWODTF * self.Bl * np.fft.fft2(TMP) ).real
				#subplot(122);imshow(TMP[0]);colorbar();show();quit()
				SIMMAPS_FG[cnt, 0] = TMP[0, ex1:ex2,ey1:ey2]
			else:
				SIMMAPS_FG[cnt, 0] = TMP

		if not return_cov: #return maps
			return SIMMAPS_FG
		else:
			FG_MAPFORCOV = SIMMAPS_FG[:,0] * 1e6
			npixels = FG_MAPFORCOV.shape[1] * FG_MAPFORCOV.shape[2]
			COV_FG = self.calcCov(FG_MAPFORCOV, noofsims, npixels)
			#imshow(COV_FG);colorbar();show();quit()
			return COV_FG

		##################################################################
		##################################################################

	def fn_add_ILC_residual(self, ilc_residual_filename, mapparams, noofsims, reqdbox = None, random_seed = None, return_cov = 0, perform_lensing = 0):

		nx, ny, dx, dy = mapparams
		dx *= self.arcmins2radians
		dy *= self.arcmins2radians
		npixels = nx * ny
		self.scalefac = np.sqrt(1./ (dx * dy))

		if noofsims == 1:
		        logline = '\t\t\tAdding ILC residual to maps now'
		        logfile = open(self.log_file,'a');logfile.writelines('%s\n' %(logline));logfile.close()

		#Temp.
		ilc_els, ilc_Cls_res = pickle.load(gzip.open(ilc_residual_filename,'rb'))['Cl_residual']
		self.ilc_els = ilc_els
		self.ilc_Cls_res = ilc_Cls_res


		#deconvolve beam now
		#Bl = pickle.load(gzip.open(ilc_residual_filename, 'rb'))['Bl_effective']
		#ilc_Cls_res = ilc_Cls_res/Bl**2.

		#Pol.
		##ilc_residual_filename_pol = ilc_residual_filename.replace('.pkl.gz','_pol.pkl.gz')
		ilc_residual_filename_pol = ilc_residual_filename.replace('with_1_f_noise.pkl.gz','_polwith_1_f_noise.pkl.gz')
		ilc_els, ilc_Cls_res_pol = pickle.load(gzip.open(ilc_residual_filename_pol, 'rb'))['Cl_residual']
		self.ilc_Cls_res_pol = ilc_Cls_res_pol

		#deconvolve beam now
		#Bl_pol = pickle.load(gzip.open(ilc_residual_filename_pol, 'rb'))['Bl_effective']
		#ilc_Cls_res_pol = ilc_Cls_res_pol/Bl_pol**2.

		ilc_Cls_res *= 1e-12
		ilc_Cls_res_pol *= 1e-12

		Dls_fac = (ilc_els) * (ilc_els+1) / 2/ np.pi
		#loglog(ilc_els, Dls_fac * ilc_Cls_res * 1e12,'k');loglog(ilc_els, Dls_fac * ilc_Cls_res_pol * 1e12,'g');show();quit()

		ilc_Cls = np.asarray( [ilc_els, ilc_Cls_res, ilc_Cls_res_pol, ilc_Cls_res_pol] ).T

		ilc_CLS = self.Cls2CLS(ilc_Cls,mapparams)
		ilc_CLS = np.sqrt(ilc_CLS)
		for ccc in range(len(ilc_CLS)):
			ilc_CLS[ccc][np.isnan(ilc_CLS[ccc])] = 0.

		ILC_MAP_FFT = ilc_CLS * self.scalefac

		'''
		if random_seed == None:
			random_seed = self.inidic['random_seed_val_for_cov']

		if noofsims>1:
			np.random.seed(random_seed) #for covariance calculation
			print 'random seed step in fn_get_ILC_power... aborting';sys.exit()
		'''

		if perform_lensing:
			#pull initialized variables here
			theta_x, theta_y = self.theta_x, self.theta_y
			to_evaluate_unlensed_theta_x, to_evaluate_unlensed_theta_y = self.to_evaluate_unlensed_theta_x, self.to_evaluate_unlensed_theta_y		

		if self.is_seq(reqdbox):
			ex1,ex2,ey1,ey2 = reqdbox
			SIMMAPS_ILC = np.zeros( (noofsims, len(ILC_MAP_FFT), ex2-ex1, ey2-ey1) )
		else:			
			SIMMAPS_ILC = np.zeros( (noofsims, len(ILC_MAP_FFT), ny, nx) )

		for cnt in range(noofsims):

			if cnt % 5000 == 0 and cnt>0:
				logline = '\t\t\tGetting COV for ILC residuals; simno: %s' %(cnt)
				logfile = open(self.log_file,'a');logfile.writelines('%s\n' %(logline));logfile.close()
				#print logline

			ILC_DUMMY = np.asarray( [np.fft.fft2( np.random.randn(nx,ny) ) for tqu in range(self.tqulen)] )
			TMP = np.fft.ifft2( ILC_MAP_FFT * ILC_DUMMY ).real

			for tqucnt in range(len(TMP)):
				TMP[tqucnt] -= np.mean(TMP[tqucnt])

			if perform_lensing: #then lens the foregrounds as well
				for tqucnt in range(len(TMP)):
					TMP[tqucnt] = intrp.RectBivariateSpline( theta_y[:,0], theta_x[0,:], TMP[tqucnt], kx = 5, ky = 5).ev(to_evaluate_unlensed_theta_y, to_evaluate_unlensed_theta_x).reshape([ny,nx])

			if noofsims == 1: SIMMAPS_ILC[cnt,:] = TMP; continue

			#for tqucnt in range(self.tqulen):
			#	subplot(2,3,tqucnt+1);imshow(TMP[tqucnt]);colorbar()

			TMP = np.fft.ifft2( self.Bl[0] * np.fft.fft2(TMP) ).real #assuming T and P beams are the same

			#for tqucnt in range(self.tqulen):
			#	subplot(2,3,self.tqulen + (tqucnt+1));imshow(TMP[tqucnt]);colorbar()
			#show();quit()

			if self.inidic['degrade_resol_fac']>1:
				TMP_DS = np.zeros( (len(TMP), self.simmapparams[0], self.simmapparams[1]) )
				for tqucnt in range(len(TMP)):
					TMP_DS[tqucnt] = self.downsample_map(TMP[tqucnt],self.inidic['degrade_resol_fac'])
				TMP = TMP_DS

			TMP = np.fft.ifft2( self.TWODTF[0] * np.fft.fft2(TMP) ).real
			if self.is_seq(reqdbox): #if generating sims to get COV then use Beam / TF as well
				ex1,ex2,ey1,ey2 = reqdbox
				SIMMAPS_ILC[cnt,:] = TMP[:,ex1:ex2,ey1:ey2]
			else:
				try:
					SIMMAPS_ILC[cnt] = TMP
				except:
					SIMMAPS_ILC = np.zeros( (noofsims, len(ILC_MAP_FFT), nx/self.inidic['degrade_resol_fac'], ny/self.inidic['degrade_resol_fac']) ) 
					SIMMAPS_ILC[cnt] = TMP

		if not return_cov: #return maps
			return SIMMAPS_ILC
		else:

			locdic = self.locdic
			identifiers = self.inidic['cov_identifiers']

			COV_ILC_RES = {}
			ILC_MAPFORCOV = SIMMAPS_ILC * 1e6
			npixels = ILC_MAPFORCOV.shape[2] * ILC_MAPFORCOV.shape[3]
			for iden in identifiers:
				if iden not in self.reqd_idens: continue
				if len(iden) == 1:
					loc = locdic[iden]
					MAPFORCOV = ILC_MAPFORCOV[:,loc]
					COV = self.calcCov(MAPFORCOV, noofsims, npixels)
				elif len(iden) == 2:
					loc1,loc2 = locdic[iden[0]],locdic[iden[1]]
					MAPFORCOV1, MAPFORCOV2 = ILC_MAPFORCOV[:,loc1], ILC_MAPFORCOV[:,loc2]
					MAPFORCOV = np.concatenate( (MAPFORCOV1,MAPFORCOV2), axis = 1 )
					COV = np.asarray(self.calcCov(MAPFORCOV, noofsims, 2*npixels))
				elif len(iden) == 3:
					loc1,loc2,loc3 = locdic[iden[0]],locdic[iden[1]],locdic[iden[2]]
					MAPFORCOV1, MAPFORCOV2, MAPFORCOV3 = ILC_MAPFORCOV[:,loc1], ILC_MAPFORCOV[:,loc2], ILC_MAPFORCOV[:,loc3]
					MAPFORCOV = np.concatenate( (MAPFORCOV1,MAPFORCOV2,MAPFORCOV3), axis = 1 )
					COV = np.asarray(self.calcCov(MAPFORCOV, noofsims, 3*npixels))


				COV_ILC_RES[iden] = COV
			#COV_FG['T'] = self.calcCov(ILC_MAPFORCOV, noofsims, npixels)
			#subplot(121);imshow(COV_ILC_RES['T']);colorbar()
			#subplot(122);imshow(COV_ILC_RES['QU']);colorbar();show();quit()
			return COV_ILC_RES
		##################################################################
		##################################################################

	def fn_plot_pow_spec(self,mapparams,MAP, binsize = None, also_plot = 0, show_theory_curves = 0, K_to_uK = 1, cross_power = 0):


		if K_to_uK: 
			K_to_uK_power = 1e12
		else:
			K_to_uK_power = 1e-12

                nx, ny, dx, dy = mapparams
                dx_rad = dx * self.arcmins2radians

                lx, ly = self.get_lxly(mapparams)

		if len(MAP) == 3:
			Q_FFT, U_FFT = np.fft.fft2(MAP[1]),np.fft.fft2(MAP[2])
			angle = self.fn_get_lxly_az_angle(lx,ly)
			E = np.fft.ifft2( np.cos(angle) * Q_FFT + np.sin(angle) * U_FFT ).real #lensed E
			B = np.fft.ifft2( -np.sin(angle) * Q_FFT + np.cos(angle) * U_FFT ).real #lensed B

			MAP = np.asarray( [MAP[0],E,B,MAP[1],MAP[2]] )

		if binsize == None:
			binsize = lx.ravel()[1] -lx.ravel()[0]

		if also_plot:
			clf();
			subplot(111, xscale ='log', yscale ='log')
			colorarr = ['k', 'r', 'g', 'b', 'm']
			lab = ['TT','EE','BB','QQ','UU']
		
		if not cross_power:
			all_cls = []
			for tqucnt in range(len(MAP)):

				MAP_PSD = abs( np.fft.fft2(MAP[tqucnt]) * dx_rad)** 2 / (nx * ny)
				cls = self.fn_radial_profile(MAP_PSD, (lx,ly), bin_size = binsize, minbin = 100, maxbin = 10000, to_arcmins = 0)
				all_cls.append([cls])

				#errorbar(cls[:,0], cls[:,1]*cls[:,0]*(cls[:,0]+1)/2/np.pi, yerr = [cls[:,2], cls[:,2]], marker = '.', color = colorarr[tqucnt], label = lab[tqucnt])
				#if tqucnt == 1 or tqucnt==2:
				#	cls[:,1] *= 2

				if also_plot:
					if show_theory_curves: Dls_fac = cls[:,0]*(cls[:,0]+1)/2/np.pi
					else: Dls_fac = 1
					plot(cls[:,0], cls[:,1]* Dls_fac* K_to_uK_power, marker = '.', color = colorarr[tqucnt], label = lab[tqucnt])
					if tqucnt <3 and show_theory_curves:
						plot(self.Cls[:,0], self.Cls[:,tqucnt+1]*self.Cls[:,0]*(self.Cls[:,0]+1)/2/np.pi * K_to_uK_power, color = colorarr[tqucnt], lw = 0.2)

			return all_cls

		else:
			#print MAP[0], MAP[1]
			#MAP_PSD = np.fft.fft2(MAP[0]) * dx_rad * np.fft.fft2(MAP[1]) * dx_rad / (nx * ny)
			MAP_PSD = np.fft.fft2(MAP[0]) * dx_rad * np.conj( np.fft.fft2(MAP[1]) ) * dx_rad / (nx * ny)
			cls = self.fn_radial_profile(MAP_PSD, (lx,ly), bin_size = binsize, minbin = 100, maxbin = 10000, to_arcmins = 0)

			return [cls]

		if also_plot and show_theory_curves and len(MAP)>1: #TE
			MAP_PSD = np.fft.fft2(MAP[0]) * dx_rad * np.conj( np.fft.fft2(MAP[1]) ) * dx_rad / (nx * ny)
			cls = self.fn_radial_profile(MAP_PSD, (lx,ly), bin_size = binsize, minbin = 100, maxbin = 10000, to_arcmins = 0)
			plot(cls[:,0], cls[:,1]*Dls_fac * K_to_uK_power, marker = '.', color = 'c', label = 'TE', ls = 'None')
			plot(self.Cls[:,0], self.Cls[:,4]*self.Cls[:,0]*(self.Cls[:,0]+1)/2/np.pi * K_to_uK_power, color = 'c', lw = 0.5)
		
		if also_plot:
			#ylim(1e-4,1e3)
			#ylim(1e-5,1.)
			xlim(100,2e4)
			if show_theory_curves: ylim(1e-2,1e4)
			legend(loc=1,fancybox=1)
			show()#;quit()
		else:
			return all_cls


        def fn_perform_lensing_nfw(self, mapparams, noofsims, param_dict, use_beam, degrade_resol_fac = 1, no_lensing  = 0, random_seed = None, passed_CMBMAP = None, Bl = None, nu = 150, add_FG = 0, pad = 1):


                from scipy import interpolate as intrp
		
		#get lx,ly
                nx, ny, dx, dy = mapparams
                dx *= self.arcmins2radians
                dy *= self.arcmins2radians
                lx, ly = self.get_lxly(mapparams)

		#things for estimators
		T_in_cov_identifiers = 0
		for idens in self.inidic['cov_identifiers']:
			if idens.find('T')>-1:
				T_in_cov_identifiers = 1
				break

		for idens in self.inidic['cov_identifiers']:
			if idens.find('E')>-1 or idens.find('B')>-1:
				self.tqulen = 5
				self.perform_EB = 1
				break

		## check the degrade_resol_fac factor 
		#new map shape would be
                nx, ny, dx, dy = mapparams
		mapparams_deg = [nx/degrade_resol_fac, ny/degrade_resol_fac, dx * degrade_resol_fac, dy * degrade_resol_fac]
		nx_deg,ny_deg,dx_deg,dy_deg = mapparams_deg
                tqulen = self.tqulen

                if not self.quiet:
	                logline = '\t\tlensing might be performed now - interpolation actually'
	                logfile = open(self.log_file,'a')
	                logfile.writelines('%s\n' %(logline))
	                logfile.close()
	                print logline		

		#starting CMB sims + lensing now
		if random_seed == None:
			random_seed = self.inidic['random_seed_val_for_cov']
		
		#from IPython import embed;embed()
		#from pylab import *
		np.random.seed(random_seed)
		self.cmb_noise_randomseedval = random_seed

		#print random_seed

		#some variables for lensing
		if self.perform_EB:
			end = tqulen - 2
		else:
			end = tqulen

		if T_in_cov_identifiers:
			start = 0
		else:
			start = 1

                start_time = time.time()

		#pull initialized variables here
		if not no_lensing:
			theta_x, theta_y = self.theta_x, self.theta_y
			to_evaluate_unlensed_theta_x, to_evaluate_unlensed_theta_y = self.to_evaluate_unlensed_theta_x, self.to_evaluate_unlensed_theta_y


		#beam stuff
		if not self.is_seq(Bl):
			Bl = self.Bl

                SIMMAPS_L = np.zeros((noofsims,tqulen,ny_deg,nx_deg))
                for cnt in range(noofsims):
                        if cnt % 5000 == 0 and cnt>0:
                                logline = '\t\t\t\tsimno: %s; time: %s' %(cnt, time.time() - start_time)
                                logfile = open(self.log_file,'a')
                                logfile.writelines('%s\n' %(logline))
                                logfile.close()
                                print logline

			if self.is_seq(passed_CMBMAP):
				CMB = np.copy(passed_CMBMAP)
				CMB_ori = np.copy(CMB)

			else:
				##################################################################
				#make CMB sims
				CAMBMAP_FFT = self.CAMBMAP_FFT
				#CMB = np.fft.ifft2( np.copy( CAMBMAP_FFT ) * np.fft.fft2( np.random.randn(nx,ny) ) ).real

				#20170403 - to include TE correaltion and different randomralisations for each component (T,E,B)
				#now we make TEB maps - convert them to TQU immediately for lensing and similar procedure as before
				#print np.random.randn(nx)[0], 'hi'
				##from IPython import embed; embed()
				if tqulen>1:
					GAUSS_REALS = np.asarray( [np.fft.fft2( np.random.randn(ny,nx) ) for teb in range(len(CAMBMAP_FFT) - 1)] )
				else:
					GAUSS_REALS = np.asarray( [np.fft.fft2( np.random.randn(ny,nx) ) for teb in range(len(CAMBMAP_FFT))] )
				CMB = np.zeros( (len(GAUSS_REALS), ny, nx) )

				for teb in range( len(GAUSS_REALS) ): #T,E,B
					if teb == 1: #for E include correlation between T
						t1 = GAUSS_REALS[0] * self.CLS[3] / self.CLS[0]**0.5
						t2 = GAUSS_REALS[1] * ( self.CLS[1] - (self.CLS[3]**2. /self.CLS[0]) )**0.5
						Emap_fft = (t1 + t2) * self.scalefac
						Emap_fft[np.isnan(Emap_fft)] = 0.
						CMB[teb] = np.fft.ifft2( Emap_fft ).real

					else:
						CMB[teb] = np.fft.ifft2( np.copy( CAMBMAP_FFT[teb] ) * GAUSS_REALS[teb] ).real

					CMB[teb] = CMB[teb] - np.mean(CMB[teb])

				if tqulen>1:
					E_FFT, B_FFT = np.fft.fft2(CMB[1]),np.fft.fft2(CMB[2])
					angle = self.fn_get_lxly_az_angle(lx,ly)
					CMB[1] = np.fft.ifft2( np.cos(angle) * E_FFT - np.sin(angle) * B_FFT ).real #Q
					CMB[2] = np.fft.ifft2( np.sin(angle) * E_FFT + np.cos(angle) * B_FFT ).real #U
					self.cambmapvars = ['T','Q','U','E','B']
					#20170403 - to include TE correaltion and different randomralisations for each component (T,E,B)
					#now we make TEB maps - convert them to TQU immediately for lensing and similar procedure as before

				"""
				make sure power specturm looks nice
				"""
				if 0==1:self.fn_plot_pow_spec(mapparams, CMB, also_plot = 1);quit()

				##################################################################
				if 0==1:#add_FG:

					#add foregrounds now to "T"
		                        logline = '\t\t\t\tadding Poisson / clustered foreground powers now'
		                        logfile = open(self.log_file,'a')
		                        logfile.writelines('%s\n' %(logline))
		                        logfile.close()
		                        print logline
					#all foregrounds will be lensed if FG are added here
					FG = self.fn_get_foreground_DFSG_power(mapparams, 1, perform_lensing = 0)
					CMB[0] = CMB[0] + FG[0, 0]
				##################################################################

				#mean removal
				for tqucnt in range(len(CMB)):
					CMB[tqucnt] = CMB[tqucnt] - np.mean(CMB[tqucnt])

				'''
				if pad:
					px, py = nx/10, ny/10
					CMB[:,0:px,:] = CMB[:,:,0:py] = CMB[:,-px:,:] = CMB[:,:,-py:] = 0.
					#imshow(CMB[0]);colorbar();show();quit()
				'''

				CMB_ori = np.copy(CMB)
				#clf();imshow(CMB[0]);colorbar();title('Mean = %s' %(np.mean(CMB[0])));show()#;quit()
				##################################################################
			
			##################################################################
			#Lensing T,Q,U
			#imshow(CMB[0]);colorbar();show();quit()
			self.CMB_unlensed = np.copy(CMB)
			for tqucnt in range(start,end): #for T, Q, U
				if not no_lensing:
					CMB[tqucnt] = intrp.RectBivariateSpline( theta_y[:,0], theta_x[0,:], CMB[tqucnt], kx = 5, ky = 5).ev(to_evaluate_unlensed_theta_y, to_evaluate_unlensed_theta_x).reshape([ny,nx])

				else: #no lensing case
					pass
			self.CMB_lensed = np.copy(CMB)
			###clf();imshow(CMB_ori[tqucnt] - CMB[tqucnt]);colorbar();show();sys.exit()

			##################################################################
			##################################################################
			"""
			adds foreground to T maps after lensing: copied from 2016_08/cluster_lensing pipeline
			"""
			#this is the diffused component
			try:
				### fg_T = self.run_s2_dict['fg_T']
				fg_T = self.inidic['fit_FG']
			except:
				fg_T = None

				####20170324 - if foreground cleaning is being performed, then add 100 per cent kSZ power spectra
				try:
					foreground_cleaning = self.run_s2_dict['foreground_cleaning']
				except:
					foreground_cleaning = 0

				if foreground_cleaning:
					fg_T = 'kSZ'

			#print fg_T
			tSZ_emission = None #temporary - fix this properly later - 20170127

			#FG_dusty_gal = self.fn_pick_DGpop_sim(mapparams, nu = 150)
			#quit()

			if fg_T <> None:#0==1:#add_FG:
				#this is the Poisson component for dusty/radio galaxies
		                #logline = '\t\t\t#%s: adding %s Gaussian / Poisson sims now to CMB T' %(cnt, fg_T);logfile = open(self.log_file,'a');logfile.writelines('%s\n' %(logline));logfile.close()
		                #print logline
						

				self.FG_added = 1
				if fg_T == 'all':
					'''
					SIM_POISSON_RADIO = self.fn_sim_radio_dusty_Poisson_pop(mapparams, noofsims = 1, nu = nu, which_pop = 'RG')
					SIM_POISSON_DUSTY = self.fn_sim_radio_dusty_Poisson_pop(mapparams, noofsims = 1, nu = nu, which_pop = 'DG-Po')
					'''
					FG_RG = self.fn_sim_radio_dusty_Poisson_pop_pick_from_sims(mapparams, noofsims = 1, nu = nu, which_pop = 'RG')
					if self.add_tSZ_emission:
						FG_DG, tSZ_emission = self.fn_sim_radio_dusty_Poisson_pop_pick_from_sims(mapparams, noofsims = 1, nu = nu, which_pop = 'DG-Po')
					else:
						FG_DG = self.fn_sim_radio_dusty_Poisson_pop_pick_from_sims(mapparams, noofsims = 1, nu = nu, which_pop = 'DG-Po')
					###FG_tSZ_kSZ_uncorr = self.fn_get_foreground_power(fg_T, mapparams, 1, perform_lensing = 0, nu1=nu, nu2=nu)[0]

					FG_tSZ_uncorr = self.fn_get_foreground_power('tSZ', mapparams, 1, perform_lensing = 0, nu1=nu, nu2=nu)[0]
					FG_kSZ_uncorr = self.fn_get_foreground_power('kSZ', mapparams, 1, perform_lensing = 0, nu1=nu, nu2=nu)[0]

					FG_tSZ_kSZ_uncorr = FG_tSZ_uncorr +  FG_kSZ_uncorr

					if len(FG_tSZ_kSZ_uncorr) < len(CMB):
						FG_tSZ_kSZ_uncorr_tmp = np.zeros ( (len(CMB), ny, nx) )
						FG_tSZ_kSZ_uncorr_tmp[0] = FG_tSZ_kSZ_uncorr
						FG_tSZ_kSZ_uncorr = FG_tSZ_kSZ_uncorr_tmp


					#clf();imshow(FG_RG[0]);colorbar();show();sys.exit(0)

					#downsample FG if necessary
					if FG_RG[0].shape[0] <> nx:
						ds_fac = FG_RG[0].shape[0] / nx
						TMP_RG = np.zeros( (len(FG_RG), nx , ny) )
						TMP_DG = np.zeros( (len(FG_RG), nx , ny) )
						for tqucnt in range(len(FG_RG)):
							TMP_RG[tqucnt] = self.downsample_map(FG_RG[tqucnt],ds_fac)
							TMP_DG[tqucnt] = self.downsample_map(FG_DG[tqucnt],ds_fac)
						FG_RG = np.copy(TMP_RG)
						FG_DG = np.copy(TMP_DG)
							
					FG = FG_tSZ_kSZ_uncorr + FG_RG + FG_DG

					#subplot(221);imshow(FG_DG[0]);colorbar();title('FG_DG');subplot(222);imshow(FG_RG[0]);colorbar();title('FG_RG')
					#subplot(223);imshow(FG_tSZ_kSZ_uncorr[0]);colorbar();title('tSZ + kSZ');show()

					
				elif fg_T == 'RG' or fg_T == 'DG-Po' or fg_T == 'DG-Cl': #Poisson realisation
					#FG = self.fn_sim_radio_dusty_Poisson_pop(mapparams, noofsims = 1, nu = nu, which_pop = fg_T, nu = nu)
					FG = self.fn_sim_radio_dusty_Poisson_pop_pick_from_sims(mapparams, noofsims = 1, nu = nu, which_pop = fg_T)
					if self.add_tSZ_emission and fg_T == 'DG-Po':
						FG, tSZ_emission = self.fn_sim_radio_dusty_Poisson_pop_pick_from_sims(mapparams, noofsims = 1, nu = nu, which_pop = 'DG-Po')
					else:
						FG = self.fn_sim_radio_dusty_Poisson_pop_pick_from_sims(mapparams, noofsims = 1, nu = nu, which_pop = 'DG-Po')
						tSZ_emission = np.zeros(FG[0].shape)

				elif fg_T == 'tSZ' or fg_T == 'kSZ': #Gaussian realisation of uncorrelated tSZ and kSZ
					FG = self.fn_get_foreground_power(fg_T, mapparams, 1, perform_lensing = 0, nu1=nu, nu2=nu)[0]
					if len(FG) < len(CMB):
						FG_tmp = np.zeros ( (len(CMB), ny, nx) )
						FG_tmp[0] = FG
						FG = FG_tmp

				#20170413
				elif fg_T == 'all_Gaussian': #Gaussian realisation of all foregrounds from George et al. 2015 data
					'''
					logline = '\t\t\t\tadding all_Gaussian foreground now to nu = %s' %(nu)
					logfile = open(self.log_file,'a');logfile.writelines('%s\n' %(logline));logfile.close()
					print logline
					'''

					#FG = self.fn_get_foreground_power(fg_T, mapparams, 1, perform_lensing = 0, nu1=nu, nu2=nu, ign_DG_RG = 0)[0]
					#20170727 - use random seed to produce same 90/150GHz foregrounds during tSZ sub. process
					### FG = self.fn_get_foreground_power(fg_T, mapparams, 1, perform_lensing = 0, nu1=nu, nu2=nu, ign_DG_RG = 0, random_seed = self.present_time_rs, set_rs = 1)[0]
					if not nu == 'tszfree':
						FG = self.fn_get_foreground_power(fg_T, mapparams, 1, perform_lensing = 0, nu1=nu, nu2=nu, ign_DG_RG = 0)[0]
						
					if nu == 'tszfree':
						clusterstuff = cluster_stuff()
						ysz_Tsz_conv_fac_150 = clusterstuff.compton_y_to_delta_Tcmb(150e9) #conversion factor for 150 GHz
						ysz_Tsz_conv_fac_90 = clusterstuff.compton_y_to_delta_Tcmb(90e9) #conversion factor for 90 GHz

						factor = ysz_Tsz_conv_fac_90/ysz_Tsz_conv_fac_150
						factor = factor **2
						
						FG_90 = self.fn_get_foreground_power(fg_T, mapparams, 1, perform_lensing = 0, nu1=90, nu2=90, ign_DG_RG = 0, just_Cls =1)[1]
						
						FG_150= self.fn_get_foreground_power(fg_T, mapparams, 1, perform_lensing = 0, nu1=150, nu2=150, ign_DG_RG = 0, just_Cls =1)[1]
						fg_Cls= (FG_150 * factor- FG_90)/(factor - 1)

						fg_CLS = self.Cls2CLS(fg_Cls,mapparams)
						FOREGROUND_MAP_FFT = np.sqrt(fg_CLS) * self.scalefac
						FG_DUMMY = np.fft.fft2( np.random.randn(ny,nx) )
						TMP = np.fft.ifft2( FOREGROUND_MAP_FFT * FG_DUMMY ).real
						
						for tqucnt in range(len(TMP)):
							TMP[tqucnt] -= np.mean(TMP[tqucnt])
						FG = TMP
						#from IPython import embed;embed()
						#sys.exit()
					'''
					if nu == 150:
						#print 'hi', nu
						self.FG_150 = FG
					elif nu == 90:
						#print 'hi', nu
						self.FG_90 = FG
						#from IPython import embed; embed()
						subplot(131);imshow(self.FG_150[0]*1e6, interpolation = 'None');colorbar()
						subplot(132);imshow(self.FG_90[0]*1e6, interpolation = 'None');colorbar()
						subplot(133);imshow((self.FG_150[0] - self.FG_90[0])*1e6, interpolation = 'None');colorbar();show()#;quit()
					'''

					if len(FG) < len(CMB):
						FG_tmp = np.zeros ( (len(CMB), ny, nx) )
						FG_tmp[0] = FG
						FG = FG_tmp

				else: #I do not know what foregrounds you are asking
					FG = None

				#imshow(tSZ_emission * 1e6);colorbar();show();quit()
				#subplot(224);imshow(FG[0]);colorbar();title(fg_T)
				#show();quit()

				#imshow(FG[0]);colorbar();show();quit()	
				
				if self.is_seq(FG): CMB += FG

				#print FG.shape, CMB.shape
				
				if self.is_seq(tSZ_emission): #DG-Po corresponding tSZ emission here
					CMB[0] += tSZ_emission

				#if 1==1:self.fn_plot_pow_spec(mapparams, CMB,also_plot = 1);quit()

				if 0 == 1:
					for tqu in range(tqulen):
						if len(FG) > tqu:
							subplot(2,3,tqu+1);imshow(CMB[tqu] - FG[tqu]);colorbar()
						else:
							subplot(2,3,tqu+1);imshow(CMB[tqu]);colorbar()
					for tqu in range(tqulen):
						subplot(2,3,tqu+1+tqulen);imshow(CMB[tqu]);colorbar()
					show();quit()

			else:
				self.FG_added = 0

			"""
			for tqucnt in range(len(CMB)):
				subplot(1,len(CMB),tqucnt+1);imshow(CMB[tqucnt]);colorbar()
			show()
			"""
			##################################################################
			
			if (1): #adding a point source here
				if nu ==150 and self.add_psource:
					s1 = nx/2
					s2 = s1+1
					self.howmuchpointsourceT = 2.5e-4
					CMB[0,s1:s2,s1:s2] += self.howmuchpointsourceT
					##imshow(CMB[0]);colorbar();show();sys.exit()
			
			##################################################################
			##################################################################
			#20170801 - ILC residual - Gaussian relaisation of the residual power
			try:
				#self.ilc_residual_file = self.run_s2_dict['ilc_residual_file']
				self.ilc_residual_file = self.inidic['ilc_residual_file']
				self.add_ILC_residual = 1
			except:
				self.add_ILC_residual = 0

			if self.add_ILC_residual:
				ILC_RESIDUAL = self.fn_add_ILC_residual(self.ilc_residual_file, mapparams, 1)[0]
				##print ILC_RESIDUAL.shape, CMB.shape
				#CMB += ILC_RESIDUAL
				for ii in range(len(CMB)):
					CMB[ii] += ILC_RESIDUAL[ii]

				'''
				#imshow(CMB[0]);colorbar();show();quit()
				for tqucnt in range(len(CMB)):
					subplot(1,len(CMB),tqucnt+1);imshow(CMB[tqucnt]);colorbar()
				show();quit()
				'''

			#if 1==1:self.fn_plot_pow_spec(mapparams, CMB,also_plot = 1);quit()
			##################################################################
			#add tSZ component to T here
			
			if nu == 90:
				tsz_comp = self.beta_model_90ghz
			elif nu == 150:
				tsz_comp = self.beta_model_150ghz
			elif nu == 'tszfree':
				tsz_comp = None

			if self.is_seq(tsz_comp):
				CMB[0] += tsz_comp
				'''
				subplot(131);imshow((CMB[0] - tsz_comp) * 1e6);colorbar()
				subplot(132);imshow(tsz_comp*1e6);colorbar()
				subplot(133);imshow(CMB[0] * 1e6);colorbar();show();quit()
				'''
			##################################################################
			#20180813 - storing this for alm calculation for LGMCA sims before doing the beam smoothing
			self.CMB_currsim = CMB
			##################################################################
			#convolve with beam now because E and B maps are made from Q, U
			#use_beam = 0
			#subplot(121);imshow(CMB[0]);colorbar();title(nu);
			if use_beam>0:#beam is compulsory
				CMB = np.fft.ifft2( np.fft.fft2( CMB ) * Bl ).real
			#subplot(122);imshow(CMB[0]);colorbar();title(nu);show();sys.exit()
			#clf();subplot(121);imshow(CMB_ori[0]);colorbar();subplot(122);imshow(CMB[0]);colorbar();show();quit()
			##################################################################

			##################################################################
			#new addition: 20160912
			#get lensed E, and B now from Q, and U
			if self.perform_EB: #for E, B
				Q_FFT, U_FFT = np.fft.fft2(CMB[1]),np.fft.fft2(CMB[2])
				angle = 2.*np.arctan2(lx, -ly)
				E = np.fft.ifft2( np.cos(angle) * Q_FFT + np.sin(angle) * U_FFT ).real #lensed E
				B = np.fft.ifft2( -np.sin(angle) * Q_FFT + np.cos(angle) * U_FFT ).real #lensed B

				CMB = np.append(CMB,[E,B],axis=0)
			##################################################################

			##################################################################
			#degrading resoltuion if necessary
			if degrade_resol_fac == 1:
				SIMMAPS_L[cnt] = CMB
			else:
				for tquebcnt in range(tqulen):
					SIMMAPS_L[cnt, tquebcnt] = self.downsample_map(CMB[tquebcnt],degrade_resol_fac)

			#imshow(CMB_ori[1] - SIMMAPS_L[cnt,1]);colorbar();show();quit()
			##################################################################

			#debug stuff - see if the CMB maps look okay here
			show_cmb_plots = 0
			if show_cmb_plots:
				TWODTF = self.fn_get_HPF(mapparams, minel = self.inidic['min_el'], maxel = self.inidic['max_el'])
				if self.tqulen == 1: #only T map
					subplot(211);css=imshow(CMB[0]);colorbar();
					subplot(212);css=imshow(CMB[0] - SIMMAPS_L[cnt,0]);colorbar();show();quit()
				else:
					clf()
					if self.tqulen == 5:
						nrows,ncols = self.tqulen - 2,self.tqulen
					else:
						nrows,ncols = self.tqulen, self.tqulen

					Q_UNLEN, U_UNLEN = np.fft.fft2(CMB_ori[1]),np.fft.fft2(CMB_ori[2])
					angle = 2.*np.arctan2(lx, -ly)
					E_UNLEN = np.fft.ifft2( np.cos(angle) * Q_UNLEN + np.sin(angle) * U_UNLEN ).real #unlensed E
					B_UNLEN = np.fft.ifft2( -np.sin(angle) * Q_UNLEN + np.cos(angle) * U_UNLEN ).real #unlensed B
					#B_UNLEN = B_UNLEN * 0.

					for ppp in range(ncols): #T,Q,U,E,B maps
						if ppp<=2:
							MAP = CMB_ori[ppp]
						elif ppp == 3:
							MAP = E_UNLEN
						elif ppp == 4:
							MAP = B_UNLEN

						if ppp<3:
							tit = self.cambmapvars[ppp]
						elif ppp == 3:
							tit = 'E modes'
						elif ppp == 4:
							tit = 'B modes'

						MAP_L = SIMMAPS_L[cnt,ppp]

						#MAP = np.fft.ifft2( np.fft.fft2(MAP) * TWODTF[ppp] ).real
						#MAP_L = np.fft.ifft2( np.fft.fft2(MAP_L) * TWODTF[ppp] ).real

					
						MAP_DIFF = MAP_L - MAP
						#MAP_DIFF = np.fft.ifft2( np.fft.fft2(MAP_DIFF) * TWODTF[0] ).real
						#MAP_DIFF = np.fft.ifft2( np.fft.fft2(MAP_DIFF) * self.WIENER_FILTER_WHU[0] ).real
						ax = subplot(nrows,ncols,ppp+1)
						css=imshow(MAP*1e6);colorbar()

						title(tit, fontsize = 14)

						ax = subplot(nrows,ncols,ppp+1+ncols)
						css=imshow(MAP_L*1e6)
						#if ppp+1+ncols == 10:
						#	css.set_clim(-0.12,.12)
						colorbar()
						title('Lensed %s' %(tit), fontsize = 10)
						ax = subplot(nrows,ncols,ppp+1+(2*ncols))
						css=imshow(MAP_DIFF*1e6);colorbar()#css.set_clim(-1,1);colorbar()
						title('LEN - UNLEN: %s' %(tit), fontsize = 10)
					show();quit()

			##################################################################
			#delete all the used objects - does not matter though
			del CMB
			if self.perform_EB: #for E, B
				del Q_FFT, U_FFT
			##################################################################

                end_time = time.time()
		
		if self.inidic['CMB_outputscale'] == 1:
			SIMMAPS_L *= 1e6
                if not self.quiet:
	                logline = '\t\tlensing complete. time take = %s seconds' %(end_time-start_time)
	                logfile = open(self.log_file,'a')
	                logfile.writelines('%s\n' %(logline))
	                logfile.close()
	                print logline
		#self.SIMMAPS_L = SIMMAPS_L

          
               
                return SIMMAPS_L

	def fn_rich_mass_M17(self, rich, z_val, A = 2.35e14, lambda_pivot = 30., z_pivot = 0.5, alpha_fit = 1.12, beta_fit = 0.18):
		#definitions from M17 (https://arxiv.org/abs/1610.06890)
		return A * (rich/lambda_pivot)**alpha_fit * ( (1+z_val) / (1+z_pivot) )**beta_fit


	def fn_flatsky_to_healpix(self, flat_sky_map, RA, DEC, nside = 2048, fill_value=None, also_hit_map = 0):##-1.6375e+30):

		import healpy as H
		nopixels = H.nside2npix(nside)
		MAP=np.zeros(nopixels)
		HIT=np.zeros(nopixels)
		for r in range(np.shape(RA)[0]):
			PP=H.pixelfunc.ang2pix(nside,np.radians(90.-DEC[r,:]),np.radians(RA[r,:]))
			MAP[PP]+=flat_sky_map[r,:]
			HIT[PP]+=1
		##from IPython import embed; embed()
		MAP[HIT>0.]/=HIT[HIT>0.]
		##print PP, MAP[PP]
		if fill_value<>None:
			MAP[HIT == 0.] = fill_valu

		if also_hit_map:
			return MAP, HIT
		else:
			return MAP

	def fn_ngrad_from_rich(self, rich, sptpol_clusters = 0, des_rich_threshold_ngrad = 60.):

		if not sptpol_clusters:
			if rich<=des_rich_threshold_ngrad:
				ngrad = 2000#1500
			else:
				ngrad = 1000

			#20180124 - fix it to 2000
			### ngrad = 2000
		else:
			if rich<=5.:
				ngrad = 2000#1500
			else:
				ngrad = 1000

		return ngrad

	def fn_TF_tszfree(self, mapparams, which_FG = 'all_Gaussian'):
		
		nu = 90.
		FG_els_90, FG_cls_90 = self.fn_get_foreground_power(which_FG, mapparams, 1, perform_lensing = 0, nu1=nu, nu2=nu, just_Cls = 1) #in K^2 already
		nu = 150.
		FG_els, FG_cls = self.fn_get_foreground_power(which_FG, mapparams, 1, perform_lensing = 0, nu1=nu, nu2=nu, just_Cls = 1) #in K^2 already

		clusterstuff = cluster_stuff()
		ysz_Tsz_conv_fac_150 = clusterstuff.compton_y_to_delta_Tcmb(150e9) #conversion factor for 150 GHz
		ysz_Tsz_conv_fac_90 = clusterstuff.compton_y_to_delta_Tcmb(90e9) #conversion factor for 90 GHz

		factor = ysz_Tsz_conv_fac_90/ysz_Tsz_conv_fac_150
		factor = factor**2. #power spectra

		#calib_factors
		try:
			FG_cls_90 = FG_cls_90 * self.inidic['calib_fac_90']**2.
			FG_cls = FG_cls * self.inidic['calib_fac_150']**2.
		except:
			pass

		#get beam
		Bl_90_1d = self.fn_read_SPTpol_el_space_beam(mapparams, nu = 90, return_1d = 1)
		Bl_90_1d = np.interp(FG_els, Bl_90_1d[0], Bl_90_1d[1]) #interpolating

		Bl_150_1d = self.fn_read_SPTpol_el_space_beam(mapparams, nu = 150, return_1d = 1)
		Bl_150_1d = np.interp(FG_els, Bl_150_1d[0], Bl_150_1d[1]) #interpolating

		Bl_gradient_map_1d = Bl_90_1d / Bl_150_1d
		Bl_gradient_map_1d[Bl_gradient_map_1d == np.inf] = 0.
		Bl_gradient_map_1d[np.where(np.isnan(Bl_gradient_map_1d))] = 0.

		FG_cls = FG_cls * Bl_gradient_map_1d**2.
		FG_cls_tSZfree = ( (FG_cls * factor) - FG_cls_90 ) / (factor - 1.)
		FG_tSZfree = np.asarray([FG_els, FG_cls_tSZfree]).T
		FG_tSZfree_lmat = self.Cls2CLS(FG_tSZfree,mapparams)[0]

		return FG_tSZfree_lmat

	def fn_get_sehgal_FG_dic(self, nu, which_comp = 'all', sehgal_folder = 'data/sehgal_2009_sims/tsz_cib_radio_ksz_extracts/'):
		clsdicfname = '%s/Cls_%s.pkl.gz' %(sehgal_folder, nu)
		selgal_fg_dic = pickle.load(gzip.open(clsdicfname, 'rb'))

		assert which_comp in [selgal_fg_dic.keys(), 'all']

		if which_comp == 'all':
			for compcnt, compkey in enumerate( selgal_fg_dic ):
				if compcnt == 0:
					ells, Cls = selgal_fg_dic[compkey]
				else:
					Cls += selgal_fg_dic[compkey][1]
		else:
			ells, Cls = selgal_fg_dic[which_comp]
		ells = ells[:-1]

		return ells, Cls

	def fn_get_kappa_QE(self, OBSMAP, mapparams, Dls_len, Dls_unlen, tszfree = 0, OBSMAP2 = None, richval = None, use_data = 1, curl_test = 0, noise_weight = 1.):

		import modules.qe_funcs as qe
		###Parameters of QE### from EBX
		#l cuts
		nlens = 60000 #babbloo changing to 27535
		nbeam = 60000
		#what is this? -related to pixel size????
		#ngaus = 27535
		#ngaus = self.inidic['ngaus']
		try:
			ngaus = int( np.sqrt(8. * np.log(2.)) / np.radians(self.exp_beam/60.) )
		except:
			ngaus = int( np.sqrt(8. * np.log(2.)) / np.radians(self.inidic['exp_beam']/60.) )
		
		if not self.inidic['cut_below_beam_scale']:
			ngaus = 6756 #1.2 am beam

		#change ngaus is szfree maps are used for lensing reconstrcution too - 20180530
		if self.just_SZfreemap:
			ngaus = int( np.sqrt(8. * np.log(2.)) / np.radians(self.exp_beam_90ghz/60.) )
			if self.use_lgmca:
				ngaus = 2048

		#don't cut at beam scale - make it an arbitarily large value
		try:
			ngaus = self.inidic['ngaus']
		except:
			pass
						
		nlensmin = 2
		#Filter for gradient
		#7777777777777777777777777777777777777777777777777777777
		#20171227 - select ngrad based on richness - otherwise gradient will be incorrectly estimated
		if richval == None or self.inidic['cross_maps'] == 0:
			ngrad = self.inidic['ngrad']#99999999999.#1500
			### ngrad = 2000.
		else:
			ngrad = self.fn_ngrad_from_rich(richval, sptpol_clusters = self.sptpol_clusters)
		### ngrad = 1000

		print richval, ngrad

		#what is this?
		#nclfilt = 5400

		C_l_len_lmat = self.Dls2map(Dls_len, mapparams, return_Clmat = 1)
		### modify Cls amplitude
		try:
			C_l_len_lmat = C_l_len_lmat * self.inidic['mod_lensed_cls_amp']
		except:
			pass

		C_l_unl_lmat = self.Dls2map(Dls_unlen, mapparams, return_Clmat = 1)

		#M1, M2 = np.fft.fftshift(C_l_len_lmat * 1e12).real, np.fft.fftshift(C_l_unl_lmat * 1e12).real
		#subplot(121);imshow(M1 - M2, vmin = 0., vmax = 1e-4);colorbar()
		#subplot(122);imshow(M2, vmin = 0., vmax = 1.);colorbar();show();quit()
		#loglog(Dls_len[:,0], Dls_len[:,1]);loglog(Dls_unlen[:,0], Dls_unlen[:,1]);show();quit()

		lx, ly = self.get_lxly(mapparams)
		l2d = np.sqrt(lx**2. + ly**2.)
		nx, ny, dx, dy = mapparams
		dxdy = dx * self.arcmins2radians * dy  * self.arcmins2radians

		#20180315 - modify TF here
		try:
			self.lx_for_qe
		except:
			self.lx_for_qe = 0
		if self.lx_for_qe > 0:
			transfer_lmat = self.fn_get_EBAX_2016_anal(mapparams, l1=self.inidic['l1'], l2=self.lx_for_qe, l3=self.inidic['l3'])[0]
			###transfer_lmat_ori = self.fn_get_EBAX_2016_anal(mapparams, l1=self.inidic['l1'], l2=self.inidic['l2'], l3=self.inidic['l3'])[0]
			self.transfer_lmat = transfer_lmat

		try:
			transfer_lmat = self.transfer_lmat		
		except:
			if self.inidic['add_TF'] == 0:
				transfer_lmat = self.fn_get_HPF(mapparams, ideal = 1)[0]
			elif self.inidic['add_TF'] == 1:			
				transfer_lmat = self.fn_get_HPF(mapparams, ideal = 0, minel=self.inidic['min_el'], maxel=self.inidic['max_el'])[0]
			elif self.inidic['add_TF'] == 2:
				transfer_lmat = self.TWODTF[0]
			elif self.inidic['add_TF'] == 3:
				#transfer_lmat = self.fn_get_EBAX_2016_anal(mapparams, minel=self.inidic['min_el'], maxel=self.inidic['max_el'])[0]
				transfer_lmat = self.fn_get_EBAX_2016_anal(mapparams, l1=self.inidic['l1'], l2=self.inidic['l2'], l3=self.inidic['l3'])[0]

			self.transfer_lmat = transfer_lmat
		

		transfer_lmat_grad = self.transfer_lmat
		try:
			if self.inidic['add_TF_grad'] == 1:
				transfer_lmat_grad= self.fn_get_HPF(mapparams, minel=self.inidic['min_el_grad'], maxel=self.inidic['max_el_grad'], all_iso = 0)[0]
				if self.use_lgmca:
					transfer_lmat_grad = self.fn_get_HPF(mapparams, ideal = 0, minel=self.inidic['min_el_grad'], maxel=self.inidic['max_el_grad'])[0]
		except:
			pass
		
		if self.inidic['use_beam'] == 1:
			beam_fwhm = self.exp_beam*(1./60.)*np.pi/180.     
			beam_lmat = np.exp(-(beam_fwhm**2.)*(l2d**2.)/(16.*np.log(2.)))
		elif self.inidic['use_beam'] == 2 or self.inidic['use_beam'] == 4:
			beam_lmat = self.Bl

		if self.nx_large<>None and self.ny_large<>None:
			beam_lmat = self.fn_beam_stuff(mapparams, exp_beam = self.inidic['use_beam'], nu  =150, return_beam = 1)

		try:
			N_lmat = self.N_lmat
			if self.is_seq(OBSMAP2):
				N_lmat_1 = self.N_lmat_1
		except:
			if self.noise_present == 1:

				if 0 ==1:#tszfree: !!! check why??
					print '\n\n\t\tFatal error: Not implement white noise for tSZfree yet. Aborting here ...'
					sys.exit(0)
				
				
				noise_level = self.expnoiselevel[0]
				DeltaT = noise_level*1e-6*(1./60.)*(np.pi/180.) #K-radian
				N_lmat = (DeltaT**2.) + np.zeros(l2d.shape)
				self.N_lmat = N_lmat
			

				if self.is_seq(OBSMAP2):
					DeltaT_1 = self.inidic['noise_level_grad_map']*1e-6*(1./60.)*(np.pi/180.) #K-radian
					N_lmat_1 = (DeltaT_1**2.) + np.zeros(l2d.shape)
					self.N_lmat_1 = N_lmat_1

			elif self.noise_present == 2:

				try:
					noisesims = self.noisesims
				except:
					noisesims = noise_sims()
					self.noisesims = noisesims

				if not tszfree:
					#N_lmat = noisesims.SPT_PSD_150[0] / 1e12
					#N_lmat[N_lmat == np.inf] = 0.

					try:
						N_lmat = self.N_lmat
						if self.is_seq(OBSMAP2):
							N_lmat_1 = self.N_lmat_1
					except:
						if not self.CR_maps:
							noise_cls_file = 'data/sanjay_maps_201705xx/noise_Cls/150ghz/combined_noise_cls_42bundlesused.pkl.gz'
						else:
							noise_cls_file = 'data/CR_maps_20170910/noise_Cls/150ghz/combined_noise_cls_89bundlesused.pkl.gz'

						noise_cls = pickle.load(gzip.open(noise_cls_file,'rb'))
						noise_cls[:,1] /= 1e12 #20180119 - missed this earlier? shit!

						N_lmat = self.Cls2CLS(noise_cls,mapparams)[0]
						self.N_lmat = N_lmat

					if self.is_seq(OBSMAP2): #20171103
						if not self.CR_maps:
							noise_cls_grad_file = 'data/sanjay_maps_201705xx/noise_Cls/tszfree/combined_noise_cls_42bundlesused.pkl.gz'
						else:
							noise_cls_grad_file = 'data/CR_maps_20170910/noise_Cls/tszfree/combined_noise_cls_89bundlesused.pkl.gz'

						if self.use_lgmca:
							noise_cls_grad_file = 'data/CR_maps_20170910/noise_Cls//lgmca/lgmca_noise_cls.pkl.gz'

						noise_cls_grad = pickle.load(gzip.open(noise_cls_grad_file,'rb'))
						noise_cls_grad[:,1] /= 1e12 #20180119 - missed this earlier? shit!
						N_lmat_1 = self.Cls2CLS(noise_cls_grad,mapparams)[0]

						self.N_lmat_1 = N_lmat_1
				else:

					if not self.CR_maps:
						noise_cls_file = 'data/sanjay_maps_201705xx/noise_Cls/tszfree/combined_noise_cls_42bundlesused.pkl.gz'
					else:
						noise_cls_file = 'data/CR_maps_20170910/noise_Cls/tszfree/combined_noise_cls_89bundlesused.pkl.gz'
					noise_cls = pickle.load(gzip.open(noise_cls_file,'rb'))
					noise_cls[:,1] /= 1e12
					N_lmat = self.Cls2CLS(noise_cls,mapparams)[0]
					#imshow(np.fft.fftshift(N_lmat[0] ));colorbar();show();quit()

					self.N_lmat = N_lmat

			### from IPython import embed; embed()
			'''
			subplot(121);imshow(N_lmat);colorbar()
			subplot(122);imshow(N_lmat_1);colorbar();show();quit()
			'''

			'''
			#now add FG elmat
			try:
				FG_lmat = self.FG_lmat
			except:
				ilc_res_file = 'data/foreground/cleaning/Cls_ILC_sptpol.pkl.gz'
				ilc_res_els, ilc_res_cls = pickle.load(gzip.open(ilc_res_file,'rb'))['Cl_residual']
				ilc_res_cls/=1e12
				ilc_res = np.asarray([ilc_res_els, ilc_res_cls]).T
				FG_lmat = self.Cls2CLS(ilc_res,mapparams)[0]

				self.FG_lmat = FG_lmat
			'''

			try:
				fit_FG = self.inidic['fit_FG']
			except:
				fit_FG = None

			'''#20180430  - turn this off and check
			#20180427 - if Sehgal sims added - use the FG in W. filter
			### from IPython import embed; embed()
			try:
				if self.sehgal_fg_added:
					fit_FG = 'all_Gaussian'
			except:
				pass
			'''

			self.fit_FG = fit_FG

			if self.fit_FG <> None or self.sehgal_fg_added:

				if self.fit_FG=='all_Gaussian':
					nu = 150.
					FG_els, FG_cls = self.fn_get_foreground_power('all_Gaussian', mapparams, 1, perform_lensing = 0, nu1=nu, nu2=nu, just_Cls = 1) #in K^2 already
				elif self.sehgal_fg_added:
					nu = 148
					FG_els, FG_cls = self.fn_get_sehgal_FG_dic(nu)

				try:
					FG_cls = FG_cls * self.inidic['calib_fac_150']**2.
				except:
					pass
				FG = np.asarray([FG_els, FG_cls]).T
				self.FG_lmat = self.Cls2CLS(FG,mapparams)[0]

				### #N_lmat = N_lmat + FG_lmat
				### self.N_lmat = self.N_lmat + FG_lmat #modified on 20171222

				if self.is_seq(OBSMAP2): #20171222
					if self.fit_FG=='all_Gaussian':
						nu = 90.
						FG_els_90, FG_cls_90 = self.fn_get_foreground_power('all_Gaussian', mapparams, 1, perform_lensing = 0, nu1=nu, nu2=nu, just_Cls = 1) #in K^2 already
					elif self.sehgal_fg_added:
						nu = 90
						FG_els_90, FG_cls_90 = self.fn_get_sehgal_FG_dic(nu)
					
					if tszfree or use_data or self.inidic['perform_tSZ_removal']:
						if self.fit_FG=='all_Gaussian':
							nu = 150.
							FG_els, FG_cls = self.fn_get_foreground_power('all_Gaussian', mapparams, 1, perform_lensing = 0, nu1=nu, nu2=nu, just_Cls = 1) #in K^2 already
						elif self.sehgal_fg_added:
							nu = 148
							FG_els, FG_cls = self.fn_get_sehgal_FG_dic(nu)							
						clusterstuff = cluster_stuff()
						ysz_Tsz_conv_fac_150 = clusterstuff.compton_y_to_delta_Tcmb(150e9) #conversion factor for 150 GHz
						ysz_Tsz_conv_fac_90 = clusterstuff.compton_y_to_delta_Tcmb(90e9) #conversion factor for 90 GHz

						factor = ysz_Tsz_conv_fac_90/ysz_Tsz_conv_fac_150
						factor = factor**2. #power spectra

						#calib_factors
						try:
							FG_cls_90 = FG_cls_90 * self.inidic['calib_fac_90']**2.
							FG_cls = FG_cls * self.inidic['calib_fac_150']**2.
						except:
							pass

						#get beam
						Bl_90_1d = self.fn_read_SPTpol_el_space_beam(mapparams, nu = 90, return_1d = 1)
						Bl_90_1d = np.interp(FG_els, Bl_90_1d[0], Bl_90_1d[1]) #interpolating

						Bl_150_1d = self.fn_read_SPTpol_el_space_beam(mapparams, nu = 150, return_1d = 1)
						Bl_150_1d = np.interp(FG_els, Bl_150_1d[0], Bl_150_1d[1]) #interpolating

						Bl_gradient_map_1d = Bl_90_1d / Bl_150_1d
						Bl_gradient_map_1d[Bl_gradient_map_1d == np.inf] = 0.
						Bl_gradient_map_1d[np.where(np.isnan(Bl_gradient_map_1d))] = 0.

						FG_cls = FG_cls * Bl_gradient_map_1d**2.
						
						FG_cls = FG_cls ## !!!hardcoded or change it accordingly
						FG_cls_tSZfree = ( (FG_cls * factor) - FG_cls_90 ) / (factor - 1.)
						FG_tSZfree = np.asarray([FG_els, FG_cls_tSZfree]).T
						self.FG_tSZfree_lmat = self.Cls2CLS(FG_tSZfree,mapparams)[0]

						### self.N_lmat_1 = self.N_lmat_1 + FG_tSZfree_lmat #modified on 20171222
					else:
						
						try:
							FG_cls_90 = FG_cls_90 * self.inidic['calib_fac_90']**2.
						except:
							pass
						FG_90 = np.asarray([FG_els_90, FG_cls_90]).T
						self.FG_90_lmat = self.Cls2CLS(FG_90,mapparams)[0]

						#N_lmat = N_lmat + FG_lmat
						### self.N_lmat_1 = self.N_lmat_1 + FG_90_lmat #included on 20180123

			#loglog(ilc_res[:,0], ilc_res[:,1]);loglog(noise_cls[:,0], noise_cls[:,1]);show()

			#subplot(121);imshow(np.fft.fftshift(N_lmat));colorbar();subplot(122);imshow(np.fft.fftshift(FG_lmat));colorbar();show();sys.exit()					FG_cls_tSZfree = ( (FG_cls * factor) - FG_cls_90 ) / (factor - 1.)
			#elif self.whichmap == '150ghz_500sqdeg_LBleem':

		#imshow(N_lmat);colorbar();show();quit()

		### from IPython import embed; embed()
		if not self.quiet:
			logline = '\t\t noise weight is %s' %(noise_weight);logfile = open(self.log_file,'a');logfile.writelines('%s\n' %(logline))
			logfile.close()
			######print logline

		if self.noise_present == 0 and self.add_ILC_residual:
			noise_ilc_cls = np.asarray( [self.ilc_els, self.ilc_Cls_res] ).T
			N_lmat_ilc = self.Cls2CLS(noise_ilc_cls,mapparams)[0]
			N_lmat = N_lmat_ilc
			self.N_lmat = N_lmat

		#correct N_lmat using weights
		N_lmat_with_weights = N_lmat/ noise_weight
		if self.is_seq(OBSMAP2):
			N_lmat_1_with_weights = N_lmat_1/ noise_weight

		#moving this here from above. Adding FGs after correct noise-weighting
		if self.fit_FG <> None:
			N_lmat_with_weights = N_lmat_with_weights + self.FG_lmat #modified on 20171222
			if self.is_seq(OBSMAP2): #20171222
				N_lmat_1_with_weights = N_lmat_1_with_weights #+ self.FG_tSZfree_lmat #modified on 20171222
			else:
				##N_lmat_1_with_weights = N_lmat_1_with_weights + self.FG_90_lmat #included on 20180123
				N_lmat_with_weights = N_lmat_with_weights + self.FG_lmat #changed on 20180731

		#beam_lmat_transfer_lmat = self.fn_convolve_using_QL(beam_lmat, transfer_lmat)
		beam_lmat_transfer_lmat = transfer_lmat * beam_lmat
		N_beam_lmat = N_lmat_with_weights/((beam_lmat_transfer_lmat)**2.)

		if not self.is_seq(OBSMAP2):
			N_beam_lmat_1 = N_beam_lmat_2 = N_beam_lmat
			beam_lmat_transfer_lmat_1 = None
		else:
			N_beam_lmat_2 = N_beam_lmat
			#get this for the gradient map
			beam_lmat_1 = self.Bl_gradient_map
			if self.nx_large<>None and self.ny_large<>None:
				beam_lmat_1 = self.fn_beam_stuff(mapparams, use_beam = self.inidic['use_beam'], nu = 90, return_beam = 1)
			beam_lmat_transfer_lmat_1 = transfer_lmat_grad * beam_lmat_1
			N_beam_lmat_1 = N_lmat_1_with_weights/((beam_lmat_transfer_lmat_1)**2.)

		##from IPython import embed; embed()

		if self.just_SZfreemap:
			N_beam_lmat_2 = np.copy( N_beam_lmat_1 )

		#subplot(121);imshow(np.fft.fftshift(N_beam_lmat_1));colorbar()
		#subplot(122);imshow(np.fft.fftshift(N_beam_lmat_2));colorbar();show()##;sys.exit()


		try:
			try_mod_code = self.inidic['try_mod_code']
		except:
			try_mod_code = 0
		l_G = ngrad
		lmax_inv_var = self.inidic['lmax_inv_var']

		C_l_unl_lmat_dic = {}
		C_l_len_lmat_dic = {}

		testing = 0
		if testing:
			OBSMAP2 = None
			N_beam_lmat_1 = N_beam_lmat_2 = N_beam_lmat
		if not try_mod_code:

			C_l_unl_lmat = C_l_unl_lmat[0]
			C_l_len_lmat = C_l_len_lmat[0]
			OBSMAP = OBSMAP[0]
			if self.is_seq(OBSMAP2):
				OBSMAP2 = OBSMAP2[0]

			#Gradient weighting and inverse variance weighting
			#print '\t\t\tGenerating gradient and inv var weights',
			#Gradient weights (W^TT)
			#weight_gradient_lmat = C_l_unl_lmat/(C_l_len_lmat + N_beam_lmat) #EB
			weight_gradient_lmat = C_l_unl_lmat/(C_l_len_lmat + N_beam_lmat_1) #EB
			
			above_lg = np.where(l2d > l_G)
			weight_gradient_lmat[above_lg] = 0.0


			#imshow(np.fft.fftshift(weight_gradient_lmat), extent = [np.min(lx), np.max(lx), np.min(ly), np.max(ly)]);colorbar();show();quit()

			#Inverse variance weights (W^T)
			#weight_inv_var_lmat = l2d /(C_l_len_lmat + N_beam_lmat) #Maturi et al. filter
			weight_inv_var_lmat = 1. /(C_l_len_lmat + N_beam_lmat_2)
			weight_inv_var_lmat[np.isnan(weight_inv_var_lmat)] = 0.
	 		

 			if self.matched_filter:# self.add_psource or self.is_seq( self.beta_model_150ghz ): #working on this: 20180522

				RADIUS = np.sqrt(self.RA**2 + self.DEC**2.)
				RADIUS *= 60.
				theta_c = 1.0


				##Bl_for_filter = self.fn_beam_stuff(mapparams, exp_beam = 1.2, return_beam = 1)
				Bl_for_filter = self.fn_beam_stuff(mapparams, exp_beam = 1.5, return_beam = 1)

				'''
				beta_value = -1.
				tszmodel = (1.0 + (RADIUS/theta_c) ** 2. ) ** beta_value
				'''

				##from IPython import embed; embed()
				aperture_1 = theta_c ###1.0 #arcmins
				aperture_2 = 10.0 #arcmins
				#if self.inidic['amp_frm_input_tsz']: 
				#	full_fit = self.tsz_amp_from_input(mapparams, OBSMAP)
				#elif not self.inidic['amp_frm_input_tsz'] :
				#full_fit= self.fit_gaussian_template(mapparams,OBSMAP)
				#filtered_map = self.filtered_cmb(mapparams, OBSMAP)
				#from IPython import embed;embed()
				
				#r_full_fit = self.robust_gaussian_fit(mapparams, OBSMAP)
				
				#input_sz = np.fft.ifft2(np.fft.fft2(self.beta_model_150ghz)*self.Bl).real
				#full_fit = r_full_fit
				if self.inidic['from_cmb_map']:

					c_x,c_y,full_fit = self.filtered_cmb( mapparams, OBSMAP)
				else:

					c_x_r,c_y_r,full_fit = self.robust_gaussian_fit(mapparams, OBSMAP)
				
				#from IPython import embed;embed()
				
				#full_fit = data_fitted.reshape(nx,ny)
				#full_fit = self.fn_fit_Arnaud_template(mapparams,OBSMAP)
				"""
				amplitude_aperture_1 = np.mean(OBSMAP[RADIUS<=aperture_1])
				amplitude_aperture_2 = np.mean(OBSMAP[(RADIUS>aperture_1) & (RADIUS<=aperture_2)])##np.mean(OBSMAP[RADIUS<=aperture_2])
				amplitude = amplitude_aperture_1 - amplitude_aperture_2
				
				what_to_remove = np.fft.ifftshift( np.fft.ifft2( Bl_for_filter * self.TWODTF ) ).real
				what_to_remove = what_to_remove/np.max(what_to_remove) * amplitude
				what_to_remove[0] = self.beta_model_150ghz
				what_to_remove[0]= np.fft.ifft2( np.fft.fft2(self.beta_model_150ghz ) * self.Bl ).real
				

				"""
				#xx =  OBSMAP - full_fit
				#from IPython import embed;embed()
				
				#sys.exit()
				#from pylab import *
				#print np.mean(OBSMAP)
				#OBSMAPMAP -= what_to_remove1[0]
				
				OBSMAP -= full_fit

				if (0):#self.add_psource:

					subplot(221);imshow(OBSMAP*1e6, vmin = -100, vmax =100.);colorbar();
					###if use_source_template: subplot(222);imshow(np.fft.ifft2(source_template).real);colorbar();
					if self.is_seq( self.beta_model_150ghz ):
						tszsignal = np.fft.ifft2( np.fft.fft2(self.beta_model_150ghz) * self.Bl).real * 1e6
						subplot(222);imshow(tszsignal);colorbar();
						axhline(nx/2., color ='k', lw = 0.1);axvline(nx/2., color ='k', lw = 0.1)
					subplot(223);imshow(what_to_remove*1e6);colorbar();
					axhline(nx/2., color ='k', lw = 0.1);axvline(nx/2., color ='k', lw = 0.1)
					dummy =np.copy(OBSMAP) - what_to_remove
					subplot(224);imshow(dummy*1e6, vmin = -100, vmax =100.);colorbar();show();sys.exit()



			#sys.exit()
			#lmax_inv_var = 15000#nlens+ngrad# - ensures that we're not using ells beyond range of input C_ls
			self.stacked_fitted_cmb = self.stacked_fitted_cmb + OBSMAP
			above_inv_var_lmax = np.where(l2d > lmax_inv_var)
			weight_inv_var_lmat[above_inv_var_lmax] = 0.0

			#imshow(np.fft.fftshift(weight_inv_var_lmat), extent = [np.min(lx), np.max(lx), np.min(ly), np.max(ly)]);colorbar();show();quit()

			#normalization of quadratic estimator
			#print '\t\t\t\tGenerating QE normalization',
			##from IPython import embed; embed()
			norm_lmat = qe.get_norm(lx, ly, C_l_unl_lmat, weight_gradient_lmat, weight_inv_var_lmat, dxdy)
			
			#subplot(121);imshow(OBSMAP);colorbar();subplot(122);imshow(OBSMAP2);colorbar();show();sys.exit()

			#kappa_qe = qe.get_kappa(lx, ly, weight_gradient_lmat, weight_inv_var_lmat, OBSMAP, norm_lmat, dxdy, ngaus)
			#kappa_qe = qe.get_kappa(lx, ly, beam_lmat_transfer_lmat, weight_gradient_lmat, weight_inv_var_lmat, OBSMAP, norm_lmat, dxdy, ngaus)
			#print "highly_filtered"
			kappa_qe = qe.get_kappa(lx, ly, beam_lmat_transfer_lmat, weight_gradient_lmat, weight_inv_var_lmat, OBSMAP, norm_lmat, dxdy, ngaus, OBSMAP2, curl_test = curl_test, beam_lmat_transfer_lmat_grad = beam_lmat_transfer_lmat_1)

			if testing: print kappa_qe, 'Original code'
			#C_l_unl_lmat_dic['TT'] = C_l_unl_lmat
			#C_l_len_lmat_dic['TT'] = C_l_len_lmat

		else:
			#subplot(131);imshow(OBSMAP[0]);colorbar();subplot(132);imshow(OBSMAP[1]);colorbar();subplot(133);imshow(OBSMAP[2]);colorbar();show();sys.exit()
			import modules.qe_funcs_temp_pol as qe_temp_pol
			estimator = 'TT'#'TT'
			N_beam_lmat_dic = {}
			if len(C_l_unl_lmat)>1:
				C_l_unl_lmat_dic['TT'], C_l_unl_lmat_dic['EE'], C_l_unl_lmat_dic['BB'], C_l_unl_lmat_dic['TE'] = C_l_unl_lmat
				C_l_len_lmat_dic['TT'], C_l_len_lmat_dic['EE'], C_l_len_lmat_dic['BB'], C_l_len_lmat_dic['TE'] = C_l_len_lmat
			else:
				C_l_unl_lmat_dic['TT'] = C_l_unl_lmat[0]
				C_l_len_lmat_dic['TT'] = C_l_len_lmat[0]

			N_beam_lmat_dic['TT'] = N_beam_lmat_dic['EE'] = N_beam_lmat_dic['BB'] = N_beam_lmat_dic['TE'] = N_beam_lmat_1
			#N_beam_lmat_dic['TT'] = N_beam_lmat_2

			weight_gradient_lmat, weight_inv_var_lmat = qe_temp_pol.fn_get_filters(l2d, C_l_len_lmat_dic, C_l_unl_lmat_dic, N_beam_lmat_dic, estimator, l_G, lmax_inv_var)

			#imshow(np.fft.fftshift(weight_inv_var_lmat), extent = [np.min(lx), np.max(lx), np.min(ly), np.max(ly)]);colorbar();show();quit()

			norm_lmat = qe_temp_pol.get_norm(lx, ly, C_l_unl_lmat_dic, weight_gradient_lmat, weight_inv_var_lmat, dxdy, estimator, mapparams)

			#imshow(np.fft.fftshift(norm_lmat.real), extent = [np.min(lx), np.max(lx), np.min(ly), np.max(ly)]);colorbar();show();quit()			

			kappa_qe = qe_temp_pol.get_kappa(lx, ly, beam_lmat_transfer_lmat, weight_gradient_lmat, weight_inv_var_lmat, OBSMAP, norm_lmat, dxdy, ngaus, estimator)

			if testing: print kappa_qe, 'Modified code'

		##if testing: 
		#imshow(kappa_qe.real);colorbar();show();quit()
		'''
		from IPython import embed; embed()
		cls_kappa = self.fn_plot_pow_spec(mapparams,[kappa_qe.real])[0][0]
		cls_tf = self.fn_plot_pow_spec(mapparams,[to_deconvolve_lmat.real])[0][0]
		'''

		'''
		subplot(221);imshow(np.fft.fftshift(weight_inv_var_lmat));colorbar();title(try_mod_code)
		subplot(222);imshow(np.fft.fftshift(weight_gradient_lmat));colorbar()
		subplot(223);imshow(kappa_qe.real);colorbar()
		subplot(224);imshow(np.fft.fftshift(norm_lmat.real));colorbar();show();quit()
		'''
		'''
		cenx, ceny = kappa_qe.real.shape
		#subplot(121);imshow(OBSMAP*1e6);colorbar();grid(1,ls='solid');title('CMB lensed')
		subplot(111);imshow(kappa_qe.real);colorbar();grid(1,ls='solid');title('Kappa_QE');
		plot(cenx/2., ceny/2., 'ko', ms = 20, mec = 'k', color = 'None', mew = 2.)
		show();quit()

		#dump these guys for Eric B.
		pickle.dump(OBSMAP, gzip.open('P1c_lensed_cmb.pkl.gz', 'wb'), protocol = 2)
		pickle.dump(kappa_qe.real, gzip.open('P1c_kappa_qe.pkl.gz', 'wb'), protocol = 2)
		subplot(121);imshow(OBSMAP*1e6);colorbar();title('CMB lensed')
		subplot(122);imshow(kappa_qe.real, vmin = -5., vmax = 5.);colorbar();title('Kappa_QE');show();quit()

		subplot(121);imshow(OBSMAP[52:76,52:76]*1e6, vmin = -48., vmax = 1.);colorbar();title('CMB lensed')
		subplot(122);imshow(kappa_qe.real[52:76,52:76], vmin = -1., vmax = 1.);colorbar();title('Kappa_QE');show();quit()
		'''

		return kappa_qe
	"""
	def tsz_amp_from_input(self,mapparams, OBSMAP):
		
		tsz_map = self.beta_model_150ghz
		# smooth it with beam
		tsz_map = np.fft.ifft2(np.fft.fft2(self.beta_model_150ghz)*self.Bl).real

		RADIUS = np.sqrt(self.RA**2 + self.DEC**2.)
		RADIUS *= 60.
		amp = np.mean(tsz_map[RADIUS<=0.5])
		def twoD_Gaussian((x, y), amplitude):
			sigma_x = sigma_y =2.0
			theta =offset = 0
			xo = float(len(x)/2)
			yo = float(len(y)/2)
			a = (np.cos(theta)**2)/(2*sigma_x**2) + (np.sin(theta)**2)/(2*sigma_y**2)
			b = -(np.sin(2*theta))/(4*sigma_x**2) + (np.sin(2*theta))/(4*sigma_y**2)
			c = (np.sin(theta)**2)/(2*sigma_x**2) + (np.cos(theta)**2)/(2*sigma_y**2)
			g = offset + amplitude*np.exp( - (a*((x-xo)**2) + 2*b*(x-xo)*(y-yo) + c*((y-yo)**2)))
			return g.ravel()
		nx,ny = OBSMAP.shape
		X,Y = np.arange(nx), np.arange(ny)
		X,Y = np.meshgrid(X, Y)
		data_fitted = twoD_Gaussian((X, Y),amp)
		data_fitted = data_fitted.reshape(nx,ny)	
		

		return data_fitted
	"""
	def twoD_Gaussian(self, (x, y), amplitude,xo, yo, sigma_x, sigma_y, offset):
		"""
		convert fwhm to sigma
		"""

		a = 1./(2*sigma_x**2) 
		c = 1./(2*sigma_y**2)
		
		g = offset + amplitude*np.exp( - (a*((x-xo)**2)  + c*((y-yo)**2)))
		
		return g.ravel()




	def filtered_cmb(self, mapparams, OBSMAP):
		"""
		Filter the 90GHz CMB map instead of 150 as tSZ is higher in 90GHz maps

		"""
		RADIUS = np.sqrt(self.RA**2 + self.DEC**2.)
		RADIUS *= 60.
		lx, ly = self.get_lxly(mapparams)
		l2d = np.sqrt(lx**2. + ly**2.)
		below_l_g_cut = np.where(l2d < 2000)
		OBSMAP_fft = np.fft.fft2(self.cmb_90[0])
		OBSMAP_fft[below_l_g_cut] = 0
		recnst_obsmap = np.fft.ifft2(OBSMAP_fft).real
		

		tsz_map = recnst_obsmap
		yy = self.cmb_90[0] - np.fft.ifft2(np.fft.fft2(self.beta_model_90ghz)*self.Bl).real*1e6
		yy_fft = np.fft.fft2(yy)
		yy_fft[below_l_g_cut] = 0
		yy = np.fft.ifft2(yy_fft).real
		self.yy = yy
		nx,ny = OBSMAP.shape
		
		fit_size = self.inidic['fit_size'] #arcmin
		i,j = nx/2 - fit_size, nx/2 + fit_size
		map_fit = tsz_map[i:j,i:j]
		amp = np.mean(tsz_map[RADIUS<=1.])
		x0, y0 = fit_size, fit_size
		x,y = np.arange(2*fit_size), np.arange(2*fit_size)
		x, y = np.meshgrid(x, y)
		X,Y = np.arange(nx), np.arange(ny)
		theta, offset  = 0,0
		bigger_box = 1
		import scipy.optimize as opt
		offset = 0
		
		X,Y = np.arange(nx), np.arange(ny)
		X,Y = np.meshgrid(X, Y)
		
		fwhm_x,fwhm_y = 2,2 # depends on mass  and is set to 1 arcminute
		sigma_x = fwhm_x / np.sqrt(8. * np.log(2.))
		sigma_y = fwhm_y / np.sqrt(8. * np.log(2.))
		x0, y0 = fit_size-1, fit_size-1
		initial_guess = amp, x0,y0, sigma_x,sigma_y, offset
		cen_x0,cen_y0 = nx/2 - (fit_size - x0), ny/2 - (fit_size - y0) 
		# bounds on parameterse
		# the bounds and initial guesses are in uK

		epsilon = 0.00001
		# put bounds for 90GHz map

		
		#param_bounds = ([-400,fit_size -1-epsilon,fit_size -1-epsilon, sigma_x,sigma_y, -10],[0,fit_size -1+epsilon,fit_size-1+epsilon, 2,2, 10]) # limits on sigma_x and sigma_y should be base on mass and redshift range of clusters 
		param_bounds = ([-400,fit_size -1,fit_size -1, sigma_x,sigma_y, -10],[0,fit_size +1,fit_size+1, 2,2, 10]) 
		try:
			popt, pcov = opt.curve_fit(self.twoD_Gaussian, (x, y), map_fit.flatten(), p0=initial_guess, bounds =param_bounds )
		except:
			p0 =initial_guess
			popt = p0
			self.no_times_sigma_blown = self.no_times_sigma_blown + 1
		amp_fit, x0_fit,y0_fit, sigma_x_fit,sigma_y_fit,offset_fit  = popt	
		cen_x_fit,cen_y_fit = nx/2 - (fit_size - x0_fit), ny/2 - (fit_size - y0_fit)
		data_fitted = self.twoD_Gaussian((X, Y),amp_fit, cen_x_fit,cen_y_fit, sigma_x_fit,sigma_y_fit,0)

		

		full_fit = data_fitted.reshape(nx,ny)

		#input_sz = np.fft.ifft2(np.fft.fft2(self.beta_model_150ghz)*self.Bl).real
		full_fit = full_fit/(1.67)
		full_fit = full_fit/(1e6)
		self.sigma_x.append(sigma_x_fit)
		self.sigma_y.append(sigma_y_fit)
		self.offset_x.append(cen_x_fit)
		self.offset_y.append(cen_y_fit)
		self.amplitude.append(amp_fit/1.67)
		self.amp_offset.append(offset_fit)
		return cen_x_fit,cen_y_fit,full_fit		
		
	def robust_gaussian_fit(self,mapparams, OBSMAP):
		
		
		RADIUS = np.sqrt(self.RA**2 + self.DEC**2.)
		RADIUS *= 60.
		theta_c = 1.0
		
		nx,ny = OBSMAP.shape
		cmb_90 = self.cmb_90[0]/1e6
		tsz_map = (cmb_90 - OBSMAP)*1e6 # converting to uK
		self.sz_map.append(tsz_map)
		fit_size = self.inidic['fit_size'] #arcmin
		i,j = nx/2 - fit_size, nx/2 + fit_size
		map_fit = tsz_map[i:j,i:j]
		amp = np.mean(tsz_map[RADIUS<=1.])		
		if self.inidic['filter_tsz_map']:
			lx, ly = self.get_lxly(mapparams)
			l2d = np.sqrt(lx**2. + ly**2.)
			#below_l_g_cut = np.where(l2d < 1000)
			above_l_g_cut = np.where(l2d >14000)
			OBSMAP_fft = np.fft.fft2(tsz_map)
			#OBSMAP_fft[below_l_g_cut] = 0
			OBSMAP_fft[above_l_g_cut] = 0 
			recnst_obsmap = np.fft.ifft2(OBSMAP_fft).real
			map_fit = recnst_obsmap[i:j,i:j]
			amp = np.mean(recnst_obsmap[RADIUS<=1.])
	
		self.amp_initial.append(amp)
		x0, y0 = fit_size, fit_size
		x,y = np.arange(2*fit_size), np.arange(2*fit_size)
		x, y = np.meshgrid(x, y)
		X,Y = np.arange(nx), np.arange(ny)
		X,Y = np.meshgrid(X, Y)
		
		fwhm_x,fwhm_y = 2,2 # depends on mass  and is set to 1 arcminute
		sigma_x = fwhm_x / np.sqrt(8. * np.log(2.))
		sigma_y = fwhm_y / np.sqrt(8. * np.log(2.))
		x0, y0 = fit_size, fit_size
		cen_x0,cen_y0 = nx/2 - (fit_size - x0), ny/2 - (fit_size - y0) 
		
		# bounds on parameterse
		# the bounds and initial guesses are in uK
		epsilon = 0.00001
		offset  = 0
		initial_guess = amp, x0,y0, sigma_x,sigma_y,offset
		
		import scipy.optimize as opt
		#param_bounds = ([-400,fit_size -1-epsilon,fit_size -1-epsilon, sigma_x,sigma_y, -10],[0,fit_size -1+epsilon,fit_size-1+epsilon, 2,2, 10])

		param_bounds = ([-400,fit_size -1,fit_size -1, sigma_x,sigma_y, -10],[0,fit_size +1,fit_size+1, 2,2, 10]) 
		

		#param_bounds = ([amp-epsilon,fit_size -1,fit_size -1, sigma_x-epsilon,sigma_y-epsilon, -10],[amp-epsilon,fit_size +1,fit_size+1, sigma_x +epsilon,sigma_y +epsilon, 10])
		#param_bounds = ([amp - epsilon,fit_size -1,fit_size -1, sigma_x-epsilon,sigma_y-epsilon, -10],[amp + epsilon,fit_size +1,fit_size+1, sigma_x +epsilon,sigma_y +epsilon, 10]) 
		try:
		
			popt, pcov = opt.curve_fit(self.twoD_Gaussian, (x, y), map_fit.flatten(), p0=initial_guess, bounds =param_bounds )	
		except:
			p0 =initial_guess
			popt = p0
			self.no_times_sigma_blown = self.no_times_sigma_blown + 1

		amp_fit, x0_fit,y0_fit, sigma_x_fit,sigma_y_fit, offset_fit  = popt	
		cen_x_fit,cen_y_fit = nx/2 - (fit_size - x0_fit), ny/2 - (fit_size - y0_fit)
		data_fitted = self.twoD_Gaussian((X, Y),amp_fit, cen_x_fit,cen_y_fit, sigma_x_fit,sigma_y_fit,0)
		full_fit = data_fitted.reshape(nx,ny)
		full_fit = full_fit*1.48
		#if self.inidic['from_cmb_amp'] ==1:
			
			
		#	amp_fit, x0_fit,y0_fit, sigma_x_fit,sigma_y_fit, offset_fit  = popt
		#	amp_fit = -25.0
		#	sigma_x_fit = sigma_x
		#	sigma_y_fit = sigma_y
		#	cen_x_fit,cen_y_fit = nx/2 - (fit_size-x0_fit), ny/2 - (fit_size- y0_fit)
		#	data_fitted = self.twoD_Gaussian((X, Y),amp_fit, cen_x_fit,cen_y_fit, sigma_x_fit,sigma_y_fit,0)
			

		full_fit = full_fit/(1e6)

#		self.cmb_150ghz.append(OBSMAP*1e6)
#		self.cmb_90ghz.append(self.cmb_90[0])
#		self.cmb_150ghz_after_fit.append((OBSMAP*1e6 - full_fit*1e6))
#		self.fits.append(full_fit)
		self.sigma_x.append(sigma_x_fit)
		self.sigma_y.append(sigma_y_fit)
		self.offset_x.append(cen_x_fit)
		self.offset_y.append(cen_y_fit)
		self.amplitude.append(amp_fit*1.48)
		self.amp_offset.append(offset_fit)
		self.amp_sz.append(full_fit.min())
		return cen_x_fit,cen_y_fit,full_fit		


	def fit_ds_map(self, mapparams, tsz_map, OBSMAP):
		nx,ny,dx,dy = mapparams
		nx,ny = nx/2, ny/2
		dx,dy = dx*2, dy *2
		from IPython import embed;embed()
		tsz_map = self.downsample_map(tsz_map)
		fit_size = self.inidic['fit_size'] #arcmin
		i,j = nx/2 - fit_size, nx/2 + fit_size
		map_fit = tsz_map[i:j,i:j]

		amp = np.mean(tsz_map[48:52,48:52])
		x0, y0 = fit_size, fit_size
		x,y = np.arange(2*fit_size), np.arange(2*fit_size)
		x, y = np.meshgrid(x, y)
		X,Y = np.arange(nx), np.arange(ny)

		X,Y = np.arange(nx), np.arange(ny)
		bigger_box = 1
		import scipy.optimize as opt
		
		
		X,Y = np.arange(2*nx), np.arange(2*ny)
		X,Y = np.meshgrid(X, Y)
		
		fwhm_x,fwhm_y = 1,1 # depends on mass  and is set to 1 arcminute
		sigma_x = fwhm_x / np.sqrt(8. * np.log(2.))
		sigma_y = fwhm_y / np.sqrt(8. * np.log(2.))
		#x0, y0 = fit_size, fit_size
		x0, y0 = fit_size-1, fit_size-1 # fixing center
		cen_x0,cen_y0 = nx/2 - (fit_size - x0), ny/2 - (fit_size - y0) 
		
		# bounds on parameterse
		# the bounds and initial guesses are in uK
		epsilon = 0.00001
		offset  = 0
		initial_guess = amp, x0,y0, sigma_x,sigma_y,offset
		

		#param_bounds = ([-400,fit_size -1-epsilon,fit_size -1-epsilon, sigma_x,sigma_y, -10],[0,fit_size -1+epsilon,fit_size-1+epsilon, 2,2, 10])
		param_bounds = ([-400,fit_size -1,fit_size -1, sigma_x-epsilon,sigma_y-epsilon, -10],[0,fit_size +1,fit_size+1, sigma_x +epsilon,sigma_y +epsilon, 10]) 
		try:
		
			popt, pcov = opt.curve_fit(self.twoD_Gaussian, (x, y), map_fit.flatten(), p0=initial_guess, bounds =param_bounds )	
		except:
			p0 =initial_guess
			popt = p0
			self.no_times_sigma_blown = self.no_times_sigma_blown + 1
		amp_fit, x0_fit,y0_fit, sigma_x_fit,sigma_y_fit, offset_fit  = popt	
		cen_x_fit,cen_y_fit = nx/2 - (fit_size - x0_fit), ny/2 - (fit_size - y0_fit)
		data_fitted = self.twoD_Gaussian((X, Y),amp_fit, cen_x_fit,cen_y_fit, 2*sigma_x_fit,2*sigma_y_fit,0)
		
		
		full_fit = data_fitted.reshape(2*nx,2*ny)
		full_fit = full_fit*1.48

		full_fit = full_fit/(1e6)
		input_sz = np.fft.ifft2(np.fft.fft2(self.beta_model_150ghz)*self.Bl).real
		

		self.cmb_150ghz.append(OBSMAP*1e6)
		self.cmb_90ghz.append(self.cmb_90[0])
		self.cmb_150ghz_after_fit.append((OBSMAP*1e6 - full_fit*1e6))
		self.fits.append(full_fit)
		self.sigma_x.append(sigma_x_fit)
		self.sigma_y.append(sigma_y_fit)
		self.offset_x.append(cen_x_fit)
		self.offset_y.append(cen_y_fit)
		self.amplitude.append(amp_fit*1.48)
		self.amp_offset.append(offset_fit)
		return cen_x_fit,cen_y_fit,full_fit		



	def fit_gaussian_template(self,mapparams, OBSMAP):
		
		
		RADIUS = np.sqrt(self.RA**2 + self.DEC**2.)
		RADIUS *= 60.
		theta_c = 1.0
		

		'''
		beta_value = -1.
		tszmodel = (1.0 + (RADIUS/theta_c) ** 2. ) ** beta_value
		'''


		nx,ny = OBSMAP.shape
		cmb_90 = self.cmb_90[0]/1e6
		tsz_map = cmb_90 - OBSMAP
		#return tsz_map
		fit_size = self.inidic['fit_size'] #arcmin
		i,j = nx/2 - fit_size, nx/2 + fit_size
		map_fit = tsz_map[i:j,i:j]
		amp = np.mean(tsz_map[RADIUS<=2.])
		x0, y0 = fit_size, fit_size
		x,y = np.arange(2*fit_size), np.arange(2*fit_size)
		x, y = np.meshgrid(x, y)
		X,Y = np.arange(nx), np.arange(ny)
		theta, offset  = 0,0
		
		import scipy.optimize as opt
		
		
		X,Y = np.arange(nx), np.arange(ny)
		X,Y = np.meshgrid(X, Y)
		if self.inidic['cutout_size'] == 'bigger':
			bigger_box = 1
		if self.inidic['cutout_size'] == 'smaller':
			bigger_box = 0

		if 1==1:
			#initial guesses... all these are in pixel sizes
			sigma_x,sigma_y = 4,4
			x0, y0 = fit_size, fit_size
			initial_guess = amp, x0,y0, sigma_x,sigma_y,theta, offset
			cen_x0,cen_y0 = nx/2 - (fit_size - x0), ny/2 - (fit_size - y0) 
			param_bounds = ([-20e-4,fit_size -4,fit_size -4, 0,0,-np.pi, -np.inf],[-20e-7,fit_size +4,fit_size+4, 4,4,np.pi, np.inf])
			
			try:
				if self.inidic['bounded_params']:
					popt, pcov = opt.curve_fit(self.twoD_Gaussian, (x, y), map_fit.flatten(), p0=initial_guess, bounds =param_bounds )
					#popt, pcov = opt.curve_fit(self.twoD_Gaussian, (x, y), map_fit.flatten(), p0=initial_guess, bounds = (( -20e-4,-20e-7),(fit_size -4, fit_size +4),(fit_size -4, fit_size +4),(-np.inf,np.inf),(0,4),(-np.inf, np.inf),(-np.inf, np.inf)))
				else:
					popt, pcov = opt.curve_fit(self.twoD_Gaussian, (x, y), map_fit.flatten(), p0=initial_guess)
			except:
				p0 =initial_guess
				popt = p0
				#if self.inidic['throw_away_cuts'] == 1.:
				#	popt =0, x0,y0, sigma_x,sigma_y,theta, 0
				self.no_times_sigma_blown = self.no_times_sigma_blown + 1
			
			
			if self.inidic['amp_frm_input_tsz'] == 1.:
				
				sigma_y,sigma_x = 2,2
				epsilon = 0.00001
				initial_guess = amp, x0,y0, sigma_x,sigma_y,theta, offset
				param_bounds = ([-20e-4,fit_size-epsilon,fit_size-epsilon, 0,0,-np.pi, -np.inf],[-20e-7,fit_size+epsilon,fit_size+epsilon, sigma_x+epsilon,sigma_y+epsilon,np.pi, np.inf])
				popt, pcov = opt.curve_fit(self.twoD_Gaussian, (x, y), map_fit.flatten(), p0=initial_guess, bounds =param_bounds )
			
			amp_fit, x0_fit,y0_fit, sigma_x_fit,sigma_y_fit,theta_fit, offset_fit  = popt	
			
			#data_fitted1 = self.twoD_Gaussian((x, y),amp, x0,y0, sigma_x,sigma_y,theta,0)
			cen_x_fit,cen_y_fit = nx/2 - (fit_size - x0_fit), ny/2 - (fit_size - y0_fit)
			if bigger_box:
				data_fitted = self.twoD_Gaussian((X, Y),amp_fit, cen_x_fit,cen_y_fit, sigma_x_fit,sigma_y_fit,theta_fit,0)
			else:
				data_fitted = self.twoD_Gaussian((x, y),amp_fit, x0_fit,y0_fit, sigma_x_fit,sigma_y_fit,theta_fit,0)
			"""
			if sigma_y > 4 or sigma_x >4: # if the fitting is blow.. similar thing for amplitude???
				
				sigma_y, sigma_x = 2., 2.
				self.no_times_sigma_blown = self.no_times_sigma_blown + 1
				

				if self.inidic['throw_away_cuts'] == 1:
					data_fitted = data_fitted*0
				elif not self.inidic['throw_away_cuts']  and bigger_box:
					data_fitted = self.twoD_Gaussian((X, Y),amp, cen_x0,cen_y0, sigma_x,sigma_y,theta,0)
				elif not self.inidic['throw_away_cuts']  and not bigger_box:
					data_fitted = self.twoD_Gaussian((x, y),amp, x0,y0, sigma_x,sigma_y,theta,0)
				
				#from Idata_fittedPython import embed;embed()
		
			elif sigma_x >4.:
				
				sigma_x = 2.
				self.no_times_sigma_blown = self.no_times_sigma_blown + 1
				#data_fitted = self.twoD_Gaussian((X, Y),amp, x0,y0, sigma_x,sigma_y,theta,0)
				if self.inidic['throw_away_cuts'] == 1.:
					data_fitted = data_fitted*0
				elif not self.inidic['throw_away_cuts']  and bigger_box:
					data_fitted = self.twoD_Gaussian((X, Y),amp, cen_x0,cen_y0, sigma_x,sigma_y,theta,0)
				elif not self.inidic['throw_away_cuts']  and not bigger_box:
					data_fitted = self.twoD_Gaussian((x, y),amp, x0,y0, sigma_x,sigma_y,theta,0)

			elif sigma_y>4.:
				
				sigma_y = 2.
				self.no_times_sigma_blown = self.no_times_sigma_blown + 1
				if self.inidic['throw_away_cuts'] == 1.:
					data_fitted = data_fitted*0
				elif not self.inidic['throw_away_cuts']  and bigger_box:
					data_fitted = self.twoD_Gaussian((X, Y),amp, cen_x0,cen_y0, sigma_x,sigma_y,theta,0)
				elif not self.inidic['throw_away_cuts']  and not bigger_box:
					data_fitted = self.twoD_Gaussian((x, y),amp, x0,y0, sigma_x,sigma_y,theta,0)

			"""
			#data_fitted1 = self.twoD_Gaussian((x, y),amp, x0,x0, sigma_x,sigma_y,theta, offset)
			#cen_x,cen_y = nx/2 - (fit_size - x0), ny/2 - (fit_size - y0)
			
			#data_fitted = self.twoD_Gaussian((X, Y),amp, cen_x,cen_y, sigma_x,sigma_y,theta,0)
			#fit = data_fitted.reshape(2*fit_size,2*fit_size)


			#data_fitted = twoD_Gaussian((X, Y),amp,theta, offset)
			if bigger_box:
				full_fit = data_fitted.reshape(nx,ny)
			else:
				fit = data_fitted.reshape(2*fit_size,2*fit_size)
				full_fit = np.zeros(OBSMAP.shape) 
				full_fit[i:j,i:j] = fit 
		
		#input_sz = np.fft.ifft2(np.fft.fft2(self.beta_model_150ghz)*self.Bl).real
		full_fit = full_fit*1.48
		self.sigma_x.append(sigma_x_fit)
		self.sigma_y.append(sigma_y_fit)
		self.offset_x.append(cen_x_fit)
		self.offset_y.append(cen_y_fit)
		self.theta.append(theta_fit)
		self.amplitude.append(amp_fit)
		self.amp_offset.append(offset_fit)
		return full_fit


	def fn_get_noise(self, mapparams, expnoiselevel = None):

		nx, ny, dx, dy = mapparams
		dx *= self.arcmins2radians
		dy *= self.arcmins2radians
		npixels = nx * ny

		if not self.is_seq(expnoiselevel):
			expnoiselevel = self.expnoiselevel
		else:
			expnoiselevel = np.asarray(expnoiselevel)
			self.tqulen = len(expnoiselevel)

		DeltaT = expnoiselevel * (np.radians(1/60.)) #in uK now
		if self.inidic['perform_extra_mul']:
			extra_mul = self.inidic['extra_mul']
			DeltaT = DeltaT * extra_mul

		NOISE = np.zeros( (self.tqulen,nx,ny) )
		for nn in range(self.tqulen): 

			DUMMY = np.random.standard_normal([nx,ny])
			
			NOISE[nn] = DUMMY * DeltaT[nn]/dx
			#NOISE[nn] = np.random.normal(loc=0.0, scale=DeltaT[nn]/dx, size=[nx,ny])

		#imshow(NOISE[0]);colorbar();show()
			#show();quit()

		return NOISE


	def fn_noise_cov(self,mapparams, noofsims, reqdbox=None, expnoiselevel = None, simple_cov = 0): #CHECK
		
		tqulen = self.tqulen
		nx, ny, dx, dy = mapparams

		NOISE = np.zeros( (tqulen, noofsims, ny, nx) )
		np.random.seed(self.cmb_noise_randomseedval)
			
		#np.random.seed(self.noise_random_seed)
		for cnt in range(noofsims):

			TMP = self.fn_get_noise(mapparams,expnoiselevel = expnoiselevel)
			#TMP = self.fn_get_noise_from_cl_noise(mapparams,expnoiselevel)
			#imshow(TMP[0]);colorbar();show();quit()

			#### 20161014 - noise must be convolved with TWODTF to produce the correct NOISE_COV
			#TMP = np.fft.ifft2( np.fft.fft2(TMP) * self.TWODTF ).real

			NOISE[:, cnt, :, :] = TMP

		npixels = nx * ny

		if simple_cov:
			COV_N = np.zeros( (tqulen, npixels, npixels) )

			for ii in range(tqulen):
				COV_N[ii, :, :] = self.calcCov(NOISE[ii], noofsims, npixels)
				#imshow(COV_N[ii, :, :]);colorbar();show();quit()

			return COV_N

		locdic = self.locdic
		identifiers = self.inidic['cov_identifiers']
		COV_N = {}

		for iden in identifiers:
			if len(iden) == 1:
				loc = locdic[iden]
				MAPFORCOV = NOISE[loc]
				COV = self.calcCov(MAPFORCOV, noofsims, npixels)
			elif len(iden) == 2:
				loc1,loc2 = locdic[iden[0]],locdic[iden[1]]
				MAPFORCOV1, MAPFORCOV2 = NOISE[loc1], NOISE[loc2]
				MAPFORCOV = np.concatenate( (MAPFORCOV1,MAPFORCOV2), axis = 1 )
				COV = np.asarray(self.calcCov(MAPFORCOV, noofsims, 2*npixels))

			COV_N[iden] = COV

		#imshow(COV_N[iden],origin='lower');colorbar();show();quit()
	
		return COV_N

	def el1d_to_EL2D(self,el1d,mapparams, interp_type = 'cubic'):

		assert  interp_type == 'cubic' or  interp_type == 'linear'

		nx, ny, dx, dy = mapparams
		dx *= self.arcmins2radians
		dy *= self.arcmins2radians

		lx, ly = self.get_lxly(mapparams)
		L = np.sqrt(lx**2. + ly**2.)

		EL2D = np.interp(L.flatten(), el1d[0], el1d[1], right = 0.).reshape(L.shape)

		#fninterp = interp1d(el1d[0], el1d[1], kind = interp_type, bounds_error = 0, fill_value = 0.)
		#EL2D = fninterp(L.flatten()).reshape(L.shape)

		return EL2D

	def fn_read_SPTpol_el_space_beam_v1(self, mapparams, fpath = 'data/sptpol_beam_20160515.pkl'):

		import gzip, pickle, numpy as np
		import scipy.interpolate as intrp

		#SPTpol el space beam
		beamdic = pickle.load(open(fpath,'r'))
		beam, beamerr = beamdic['B_ell_hybrid'], beamdic['B_ell_hybrid_error']

		#plot(beam[0],beam[1]);show();quit()

		Bl = self.el1d_to_EL2D(beam,mapparams)
		#clf();imshow(np.fft.fftshift(Bl));colorbar();show();quit()
		
		return Bl

	def fn_interp_beam(self, els_binned, bls_binned, interp_type = 'cubic', norm_el = 800.5):
		#remove nans
		inds = np.where( np.isnan(bls_binned) == 0)[0]
		els_binned, bls_binned = els_binned[inds], bls_binned[inds]


		fninterp = interp1d(els_binned, bls_binned, kind = interp_type, bounds_error = 0, fill_value = 0.)
		els = np.arange( min(els_binned), max(els_binned))
		bls = fninterp(els)

		if norm_el<>None:
			norm_ind = np.where( els == norm_el )[0]
			norm_val = bls[norm_ind]
			bls /= norm_val

		return els, bls		

	#def fn_read_SPTpol_el_space_beam(self, mapparams, fpath = 'data/500sqdeg_henning/bb500d_beams_temp_edited.pkl', nu = 150):
	def fn_read_SPTpol_el_space_beam(self, mapparams, fpath = 'data/sptpol/20171107/500d_venus_beams_for_Christian_extracted.pkl.gz', nu = 150, return_1d = 0, beam_error = 0):

		#SPTpol el space beam
		if fpath == 'data/500sqdeg_henning/bb500d_beams_temp_edited.pkl':
			beamdic = pickle.load(open(fpath,'r'))[nu]
			beam, beamerr = beamdic['B_ell_hybrid'], beamdic['B_ell_hybrid_error']
		elif fpath == 'data/sptpol/20171107/500d_venus_beams_for_Christian_extracted.pkl.gz':
			#20180319 - incude beam error
			if abs(beam_error)>0:
				fpath = 'data/sptpol/20171107/500d_venus_beams_for_Christian_extracted_with_errors.pkl.gz'

			beamdic = pickle.load(gzip.open(fpath,'r'))
			Bl_els = beamdic['els']
			if nu == 90:
				### beam = np.asarray( [Bl_els, beamdic['Bl_90']] )
				beam_val = beamdic['Bl_90']**0.5
				if abs(beam_error)>0:
					beam_1sigma = beamdic['cov_B_ell_90']**0.5

			elif nu == 150:
				### beam = np.asarray( [Bl_els, beamdic['Bl_150']] )
				beam_val = beamdic['Bl_150']**0.5
				if abs(beam_error)>0:
					beam_1sigma = beamdic['cov_B_ell_150']**0.5

			if abs(beam_error)>0:
				beam_val = beam_val + (beam_1sigma * beam_error)
				beam_val /= max(beam_val) #normalise this

			#normalise the beam at ell = 800 consistent with SPTpol beam #added on 20180809
			try:
				norm_beam_el = self.inidic['norm_beam_el']
			except:
				norm_beam_el = None
			if norm_beam_el<>None:
				Bl_els, beam_val = self.fn_interp_beam(Bl_els, beam_val, interp_type = 'cubic', norm_el = norm_beam_el)

			'''
			print nu, norm_beam_el, beam_val[0]
			plot(Bl_els, beam_val);show();sys.exit()
			'''

			beam = np.asarray( [Bl_els, beam_val] ) #20180123 - JT provides Bl**2.

		if return_1d:
			return beam

		'''
		plot(beam[0],beam[1])
		if beam_error>0:
			clf()
			errval = beam_1sigma * beam_error
			dummy = beamdic['Bl_150']**0.5
			dummy1 = (dummy - errval)
			dummy2 = (dummy + errval)
			semilogy(beam[0],dummy/max(dummy), 'k--');
			semilogy(beam[0],dummy1/max(dummy1));
			semilogy(beam[0],dummy2/max(dummy2));
		show();quit()
		'''
		#print beam.shape
		#plot(beam[0],beam[1]);show()

		Bl = self.el1d_to_EL2D(beam,mapparams)

		Bl[Bl == np.inf] = 0.
		Bl[np.where(np.isnan(Bl))] = 0.

		return Bl


	def fn_read_ACTPol_el_space_beam(self, mapparams, fpath = 'data/ACTPol/ACTPol_consolidated_beam.pkl.gz', keyname = ('2014','pa1','deep56')):

		#https://lambda.gsfc.nasa.gov/product/act/actpol_beams_info.cfm
		possible_keynames = [('2013','pa1','deep5'), ('2013','pa1','deep6'), ('2014','pa1','deep56'), ('2014','pa2','deep56')]
		assert keyname in possible_keynames

		import gzip, pickle, numpy as np
		import scipy.interpolate as intrp

		#ACTPol el space beam
		actpol_beam = pickle.load(gzip.open(fpath,'rb'))[keyname]
		#els, Bls = actpol_beam
		#plot(els, Bls);show();quit()

		Bl = self.el1d_to_EL2D(actpol_beam,mapparams)
		#clf();imshow(np.fft.fftshift(Bl));colorbar();show();quit()
		
		return Bl


	def fn_tf_noise(self, SIMMAPS, mapparams, noofsims = 10, reqdbox = None, add_TF=0, add_noise=0, expnoiselevel = None, noise_random_seed = None, cluscnt = None, nu = 150,ra_for_rot = None, dec_for_rot = None, return_full_T_map = 0, special_map = 0):

		########################################################################
		########################################################################

		#MAP RESOL, EL stuffs
		nx, ny, dx, dy = mapparams
		dx *= self.arcmins2radians
		dy *= self.arcmins2radians
		lx, ly = self.get_lxly(mapparams)
		L = np.sqrt(lx**2. + ly**2.)

		#imshow(Bl);show();quit()

		if not self.quiet:
			logline = '\t\tperform TF convolution to sims (and might add noise here for sims)'
			logfile = open(self.log_file,'a')
			logfile.writelines('%s\n' %(logline))
			logfile.close()
			print logline
		####################################################################
		####################################################################
		#transfer function stuff
		#### 
		###TWODTF = self.fn_get_TWODTF(mapparams, RA, DEC, ra_val, dec_val)
		#TWODTF = self.fn_get_spt_centre_TF(mapparams) #rotation not implemented yet
		####
		
		if not self.is_seq(self.TWODTF):
			if add_TF == 0:
				TWODTF = self.fn_get_HPF(mapparams, ideal = 1) #still analytic one
				logline = '\t\t\tgetting ideal TF - all ones'
			elif add_TF == 1:
				#analytic TF
				TWODTF = self.fn_get_HPF(mapparams, minel = self.inidic['min_el'], maxel = self.inidic['max_el']) #still analytic one
				logline = '\t\t\tanalytic TF - minel = %s; maxel = %s' %(self.inidic['min_el'], self.inidic['max_el'])
			elif add_TF == 2:
				#sptpol TF
				TWODTF = self.fn_get_SPTpol_TWODTF(mapparams)
				logline = '\t\t\tSPTpol TF'
			elif self.inidic['add_TF'] == 3:
				#transfer_lmat = self.fn_get_EBAX_2016_anal(mapparams, minel=self.inidic['min_el'], maxel=self.inidic['max_el'])[0]
				TWODTF = self.fn_get_EBAX_2016_anal(mapparams, l1=self.inidic['l1'], l2=self.inidic['l2'], l3=self.inidic['l3'])[0]
			
			try:
				if self.inidic['add_TF_grad']:
					TWODTF_GRAD = self.fn_get_HPF(mapparams, minel=self.inidic['min_el_grad'], maxel=self.inidic['max_el_grad'], all_iso = 0)[0]
					self.TWODTF_GRAD = TWODTF_GRAD
				if not self.inidic['add_TF_grad']:
					self.TWODTF_GRAD = TWODTF
			except:
				self.TWODTF_GRAD = TWODTF
				pass
			

			self.TWODTF = TWODTF

			if 0 == 1:#not self.quiet:
				logfile = open(self.log_file,'a')
				logfile.writelines('%s\n' %(logline))
				logfile.close()
				print logline

		#now perform rotation of TF if needed
		#TWODTF_ORI = np.copy(TWODTF)
		if ra_for_rot <> None and dec_for_rot <> None:
			TWODTF[0] = self.rotate_tf(TWODTF[0], ra_for_rot, dec_for_rot, in_elspace = 1)

		#subplot(121);imshow(np.fft.fftshift(TWODTF_ORI[0]), extent = [np.min(lx), np.max(lx), np.min(ly), np.max(ly)]);colorbar()
		#subplot(122);imshow(np.fft.fftshift(TWODTF[0]), extent = [np.min(lx), np.max(lx), np.min(ly), np.max(ly)]);colorbar();show();quit()

		TWODTF = self.TWODTF
		if self.is_seq(self.TWODTF_GRAD):
			TWODTF_GRAD = self.TWODTF_GRAD 

		####################################################################
		####################################################################


		self.npixels = nx * ny

		#if not reqdbox == None:
		if hasattr(reqdbox, "__len__") == 1:
			ex1, ex2, ey1, ey2 = reqdbox
			nx, ny = ex2-ex1, ey2-ey1
			self.npixels = nx * ny
			self.mapparams = [nx, ny, np.degrees(dx) * 60.,np.degrees(dy) * 60.]

		tqulen = self.tqulen

		#self.noise_random_seed = noise_random_seed
		#np.random.seed(noise_random_seed)

		#SIMS --> #([noofsims], [3 -> T,Q,U], [nx, ny])
		start = time.time()
		SIMS = np.zeros( (noofsims, tqulen, ny, nx) )
		FULL_T_MAP = []
		for cnt in range(noofsims):

			if cnt % 5000 == 0 and cnt>0:
				#logline = '\t\t\t\tsimno: %s for %s' %(cnt, typestr[ii])
				logline = '\t\t\t\tsimno: %s' %(cnt)
				logfile = open(self.log_file,'a')
				logfile.writelines('%s\n' %(logline))
				logfile.close()
				print logline

			MAP = SIMMAPS[cnt,:]

			if not special_map:
				self.expnoiselevel = expnoiselevel
				self.noise_present = 0

			if not cluscnt == None:
				rs = int( (cluscnt+1) * self.cmb_noise_randomseedval * nu)
			else:
				
				rs = int( time.time()/1e5 * self.cmb_noise_randomseedval * nu)

			while rs>1e8:
				rs = int(rs/1e3)

			if noise_random_seed == 0: #then do not set it
				rs = None
			#print rs, self.cmb_noise_randomseedval, nu
			
			MAP_BEFORE_NOISE = np.copy(MAP)
			if add_noise == 1:
				'''
				logline = '\t\tAdding white noise now; Level = %s' %(expnoiselevel)
				logfile = open(self.log_file,'a')
				logfile.writelines('%s\n' %(logline))
				logfile.close()
				#print logline
				'''

				self.noise_present = 1
				#this makes sure that the result is not affected by number of clsuters used - same noise for a particular cluster
				#not that the previous result was bad

				'''#20171201 - moving this up and making it right for all frequencies and cluscnt
				if not cluscnt == None:
					rs = cluscnt * nu
				else:
					rs = int((cnt+1) * time.time()/1e4)

				while rs>1e8:
					rs = int(rs/1e5)
				'''
				np.random.seed(rs)#;
				#from IPython import embed;embed()
				#print 'stopping here ... ';quit()
				self.NOISE = self.fn_get_noise(mapparams, expnoiselevel) #different noise for each sim
				#self.NOISE = self.fn_get_noise_from_cl_noise(mapparams, expnoiselevel) #different noise for each sim

				'''
				if cluscnt == 0:
					if nu == 150:
						self.randomnoise_150 = np.copy(self.NOISE)[0]
					else:
						self.randomnoise_tszfree = np.copy(self.NOISE)[0]

						clf()
						subplot(131);imshow(self.randomnoise_150);colorbar()
						subplot(132);imshow(self.randomnoise_tszfree);colorbar()
						subplot(133);imshow(self.randomnoise_tszfree - self.randomnoise_150);colorbar();show();sys.exit()
				'''

				#print np.std(self.NOISE) * 0.5, np.std(MAP[0]) * 0.5, 'just map'
				MAP = MAP + self.NOISE

				#imshow(self.NOISE[0]);colorbar();show();quit()
				#MAP = np.copy(self.NOISE)
				#print np.std(self.NOISE[0]) * 0.5, np.std(MAP[0]) * 0.5, 'map = map + noise'#;quit()

			elif add_noise == 2:
				'''
				logline = '\t\tAdding SPTpol like noise now'
				logfile = open(self.log_file,'a')
				logfile.writelines('%s\n' %(logline))
				logfile.close()
				#print logline
				'''

				self.noise_present = 2
				try:
					noisesims = self.noisesims
					tSZfitting = self.tSZfitting
				except:
					noisesims = noise_sims()
					tSZfitting = tSZ_fitting()
					self.noisesims = noisesims
					self.tSZfitting = tSZfitting

				try:
					tszfree_map = self.inidic['tszfree_map']
				except:
					tszfree_map = 0

				if not special_map:
					whichmap = self.whichmap
				else:
					whichmap = self.inidic['whichmap_grad']

				if tszfree_map:
				        nu = 90
				        #whichmap = '90ghz_500sqdeg'
				        NOISE_90 = noisesims.fn_get_noise_sims_from_noise_cls(mapparams, 1, random_seed = int(time.time()), nu = nu, whichmap = whichmap, quiet = self.quiet)[0,0]

				        nu = 150
				        #whichmap = '150ghz_500sqdeg'
				        NOISE_150 = noisesims.fn_get_noise_sims_from_noise_cls(mapparams, 1, random_seed = int(time.time() / nu), nu = nu, whichmap = whichmap)[0,0]
			                use_beam_90 = use_beam_150 = self.inidic['use_beam']
			                NOISE_tSZ_free = tSZfitting.fn_subtract_tSZ(mapparams, NOISE_150, NOISE_90, use_beam_90 = use_beam_90, use_beam_150 = use_beam_150)

					self.NOISE = NOISE_tSZ_free

				else:
					self.NOISE = noisesims.fn_get_noise_sims_from_noise_cls(mapparams, 1, random_seed = rs, nu = nu, whichmap = whichmap, quiet = self.quiet)[0]
					self.NOISE = self.NOISE[0:self.tqulen]

				MAP = MAP + np.copy(self.NOISE)



			if 0==1:#self.debug:
				for tqucnt in range(self.tqulen):
					subplot(3, self.tqulen, tqucnt+1);imshow(MAP_BEFORE_NOISE[tqucnt]);colorbar()
					subplot(3, self.tqulen, self.tqulen + tqucnt+1);imshow(self.NOISE[tqucnt]);colorbar()
					subplot(3, self.tqulen, 2 * self.tqulen + tqucnt+1);imshow(MAP[tqucnt]);colorbar()
				show();quit()

			#if add_TF:
			#subplot(121);imshow(MAP[0]);colorbar()
			
			if self.inidic['add_TF']: #adding transfer function
				#MAP = self.fn_apodize(MAP, mapparams)
				if self.is_seq(self.TWODTF_GRAD) and special_map:
					MAP = np.fft.ifft2( np.fft.fft2(MAP) * self.TWODTF_GRAD ).real
				else:
					MAP = np.fft.ifft2( np.fft.fft2(MAP) * TWODTF ).real
				#for mm in range(self.tqulen):
				#	MAP[mm] = self.Convolve(MAP[mm], np.fft.ifft2(TWODTF[mm]).real)[0:len(MAP[mm]),0:len(MAP[mm])] #Check
			#subplot(122);imshow(MAP[0]);colorbar();show();quit()
			
			if 0==1:#self.debug:
				for tqucnt in range(1):#self.tqulen):
					MMM = MAP[tqucnt]
					subplot(2, self.tqulen, tqucnt+1);imshow(MMM);colorbar()
					PSD = np.fft.fftshift( abs( np.fft.fft2(MMM) ) )
					subplot(2, self.tqulen, self.tqulen+tqucnt+1);imshow(PSD,extent = [np.min(lx),np.max(lx),np.min(ly),np.max(ly)]);colorbar()
					#xlim(-500,500);ylim(-500,500.)

					#subplot(2, self.tqulen, tqucnt+2);imshow(MMM[ex1:ex2, ey1:ey2]);colorbar()
					#PSD = np.fft.fftshift( abs( np.fft.fft2(MMM[ex1:ex2, ey1:ey2]) ) )
					#subplot(2, self.tqulen, self.tqulen+tqucnt+2);imshow(PSD,extent = [np.min(lx),np.max(lx),np.min(ly),np.max(ly)]);colorbar()
				show();quit()


			#print np.std(self.NOISE[0]) * 0.5, np.std(MAP[0]) * 0.5, 'map = map + noise + tf'#;quit() #looks alright here

			#################################################################################
			#################################################################################
			#one more filtering for SPTPol noise
			if 0==1:
				logline = '\t\t\t\tperforming one more noise filtering'
				logfile = open(self.log_file,'a')
				logfile.writelines('%s\n' %(logline))
				logfile.close()
				print logline

				MAP = self.fn_apodize(MAP, mapparams)
				MAP_BEFORE_NOISE_FILTER = np.copy(MAP)
				NOISE_FILTER = self.fn_get_HPF(mapparams, minel = 500, maxel = self.inidic['max_el']) #still analytic one
				MAP = np.fft.ifft2( np.fft.fft2(MAP) * NOISE_FILTER ).real
				'''
				for tqucnt in range(self.tqulen):
					subplot(3, self.tqulen, tqucnt+1);imshow(MAP_BEFORE_NOISE_FILTER[tqucnt]);colorbar()
					subplot(3, self.tqulen, self.tqulen + tqucnt+1);imshow(MAP[tqucnt]);colorbar()
					subplot(3, self.tqulen, 2 * self.tqulen + tqucnt+1);imshow(MAP_BEFORE_NOISE_FILTER[tqucnt] - MAP[tqucnt]);colorbar()
				show();quit()
				'''

			#################################################################################
			#################################################################################

			#subplot(122);imshow(MAP[0]);colorbar();show();quit()

			if return_full_T_map == 1:
				FULL_T_MAP.append(MAP[0])
			
			if hasattr(reqdbox, "__len__") == 1:
				MAP = MAP[:, ex1:ex2, ey1:ey2]

			SIMS[cnt,:,:,:] = MAP

			'''
			for i in range(len(self.SIMMAPS[0])):
				subplot(2,5,i+1);css=imshow(self.SIMMAPS[0,i][ex1:ex2, ey1:ey2],cmap=cm.gnuplot);colorbar()
			for i in range(len(self.SIMMAPS_L[0])):
				subplot(2,5,i+5+1);css=imshow(self.SIMMAPS_L[0,i][ex1:ex2, ey1:ey2],cmap=cm.gnuplot);colorbar()
			show();quit()#looks okay here
			'''

			"""
			#to check with Fig. 1 of http://arxiv.org/abs/astro-ph/0512104v2
			MAP, MAP_L = self.SIMMAPS[0,0][ex1:ex2, ey1:ey2], self.SIMMAPS_L[0,0][ex1:ex2, ey1:ey2]
			MAP_DIFF = MAP_L - MAP
			subplot(131);css=imshow(MAP,cmap=cm.gnuplot);css.set_clim(-100,100);colorbar()
			subplot(132);css=imshow(MAP_L,cmap=cm.gnuplot);css.set_clim(-100,100);colorbar()
			subplot(133);css=imshow(MAP_DIFF,cmap=cm.gnuplot);css.set_clim(-10,10);colorbar()
			show();quit()
			"""

		end = time.time()

		if not self.quiet:
			logline = '\t\tmap cutouts extracted. time taken = %s seconds. next is COV calc.' %(end-start)
			logfile = open(self.log_file,'a')
			logfile.writelines('%s\n' %(logline))
			logfile.close()
			print logline

		if return_full_T_map:

			return SIMS, FULL_T_MAP
		else:
			return SIMS

	def fn_apodize(self, MAP, mapparams, mask = 'circle', just_return_mask = 0):

		import numpy as np, scipy.ndimage as ndimage, scipy.signal as signal

		nx, ny, dx, dy = mapparams

		start =  time.time()
		pix = dx
		radius = (nx * pix)/10. #(nx * pix)/3. #changed on 20171206
		npix_cos = int(radius/pix)
		ker=np.hanning(npix_cos)
		ker2d=np.asarray( np.sqrt(np.outer(ker,ker)) )

		MASKf=np.zeros((nx,ny))
		minval, maxval = -(nx*pix)/2,  (nx*pix)/2
		x = y = np.linspace(minval, maxval, nx)
		X, Y = np.meshgrid(x,y)
		xc, yc = 0., 0.
		if mask == 'circle':
				radius = (nx * dx/2) - (nx/50) #changed on 20171206
				inds=np.where((X-xc)**2. + (Y-yc)**2. <= radius**2.) #all in arcmins
		elif mask == 'square':
				radius = (nx * dx/2) - 1.
				inds=np.where((abs(X)<=radius) & (abs(Y)<=radius)) #all in arcmins
		MASKf[inds]=1.

		apodMASKf=ndimage.convolve(MASKf, ker2d)#, mode='wrap')
		apodMASKf/=apodMASKf.max()

		if just_return_mask: return apodMASKf

		self.apodMASKf = apodMASKf

		return apodMASKf * MAP
		
	def c_Duffy(self, M, z, h, kind='200', profile_name = 'NFW'):
		"""
		Concentration from c(M) relation published in Duffy et al. (2008).
		"""

		M_pivot = 2.e12/h # [M_sun]
		
		if profile_name == 'NFW':
			if kind == '200':
				A = 5.71
				B = -0.084
				C = -0.47
			elif kind == 'vir':
				A = 7.85 
				B = -0.081
				C = -0.71
		elif profile_name == 'Einasto':
			if kind == '200':
				A = 6.4
				B = -0.108
				C = -0.62
			elif kind == 'vir':
				A = 8.82 
				B = -0.106
				C = -0.87


		if M>0:
			concentration = A * ((M / M_pivot)**B) * (1+z)**C
		else:
			concentration = 0.

		return concentration

	def fn_select_clusters(self, kappa_qe_arr_full, CLUS_IDENTIFIER_full, rich1, rich2, z1 = None, z2 = None, weights_full = None, quiet = 1, passed_ra_dec_z = None):

		passed_inds = []
		for kcnt, cluskey in enumerate( CLUS_IDENTIFIER_full ):
			try:
				ra, dec, z_val, rich, weight, weights_norm = cluskey
			except:
				ra, dec, z_val, rich, weight = cluskey

			if passed_ra_dec_z<>None:
				ind = np.where( (passed_ra_dec_z[0] == ra) & (passed_ra_dec_z[1] == dec) & (passed_ra_dec_z[2] == z_val) )[0]
				if len(ind) == 0: continue

			passed = 0
			if rich >= rich1 and rich<rich2:
				passed = 1

			if z1<>None:
				if z_val<z1: 
					passed = 0

			if z2>None:
				if z_val>z2: 
					passed = 0

			if passed: passed_inds.append(kcnt)

		if not quiet:
			print len(passed_inds), len(kappa_qe_arr_full), len(CLUS_IDENTIFIER_full)

		passed_inds = np.asarray( passed_inds )

		if self.is_seq(weights_full):
			return kappa_qe_arr_full[passed_inds], CLUS_IDENTIFIER_full[passed_inds], weights_full[passed_inds]
		else:
			return kappa_qe_arr_full[passed_inds], CLUS_IDENTIFIER_full[passed_inds]

	def fn_simple_JK(self, raarr, noofsims):

		total = len(raarr)
		each_split_should_contain = int(total * 1./noofsims)

		fullarr = np.arange(total)
		inds_to_pick = np.arange(len(fullarr))
		already_picked_inds = []

		JK_SAMPLES = []
		for n in range(noofsims):

			inds = np.random.choice(inds_to_pick, size = each_split_should_contain, replace = 0)
			inds_to_delete = np.where (np.in1d(inds_to_pick, inds) == True)[0]
			inds_to_pick = np.delete(inds_to_pick, inds_to_delete)

			#push all on the non inds dic into - because for each JK we will ignore the files for this respective sim
			tmp = np.in1d(fullarr, inds)
			non_inds = np.where(tmp == False)[0]

			JK_SAMPLES.append( (non_inds) )

		##from IPython import embed; embed()

		return JK_SAMPLES


	def fn_healpix_JK(self, RA, DEC, noofsims, nside = 1024, dx = 0.5, testing = 0):

		import healpy as H

		'''
		#create a HP mask based on SPTpol mask
		npix = H.nside2npix(nside)
		pixel_coords = np.indices(totalmask.shape)
		map_centre = np.asarray( [0.0, -57.5] )
		map_resol = dx #arcmins
		mapshape = np.asarray( totalmask.shape )
		map_proj = 0
		sptpol_ra, sptpol_dec = sky_local.pix2Ang(pixel_coords, map_centre, map_resol, mapshape, proj = map_proj)

		sptpol_mask = np.copy(totalmask)
		sptpol_mask[sptpol_mask>0] = 1.
		HMASK=np.zeros(npix)
		SPTPOLPIXELS = []
		for r in range(np.shape(sptpol_ra)[0]):
			PP=H.pixelfunc.ang2pix(nside,np.radians(90.-sptpol_dec[r,:]),np.radians(sptpol_ra[r,:]), nest = 1)
			SPTPOLPIXELS.extend(PP)
		SPTPOLPIXELS = np.unique(SPTPOLPIXELS)
		'''

		npix = H.nside2npix(nside)
		SPTPOLPIXELS=H.pixelfunc.ang2pix(nside,np.radians(90.-DEC),np.radians(RA), nest = 1)

		groups = np.array_split(SPTPOLPIXELS, noofsims)

		JK_SAMPLES = []
		if testing: PICKED = np.zeros( npix )
		for gcnt, g in enumerate( groups ):
			non_inds = np.where( np.in1d( SPTPOLPIXELS, g ) == False )[0]
			JK_SAMPLES.append( (non_inds) )
			if testing: PICKED[g] = gcnt

		if testing: H.mollview(PICKED, nest = 1);show()

		return JK_SAMPLES

	def fn_get_COV_from_JK(self, kappa_qe_arr_full, CLUS_IDENTIFIER_full, noofsims, richness_bins, nside = 1024, maxbin = 10.0, binsize = 1.0, multiple_richness_bins = 0, weights = None, kappa_inv_var_weights = 0, passed_ra_dec_z = None, quiet = 0):

		totbins = int(maxbin/binsize)
		totrichs = len(richness_bins)

		JK_DIC = {}
		RICH_DIC = {}
		for rcnt, r1r2 in enumerate(richness_bins):

			print r1r2

			minrich, maxrich = r1r2
			'''
			if multiple_richness_bins:
				kappa_qe_arr_this_bin, CLUS_IDENTIFIER_this_bin = fn_select_clusters(kappa_qe_arr_full, CLUS_IDENTIFIER_full, minrich, maxrich)
			else:
				kappa_qe_arr_this_bin, CLUS_IDENTIFIER_this_bin = kappa_qe_arr_full, CLUS_IDENTIFIER_full
			'''
			if self.is_seq(weights):
				kappa_qe_arr_this_bin, CLUS_IDENTIFIER_this_bin, weights_this_bin = self.fn_select_clusters(kappa_qe_arr_full, CLUS_IDENTIFIER_full, minrich, maxrich, weights_full = weights, passed_ra_dec_z = passed_ra_dec_z)
			else:
				kappa_qe_arr_this_bin, CLUS_IDENTIFIER_this_bin = self.fn_select_clusters(kappa_qe_arr_full, CLUS_IDENTIFIER_full, minrich, maxrich, passed_ra_dec_z = passed_ra_dec_z)
			raarr, decarr = CLUS_IDENTIFIER_this_bin[:,0], CLUS_IDENTIFIER_this_bin[:,1]
			#make sure raarr and decarr have floats in them and not lists
			if isinstance(raarr[0], list):
				raarr = np.asarray( map(lambda x: x[0], raarr) )
				decarr = np.asarray( map(lambda x: x[0], decarr) )

			if len(np.unique(raarr)) == 1:
				JK_DIC[r1r2] = self.fn_simple_JK(raarr, noofsims)
			else:
				JK_DIC[r1r2] = self.fn_healpix_JK(raarr, decarr, noofsims, nside = nside)
			if not kappa_inv_var_weights and not self.is_seq(weights):
				RICH_DIC[r1r2] = kappa_qe_arr_this_bin
			else:
				if not self.is_seq(weights):
					kappa_qe_arr_for_std = np.copy(kappa_qe_arr_this_bin)
					kappa_qe_arr_for_std[:,80:120,80:120] = 0.
					weights_this_bin = 1./np.std(kappa_qe_arr_for_std, axis =(1,2))**2.
					RICH_DIC[r1r2] = [kappa_qe_arr_this_bin, weights_this_bin]
				else:
					RICH_DIC[r1r2] = [kappa_qe_arr_this_bin, weights_this_bin]

		'''
		#get the number of sims based on the JK that exists in all z-bins
		simarr = []
		JK_SAMPLE_DIC = JK_DIC[richness_bins[0]]
		for n in JK_SAMPLE_DIC:
			if n in JK_DIC[richness_bins[1]] and n in JK_DIC[richness_bins[2]]:
				simarr.append(n)
		'''

		## noofsims = len(simarr)
		simarr = np.arange(noofsims)
		VEC_MAT = np.zeros( (totrichs, totbins, noofsims) )

		nx, ny = kappa_qe_arr_this_bin[0].shape
		dx = dy = 0.5
		boxsize = nx * dx
		mapparams = [nx, ny, dx, dy]
		clra, cldec = 0., 0.
		minval, maxval = clra-boxsize/2/60.,  clra+boxsize/2/60.
		ra = dec = np.linspace(minval, maxval, nx)
		RA, DEC = np.meshgrid(ra,dec)
		RADEC = [RA, DEC]

		for rcnt, r1r2 in enumerate(richness_bins):

			print '\n\tBin=', r1r2
			JK_SAMPLE_DIC = JK_DIC[r1r2] #JK samples
			if not kappa_inv_var_weights:
				### from IPython import embed; embed()
				kappa_qe_arr = RICH_DIC[r1r2] #kappa_qe arr for this z-bin
				##STACKED_KAPPA_FULL = np.mean( kappa_qe_arr, axis = 0 )
			else:
				kappa_qe_arr, weights = RICH_DIC[r1r2] #kappa_qe arr for this z-bin
				##STACKED_KAPPA_FULL = np.sum( np.asarray( map(lambda x, y: x * y, kappa_qe_arr, weights) ), axis = 0) / np.sum( weights )

			fullarr = np.arange(len(kappa_qe_arr))

			for simcnt, n in enumerate( simarr ):

				non_inds = JK_SAMPLE_DIC[n]

				if not quiet:
					print '\t\tJK #%s of %s: Total inds here = %s' %(simcnt+1, noofsims, len(non_inds))

				if not kappa_inv_var_weights:
					STACKED_KAPPA = np.mean( kappa_qe_arr[non_inds], axis = 0 ) #this is faster than summing non inds
				else:
					STACKED_KAPPA = np.sum( np.asarray( map(lambda x, y: x * y, kappa_qe_arr[non_inds], weights[non_inds]) ), axis = 0) / np.sum( weights[non_inds] )

				### if self.is_seq(MEAN_FIELD): STACKED_KAPPA = STACKED_KAPPA - MEAN_FIELD
				RADPROFILES = self.fn_radial_profile(STACKED_KAPPA, RADEC, bin_size=binsize, minbin=0.0, maxbin=maxbin)[:,1]

				VEC_MAT[rcnt, :, simcnt] = RADPROFILES

		FULL_COV = np.zeros( (totrichs * totbins, totrichs * totbins) )
		prevrow = 0
		for i in range(totrichs):
			V1 = VEC_MAT[i]
			prevcol = 0
			for j in range(totrichs):
				#if i == j: continue
				V2 = VEC_MAT[j]
				#FULL_COV[]
				r1, r2 = prevrow, prevrow + totbins 
				c1, c2 = prevcol, prevcol + totbins
				FULL_COV[r1 :r2, c1:c2] = np.cov(V1, V2)[totbins:, 0:totbins]
				prevcol += totbins 
			prevrow += totbins

		factor = (noofsims - 1)
		FULL_COV = FULL_COV * factor

		return FULL_COV

	def calcCov(self, MAPMAT, noofsims, npixels, perform_mode_ref=0, perform_mean_sub = 0):
		
		start = time.time()

		"""
		M = MAPMAT.flatten().reshape(noofsims, npixels)
		#Mmean = np.mean(M, axis = 1).squeeze()
		M = np.mat( M )# - Mmean[:,None] )
		Mt = M.T
		COV = (Mt * M) / (noofsims)
		"""
		M = MAPMAT.flatten().reshape(noofsims,npixels)
		#subplot(121);pcolor(MAPMAT[0]);colorbar()
		#subplot(122);pcolor(M[0:1,:]);colorbar();show();quit()
		if perform_mean_sub:
			print '\n\t\tAlert: Performing mean subtraction before COV calc.'

			"""
			Mmean = np.mean(M, axis = 0).squeeze()
			M = M - Mmean
			"""
			Mmean = np.mean(M, axis = 1).squeeze()
			M = M - Mmean[:,None]

		perform_var_norm = 0
		if perform_var_norm:
			print '\n\t\tAlert: Performing variance normalisation before COV calc.'
			Mvar = np.var(M, axis = 0).squeeze()
			M = M/np.sqrt(Mvar)

		M = np.mat( M ).T
		Mt = M.T

		COV = (M * Mt) / (noofsims)# - 1)

		#imshow(COV);colorbar();show();quit()

		return COV
		'''
		print COV.shape

		#subplot(121);imshow(COV);colorbar()

		COV_TMP = np.zeros((npixels,npixels))
		for n in range(noofsims):

			M1 =  np.mat( MAPMAT[n].flatten() )
			Mt = M.T

			COV_CURR = Mt * M

			COV_TMP += COV_CURR
		COV_TMP /= noofsims
		#subplot(122);imshow(COV_TMP);colorbar()
		#show();quit()

		return COV
		'''

	def fn_calc_cov(self, SIMS, identifiers = None):

		#calculate covariance matrix now
		if self.npixels == None:
			self.npixels = np.prod(np.shape(SIMS)[-2:])
		npixels = self.npixels
		tqulen = self.tqulen#SIMS.shape[2]
		noofsims = len(SIMS)

		locdic = self.locdic
		if identifiers == None:
			identifiers = self.inidic['cov_identifiers']
		PARENT_C = {}
		for iden in identifiers:
			if len(iden) == 1: #just sample COV for a single field
				loc = locdic[iden]
				MAPFORCOV = SIMS[:,loc,:,:]
				COV = np.asarray(self.calcCov(MAPFORCOV, noofsims, npixels))
			elif len(iden) == 2: #sample covariance for 2 fields
				loc1, loc2 = locdic[iden[0]], locdic[iden[1]]
				MAPFORCOV1, MAPFORCOV2 = SIMS[:,loc1,:,:], SIMS[:,loc2,:,:]
				MAPFORCOV = np.concatenate( (MAPFORCOV1,MAPFORCOV2), axis=1 )
				COV = np.asarray(self.calcCov(MAPFORCOV, noofsims, 2*npixels))
				#imshow(COV);colorbar();show();quit()
			elif len(iden) == 3: #sample covariance for 2 fields
				loc1, loc2, loc3 = locdic[iden[0]], locdic[iden[1]], locdic[iden[2]]
				MAPFORCOV1, MAPFORCOV2, MAPFORCOV3 = SIMS[:,loc1,:,:], SIMS[:,loc2,:,:], SIMS[:,loc3,:,:]
				MAPFORCOV = np.concatenate( (MAPFORCOV1,MAPFORCOV2,MAPFORCOV3), axis=1 )
				COV = np.asarray(self.calcCov(MAPFORCOV, noofsims, 3*npixels))

			#title(iden);imshow(COV,origin='lower');colorbar();show();quit()

			PARENT_C[iden] = COV

		return PARENT_C


	def fn_Cov_precomputed(self, iden, M_cl, z_cl, noofsims = 25000, mode_refinement = 0, check_for_pos_def = 0, ra_for_rot = None, dec_for_rot = None):

		if hasattr(self.covfiles, "__len__") == 0:
			covfolder = self.covfolder
			#covfiles = sorted(glob.glob('%s/SIMS_COV_LEN*%s*' %(covfolder,noofsims)))
			searchstr = '%s/SIMS_COV_LEN*pkl*%s*_%s_%s' %(covfolder,noofsims,ra_for_rot,dec_for_rot)
			covfiles = sorted(glob.glob(searchstr))
		else:
			covfiles = self.covfiles

		if len( self.covdic[iden].keys() ) == 0:

			covdic_tmp = {}
			for ff in covfiles:
				dic = pickle.load(gzip.open(ff,'r'))
				#print 'file read'
				for keys in sorted(dic.keys()):

					if isinstance(keys,tuple):
						currM, currz = keys
						currM /= 1e14
					else:
						currM, currz = keys.split(',')
						currM, currz = float(currM), float(currz)

					if self.all_clusters_at_same_redshift:
						if currz <> z_cl:
							continue

					covdic_tmp[('%.03f, %.03f' %(currM, currz))] = dic[keys][iden]

			#for a particular ra, dec #not yet implemented
			

			self.covdic[iden] = covdic_tmp

		keyname = ('%.03f, %.03f' %(M_cl, z_cl))
		if keyname in self.covdic[iden].keys():
			COV_precompu = self.covdic[iden][('%.03f, %.03f' %(M_cl, z_cl))]

			if mode_refinement:
				COV_precompu = self.fn_mode_refinement(COV_precompu)#, totiter = 10)

			#imshow(COV_precompu, origin = 'lower'); colorbar(); show();quit()

			"""
			u, t, v = np.linalg.svd(COV_precompu)
			detval = np.prod( t )
			print detval, np.linalg.det(COV_precompu), np.linalg.slogdet(COV_precompu)
			
			imshow(COV_precompu, origin = 'lower'); colorbar(); show()

			L = np.linalg.cholesky(COV_precompu)
			#print L
			quit()
			"""
			
			if check_for_pos_def:
				try:
					L = np.linalg.cholesky(COV_precompu)
				except np.linalg.LinAlgError as errvals:
					print errvals, keyname
					return None

			#COV_precompu = np.dot(L, L.T.conj())
			#COV_precompu = np.dot(COV_precompu, COV_precompu.T.conj())
			
		else:
			COV_precompu = None
		return COV_precompu

	'''
	def fn_get_tSZ_COV(self,mapparams, reqdbox = None): #will use beam, TF, etc.

		logline = '\t\tCOV for tSZ emission now for this particulat mass limit - Hacked from Sanjay\'s suggestion instead of complicated tSZ fitting'
		logfile = open(self.log_file,'a')
		logfile.writelines('%s\n' %(logline))
		logfile.close()
		print logline

		tSZ_fitting = self.tSZ_fitting
		nsims_CtSZ = 10

		C_tSZ = {}
		#M_arr = np.asarray([0.001, 2.5, 5., 7.5]) * 1e14
		M_arr = self.M_arr * 1e14
		for masscnt,MM in enumerate(M_arr):
			T_0 = tSZ_fitting.fn_guess_Tsz_fromega_mass( MM / 2.) #2.5 factor to roughly get M500 from M200
			#print T_0, MM
			p0 = np.asarray([T_0,tSZ_fitting.theta_core,tSZ_fitting.tSZ_beta])

			#get the beta model
			beta_model = tSZ_fitting.fitting_func(p0,p0,self.RADIUS_for_beta_model,return_fit=1)
			#subplot(121);imshow(beta_model);colorbar()
			#add beam, TF to this guy
			beta_model = np.fft.ifft2( np.fft.fft2(beta_model) * self.Bl * self.TWODTF[0]).real
			#subplot(122);imshow(beta_model);colorbar();show();quit()

			#ideally we want to perform several realisations of this; but for now restricting to just 1##CHECK THIS
			ex1,ex2,ey1,ey2 = reqdbox
			beta_model = beta_model[ex1:ex2, ey1:ey2]
			nx, ny = beta_model.shape
			npixels = nx * ny
			
			BETAMODELFORCOV = np.asarray( [beta_model for iii in range(nsims_CtSZ) ] )
			keyname = '%.3f' %(MM/1e14)
			C_tSZ[keyname] = self.calcCov(BETAMODELFORCOV, len(BETAMODELFORCOV), npixels)

			#subplot(2,2,masscnt+1);imshow(C_tZ[MM]);colorbar();title(MM);
			"""
			#add beam to this guy
			beta_model = np.fft.ifft2( np.fft.fft2(beta_model) * self.Bl ).real
			subplot(132);imshow(beta_model);colorbar()

			#add TF to this guy
			beta_model = np.fft.ifft2( np.fft.fft2(beta_model) * self.TWODTF[0] ).real
			subplot(133);imshow(beta_model);colorbar();show();quit()
			"""
		#show();quit()

		return C_tSZ

	'''
	def fn_get_tSZ_COV(self, MAP_to_fit, mapparams, M_arr = None, reqdbox = None): #will use beam, TF, etc.

		logline = '\t\tCOV for tSZ emission now for this particulat mass limit - Hacked from Sanjay\'s suggestion instead of complicated tSZ fitting'
		logfile = open(self.log_file,'a')
		logfile.writelines('%s\n' %(logline))
		logfile.close()
		print logline

		tSZ_fitting = self.tSZ_fitting
		nsims_CtSZ = 5000

		clra, cldec = 0., 0.
		nx, ny, dx, dy = mapparams
		minval, maxval = clra - (nx * dx),  clra + (nx * dx)
		ra = dec = np.linspace(minval, maxval, ny)
		RA, DEC = np.meshgrid(ra,dec)
		RADIUS = (RA ** 2. + DEC ** 2.)**0.5

		C_tSZ = {}
		#M_arr = np.asarray([5.6, 7.5]) * 1e14
		if M_arr == None:
			M_arr = self.M_arr * 1e14

		for masscnt,MM in enumerate(M_arr):

			T_0 = tSZ_fitting.fn_guess_Tsz_fromega_mass( MM / 2.5) #2.5 factor to roughly get M500 from M200
			#print T_0, MM
			p0 = np.asarray([T_0,tSZ_fitting.theta_core,tSZ_fitting.tSZ_beta])

			fixed = [False,False,False]
			lbounds = np.asarray([0.9*T_0,0.25,0.6])
			ubounds = np.asarray([1.1*T_0,0.75,1.2])
			#fit the beta model
			p1, success = optimize.leastsq(tSZ_fitting.fitting_func, p0[:], args=(p0, RADIUS, MAP_to_fit, lbounds, ubounds,fixed))
			
			#get the beta model
			beta_model = tSZ_fitting.fitting_func(p1,p0,RADIUS,return_fit=1)
			ex1,ex2,ey1,ey2 = reqdbox
			BETAMODELFORCOV = np.zeros( (nsims_CtSZ, ex2-ex1, ey2-ey1) )
			for simcnt in range(nsims_CtSZ):

				BETA_MAP = beta_model + np.random.randn(ny, nx)
				#add beam, TF to this guy
				BETA_MAP = np.fft.ifft2( np.fft.fft2(BETA_MAP) * self.Bl * self.TWODTF[0]).real

				#ideally we want to perform several realisations of this; but for now restricting to just 1##CHECK THIS
				BETA_MAP = beta_model[ex1:ex2, ey1:ey2]
				BETAMODELFORCOV[simcnt] = BETA_MAP

			npixels = (ex2-ex1) * (ey2-ey1)
			
			keyname = '%.3f' %(MM/1e14)
			C_tSZ[keyname] = self.calcCov(BETAMODELFORCOV, len(BETAMODELFORCOV), npixels)


			imshow(C_tSZ[keyname]);colorbar();title(keyname);show();quit()
		#show();quit()

		return C_tSZ

	'''
	def fn_get_tSZ_COV(self, MAP_to_fit, MM, mapparams, reqdbox = None): #will use beam, TF, etc.

		"""
		use the map, guess the tSZ for a specific mass, fit beta / theta_core
		"""
		tSZ_fitting = self.tSZ_fitting
		nsims_CtSZ = 10

		if MM<1e2:
			MM *= 1e14

		T_0 = tSZ_fitting.fn_guess_Tsz_fromega_mass( MM / 2.5) #2.5 factor to roughly get M500 from M200

		clra, cldec = 0., 0.
		nx, ny, dx, dy = mapparams
		minval, maxval = clra - (nx * dx),  clra + (nx * dx)
		ra = dec = np.linspace(minval, maxval, ny)
		RA, DEC = np.meshgrid(ra,dec)
		RADIUS = (RA ** 2. + DEC ** 2.)**0.5

		p0 = np.asarray([T_0,tSZ_fitting.theta_core,tSZ_fitting.tSZ_beta])
		fixed = [1,0,0]
		lbounds, ubounds = None, None
		#fit the beta model
		p1, success = optimize.leastsq(tSZ_fitting.fitting_func, p0[:], args=(p0, RADIUS,MAP_to_fit, lbounds, ubounds,fixed))

		#get the beta model
		beta_model = tSZ_fitting.fitting_func(p1,p0,RADIUS,return_fit=1)

		#add beam, TF to this guy
		###### beta_model = np.fft.ifft2( np.fft.fft2(beta_model) * self.Bl * self.TWODTF[0]).real
		#subplot(122);imshow(beta_model);colorbar();show();quit()

		#ideally we want to perform several realisations of this; but for now restricting to just 1##CHECK THIS
		if not reqdbox == None:
			ex1,ex2,ey1,ey2 = reqdbox
			beta_model = beta_model[ex1:ex2, ey1:ey2]
			nx, ny = beta_model.shape

		npixels = nx * ny
		
		BETAMODELFORCOV = np.asarray( [beta_model for iii in range(nsims_CtSZ) ] )
		C_tSZ = self.calcCov(BETAMODELFORCOV, len(BETAMODELFORCOV), npixels)

		#imshow(C_tSZ);colorbar();title(MM);show();quit()

		return C_tSZ
	'''

	def lnlike(self, iden, M_cl, D_LEN = None, z_cl_arr = None, noofsims = 20000, use_cov_interp = 0, expnoiselevel = None, cls_dic_keys = None):

		if self.noise_present == 1:
			nx, ny, dx, dy = self.mapparams
			degrade_resol_fac = self.inidic['degrade_resol_fac']
			mapparams_deg = [nx/degrade_resol_fac, ny/degrade_resol_fac, dx * degrade_resol_fac, dy * degrade_resol_fac]

			if not self.is_seq(self.C_NOISE): #will be calculated only once for a run
				if not self.is_seq(expnoiselevel):
					expnoiselevel = self.expnoiselevel

				logline = '\t\t"White" Noise present. Level = %s. Calculating "white" noise covariance matrix now and storing for other iterations' %(expnoiselevel)
				logfile = open(self.log_file,'a')
				logfile.writelines('%s\n' %(logline))
				logfile.close()
				print logline

				self.C_NOISE = self.fn_noise_cov(self.mapparams, noofsims, reqdbox = None, expnoiselevel = expnoiselevel)
				#C_NOISE = self.fn_noise_cov(mapparams_deg, noofsims, reqdbox = None, expnoiselevel = expnoiselevel)
				#imshow(C_NOISE[iden]);colorbar();show();quit()


		elif self.noise_present == 2:

			if not self.is_seq(self.C_NOISE): #will be calculated only once for a run

				logline = '\t\t"SPTpol like" Noise present. Calculating noise covariance matrix using SPTpol noise model and storing for other iterations'
				logfile = open(self.log_file,'a')
				logfile.writelines('%s\n' %(logline))
				logfile.close()
				print logline

				try:
					self.noisesims
				except:
					noisesims = noise_sims()
					self.noisesims = noisesims
				noisesims = self.noisesims
				#print self.cmb_noise_randomseedval, self.whichmap;quit()
				#NOISE_SIMS = noisesims.fn_get_noise_sims_from_noise_cls(self.mapparams, noofsims, random_seed = self.cmb_noise_randomseedval, whichmap = self.whichmap)
				#20161027 - make NOISE sims on a bigger "original" boxsize and then extract
				NOISE_SIMS = noisesims.fn_get_noise_sims_from_noise_cls(self.simmapparams, noofsims, reqdbox = self.reqdbox, random_seed = self.cmb_noise_randomseedval, whichmap = self.whichmap)
				#print 'implement TF multiplication here - quitting'
				self.C_NOISE = self.fn_calc_cov(NOISE_SIMS)

		if not self.is_seq(self.C_for_regular): #will be calculated only once for a run
			self.C_for_regular = self.fn_noise_cov(self.mapparams, noofsims, reqdbox = None, expnoiselevel = np.tile(1e-1,self.tqulen))

		C_NOISE = self.C_NOISE

		#imshow(C_NOISE[iden]);colorbar();show()

		#get the FG cov matrix now
		if self.FG_added:
			if not self.is_seq(self.C_FG):
				self.C_FG = self.fn_get_foreground_DFSG_power(self.simmapparams, noofsims, reqdbox = self.reqdbox, random_seed = self.cmb_noise_randomseedval, perform_lensing = 0, return_cov = 1)
		#imshow(self.C_FG);colorbar();show();quit()

		Cinv_dic = {}
		logL_arr = []
		start = time.time()
		for cnt in range(len(D_LEN)):

			#get interpolated COV matrices
			if not use_cov_interp:
				C = self.fn_Cov_precomputed(iden, M_cl, z_cl_arr[cnt])#, ra_for_rot = self.ra_arr[cnt], dec_for_rot = self.dec_arr[cnt])
			else:
				#C = self.fn_Cov_interp(M_cl, z_cl_arr[cnt])
				C = self.fn_Cov_interp2(iden, M_cl, z_cl_arr[cnt])

			if hasattr(C, "__len__") == 0:
				#print cnt, M_cl, z_cl_arr[cnt], C#;quit()
				logline = '\tAlert *** Cluster number = %s; No COV obtained for (%s, %s),' %(cnt, M_cl, z_cl_arr[cnt])
				logfile = open(self.log_file,'a')
				logfile.writelines('%s\n' %(logline))
				logfile.close()
				print logline
				logL_arr.append(np.nan)
				continue

			#do it here as COV calc. is just C_CMB
			if self.noise_present:
				C = C + C_NOISE[iden]

			"""
			if self.add_tSZ_cov and iden == 'T':
				C_tSZ = self.tSZ_cov_dic[cls_dic_keys[cnt]][0] #contains COV_tSZ, beta_model_params for the cluster, fit_flag
				C = C + C_tSZ
				#imshow(C);colorbar();show();quit()
			"""
			'''
			if self.perform_tSZ_removal == 0 and iden == 'T' and M_cl>0.:
				if not self.is_seq(self.C_tSZ):
					C_tSZ = self.fn_get_tSZ_COV(self.simmapparams, reqdbox = self.reqdbox)
					self.C_tSZ = C_tSZ

				C = C + self.C_tSZ['%.3f' %(M_cl)]
				#imshow(C);colorbar();show();quit()
			if self.perform_tSZ_removal == 0 and iden == 'T' and M_cl>0.:
				nx, ny = self.mapparams[0], self.mapparams[1]
				npixels = nx * ny
				#CURRMAP = D_LEN[cnt].reshape(nx, ny)
				#C_tSZ = self.fn_get_tSZ_COV(CURRMAP, M_cl, self.mapparams)
				C = C + C_tSZ
			'''

			#imshow(C);colorbar();show();quit()

			#add COV for FOREGROUND
			if self.is_seq(self.C_FG):
				C = C + self.C_FG

			#for matrix regularisation
			if self.is_seq(self.C_for_regular): 
				C = C + self.C_for_regular[iden]

			#imshow(C_NOISE[iden] + self.C_FG);colorbar();show();quit()

			logdetval = None
			C = np.mat(C)
			sign, logdetval = np.linalg.slogdet(C)
			logdetval = logdetval * sign

			#subplot(121);imshow(C, origin='lower');colorbar();subplot(122);imshow(C.T, origin='lower');colorbar();show();quit()

			d_len = D_LEN[cnt].squeeze() #this count runs for different clusters
			d_len = d_len - np.mean(d_len) #subtract the mean? #CHECK
			#hist(d_len,bins=10,histtype='step')
			d_len = np.mat(d_len).T

			#start = time.time()
			currkey = ('%.03f, %.03f' %(M_cl, z_cl_arr[cnt]))
			if not currkey in Cinv_dic.keys():
				Cinv = sc.linalg.pinv2(C)
				#Cinv2 = sc.linalg.pinv(C)
				#Cinv3 = sc.linalg.inv(C)
				Cinv_dic[currkey] = Cinv
			else:
				Cinv = Cinv_dic[currkey]

			'''
			if M_cl>0.:
				subplot(131);imshow(sc.linalg.pinv(C));colorbar()
				subplot(132);imshow(sc.linalg.pinv2(C));colorbar()
				subplot(133);imshow(sc.linalg.inv(C));colorbar();show();quit()
			'''

			npixels = len(C)
			t1 = -.5 * npixels * np.log( 2 * np.pi )

			if logdetval == None:
				#t2 = -0.5 * np.log( detval)
				logline = '\tAlert *** Cluster number = %s; determinant is bad. (%s, %s),' %(cnt, M_cl, z_cl_arr[cnt])
				logfile = open(self.log_file,'a')
				logfile.writelines('%s\n' %(logline))
				logfile.close()
				print logline

				logL_arr.append(np.nan)
			else:
				t2 = -0.5 * logdetval

			#print C.shape;quit()
			#dummy = np.asarray(np.dot( Cinv, d_len )).squeeze()
			t3 = -0.5 * np.asarray( np.dot(d_len.T, np.dot( Cinv, d_len ))).squeeze()
			logL = (t1 + t2 + t3)
			#logL = (t2 + t3)

			#logL = np.exp(t3) / (detval**0.5)
			#logL = t1 - (t2 + t3)

			#t3 = -0.5 * np.asarray( np.dot(d_len.T, np.dot( Cinv, d_len ))).squeeze()
			#logL = np.exp(t3)  / detval**0.5
			#print t3, np.exp(t3), detval, detval**0.5, logL;quit()
			#print '\t\t', M_cl, z_cl_arr[cnt], logdetval, t1, t2, t3, np.linalg.slogdet(np.mat(Cinv))[1]
			print '\t\t', cnt, M_cl, z_cl_arr[cnt], logdetval, t2, t3
			#quit()

			val = np.asarray(logL).squeeze()
			logL_arr.append(val)

		retval = np.sum(logL_arr)
		end = time.time()

		logline = '\t(%s, %s),' %(M_cl, retval)
		logfile = open(self.log_file,'a')
		logfile.writelines('%s\n' %(logline))
		logfile.close()
		print logline
		
		return retval

	def fn_normalize(self,data,Nmin=0.,Nmax=1.):

		arrmin,arrmax=min(data),max(data)
		normed_data=(data-arrmin)/(arrmax-arrmin)
		normed_data=normed_data*Nmax+Nmin

		return normed_data

	def fn_max_likelihood(self, D_UNLEN, D_LEN, z_cl_arr, actual_masses, covfolder = 'data/covs', opfolder = 'pkls/', noofsims = 20000, use_cov_interp = 1, no_lensing = 0, expnoiselevel = None, C_NOISE = 0, cls_dic_keys = None):

		import numpy as np, glob, time, pickle, gzip#multiprocessing as mp
		self.tSZ_fitting = tSZ_fitting()

		'''#### 20161026 ####
		#M_arr = np.arange(4.,6.1,0.1)
		self.covfolder = str(covfolder)
		covfiles_all = sorted(glob.glob('%s/SIMS_COV_LEN*pkl*%s*' %(covfolder,noofsims)))

		covfiles = []
		for f in covfiles_all:
			covfiles.append(f)

		self.covfiles = covfiles
		#### 20161026 ####'''

		if self.is_seq(C_NOISE):
			self.C_NOISE = C_NOISE

		###
		#initialise self.covdic here
		for iden in self.inidic['cov_identifiers']:
			self.covdic[iden] = {}
		###

		'''
		M_arr = np.asarray( [float(x[x.rfind('_')+1:]) for x in covfiles] )
		if use_cov_interp:
			delM = 0.2
			M_arr = np.arange(minM,maxM+delM, delM)
		'''

		M_arr = np.arange(0.,15.5,0.5)

		M_arr = np.asarray( sorted(M_arr) )
		self.M_arr = M_arr
		logline = '\n\tReal likelihood calc. starts now'
		logfile = open(self.log_file,'a')
		logfile.writelines('%s\n' %(logline))
		logfile.close()
		print logline

		########################################################
		########################################################
		start = time.time()
		opdic = {}

		timedatestr='%s_%s' %(time.strftime('%Y%m%d_%H%M%S', time.localtime()),self.cmb_noise_randomseedval)
		compdic = {0:'T, r', 1:'Q, b', 2:'U, g'}

		## modify D_LEN to handle xcovs of TQ,QU,UU
		identifiers = self.inidic['cov_identifiers']
		locdic = self.locdic

		if self.add_tSZ_cov: #pick the tSZ COVS corresponding to each cluster for "T" estimator
			#tSZ_cov_file = str(self.covfolder) + '/tSZ_COV.pkl.gz_%ssims' %(noofsims)
			self.tSZ_cov_dic = pickle.load(gzip.open(self.tSZ_cov_file,'r'))

		#for tqucnt in range(self.tqulen):
		for iden in identifiers:
			resarr = []
			M_arr_ref = []
			for cnt, MM in enumerate(M_arr):

				if len(iden) == 1:
					loc = locdic[iden]
					MOD_D_LEN = D_LEN[loc]
					tqucnt = 0 
				elif len(iden) == 2:
					loc1, loc2 = locdic[iden[0]], locdic[iden[1]]
					NEW_D_LEN_1,NEW_D_LEN_2 = D_LEN[loc1],D_LEN[loc2]
					MOD_D_LEN = np.concatenate( (NEW_D_LEN_1, NEW_D_LEN_2), axis = 1)
					tqucnt = 0 

				#check if all clusters are at same redshift (for faster operation - store only those COVs)
				z_cl_arr = np.asarray(z_cl_arr)
				if len(np.unique(z_cl_arr))==1:
					logline = '\tAlert *** All clusters are at same redshift. I will not read COVs for all redshifts. z is %s' %(np.unique(z_cl_arr)[0])

					self.all_clusters_at_same_redshift = 1

				else:
					logline = '\tClusters are not at same redshift. I will read COVs for all redshifts'

				logfile = open(self.log_file,'a')
				logfile.writelines('%s\n' %(logline))
				logfile.close()
				print logline

				logL_val = self.lnlike(iden, MM, MOD_D_LEN, z_cl_arr, noofsims = noofsims, use_cov_interp = use_cov_interp, expnoiselevel = expnoiselevel, cls_dic_keys = cls_dic_keys)
				if np.isnan(logL_val):
					continue
				#plt.plot(MM,logL_val,'ro')
				M_arr_ref.append(MM)
				resarr.append(logL_val)
				opdic['%.3f' %(MM)] = logL_val

			tqucomp = compdic[tqucnt].split(',')[0].strip()
			colorval = compdic[tqucnt].split(',')[1].strip()

			'''
			#show();quit()
			M_arr = np.asarray(M_arr_ref)
			resarr = np.asarray(resarr)
			resarr = self.fn_normalize(resarr)
			#resarr = np.exp(resarr)

			plot(M_arr,resarr, lw=3., label = iden)#, color = colorval)

			xlabel('$Mass\ [10^{14}\ M_{\odot}]$ -->')
			#ylabel('$Normalised\ L$ -->')
			#vlines(np.median(actual_masses),0,1,linestyle = 'dashed', color = 'r', lw = 2.)
			legend(loc = 1, fancybox = 1)
			#ylim(0.,1.1)

		show();quit()
		if 1==1:
			'''

			resdic = {}
			resdic['clustermasses']	= actual_masses
			resdic['results'] = opdic

			totalclus = len(D_UNLEN[0])
			modopfolder = '%s/%s' %(opfolder, iden)

			if not os.path.exists(modopfolder):
				os.system('mkdir %s' %(modopfolder))

			if self.noise_present == 1:
				if len(iden) == 1:
					loc = self.locdic[iden]
					noise_folder = 'noise_%.3f' %(expnoiselevel[loc])
				elif len(iden) == 2:
					loc_1,loc_2 = self.locdic[iden[0]],self.locdic[iden[1]]
					if loc_1>2:
						loc_1=2
					if loc_2>2:
						loc_2=2
					noise_folder = 'noise_%.3f_%.3f' %(expnoiselevel[loc_1],expnoiselevel[loc_2])
			elif self.noise_present == 2:
				noise_folder = 'sptpol_noise_model'

			modopfolder = '%s/%s/%s' %(opfolder, iden,noise_folder)
			os.system('mkdir %s' %(modopfolder))

			if no_lensing:
				modopfolder = '%s/no_lensing' %(modopfolder)
				
				os.system('mkdir %s' %(modopfolder))

			if self.noise_present == 0:
				opfilename = '%s/%s_%s_logL.pkl.gz'  %(modopfolder,timedatestr,totalclus)
			else:
				#opfilename = '%s/%s_%s_noise%s_logL.pkl.gz'  %(modopfolder,timedatestr,totalclus,expnoiselevel[tqucnt])
				opfilename = '%s/%s_%s_logL.pkl.gz'  %(modopfolder,timedatestr,totalclus)
			pickle.dump(resdic, gzip.open(opfilename,'w'))
			#savefig('plots/%sclusters_%s.png' %(totalclus,timedatestr))
	
		end = time.time()
		logline = '\n\tLikelihood calc. complete. time taken = %s secs' %(end-start)
		logfile = open(self.log_file,'a')
		logfile.writelines('%s\n' %(logline))
		logfile.close()
		print logline

		quit()

	def fn_get_EB_from_QU(self, Q, U, mapparams):

		lx, ly = self.get_lxly(mapparams)
		angle = 2.*np.arctan2(lx, -ly)

		QFFT, UFFT = np.fft.fft2(Q),np.fft.fft2(U)
		E = np.fft.ifft2( np.cos(angle) * QFFT + np.sin(angle) * UFFT ).real
		B = np.fft.ifft2( -np.sin(angle) * QFFT + np.cos(angle) * UFFT ).real

		return E, B

	def fn_YDELTA_MDELTA_scaling_shaw(self, M, z_cl, nu = 150e9, logAy = -5.54, alpha = 1.67, beta = 0.66, norm_M = 1e14, theta_c_physical = 0.5, theta_int_physical = 1., b = 0., sigma_logY_scatter = 0.2):

		"""
		arXiv: 0710.4555 (table 1)
		"""		

		omega_m, omega_k, omega_lambda = self.inidic['omega_m'], 0., self.inidic['omega_lambda']

		def e_z(z):
			return (1/np.sqrt(omega_m*((1+z)**3)+ omega_k*((1+z)**2) + omega_lambda))

		def fn_r_delta(M_delta,z, delta_overdensity = 200., rho_def = 'crit'): #M_200 in solar mass
			"""
			r_delta = M_delta / ( 4*pi/3 * delta * rho(z) )
			"""
			import numpy as np
			solar_mass_to_kgs = 1.988e30
			Mpc_to_metres = 3.086e22

			from astropy.cosmology import FlatLambdaCDM
			cosmo = FlatLambdaCDM(H0 = self.inidic['h']*100., Om0 = self.inidic['omega_m'])

			M_delta = M_delta * solar_mass_to_kgs

			rho_z = cosmo.critical_density(z).to('kg/m3').value
			if rho_def == 'mean':
			    rho_z = rho_z * cosmo.Om(z)

			r_delta_cubed = M_delta / ( (4*np.pi/3) * delta_overdensity * rho_z)
			r_delta = (r_delta_cubed)**(1/3.) / Mpc_to_metres #metres

			return r_delta


		def fn_distances_volume(e_z_int, z, dz = 0.01):

			#H0 = 100.0    # HUBBLE CONSTANT = 100h Km/s/Mpc

			# HUBBLE CONSTANT IN STANDARD UNITS - DIMENSIONS sec-1
			H0 = self.inidic['h'] * 100.
			c = 2.99797e5
			H0_std = (H0/(3.08568025 * 10**19))

			# HUBBLE DISTANCE IN Mpc
			d_h = c/H0

			# HUBBLE TIME IN h-1 yr
			t_h = 1/H0_std

			# TOTAL LINE-OF-SIGHT COMOVING DISTANCE	
			d_c = d_h * e_z_int

			# TRANSVERSE COMOVING DISTANCE
			if (omega_k==0.0):
				d_t = d_c
			elif (omega_k>0.0):
				d_t = d_h/np.sqrt(omega_k) * np.sinh(np.sqrt(omega_k)*d_c/d_h)
			else:
				d_t = d_h/np.sqrt(abs(omega_k)) * np.sin(np.sqrt(abs(omega_k))*d_c/d_h)

			if (omega_lambda==0.0):
				d_t = d_h * 2 *(2 - (omega_m *(1-z)) - ((2-omega_m) * (np.sqrt(1+(omega_m*z))))) / (omega_m**2 * (1+z))


			# ANGULAR DIAMETER DISTANCE
			d_a = d_t / (1+z)

			# LUMINOSITY DISTANCE
			d_l = (1+z) * d_t

			# DISTANCE MODULUS
			dm = 5.0 * np.log10(d_l*10**6/10) # 1 Mpc = 10**6 Mpc

			# COMOVING VOLUME
			v_c = (4.0 * np.pi * d_t**3) / 3.0

			#comoving volume element
			dv_c = d_h * ( (1+z)**2. * d_a**2. * dz) * 4 * np.pi / e_z_int

			return d_c, d_t, d_a, d_l, dm, v_c, dv_c


		norm_M = norm_M/self.inidic['h']

		Ay = 10**(logAy)
		E_z = 1./e_z(z_cl)

		e_z_int, e_z_int_err = integrate.quad(e_z,0.,z_cl)
		d_c, d_t, d_a, d_l, dm, v_c, dv_c = fn_distances_volume(e_z_int, z_cl)

		t1 = Ay
		t2 = E_z**beta
		t3 = ( M / norm_M ) ** alpha
		#t4 = np.random.lognormal(sigma = sigma_logY_scatter)

		yflux_integrated = t1 * t2 * t3

		if (1):
			##sigma_logY_scatter = 0.3 #https://arxiv.org/pdf/1312.3015.pdf: table 1
			t4 = np.random.normal(0., sigma_logY_scatter)
			lnyflux_integrated = np.log( yflux_integrated ) + t4
			yflux_integrated = np.exp(lnyflux_integrated)
		f_x = fn_get_f_x(nu)

		#fixing theta_int and theta_c in arcmins
		##https://arxiv.org/pdf/1312.3015.pdf eq. (6)
		t1 = yflux_integrated * f_x * self.inidic['T_cmb']  * 1e6
		t2 = np.pi * self.shaw_theta_c**2. 
		t3 = np.log10( 1. + (self.shaw_theta_int/self.shaw_theta_c)**2. )
		peak_T = t1 / t2 / t3

		return peak_T


		'''
		#divide yflux_integrated by the intergation quantity based on cluster r_200, etc. radius
		r_delta = fn_r_delta(M, z_cl, delta_overdensity = self.inidic['mass_def'], rho_def = self.inidic['rho_def'])
		theta_int_physical = r_delta
		theta_int = np.degrees(theta_int_physical/d_a) * 60. #arcmins

		#beta model
		theta_c = np.degrees(theta_c_physical / d_a) * 60.

		theta_int = self.2.
		theta_c = 0.5 #arcmins

		print theta_int, theta_c
		'''

		#fixing theta_int and theta_c in arcmins
		##https://arxiv.org/pdf/1312.3015.pdf eq. (6)
		t1 = yflux_integrated * f_x * self.inidic['T_cmb']  * 1e6
		t2 = np.pi * self.theta_c**2. 
		t3 = np.log10( 1. + (self.theta_int/self.theta_c)**2. )
		peak_T = t1 / t2 / t3

		return peak_T

		'''
		t1 = yflux_integrated * f_x * param_dict['T_cmb'] * 1e6 / np.pi / theta_c**2.
		t2 = np.log10(1 + (theta_int/theta_c)**2.)

		return t1/t2
		'''

		def Ysz_int(theta):
			return 2 * np.pi * theta * (1 + (theta/theta_c)**2.) ** -1

		small_number = self.small_number
		Ysz_int, Ysz_int_err = integrate.quad(Ysz_int, small_number, theta_int)
		
		peak_compton_y = (yflux_integrated / Ysz_int)

		peak_T = peak_compton_y * f_x *  self.inidic['T_cmb']  * 1e6

		return peak_T

		"""
		clusterstuff = cluster_stuff()
		ysz_Tsz_conv_fac = clusterstuff.compton_y_to_delta_Tcmb(nu)

		print ysz_Tsz_conv_fac, peak_T, peak_compton_y * ysz_Tsz_conv_fac * 1e6
		return peak_T
		"""
	def fn_fit_Arnaud_template(self,mapparams, OBSMAP): #Mval_arr, zval_arr, dx_finer = 0.5, cosmo_params_dict = None, mass_def = 200., rho_def = 'mean', howmany_theta_500 = 4.5, convert_to_T = 1, nu = 150e9, sigma_logY_scatter = 0.2, random_seed = None):

		#how many clusters
		totalclus = 1 #len(Mval_arr)

		nx,ny,dx,dy = mapparams
		#map details
		pixel_size_arcm = dx
		pixels = nx/2
		MAPS = np.zeros((totalclus, 2*pixels,2*pixels))

		#cosmology details
		if 1==1:#cosmo_params_dict == None:
			cosmo_params_dict = {}
			cosmo_params_dict['h'] = 0.6774
			cosmo_params_dict['omega_m'] = 0.307320
			cosmo_params_dict['omega_lambda'] = 0.692680
			cosmo_params_dict['z_lss'] = 1089.
			cosmo_params_dict['T_cmb'] = 2.725
		H0 = cosmo_params_dict['h'] * 100.
		zval_arr = [0.7] # fixing this for methods paper
		Mval_arr = [2*1e14]
		param_dict = self.inidic
		mass_def = param_dict['mass_def']
		rho_def = param_dict['rho_def']
		
		for cntr, (Mval, zval) in enumerate(zip(Mval_arr, zval_arr)):

			if mass_def == 200.: #convert to M500 since Arnaud profile is for M500
				mdef = '%s%s' %(int(mass_def), rho_def[0])
				mdefout = '500%s' %(rho_def[0])
				cval = concentration.concentration(Mval, mdef, zval)
				Mval, r200val, c200val = mass_defs.changeMassDefinition(Mval, cval, zval, mdef, mdefout, profile='nfw')


			r500 = nc.M500_to_R500(Mval, zval, cosmo_params_dict['omega_m'], cosmo_params_dict['omega_lambda'], H0)
			theta500 = nc.radius_to_theta_arcm (r500, zval, cosmo_params_dict['omega_m'], cosmo_params_dict['omega_lambda'],  H0) # variable qunatity


			ymax, xmax = MAPS[cntr].shape
			xcoord, ycoord = np.arange(0,xmax,1), np.arange(0, ymax,1)
			xcoord_mesh, ycoord_mesh = np.meshgrid(xcoord, ycoord)
			xcen = MAPS[cntr].shape[1] / 2 - 1
			ycen = MAPS[cntr].shape[0] / 2 - 1
			distance_each_pix = np.hypot(xcoord_mesh - xcen, ycoord_mesh - ycen)
			initial_guess =  Mval
			import scipy.optimize as opt
			cmb_90 = self.cmb_90[0]/1e6
			tsz_map = ((cmb_90 - OBSMAP)*1.48)[95:105,95:105]
			x_1, y_1 = np.arange(0,10,1), np.arange(0, 10,1)
			x,y = np.meshgrid(x_1, y_1)
			xcen_1 = 10 / 2 - 1
			ycen_1 = 10/2 -1 
			distance_each_pix = np.hypot(x_1 - xcen, y_1- ycen)
			howmany_theta_500  = 3.6
			MAPS = np.zeros((totalclus, 2*5,2*5))
			nx, ny = tsz_map.shape
			boxsize = 5
			dx = 0.5
			mapparams = [nx, ny, dx, dy]
			clra, cldec = 0., 0.
			minval, maxval = clra-boxsize/2/60.,  clra+boxsize/2/60.
			ra = dec = np.linspace(minval, maxval, nx)
			RA, DEC = np.meshgrid(ra,dec)
			RADEC = [RA, DEC]
			tsz_data = self.fn_radial_profile(tsz_map,RADEC, maxbin = 5., bin_size = 0.5)
			
			tsz_data = tsz_data[:,1]
			r =np.asarray([ 0.25,  0.75,  1.25,  1.75,  2.25,  2.75,  3.25,  3.75,  4.25,  4.75])
			def get_tsz_map(r, Mval):
				map_fit = np.zeros(len(r))
				theta500 = nc.radius_to_theta_arcm (r500, zval, cosmo_params_dict['omega_m'], cosmo_params_dict['omega_lambda'],  H0)
				for i,rr in enumerate(r): 
					Yout = nc1.M500_to_Yx500_arcm2 ((i*pixel_size_arcm+pixel_size_arcm)/theta500,Mass=Mval)
					Yin = nc1.M500_to_Yx500_arcm2 ((i*pixel_size_arcm)/theta500,Mass=Mval)
					pixel_area = pixel_size_arcm**2.
					#put in the values of Yout - Yin in diffferent pixels
					#number_pixels_between_radius = MAPS[cntr, (distance_each_pix >= radius/pixel_size_arcm) &  (distance_each_pix < (radius+pixel_size_arcm)/pixel_size_arcm)].shape[0]
					map_fit[i]    =       (Yout-Yin)# / number_pixels_between_radius / pixel_area
				def fn_get_f_x(nu = 150e9):

					h=6.62607004e-34 #Planck constant in m2 kg / s
					k_B=1.38064852e-23 #Boltzmann constant in m2 kg s-2 / K-1

					x = h * nu / (k_B * self.inidic['T_cmb'])
					g_nu = (np.exp(x) + 1) / (np.exp(x) - 1)
					f_x = (x * g_nu - 4) # +  rela. corrections

					return f_x

				f_x = fn_get_f_x(nu = 150e9)
				map_fit = map_fit * f_x *self.inidic['T_cmb'] 
				# getting_beam
				dx = 0.5*self.arcmins2radians #hardcoded
				freq = np.fft.fftfreq(len(r),dx)
				freq = 2*np.pi*abs(freq)
				fwhm = self.exp_beam* self.arcmins2radians
				sigma = fwhm/np.sqrt(8. * np.log(2.))
				bm = np.exp(-.5 * freq**2. * sigma**2)
				map_fit = np.fft.ifft(np.fft.fft(map_fit)*bm).real
				return map_fit
			
			popt, pcov = opt.curve_fit(get_tsz_map, r, tsz_data, p0=initial_guess)
		
			Mval_fit = popt
			fitted_template = get_tsz_map(r,Mval_fit)
			mapparams = 200,200,0.5,0.5

			full_fit = self.fn_get_Arnaud_ymap_from_Nikhel_Gupta_1(mapparams, Mval_arr, zval_arr, dx_finer = 0.1, cosmo_params_dict = None, mass_def = 200., rho_def = 'mean', howmany_theta_500 = 4.5, convert_to_T = 1, nu = 150e9, sigma_logY_scatter = None , random_seed = None,m_val_fit = Mval_fit)[0]
			
		return full_fit

	def fn_get_Arnaud_ymap_from_Nikhel_Gupta_1(self, mapparams, Mval_arr, zval_arr, dx_finer = 0.5, cosmo_params_dict = None, mass_def = 200., rho_def = 'mean', howmany_theta_500 = 4.5, convert_to_T = 1, nu = 150e9, sigma_logY_scatter = 0.2, random_seed = None,m_val_fit = None ):

		#how many clusters
		totalclus = len(Mval_arr)

		if (1): #do this on a finer grid
			nx, ny, dx, dy = mapparams
			ds_fac = int(dx/dx_finer)
			nx_finer = ds_fac * nx

		#map details
		pixel_size_arcm = dx_finer
		pixels = nx_finer/2
		MAPS = np.zeros((totalclus, 2*pixels,2*pixels))

		#cosmology details
		if cosmo_params_dict == None:
			cosmo_params_dict = {}
			cosmo_params_dict['h'] = 0.6774
			cosmo_params_dict['omega_m'] = 0.307320
			cosmo_params_dict['omega_lambda'] = 0.692680
			cosmo_params_dict['z_lss'] = 1089.
			cosmo_params_dict['T_cmb'] = 2.725
		H0 = cosmo_params_dict['h'] * 100.


		for cntr, (Mval, zval) in enumerate(zip(Mval_arr, zval_arr)):

			if mass_def == 200.: #convert to M500 since Arnaud profile is for M500
				mdef = '%s%s' %(int(mass_def), rho_def[0])
				mdefout = '500%s' %(rho_def[0])
				cval = concentration.concentration(Mval, mdef, zval)
				Mval, r200val, c200val = mass_defs.changeMassDefinition(Mval, cval, zval, mdef, mdefout, profile='nfw')


			r500 = nc.M500_to_R500(Mval, zval, cosmo_params_dict['omega_m'], cosmo_params_dict['omega_lambda'], H0)
			theta500 = nc.radius_to_theta_arcm (r500, zval, cosmo_params_dict['omega_m'], cosmo_params_dict['omega_lambda'],  H0)


			ymax, xmax = MAPS[cntr].shape
			xcoord, ycoord = np.arange(0,xmax,1), np.arange(0, ymax,1)
			xcoord_mesh, ycoord_mesh = np.meshgrid(xcoord, ycoord)
			xcen = MAPS[cntr].shape[1] / 2 - 1
			ycen = MAPS[cntr].shape[0] / 2 - 1
			distance_each_pix = np.hypot(xcoord_mesh - xcen, ycoord_mesh - ycen)
			
			for radius in nc.my_range (0., howmany_theta_500*theta500, pixel_size_arcm):
				print radius
				Yout = nc.M500_to_Yx500_arcm2 ((radius+pixel_size_arcm)/theta500,Mass = m_val_fit ,redshift=zval, Om_M= cosmo_params_dict['omega_m'], Om_L=cosmo_params_dict['omega_lambda'], H0=H0)
				Yin = nc.M500_to_Yx500_arcm2 (radius/theta500,Mass=Mval,redshift=m_val_fit , Om_M=cosmo_params_dict['omega_m'], Om_L=cosmo_params_dict['omega_lambda'], H0=H0)
				pixel_area = pixel_size_arcm**2.
				#put in the values of Yout - Yin in diffferent pixels
				number_pixels_between_radius = MAPS[cntr, (distance_each_pix >= radius/pixel_size_arcm) &  (distance_each_pix < (radius+pixel_size_arcm)/pixel_size_arcm)].shape[0]
				MAPS[cntr, (distance_each_pix >= radius/pixel_size_arcm) &  (distance_each_pix <= (radius+pixel_size_arcm)/pixel_size_arcm)]    =       (Yout-Yin) / number_pixels_between_radius / pixel_area
			
		#return downsampled map now
		if ds_fac>1:
			MAPS_ds = []
			for currMAP in MAPS:
				MAPS_ds.append( self.downsample_map(currMAP, ds_fac) )
			MAPS = np.asarray(MAPS_ds)

		#add scatter here
		if random_seed <>None:
			np.random.seed(random_seed)
		if sigma_logY_scatter<>None:
			MAPS_WITH_SCATTER = []
			for currMAP in MAPS:
				#scatterval = np.random.normal(0., sigma_logY_scatter)
				#scatterval = np.random.normal(0., sigma_logY_scatter)
				#scatterval = np.random.normal(0., sigma_logY_scatter)
				scatterval = np.random.normal(0., sigma_logY_scatter)
				#from IPython import embed; embed()
				currMAP[currMAP==0.] = self.small_number/1e5
				lnyflux_integrated = np.log( currMAP ) + scatterval
				yflux_integrated = np.exp(lnyflux_integrated)

				MAPS_WITH_SCATTER.append(yflux_integrated)

			MAPS = np.asarray( MAPS_WITH_SCATTER )

		#convert y-maps into CMB-T at the specified nu
		if convert_to_T:

			def fn_get_f_x(nu = 150e9):

				h=6.62607004e-34 #Planck constant in m2 kg / s
				k_B=1.38064852e-23 #Boltzmann constant in m2 kg s-2 / K-1

				x = h * nu / (k_B * self.inidic['T_cmb'])
				g_nu = (np.exp(x) + 1) / (np.exp(x) - 1)

				f_x = (x * g_nu - 4) # +  rela. corrections

				return f_x

			f_x = fn_get_f_x(nu)
			MAPS = MAPS * f_x *  self.inidic['T_cmb']  * 1e6 #maps in uK


		return MAPS



	def fn_get_Arnaud_ymap_from_Nikhel_Gupta(self, mapparams, Mval_arr, zval_arr, dx_finer = 0.5, cosmo_params_dict = None, mass_def = 200., rho_def = 'mean', howmany_theta_500 = 4.5, convert_to_T = 1, nu = 150e9, sigma_logY_scatter = 0.2, random_seed = None):

		#how many clusters
		totalclus = len(Mval_arr)

		if (1): #do this on a finer grid
			nx, ny, dx, dy = mapparams
			ds_fac = int(dx/dx_finer)
			nx_finer = ds_fac * nx

		#map details
		pixel_size_arcm = dx_finer
		pixels = nx_finer/2
		MAPS = np.zeros((totalclus, 2*pixels,2*pixels))

		#cosmology details
		if cosmo_params_dict == None:
			cosmo_params_dict = {}
			cosmo_params_dict['h'] = 0.6774
			cosmo_params_dict['omega_m'] = 0.307320
			cosmo_params_dict['omega_lambda'] = 0.692680
			cosmo_params_dict['z_lss'] = 1089.
			cosmo_params_dict['T_cmb'] = 2.725
		H0 = cosmo_params_dict['h'] * 100.


		for cntr, (Mval, zval) in enumerate(zip(Mval_arr, zval_arr)):

			if mass_def == 200.: #convert to M500 since Arnaud profile is for M500
				mdef = '%s%s' %(int(mass_def), rho_def[0])
				mdefout = '500%s' %(rho_def[0])
				cval = concentration.concentration(Mval, mdef, zval)
				Mval, r200val, c200val = mass_defs.changeMassDefinition(Mval, cval, zval, mdef, mdefout, profile='nfw')


			r500 = nc.M500_to_R500(Mval, zval, cosmo_params_dict['omega_m'], cosmo_params_dict['omega_lambda'], H0)
			theta500 = nc.radius_to_theta_arcm (r500, zval, cosmo_params_dict['omega_m'], cosmo_params_dict['omega_lambda'],  H0)


			ymax, xmax = MAPS[cntr].shape
			xcoord, ycoord = np.arange(0,xmax,1), np.arange(0, ymax,1)
			xcoord_mesh, ycoord_mesh = np.meshgrid(xcoord, ycoord)
			xcen = MAPS[cntr].shape[1] / 2 - 1
			ycen = MAPS[cntr].shape[0] / 2 - 1
			distance_each_pix = np.hypot(xcoord_mesh - xcen, ycoord_mesh - ycen)

			for radius in nc.my_range (0., howmany_theta_500*theta500, pixel_size_arcm):
				Yout = nc.M500_to_Yx500_arcm2 ((radius+pixel_size_arcm)/theta500,Mass=Mval,redshift=zval, Om_M= cosmo_params_dict['omega_m'], Om_L=cosmo_params_dict['omega_lambda'], H0=H0)
				Yin = nc.M500_to_Yx500_arcm2 (radius/theta500,Mass=Mval,redshift=zval, Om_M=cosmo_params_dict['omega_m'], Om_L=cosmo_params_dict['omega_lambda'], H0=H0)
				pixel_area = pixel_size_arcm**2.
				#put in the values of Yout - Yin in diffferent pixels
				number_pixels_between_radius = MAPS[cntr, (distance_each_pix >= radius/pixel_size_arcm) &  (distance_each_pix < (radius+pixel_size_arcm)/pixel_size_arcm)].shape[0]
				MAPS[cntr, (distance_each_pix >= radius/pixel_size_arcm) &  (distance_each_pix <= (radius+pixel_size_arcm)/pixel_size_arcm)]    =       (Yout-Yin) / number_pixels_between_radius / pixel_area
			
		#return downsampled map now
		if ds_fac>1:
			MAPS_ds = []
			for currMAP in MAPS:
				MAPS_ds.append( self.downsample_map(currMAP, ds_fac) )
			MAPS = np.asarray(MAPS_ds)

		#add scatter here
		if random_seed <>None:
			np.random.seed(random_seed)
		if sigma_logY_scatter<>None:
			MAPS_WITH_SCATTER = []
			for currMAP in MAPS:
				#scatterval = np.random.normal(0., sigma_logY_scatter)
				#scatterval = np.random.normal(0., sigma_logY_scatter)
				#scatterval = np.random.normal(0., sigma_logY_scatter)
				scatterval = np.random.normal(0., sigma_logY_scatter)
				#from IPython import embed; embed()
				currMAP[currMAP==0.] = self.small_number/1e5
				lnyflux_integrated = np.log( currMAP ) + scatterval
				yflux_integrated = np.exp(lnyflux_integrated)

				MAPS_WITH_SCATTER.append(yflux_integrated)

			MAPS = np.asarray( MAPS_WITH_SCATTER )

		#convert y-maps into CMB-T at the specified nu
		if convert_to_T:

			def fn_get_f_x(nu = 150e9):

				h=6.62607004e-34 #Planck constant in m2 kg / s
				k_B=1.38064852e-23 #Boltzmann constant in m2 kg s-2 / K-1

				x = h * nu / (k_B * self.inidic['T_cmb'])
				g_nu = (np.exp(x) + 1) / (np.exp(x) - 1)

				f_x = (x * g_nu - 4) # +  rela. corrections

				return f_x

			f_x = fn_get_f_x(nu)
			
			MAPS = MAPS * f_x *  self.inidic['T_cmb']  * 1e6 #maps in uK

		return MAPS

	####################################################################################################
	####################################################################################################
	####################################################################################################
	#two halo term stuff
	def fn_get_halo_bias_terms(self, delta):
		"""
		http://iopscience.iop.org/article/10.1088/0004-637X/724/2/878/pdf
		Table 2
		"""
		y = np.log10(delta)

		A = 1 + 0.24 * y * np.exp(- (4/y)**4. )
		alpha = 0.44 * y - 0.88
		B = 0.183
		beta = 1.5
		C = 0.019 + 0.107 * y + 0.19 * np.exp(- (4/y)**4. )
		gamma = 2.4

		return A, B, C, alpha, beta, gamma


	def fn_get_halo_bias(self, M, z, fb = None, delta = 200.):

		"""
		http://iopscience.iop.org/article/10.1088/0004-637X/724/2/878/pdf
		Eq. 5
		"""

		nu = fb.nu_M(M, z)
		delta_c =fb.delta_c()

		A, B, C, alpha, beta, gamma = self.fn_get_halo_bias_terms(delta)
		t1 = 1.
		t2 = A * ( (nu**alpha) / ( nu**alpha + delta_c**alpha) ) 
		t3 = B * nu**beta
		t4 = C * nu**gamma

		b_nu = t1 - t2 + t3 + t4

		return b_nu, nu

	def fn_kappa_two_halo_term_helper_function(self, theta, clus_identifier, param_dict, opfname, fb = None, cosmo = None, lmin = 0, lmax = 30000, delta_els = 10):

		theta = np.radians(theta)

		import astropy.constants as const
		import scipy.special as special
		import scipy.integrate as integrate

		if not os.path.exists(opfname):
			kappa_two_halo_int_term_dic = {}
		else:
			kappa_two_halo_int_term_dic = pickle.load(gzip.open(opfname, 'rb'))
		for cicnt, CI in enumerate(clus_identifier):

			keyname = tuple(CI)

			if keyname in kappa_two_halo_int_term_dic:
				print '\n\t', keyname, ' already comupted'
				continue

			##from IPython import embed; embed()

			print keyname

			z = CI[2]

			rho_crit_z = cosmo.critical_density(z).to('M_sun/Mpc3').value
			rho_matter_z = rho_crit_z * cosmo.Om(z)

			z_lss = param_dict['z_lss']
			D_L = cosmo.comoving_distance(z)/(1.+z)
			D_S = cosmo.comoving_distance(z_lss)/(1.+z_lss)

			kappa_two_halo_int_term = np.zeros(theta.shape)
			for tcnt, theta_val in enumerate( theta ):

				###if tcnt%1000 == 0: print tcnt

				###print theta_val, lmin, lmax
				def fn_for_int_for_simps(el, bessel_order = 0):
					t1 = el * special.jn(bessel_order, el * theta_val) /2. /np.pi
					kl = el / (1.+z) / D_L.value
					t3 = fb.cambmatterpower.P(z, kl)
					return t1 * t3

				els = np.arange(lmin, lmax)
				kappa_two_halo_int_term[tcnt] = integrate.simps( fn_for_int_for_simps(els), x=els , dx = delta_els)
				#print kappa_two_halo[tcnt];sys.exit()

			kappa_two_halo_int_term_dic[keyname] = kappa_two_halo_int_term

			pickle.dump(kappa_two_halo_int_term_dic, gzip.open(opfname, 'wb'), protocol = 2)


	def fn_kappa_two_halo_term(self, theta, M, z, param_dict, fb = None, cosmo = None, kappa_two_halo_int_term = None, lmin = 0, lmax = 30000, delta_els = 10):

		###print '\n\tlmin = %s; lmax = %s' %(lmin, lmax)

		theta = np.radians(theta)
		import astropy.constants as const
		import scipy.special as special
		import scipy.integrate as integrate

		b_M_z, nu = self.fn_get_halo_bias(M, z, fb = fb)
		#print M, z
		#print b_M_z;sys.exit()

		rho_crit_z = cosmo.critical_density(z).to('M_sun/Mpc3').value
		rho_matter_z = rho_crit_z * cosmo.Om(z)

		z_lss = param_dict['z_lss']
		D_L = cosmo.comoving_distance(z)/(1.+z)
		D_S = cosmo.comoving_distance(z_lss)/(1.+z_lss)
		'''
		#D_LS = (D_S - D_L)/(1.+z_lss)
		D_LS = (D_S - D_L)/(1.+z)
		'''

		t1 = (1. / (1+z_lss))
		d_m2, d_m1 = cosmo.comoving_distance(z_lss), cosmo.comoving_distance(z)
		t2 = d_m2# * np.sqrt( 1 +  (omega_k * d_m1**2 / d_h**2.)  )
		t3 = d_m1# * np.sqrt( 1 +  (omega_k * d_m2**2 / d_h**2.)  )

		D_LS = t1 * (t2 - t3)

		#print D_LS;sys.exit()

		sigma_c = (((const.c.cgs**2.)/(4.*np.pi*const.G.cgs))*(D_S/(D_L*D_LS))).to('M_sun/Mpc2')

		t2_nr = rho_matter_z * b_M_z
		t2_dr = (1. + z)**3. * sigma_c.value * D_L.value**2.
		t2 = t2_nr / t2_dr

		if not self.is_seq(kappa_two_halo_int_term):

			kappa_two_halo = np.zeros(theta.shape)
			for tcnt, theta_val in enumerate( theta ):

				#if tcnt%1000 == 0:
				#	print tcnt, len(theta)

				"""
				def fn_for_int(el, bessel_order = 0):
					t1 = el * special.jn(bessel_order, el * theta_val) /2. /np.pi

					#picking from the dumped file
					kl = el / (1.+z) / D_L.value
					t3 = fb.cambmatterpower.P(z, kl)

					return t1 * t2 * t3

				kappa_two_halo[tcnt] = integrate.quad(fn_for_int, lmin, lmax)[0]
				print kappa_two_halo[tcnt]
				"""

				###print theta_val, lmin, lmax
				def fn_for_int_for_simps(el, bessel_order = 0):
					t1 = el * special.jn(bessel_order, el * theta_val) /2. /np.pi
					kl = el / (1.+z) / D_L.value
					t3 = fb.cambmatterpower.P(z, kl)
					return t1 * t2 * t3
					
				els = np.arange(lmin, lmax)
				kappa_two_halo[tcnt] = integrate.simps( fn_for_int_for_simps(els), x=els , dx = delta_els)

		else:

			kappa_two_halo = kappa_two_halo_int_term * t2

		return kappa_two_halo

	####################################################################################################
	####################################################################################################
	####################################################################################################

	def gaussian_fitting_func(self, p, p0, X, Y, MAP, lbounds = None, ubounds = None, fixed = None, return_fit = 0):
		'''
		if hasattr(fixed, '__len__'):
			p[fixed] = p0[fixed]

		if hasattr(lbounds, '__len__'):
			linds = abs(p)<abs(lbounds)
			p[linds] = lbounds[linds]

		if hasattr(ubounds, '__len__'):
			uinds = abs(p)>abs(ubounds)
			p[uinds] = ubounds[uinds]
		'''

		def fn_gaussian(p, xp, yp):
			#rota = np.radians(p[6])
			#xp = xp * np.cos(rota) - yp * np.sin(rota)
			#yp = xp * np.sin(rota) + yp * np.cos(rota)
			#g = p[0]+
			g = p[1]*np.exp( -(((p[2]-xp)/p[4])**2+ ((p[3]-yp)/p[5])**2)/2.)

			return g

		if not return_fit:
			return np.ravel(fn_gaussian(p, X, Y) - MAP)
		else:
			return fn_gaussian(p, X, Y)		

	def fitting_func(self, p, p0, X, cluster_mass, DATA = None, lbounds = None, ubounds = None, fixed = None, return_fit = 0, commit_20161221_change = 0):
		if hasattr(fixed, '__len__'):
			p[fixed] = p0[fixed]

		if hasattr(lbounds, '__len__'):
			linds = abs(p)<abs(lbounds)
			p[linds] = lbounds[linds]

		if hasattr(ubounds, '__len__'):
			uinds = abs(p)>abs(ubounds)
			p[uinds] = ubounds[uinds]

		if commit_20161221_change:
			fitfunc = lambda p, X: p[2] + ((X - p[1])) ** 2. / p[0] ### added on 20161221 for better fitting
		else:
			fitfunc = lambda p, X: ((X - cluster_mass)) ** 2. / p[0]
		if not return_fit:
			return fitfunc(p, X) - DATA
		else:
			return fitfunc(p, X)

	def fn_calc_delta_chisq_one_width(self, Marr, del_chi_sq_arr, cluster_mass, perform_ip = 0, fit_parabola = 1, testing = 0):

		default_width = 1.
		######################################################
		commit_20161221_change = 0
		if 1==1:#commit_20161221_change:
			### added on 20161221 for better fitting
			reqdinds = np.where( (del_chi_sq_arr<20.) )
			Marr = Marr[reqdinds]
			del_chi_sq_arr = del_chi_sq_arr[reqdinds]
		######################################################

		if perform_ip:
			x, y = Marr, del_chi_sq_arr
			x_ip = np.arange(0.,max(x),0.1)

			poly_deg = 5
			poly_fit_eq=np.polyfit(x, y, poly_deg)
			poly_fit=np.poly1d(poly_fit_eq)

			y_ip = poly_fit(x_ip)

			#from pylab import *
			plot(x_ip, y_ip, 'k');plot(Marr, del_chi_sq_arr, 'ro');ylim(0,20.)
			show();quit()

			Marr, del_chi_sq_arr = x_ip, y_ip

		#from pylab import *
		#plot(Marr, del_chi_sq_arr);show();quit()

		#del_chi_sq_arr[del_chi_sq_arr<0.] = 0.
		#del_chi_sq_arr[del_chi_sq_arr>500.] = np.nan
		#width of L at \Delta\chi^2<=1.
		value_for_error = 1.
		interp_type = 'linear'
		linds, uinds = np.where(Marr<cluster_mass)[0], np.where(Marr>cluster_mass)[0]

		minreqdindsforfitting = 2
		if len(linds) < minreqdindsforfitting  or len(uinds) < minreqdindsforfitting: #5 is just an arbitray value. If it either hits the minimum or maximum there will not be equal number of linds, uinds
			width = default_width
			return width

		fninterp = interp1d(del_chi_sq_arr[linds], Marr[linds], kind = interp_type, bounds_error = 0, fill_value = 0.)
		l_err = fninterp(value_for_error)
		fninterp = interp1d(del_chi_sq_arr[uinds], Marr[uinds], kind = interp_type, bounds_error = 0, fill_value = 0.)
		u_err = fninterp(value_for_error)
		#width = (_err - l_err)/2.
		width = 0.001
		default_width = 0.001
		#print width
		if width == 0.:
			width = .05

		if fit_parabola:

			import scipy.optimize as optimize

			p0 = np.asarray([width,Marr[np.argmin(del_chi_sq_arr)], min(del_chi_sq_arr)])
			lbounds = np.asarray([p0[0] * 1e-4, 0., min(del_chi_sq_arr)])
			ubounds = np.asarray([p0[0] * 10., 10., min(del_chi_sq_arr)*1.1])
			fixed = np.asarray([False, False,True])
			#inds_to_fit = np.where(del_chi_sq_arr<=10.)[0]
			#p1, success = optimize.leastsq(fitting_func, p0, args=(p0, Marr[inds_to_fit], del_chi_sq_arr[inds_to_fit], lbounds, ubounds))
			if commit_20161221_change:
				p1, success = optimize.leastsq(self.fitting_func, p0, args=(p0, Marr,  cluster_mass, del_chi_sq_arr, lbounds, None, fixed))#, ubounds))
			else:
				p1, success = optimize.leastsq(self.fitting_func, p0, args=(p0, Marr,  cluster_mass, del_chi_sq_arr, lbounds))#, ubounds))

			x_ip = np.arange(min(Marr),max(Marr),cluster_mass/100.)
			y_ip = self.fitting_func(p1,p1,x_ip, cluster_mass,return_fit = 1)

			if commit_20161221_change:
				linds, uinds = np.where(x_ip<p1[1])[0], np.where(x_ip>p1[1])[0] #added on 20161221 for better fitting
			else:
				linds, uinds = np.where(x_ip<=cluster_mass)[0], np.where(x_ip>=cluster_mass)[0]

			value_for_error = 1.
			if 1==1:#sci_interp:
				interp_type = 'linear'

				fninterp = interp1d(y_ip[linds], x_ip[linds], kind = interp_type, bounds_error = 0, fill_value = 0.)
				l_err = fninterp(value_for_error)
				fninterp = interp1d(y_ip[uinds], x_ip[uinds], kind = interp_type, bounds_error = 0, fill_value = 0.)
				u_err = fninterp(value_for_error)
			else:
				l_err = np.interp(value_for_error,y_ip[linds], x_ip[linds])
				u_err = np.interp(value_for_error,y_ip[uinds], x_ip[uinds])
			width = (u_err - l_err)/2.#/cluster_mass
			if width == 0:
				width = default_width#small_number
			
			#width = (u_err - cluster_mass)#/2.

			if testing:#0==1:
				clf()
				y_ip_ori = self.fitting_func(p0,p0,x_ip, cluster_mass,return_fit = 1)

				#from pylab import *
				clf()
				plot(x_ip, y_ip, 'k', lw= 3.);plot(Marr, del_chi_sq_arr, 'ro');ylim(0,10.)
				#plot(x_ip, y_ip_ori, 'r')
				plot(l_err,1.,'go')
				plot(u_err,1.,'go')
				#print u_err - l_err
				hlines(1.,0.,max(Marr)*2.,color='k',linestyle='solid')

				show()#;sys.exit()

		return width


	def fn_get_det_significance(self, Marr, res_arr, massval_for_L_peak = None, two_lnL = 1):
		zero_ind = np.where(Marr == 0)[0]

		Marr = np.asarray(Marr)
		res_arr = np.asarray(res_arr)

		if len(zero_ind) > 0:
			if massval_for_L_peak == None:
				maxval = max(res_arr)
			else:
				L_peak_ind = np.argmin(abs(Marr - massval_for_L_peak))
				maxval = res_arr[L_peak_ind]

				#print maxval, max(res_arr)

			diffval = maxval - res_arr[zero_ind][0]
			#print maxval, res_arr[zero_ind][0]

			if diffval == 0.: 
				det_sig = 0.
			else:
				if two_lnL:
					det_sig =  np.sqrt( diffval )
				else:
					det_sig =  np.sqrt( 2 * ( diffval) ) #factor of 2 becuase we are using lnL and not 2lnL

			if np.isnan(det_sig):
				det_sig = 0.

		else:
			#in this try to interpolate likelihood curve and get the SNR
			x, y = Marr, res_arr
			delM = np.diff(Marr)[0]/10.
			x_ip = np.arange(0.,max(x),delM)

			poly_deg = 3.
			poly_fit_eq=np.polyfit(x, y, poly_deg)
			poly_fit=np.poly1d(poly_fit_eq)

			y_fit = poly_fit(x_ip)
			y_guess_val = np.interp(0., x_ip, y_fit)

			if massval_for_L_peak == None:
				maxval = max(res_arr)
			else:
				L_peak_ind = np.argmin(abs(Marr - massval_for_L_peak))[0]
				maxval = res_arr[L_peak_ind]

			det_sig =  np.sqrt( 2 * ( maxval - y_guess_val) )

			#y_fit = fn_normalize(y_fit)#, arrmin = meanarr[zero_ind])
			#plot(x_ip, y_fit, 'm', lw = 0.5, label = 'intrp for SNR')
			#legend(loc=1, fancybox = 1)

		return det_sig


class tSZ_fitting():

	def __init__(self):
		self.tSZ_beta = 1. #1.
		self.theta_core = 0.75 #arcmins
		self.sptsz_cat_file = 'data/2500d_cluster_sample_fiducial_cosmology.fits'

	def is_seq(self,o):
		return hasattr(o, '__len__')

	def fn_guess_Tsz_fromega_mass(self, x_to_guess_y, guess_Tsz = 1, guess_M500 = 0, poly_deg = 2, nu = 150e9):

		assert guess_Tsz <> guess_M500

		clusterstuff = cluster_stuff()

		#x_to_guess_y should be M500 or Tsz to guess the other based on SPT-SZ fitting

		sptsz_cat_file = self.sptsz_cat_file
		#sptsz = pyfits.getdata(sptsz_cat_file)
		sptszfits = fits.open(sptsz_cat_file)[1]
		sptsz = sptszfits.data
		ysz = np.asarray( map(lambda x: x[5], sptsz) )
		m500 = np.asarray( map(lambda x: x[11], sptsz) )

		#print min(m500), max(m500);quit()

		ysz_Tsz_conv_fac = clusterstuff.compton_y_to_delta_Tcmb(nu)
		Tsz = ysz * float( ysz_Tsz_conv_fac ) * 1e6

		#Tsz vs. M500 ploynomial fit
		if guess_Tsz:
			y, x = Tsz[abs(m500)>0.],m500[abs(m500)>0.]
			minM = 1e11
			x_ip = np.arange(minM,5e15,2*minM)

		elif guess_M500:
			x, y = Tsz[abs(m500)>0.],m500[abs(m500)>0.]
			x_ip = np.arange(min(Tsz),0.)

		poly_fit_eq=np.polyfit(x, y, poly_deg)
		poly_fit=np.poly1d(poly_fit_eq)

		y_fit = poly_fit(x_ip)

		y_guess_val = np.interp(x_to_guess_y, x_ip, y_fit)

		'''
		plot(x,y,'r.');plot(x_ip, y_fit,'k')
		plot(x_to_guess_y,y_guess_val,'go')
		show();quit()
		'''

		return y_guess_val #either M500 or Tsz guessed based on the other

	#def fn_subtract_tSZ(self, mapparams, MAP_150, MAP_90 = None, smooth_150ghz_map = 1, use_beam_150 = 1, use_beam_90 = 1, exp_beam_150 = 1.2, exp_beam_90 = 1.8, return_SZ_150ghz_map = 0, factor = None):
	def fn_subtract_tSZ(self, mapparams, MAP_150, MAP_90 = None, smooth_150ghz_map = 1, use_beam_150 = 1, use_beam_90 = 1, exp_beam_150 = 1.2, exp_beam_90 = 1.7, return_SZ_150ghz_map = 0, factor = None, Bl_90 = None, Bl_150 = None, calib_facs = None):

		clusterstuff = cluster_stuff()
		sims = simulations()

		try:
			use_beam_150 = use_beam_90 = sims.inidic['use_beam']
		except:
			pass

		if Bl_90 == None:
			if use_beam_150 == 1:
				sims.exp_beam = exp_beam_150
		if Bl_150 == None:
			if use_beam_90 == 1:
				sims.exp_beam_90ghz = exp_beam_90

		if not self.is_seq(MAP_90):
			print 'Not implemented yet'
			return None

		if factor == None:
			ysz_Tsz_conv_fac_150 = clusterstuff.compton_y_to_delta_Tcmb(150e9) #conversion factor for 150 GHz
			ysz_Tsz_conv_fac_90 = clusterstuff.compton_y_to_delta_Tcmb(90e9) #conversion factor for 90 GHz

			factor = ysz_Tsz_conv_fac_90/ysz_Tsz_conv_fac_150

		if calib_facs <> None:
			MAP_150 = MAP_150 * calib_facs[150]
			MAP_90 = MAP_90 * calib_facs[90]

		"""
		def fn_correct_deconv_beam(MAP, Bl = None, Bl_deconv = None):
			if Blinv == None:
				bad = np.where(Bl==0.)
				Bl[bad] = 1.
				Bl_deconv = 1./Bl
		    	Bl_deconv[bad] = 0.
		    else:
			    highly_filtered = np.where((Bl_deconv > 1.0e8) | (Bl_deconv < 0) )
			    Bl_deconv[highly_filtered] = 0.0

			return np.fft.ifft2( np.fft.fft2(MAP) * Bl_deconv).real
		"""

		ORI_MAP_150 = np.copy(MAP_150)
		#smooth the 150 ghz map using 90 ghz beam
		if smooth_150ghz_map:

			#using 2d beams
			if Bl_90 == None:
				Bl_90 = sims.fn_beam_stuff(mapparams, use_beam = use_beam_90, nu = 90, return_beam = 1)
			if Bl_150 == None:
				Bl_150 = sims.fn_beam_stuff(mapparams, use_beam = use_beam_150, nu = 150, return_beam = 1)

			Bl_tSZ_for_150 = Bl_90 / Bl_150

			#imshow(Bl_tSZ_for_150);colorbar();show();sys.exit()

			'''
			Bl_150 = sims.fn_beam_stuff(mapparams, use_beam = 2, nu = 150, return_beam = 1)
			Bl_tSZ_for_150 = Bl_90 / Bl_150
			Bl_tSZ_for_150[Bl_tSZ_for_150 == np.inf] = 0.
			'''

			if len(Bl_tSZ_for_150.shape)>2:
				Bl_tSZ_for_150 = Bl_tSZ_for_150[0]

			Bl_tSZ_for_150[Bl_tSZ_for_150 == np.inf] = 0.
			Bl_tSZ_for_150[np.where(np.isnan(Bl_tSZ_for_150))] = 0.

			'''
			subplot(131);imshow(np.fft.fftshift( Bl_90[0] ));colorbar()
			subplot(132);imshow(np.fft.fftshift( Bl_150[0] ));colorbar()
			subplot(133);imshow(np.fft.fftshift( Bl_tSZ_for_150 ));colorbar();show();quit()
			'''
				
			'''
			#using 1d Bls
			import healpy as H
			els_1d = np.arange(30000)
			Bl_1d_150 = H.gauss_beam(np.radians( sims.exp_beam/ 60.), lmax = max(els_1d))
			Bl_1d_90 = H.gauss_beam(np.radians( sims.exp_beam_90ghz/ 60.), lmax = max(els_1d))

			Bl_1d_tSZ = Bl_1d_90 / Bl_1d_150
			tsz_spl_Bl = np.asarray( [els_1d, Bl_1d_tSZ] )

			Bl_tSZ_for_150_2d = sims.el1d_to_EL2D(tsz_spl_Bl,mapparams)
			Bl_tSZ_for_150 = Bl_tSZ_for_150_2d
			'''

			#Bl_tSZ_for_150 = sims.fn_beam_stuff(mapparams, use_beam = 1, nu = 150, return_beam = 1, exp_beam = 1.2)[0]
			#imshow(Bl_tSZ_for_150_2d - Bl_tSZ_for_150);colorbar();show();quit()

			#subplot(121);imshow(np.fft.fftshift( Bl_tSZ_for_150 ));colorbar()
			#subplot(122);imshow(np.fft.fftshift( Bl_tSZ_for_150_1d ));colorbar();show();quit()

			#plot(els_1d, Bl_1d_90,'b');plot(els_1d, Bl_1d_150,'r');plot(els_1d, Bl_1d_tSZ,'k');show();quit()

			MAP_150 = np.fft.ifft2( np.fft.fft2(MAP_150) * Bl_tSZ_for_150 ).real
			#MAP_150 = fn_correct_deconv_beam(MAP_150, Bl_tSZ_for_150)

		tSZ_free_map = ( (MAP_150 * factor) - MAP_90 ) / (factor - 1.)


		if 0 == 1:

			ds_fac = 2
			M1 = sims.downsample_map(MAP_90,ds_fac)
			M2 = sims.downsample_map(MAP_150,ds_fac)
			M3 = sims.downsample_map(tSZ_free_map,ds_fac)

			nx, ny, dx, dy = mapparams
			nx_mod, ny_mod = M1.shape
			dx_mod, dy_mod = dx * ds_fac, dy * ds_fac
			mapparams_mod = [nx_mod, ny_mod, dx_mod, dy_mod]

			cls_90 = sims.fn_plot_pow_spec(mapparams_mod,[M1])[0][0]
			cls_150 = sims.fn_plot_pow_spec(mapparams_mod,[M2])[0][0]
			cls_tszfree = sims.fn_plot_pow_spec(mapparams_mod,[M3])[0][0]

			Dl_fac = cls_90[:,0] * (cls_90[:,0] + 1)/2/np.pi * 1e12
			clf()
			loglog(cls_90[:,0], Dl_fac * cls_90[:,1],'ro-');loglog(cls_150[:,0], Dl_fac * cls_150[:,1],'go-')
			loglog(cls_tszfree[:,0], Dl_fac * cls_tszfree[:,1],'ro-')
			xlim(1e2,3e4);show();quit()
			
			subplot(131);imshow(M1, vmin=-100, vmax=100);colorbar()
			subplot(132);imshow(M2, vmin=-100, vmax=100);colorbar()
			subplot(133);imshow(M3, vmin=-100, vmax=100);colorbar();show();quit()

		if return_SZ_150ghz_map:
			SZ_150ghz_map = (MAP_150 - tSZ_free_map)
			#SZ_150ghz_map = np.fft.ifft2( np.fft.fft2(SZ_150ghz_map)  ).real
			
		'''
		subplot(321);imshow(MAP_90);colorbar()
		subplot(322);imshow(ORI_MAP_150);colorbar()
		subplot(323);imshow(MAP_90);colorbar()
		subplot(324);imshow(MAP_150);colorbar()
		subplot(325);imshow(tSZ_free_map);colorbar()
		subplot(326);imshow(SZ_150ghz_map);colorbar();show();quit()
		'''

		if not return_SZ_150ghz_map:
			return tSZ_free_map
		else:
			return tSZ_free_map, SZ_150ghz_map

	def fn_subtract_tSZ_v1(self, mapparams, MAP_150, MAP_90 = None, smooth_150ghz_map = 1, use_beam_150 = 1, use_beam_90 = 1, exp_beam_150 = 1.2, exp_beam_90 = 1.8, return_SZ_150ghz_map = 0):


		clusterstuff = cluster_stuff()
		sims = simulations()

		try:
			use_beam_150 = use_beam_90 = sims.inidic['use_beam']
		except:
			pass

		if use_beam_150 == 1:
			sims.exp_beam = exp_beam_150
		if use_beam_90 == 1:
			sims.exp_beam_90ghz = exp_beam_90


		if not self.is_seq(MAP_90):
			print 'Not implemented yet'
			return None

		ysz_Tsz_conv_fac_150 = clusterstuff.compton_y_to_delta_Tcmb(150e9) #conversion factor for 150 GHz
		ysz_Tsz_conv_fac_90 = clusterstuff.compton_y_to_delta_Tcmb(90e9) #conversion factor for 90 GHz

		factor = ysz_Tsz_conv_fac_90/ysz_Tsz_conv_fac_150


		#embed()

		ORI_MAP_150 = np.copy(MAP_150)
		#smooth the 150 ghz map using 90 ghz beam
		if smooth_150ghz_map:

			#using 2d beams
			Bl_90 = sims.fn_beam_stuff(mapparams, use_beam = use_beam_90, nu = 90, return_beam = 1)
			Bl_150 = sims.fn_beam_stuff(mapparams, use_beam = use_beam_150, nu = 150, return_beam = 1)
			Bl_tSZ_for_150 = Bl_90 / Bl_150

			'''
			Bl_150 = sims.fn_beam_stuff(mapparams, use_beam = 2, nu = 150, return_beam = 1)
			Bl_tSZ_for_150 = Bl_90 / Bl_150
			Bl_tSZ_for_150[Bl_tSZ_for_150 == np.inf] = 0.
			'''

			if len(Bl_tSZ_for_150.shape)>2:
				Bl_tSZ_for_150 = Bl_tSZ_for_150[0]

			
			Bl_tSZ_for_150[Bl_tSZ_for_150 == np.inf] = 0.
			Bl_tSZ_for_150[np.where(np.isnan(Bl_tSZ_for_150))] = 0.

			'''
			subplot(131);imshow(np.fft.fftshift( Bl_90[0] ));colorbar()
			subplot(132);imshow(np.fft.fftshift( Bl_150[0] ));colorbar()
			subplot(133);imshow(np.fft.fftshift( Bl_tSZ_for_150 ));colorbar();show();quit()
			'''
				
			'''
			#using 1d Bls
			import healpy as H
			els_1d = np.arange(30000)
			Bl_1d_150 = H.gauss_beam(np.radians( sims.exp_beam/ 60.), lmax = max(els_1d))
			Bl_1d_90 = H.gauss_beam(np.radians( sims.exp_beam_90ghz/ 60.), lmax = max(els_1d))

			Bl_1d_tSZ = Bl_1d_90 / Bl_1d_150
			tsz_spl_Bl = np.asarray( [els_1d, Bl_1d_tSZ] )

			Bl_tSZ_for_150_2d = sims.el1d_to_EL2D(tsz_spl_Bl,mapparams)
			Bl_tSZ_for_150 = Bl_tSZ_for_150_2d
			'''

			#Bl_tSZ_for_150 = sims.fn_beam_stuff(mapparams, use_beam = 1, nu = 150, return_beam = 1, exp_beam = 1.2)[0]
			#imshow(Bl_tSZ_for_150_2d - Bl_tSZ_for_150);colorbar();show();quit()

			#subplot(121);imshow(np.fft.fftshift( Bl_tSZ_for_150 ));colorbar()
			#subplot(122);imshow(np.fft.fftshift( Bl_tSZ_for_150_1d ));colorbar();show();quit()

			#plot(els_1d, Bl_1d_90,'b');plot(els_1d, Bl_1d_150,'r');plot(els_1d, Bl_1d_tSZ,'k');show();quit()

			MAP_150 = np.fft.ifft2( np.fft.fft2(MAP_150) * Bl_tSZ_for_150 ).real

			#imshow(MAP_150);colorbar();show();quit()

		tSZ_free_map = ( (MAP_150 * factor) - MAP_90 ) / (factor - 1.)

		if return_SZ_150ghz_map:
			SZ_150ghz_map = (MAP_150 - tSZ_free_map)
			#SZ_150ghz_map = np.fft.ifft2( np.fft.fft2(SZ_150ghz_map)  ).real
			
		'''
		subplot(321);imshow(MAP_90);colorbar()
		subplot(322);imshow(ORI_MAP_150);colorbar()
		subplot(323);imshow(MAP_90);colorbar()
		subplot(324);imshow(MAP_150);colorbar()
		subplot(325);imshow(tSZ_free_map);colorbar()
		subplot(326);imshow(SZ_150ghz_map);colorbar();show();quit()
		'''

		if not return_SZ_150ghz_map:
			return tSZ_free_map
		else:
			return tSZ_free_map, SZ_150ghz_map

	#from http://www.rzuser.uni-heidelberg.de/~ge6/Programing/convolution.html
	def Convolve(self, image1, image2, MinPad=True, pad=True):
		""" Not so simple convolution """

		#The size of the images:
		r1,c1 = image1.shape
		r2,c2 = image2.shape

		#MinPad results simpler padding,smaller images:
		if MinPad:
			r = r1+r2
			c = c1+c2
		else:
			#if the Numerical Recipies says so:
			r = 2*max(r1,r2)
			c = 2*max(c1,c2)

		#For nice FFT, we need the power of 2:
		if pad:
			pr2 = int(np.log(r)/np.log(2.0) + 1.0 )
			pc2 = int(np.log(c)/np.log(2.0) + 1.0 )
			rOrig = r
			cOrig = c
			r = 2**pr2
			c = 2**pc2
			#end of if pad
	
		#numpy fft has the padding built in, which can save us some steps
		#here. The thing is the s(hape) parameter:
		fftimage = np.fft.fft2(image1,s=(r,c)) * np.fft.fft2(image2,s=(r,c))

		if pad:
			return ((np.fft.ifft2(fftimage)))[:rOrig,:cOrig].real
		else:
			return (np.fft.ifft2(fftimage)).real
		

	def fitting_func(self, p, p0, RADIUS, DATA = None, lbounds = None, ubounds = None, fixed = None, return_fit = 0):
		if hasattr(fixed, '__len__'):
			p[fixed] = p0[fixed]

		if hasattr(lbounds, '__len__'):
			linds = abs(p)<abs(lbounds)
			p[linds] = lbounds[linds]
	
		if hasattr(ubounds, '__len__'):
			uinds = abs(p)>abs(ubounds)
			p[uinds] = ubounds[uinds]

		#print p, lbounds, ubounds, linds, uinds
		fitfunc = lambda p, RADIUS: p[0] * ( 1.0 + (RADIUS/p[1]) ** 2. ) ** (0.5 - (1.5 * p[2]) ) #in uK
		### fitfunc = lambda p, RADIUS: p[0] * ( 1.0 + (RADIUS/p[1]) ** 2. ) ** (0.5 - (3. * p[2]) ) #in uK
		if not return_fit:
			return np.ravel(fitfunc(p, RADIUS) - DATA)
		else:
			return fitfunc(p, RADIUS)


	def fn_tSZ_beta_model(self, DATA, dx, X = None, Y = None, nu = 150): #cutout, dx = ang. resol in arcmins

		if X == None and Y == None:

			nx = ny = len(DATA)
			
			boxsize = nx * dx
			minval, maxval = -boxsize/2,  boxsize/2
			x = y = np.linspace(minval,maxval,int(boxsize/dx))

			X, Y = np.meshgrid(x,y)
			RADIUS = (X ** 2. + Y ** 2.)**0.5

			'''
			MASK = RADIUS**2.<=self.theta_core**2.
			DATA[RADIUS**2.<=self.theta_core**2.] = 0.
			imshow(DATA);colorbar();show();quit()
			'''
		else:
			
			RADIUS = (X ** 2. + Y ** 2.)**0.5

		ret_dic = {}
		#initial guess
		tSZ_beta = 1.
		theta_core = 0.5 #arcmins
		T_0 = np.mean(DATA[RADIUS<=self.theta_core]) #cluster core T_0 (or) y_0 depending on the map type
		#print T_0;quit()

		#some bounds
		lbounds = np.asarray([10*T_0,0.25,0.5])
		ubounds = np.asarray([.1*T_0,1.0,1.2])
		fixed = np.asarray([False,False,False])

		#initial params
		p0 = np.asarray([T_0,self.theta_core,self.tSZ_beta])

		if nu <= 220: #220 GHz channel must be null
			if T_0 > 0.:
				print '\n *** Alert: tSZ cannot be positive for frequencies below 220 GHz. Aborting fitting here and returning default params'
				ret_dic['beta_model'], ret_dic['beta_model_params'], ret_dic['beta_model_success'] = self.fitting_func(p0, p0, RADIUS,return_fit = 1), p0, 0
				return ret_dic

		p1, success = optimize.leastsq(self.fitting_func, p0[:], args=(p0, RADIUS,DATA, lbounds, ubounds,fixed))

		if success:

			ret_dic['beta_model']=self.fitting_func(p1, p0, RADIUS,return_fit = 1)
			ret_dic['beta_model_params'] = p1
			ret_dic['beta_model_success'] = 1

		else:

			print '\n *** Alert: tSZ beta model fitting did not converge. Returning default params'
			ret_dic['beta_model'], ret_dic['beta_model_params'], ret_dic['beta_model_success'] = self.fitting_func(p0, p0, RADIUS,return_fit = 1), p0, 0


		return ret_dic



class noise_sims():

	def __init__(self):

		self.noisemap_file = 'data/500sqdeg_henning/20160518_left_right_noise_map.pkl.gz'
		self.mask_file = 'data/500sqdeg_henning/total_mask_extended.pkl'
		self.map_resol = 0.5 #arcmins
		self.sims = simulations()
		self.do_pol = 1
		#self.noise_cls_file = 'data/500sqdeg_henning/noise_cls.pkl.gz' #noise Cls [ells, TT, QQ, UU]
		self.noise_cls_file = 'data/500sqdeg_henning/noise_cls.pkl.gz_v_JT_20161004' #Updated: noise Cls [ells, TT, QQ, UU] - looks correct

		#spt psd
		self.sptpol_clus_map_dx = 0.25
		self.SPT_PSD_150 = None
		self.SPT_PSD_90 = None

		#self.SPT_PSD_file_90ghz = 'data/100sqdeg_huang/090ghz_psd.sav' #100 sq.deg. NDHUANG
		#self.SPT_PSD_file_150ghz = 'data/100sqdeg_huang/150ghz_psd.sav' #100 sq.deg. NDHUANG
		self.SPT_PSD_file_90ghz = 'data/500sqdeg_bleem/noise_psd_090ghz.pkl.gz' #500 sq.deg. LBLEEM
		self.SPT_PSD_file_150ghz = 'data/500sqdeg_bleem/noise_psd_150ghz.pkl.gz' #500 sq.deg. LBLEEM

		#500sq. deg. Sanjay maps
		self.noise_cls_file_90ghz_SPatil = 'data/sanjay_maps_201705xx/noise_Cls/90ghz/combined_noise_cls_42bundlesused.pkl.gz'
		self.noise_cls_file_150ghz_SPatil = 'data/sanjay_maps_201705xx/noise_Cls/150ghz/combined_noise_cls_42bundlesused.pkl.gz'
		self.noise_cls_file_tszfree_SPatil = 'data/sanjay_maps_201705xx/noise_Cls/tszfree/combined_noise_cls_42bundlesused.pkl.gz'

		#500sq. deg. CR maps
		self.noise_cls_file_90ghz_CR = 'data/CR_maps_20170910/noise_Cls/90ghz/combined_noise_cls_89bundlesused.pkl.gz'
		self.noise_cls_file_150ghz_CR = 'data/CR_maps_20170910/noise_Cls/150ghz/combined_noise_cls_89bundlesused.pkl.gz'
		self.noise_cls_file_tszfree_CR = 'data/CR_maps_20170910/noise_Cls/tszfree/combined_noise_cls_89bundlesused.pkl.gz'

		#3g july 2018 ilc curves
		self.noise_cls_file_ilc_3g_july_2018 = 'data/CR_maps_20170910/noise_Cls/tszfree/combined_noise_cls_89bundlesused.pkl.gz'

		#Planck LGMCA maps
		self.noise_cls_file_lgmca = 'data/CR_maps_20170910/noise_Cls/lgmca/lgmca_noise_cls.pkl.gz'

		#500sq. deg. LBleem maps
		self.noise_cls_file_90ghz_LBleem = 'data/500sqdeg_bleem/noise_psd_090ghz.pkl.gz' #500 sq.deg. LBLEEM
		self.noise_cls_file_150ghz_LBleem = 'data/500sqdeg_bleem/noise_psd_150ghz.pkl.gz' #500 sq.deg. LBLEEM

		#500sq. deg. Henning maps
		self.noise_cls_file_90ghz_500sqdeg_JHenning = 'data/500sqdeg_henning/noise_cls_500sqdeg/90ghz/combined_noise_cls.pkl.gz'
		self.noise_cls_file_150ghz_500sqdeg_JHenning = 'data/500sqdeg_henning/noise_cls_500sqdeg/150ghz/combined_noise_cls.pkl.gz'

		#self.noise_cls_file_90ghz_100sqdeg = 'data/100sqdeg_huang/ndhuang_bundles/90ghz/combined_noise_cls.pkl.gz'
		#self.noise_cls_file_150ghz_100sqdeg = 'data/100sqdeg_huang/ndhuang_bundles/150ghz/combined_noise_cls.pkl.gz'
		#self.noise_cls_file_90ghz_100sqdeg = 'data/noise_cls/90ghz/combined_noise_cls.pkl.gz'
		#self.noise_cls_file_150ghz_100sqdeg = 'data/noise_cls/150ghz/combined_noise_cls.pkl.gz'
		#self.noise_cls_file_tszfree_100sqdeg = 'data/noise_cls/tszfree_own/combined_noise_cls.pkl.gz'
		#self.noise_cls_file_tszfree_100sqdeg = 'data/noise_cls/tszfree/combined_noise_cls.pkl.gz'
		#self.noise_cls_file_tszfree_100sqdeg = 'data/noise_cls/tszfree_better/combined_noise_cls.pkl.gz'


	def fn_get_EB_from_QU(self, Q, U, mapparams):

		sims =self.sims

		lx, ly = sims.get_lxly(mapparams)
		angle = 2.*np.arctan2(lx, -ly)

		QFFT, UFFT = np.fft.fft2(Q),np.fft.fft2(U)
		E = np.fft.ifft2( np.cos(angle) * QFFT + np.sin(angle) * UFFT ).real
		B = np.fft.ifft2( -np.sin(angle) * QFFT + np.cos(angle) * UFFT ).real

		return E, B

	def fn_read_noise_map(self, noisemap_file = None, mask_file = None):
	
		if noisemap_file == None:
			noisemap_file = self.noisemap_file

		if mask_file == None:
			mask_file = self.mask_file

		NOISE_MAPS = pickle.load(gzip.open(noisemap_file,'rb'))

		MASK = pickle.load(open(mask_file,'r'))

		NOISE_MAPS *= MASK

		if self.do_pol:

			NOISE_T,NOISE_Q,NOISE_U = NOISE_MAPS

			#get E and B noise maps now
			ny, nx = NOISE_T.shape
			dx = dy = self.map_resol
			mapparams = [nx, ny, dx, dy]

			NOISE_E, NOISE_B = self.fn_get_EB_from_QU(NOISE_Q,NOISE_U,mapparams)
		
			return np.asarray( [NOISE_T,NOISE_Q,NOISE_U,NOISE_E, NOISE_B] )
		else:

			return np.asarray( [NOISE_MAPS[0]] )

	def fn_get_noise_PSD(self, noisemap_file = None, mask_file = None, simmapparams = None, get_box_noise_PSD = 0):

		SPTpol_NOISE_MAPS = self.fn_read_noise_map(noisemap_file, mask_file)
		SPTpol_NOISE_FFT = np.fft.fft2(SPTpol_NOISE_MAPS)

		SPTpol_NOISE_PSD = abs(SPTpol_NOISE_FFT)**2.

		dx = dy = np.radians( self.map_resol / 60. )
		#scalefac = 1./(dx * dy)
		print '*** Alert: Check the 2pi factor in PSD normalisation'
		scalefac = 1./( 2 * np.pi * dx * dy) #CHECK

		SPTpol_NOISE_PSD/=scalefac

		if get_box_noise_PSD:

			sims = self.sims
			assert sims.is_seq(simmapparams)

			tqulen = len(SPTpol_NOISE_MAPS)
			ny_full, nx_full = SPTpol_NOISE_MAPS[0].shape
			mapparams_full = [nx_full, ny_full, dx, dy]

			lx_full, ly_full = sims.get_lxly(mapparams_full)
			lx, ly = sims.get_lxly(simmapparams)

			lx_full = np.fft.fftshift( lx_full )
			ly_full = np.fft.fftshift( ly_full )
			lx = np.fft.fftshift( lx )
			ly = np.fft.fftshift( ly )

			nx, ny, dx, dy = simmapparams

			NOISE_PSD = np.zeros( (tqulen, nx, ny ) )
			
			for tqucnt in range(tqulen):
				TMP = np.fft.fftshift( intrp.RectBivariateSpline( ly_full[:,0], lx_full[0,:], SPTpol_NOISE_PSD[tqucnt], kx = 5, ky = 5).ev(ly_full.flatten(),lx_full.flatten()).reshape([ny_full,nx_full]) )
				subplot(121);imshow(TMP, extent = [np.min(lx_full), np.max(lx_full), np.min(ly_full), np.max(ly_full) ]);colorbar()

				TMP = np.fft.fftshift( intrp.RectBivariateSpline( ly_full[:,0], lx_full[0,:], SPTpol_NOISE_PSD[tqucnt], kx = 5, ky = 5).ev(ly.flatten(),lx.flatten()).reshape([ny,nx]) )
				subplot(122);imshow(TMP, extent = [np.min(lx), np.max(lx), np.min(ly), np.max(ly) ]);colorbar();show();quit()

			quit()
			return NOISE_PSD
		else:

			return SPTpol_NOISE_PSD


	def fn_get_noise_sims_from_PSD(self, mapparams, noofsims, random_seed = None, PSD2D = None, do_pol = 1):

		if not do_pol:
			self.do_pol = 0

		if PSD2D == None:
			PSD2D = self.fn_get_noise_PSD(simmapparams = mapparams, get_box_noise_PSD = 1)

		print PSD2D

		tqulen = len(PSD2D)
		nx, ny, dx, dy = mapparams
		NOISE_SIMS = np.zeros( (tqulen, noofsims, ny, nx) )

		if random_seed == None:
			random_seed = int(time.time())
		np.random.seed(random_seed)

		for simcnt in range(noofsims):
			RANDOM_REALS = np.fft.fft2( np.random.randn(ny,nx) )
			for tqucnt in range(tqulen):

				NOISE_SIMS[tqucnt,simcnt] = np.fft.ifft2( RANDOM_REALS * np.sqrt(PSD2D[tqucnt]) ).real

				imshow(NOISE_SIMS[tqucnt,simcnt]);colorbar();show();quit()

	def fn_get_noise_sims_from_noise_cls(self, mapparams, noofsims, reqdbox = None, noise_cls_file = None, random_seed = None, min_el = None, max_el = None, nu = 150, whichmap = None, quiet = 0):

		"""
		look into reading_maps.py about calculation of Cls on sptcloud machine
		note that these are binned Cls - not theoretical
		"""
		sims = self.sims

		################################################
		#should we use SPTpol PSD or Nls for noise
		if whichmap.find('_500sqdeg_LBleem')>-1:
			use_psd = 1
		else:
			use_psd = 0

		################################################

		################################################
		#map stuff
		nx, ny, dx, dy = mapparams
		dx *= sims.arcmins2radians
		dy *= sims.arcmins2radians
		################################################

		################################################
		#self.scalefac = np.sqrt((nx * ny) / (dx * dy))
		self.scalefac = np.sqrt(1./ (dx * dy))
		################################################

		if not use_psd:
			if noise_cls_file == None:
				if whichmap == None:
					noise_cls_file = self.noise_cls_file
				else:
					cmd = 'noise_cls_file = self.noise_cls_file_%s'  %(whichmap)
					exec( cmd )
		
			if not quiet:
				logline = '\t\tNoise cls file = %s' %(noise_cls_file)
				#logfile = open(sims.log_file,'a')
				#logfile.writelines('%s\n' %(logline))
				#logfile.close()
				print logline

			noise_cls = pickle.load(gzip.open(noise_cls_file,'rb'))

			################################################
			#convert cls2CLS matrix
			if len(noise_cls)<10: #then els, Cls must be transposed
				noise_cls = noise_cls.T
			NOISE_CLS = sims.Cls2CLS(noise_cls,mapparams)
			NOISE_CLS = np.sqrt(NOISE_CLS)

			'''
		        #noise cannot be spherically symmetric
		        lx, ly = sims.get_lxly(mapparams)
		        L = np.sqrt(lx**2. + ly**2.)
		        TWODTF_FOR_NOISE =1.+ np.zeros(L.shape)
		        scan_hipass = np.where( (np.abs(lx) < 400.) )
		        TWODTF_FOR_NOISE[scan_hipass] = 0.
		        NOISE_CLS = NOISE_CLS * TWODTF_FOR_NOISE
			'''

			NOISE_FFT = NOISE_CLS * self.scalefac

			#NOISE_FFT *= 2#for 500 sq. deg. field
		        #imshow(np.fft.fftshift( NOISE_CLS[0] ));colorbar();show();quit()

		else:
			################################################
			#how does the SPT PSD look?
			deg_fac_sptpol_clus_map = int(np.degrees(dx) * 60. / self.sptpol_clus_map_dx )

			if not sims.is_seq(self.SPT_PSD_90):
				SPT_PSD_file = self.SPT_PSD_file_90ghz
				logline = '\t\tAlert!!!! replacing noise power spectra with SPT PSD - %s' %(SPT_PSD_file)
				print logline
				if SPT_PSD_file.find('.sav')>-1:
					SPT_PSD = readsav(SPT_PSD_file)['psd']
				else:
					SPT_PSD = pickle.load(gzip.open(SPT_PSD_file,'rb'))

				#subplot(121);imshow(np.fft.fftshift( SPT_PSD )  *1e6); colorbar();
				SPT_PSD = sims.downsample_map(SPT_PSD,deg_fac_sptpol_clus_map) #CHECK - is this correct?
				#subplot(122);imshow(np.fft.fftshift( SPT_PSD )  *1e6); colorbar();show();quit()

				SPT_PSD = np.fft.fftshift( sims.fn_get_SPTpol_TWODTF(mapparams, TWODTF_FULL = SPT_PSD)[0] * 1e6 )
				#SPT_PSD = sims.fn_get_SPTpol_TWODTF(mapparams, TWODTF_FULL = SPT_PSD)[0] * 1e6 #wrong
				SPT_PSD = np.asarray( [SPT_PSD] )

				self.SPT_PSD_90 = SPT_PSD

			if not sims.is_seq(self.SPT_PSD_150):
				SPT_PSD_file = self.SPT_PSD_file_150ghz
				logline = '\t\tAlert!!!! replacing noise power spectra with SPT PSD - %s' %(SPT_PSD_file)
				print logline
				if SPT_PSD_file.find('.sav')>-1:
					SPT_PSD = readsav(SPT_PSD_file)['psd']
				else:
					SPT_PSD = pickle.load(gzip.open(SPT_PSD_file,'rb'))

				#subplot(121);imshow(np.fft.fftshift( SPT_PSD )  *1e6); colorbar();
				SPT_PSD = sims.downsample_map(SPT_PSD,deg_fac_sptpol_clus_map) #CHECK - is this correct?
				#subplot(122);imshow(np.fft.fftshift( SPT_PSD )  *1e6); colorbar();show();quit()

				SPT_PSD = np.fft.fftshift( sims.fn_get_SPTpol_TWODTF(mapparams, TWODTF_FULL = SPT_PSD)[0] * 1e6 )
				#SPT_PSD = sims.fn_get_SPTpol_TWODTF(mapparams, TWODTF_FULL = SPT_PSD)[0] * 1e6 #wrong
				SPT_PSD = np.asarray( [SPT_PSD] )
				self.SPT_PSD_150 = SPT_PSD

			if nu == 90:
				SPT_PSD = self.SPT_PSD_90
			elif nu == 150:
				SPT_PSD = self.SPT_PSD_150

			#subplot(111);imshow(np.fft.fftshift( SPT_PSD[0] ));colorbar();show();quit()
			#subplot(122);imshow(np.fft.fftshift( NOISE_CLS[0] ));colorbar();show();quit()

			NOISE_FFT = SPT_PSD * self.scalefac

		################################################
		#now make sims
		if random_seed == None:
			random_seed = int(time.time())
		else: ### adding this 20171213
			np.random.seed(random_seed)

		noise_tqulen = len(NOISE_FFT) 
		if noise_tqulen>1:
			noise_tqulen += 2#+2 is for E, and B

		nx, ny, dx, dy = mapparams
		if reqdbox == None:
			NOISE_SIMS = np.zeros( (noofsims, noise_tqulen, ny, nx) ) #T, U, U, E, B
		else:
			ex1, ex2, ey1, ey2 = reqdbox
			nx_reqd, ny_reqd = ex2-ex1, ey2-ey1
			NOISE_SIMS = np.zeros( (noofsims, noise_tqulen, ny_reqd, nx_reqd) ) #T, U, U, E, B

		if min_el <> None:
			TWODTF = sims.fn_get_HPF(mapparams, minel = min_el, maxel = max_el) #filtering LSS from noise modelling

		for simcnt in range(noofsims):
			if simcnt % 5000 == 0 and simcnt > 0:
				logline = '\t\t\t\tsimno: %s' %(simcnt)
				print logline

			NOISE_DUMMY = np.fft.fft2( np.random.randn(ny,nx) )

			#print mapparams, noofsims
			if min_el <> None:
				TMP = np.fft.ifft2( NOISE_FFT * NOISE_DUMMY * TWODTF).real
			else:
				TMP = np.fft.ifft2( NOISE_FFT * NOISE_DUMMY ).real

			#remove mean just in case
			for tqucnt in range(noise_tqulen-2):
				TMP[tqucnt] = TMP[tqucnt] - np.mean(TMP[tqucnt])

			'''
			if nu == 90.:
				print '\t\t\tIncreasing noise by factor of 2. map depth is twice higher for 90 GHz channel'
				TMP *= 2. #map depth is twice higher for 90 GHz channel
			'''

			NOISE_T = TMP[0]

			#imshow(NOISE_T, interpolation = 'None');colorbar();show();quit()
			#print np.var(NOISE_Q)**.5 * dx
			#subplot(111);imshow(NOISE_T);colorbar();show();quit()
			#subplot(122);hist(NOISE_T.ravel(),bins=100);title(np.std(NOISE_T.ravel()) * dx);show();quit()
			

			if noise_tqulen > 1:
				NOISE_Q = NOISE_U = np.copy(NOISE_T) * 0.

				#get E and B noise maps now
				NOISE_E, NOISE_B = self.fn_get_EB_from_QU(NOISE_Q,NOISE_U,mapparams)

				#push all into array now
				if reqdbox == None:
					NOISE_SIMS[simcnt,:] = NOISE_T,NOISE_Q,NOISE_U,NOISE_E, NOISE_B
				else:
					NOISE_SIMS[simcnt,:] = NOISE_T[ex1:ex2, ey1:ey2],NOISE_Q[ex1:ex2, ey1:ey2],NOISE_U[ex1:ex2, ey1:ey2],NOISE_E[ex1:ex2, ey1:ey2], NOISE_B[ex1:ex2, ey1:ey2]

			else:

				if reqdbox == None:
					NOISE_SIMS[simcnt,:] = NOISE_T
				else:
					NOISE_SIMS[simcnt,:] = NOISE_T[ex1:ex2, ey1:ey2]

			'''
			for tqucnt in range(noise_tqulen):
				subplot(2,3,tqucnt+1);css=imshow(NOISE_SIMS[tqucnt,simcnt]);colorbar()

			whitenoise = np.fft.ifft2(NOISE_DUMMY).real * 7./dx
			subplot(2,3,6);css=imshow(whitenoise);colorbar()
			show();quit()

			'''
		return NOISE_SIMS		

class cluster_stuff():

	"""
	class to have models/calculations related to tSZ, kSZ, tSZ power spectra, etc.
	"""

	def __init__(self):
		self.sims = simulations()
		self.h=6.62607004e-34 #Planck constant in m2 kg / s
		self.k_B=1.38064852e-23 #Boltzmann constant in m2 kg s-2 / K-1
		self.Tcmb=2.73 #Kelvin
		self.delta_nu = 1e9 #Hz
		self.tSZ_cutouts = None
		self.tSZ_sim_file = 'data/tsz_y_cutouts_M200m_v01.pk' #contains ['M500c', 'y_cutouts', 'npix', 'reso_arcmin', 'z', 'M200m']
		self.mass_tolerance =0.25e14#1e14
		self.K_uK = 1e6

	def coth(self,x):
		return (np.exp(x) + np.exp(-x)) / (np.exp(x) - np.exp(-x))

	def compton_y_to_delta_Tcmb(self, nu_1, nu_2 = None):

		"""
		c.f:  table 1, sec. 3 of arXiv: 1303.5081; 
		table 8 of http://arxiv.org/pdf/1303.5070.pdf
		no relativistic corrections included.
		nu_1, nu_2 = frequencies in GHz to cover the bandpass
		nu_2 = None will force nu_1 to be the centre frequency
		"""
		if not nu_2 == None:
			nu = np.arange(nu_1,nu_2,self.delta_nu)
		else:
			nu = [nu_1]

		
		x = np.asarray( map( lambda n: (self.h * n) / (self.k_B * self.Tcmb), nu) )
		g_nu = np.asarray( map( lambda n: n * self.coth(n/2.) - 4., x) )

		return self.Tcmb * np.mean(g_nu)

		'''
		from scipy import integrate
		def fn_conv_y_deltaTcmb(nu):
			x = self.h * nu / (self.k_B * self.Tcmb)
                        return (self.Tcmb * (x * self.coth(x/2.) - 4.) )

		if nu_2 <> None:
	                conversion_fac, conversion_fac_err = integrate.quad(fn_conv_y_deltaTcmb,nu_1,nu_2)
		else:
			conversion_fac = fn_conv_y_deltaTcmb(nu_1)
			
		return conversion_fac
		'''

	def fn_read_cutouts(self):
		tSZ = pickle.load(open(self.tSZ_sim_file,'r'))

		return tSZ
		
	def fn_get_closest_mass_halos(self, halo_mass):

		diffmasses = self.tSZ_M200 - halo_mass
		inds = np.where(abs(diffmasses)<self.mass_tolerance)[0]
		itercnt = 1
		while len(inds) == 0:
			inds = np.where(abs(diffmasses)< (self.mass_tolerance * (itercnt+1)))[0]
			itercnt += 1
		return inds

	def fn_pick_add_sehgal_sims(self, halo_mass, which_comps = 'all', also_90 = 1, zmin = 0.25, percent_mass_tolerance = 0.05, ipfolder = 'data/sehgal_2009_sims/tsz_cib_radio_ksz_extracts/', mass_def = 'mean'):

		if also_90:
			fname = '%s/tSZ_CIB_radio_kSZ_summed_extracts_090_M200min1.3_zmin0.1_boxsize100.0am_dx0.5am.pkl.gz' %(ipfolder)
			#fname = '%s/tSZ_CIB_radio_kSZ_summed_extracts_090_M200min1.3_zmin0.35_boxsize100.5am_dx0.5am.pkl.gz' %(ipfolder)
			if not os.path.exists(fname): fname = '%s_test' %(fname)
			sehgal_dic_90 = pickle.load(gzip.open(fname, 'rb'))

		fname = '%s/tSZ_CIB_radio_kSZ_summed_extracts_148_M200min1.3_zmin0.1_boxsize100.0am_dx0.5am.pkl.gz' %(ipfolder)
		#fname = '%s/tSZ_CIB_radio_kSZ_summed_extracts_148_M200min1.3_zmin0.35_boxsize100.5am_dx0.5am.pkl.gz' %(ipfolder)
		if not os.path.exists(fname): fname = '%s_test' %(fname)
		sehgal_dic = pickle.load(gzip.open(fname, 'rb'))

		if also_90: assert sehgal_dic.keys() == sehgal_dic_90.keys()

		keynames = sehgal_dic.keys()
		M200m_sehgal_sims = np.asarray(keynames)[:,3]
		redshift_sehgal_sims = np.asarray(keynames)[:,2]

		### from IPython import embed; embed()

		if also_90:
			sel_cutouts_90 = []
		sel_cutouts = []

		fixed_M = 0
		#if len(np.unique(halo_mass)) == 1: fixed_M = 1

		for mmcnt, mm in enumerate(halo_mass):
			mm = mm/1e14
			mtoldel = mm * percent_mass_tolerance
			closestinds = np.where( (M200m_sehgal_sims>=mm-mtoldel) & (M200m_sehgal_sims<=mm+mtoldel) & (redshift_sehgal_sims>=zmin))[0]

			### print mm, len(closestinds)

			if len(closestinds) == 0:
				print 'No sehgal sims in the given mass range: %s, %s' %(mm - mtoldel, mm + mtoldel)
				from IPython import embed; embed()
				sys.exit()

			if not fixed_M:
				selkey = keynames[closestinds[np.random.randint(len(closestinds))]]
			else:
				selkey = keynames[closestinds[0]]

			iterrange = 1
			if also_90: 
				iterrange = 2

			for nn in range(iterrange):
				if nn == 0:
					tmp_sel_cutout = np.asarray( sehgal_dic[selkey] )* self.K_uK
				else:
					tmp_sel_cutout = np.asarray( sehgal_dic_90[selkey] )* self.K_uK

				if nn == 0:
					sel_cutouts.append(tmp_sel_cutout)
				else:
					sel_cutouts_90.append(tmp_sel_cutout)

			if fixed_M:
				sel_cutouts = [sel_cutouts[0] for i in range(len(halo_mass))]
				if also_90:
					sel_cutouts_90 = [sel_cutouts_90[0] for i in range(len(halo_mass))]
				break

		if also_90:
			return np.asarray( sel_cutouts ), np.asarray( sel_cutouts_90 )
		else:
			return np.asarray( sel_cutouts )

	def fn_pick_add_sehgal_sims_v1(self, halo_mass, which_comps = 'all', also_90 = 1, zmin = 0.1, percent_mass_tolerance = 0.05, ipfolder = 'data/sehgal_2009_sims/tsz_cib_radio_ksz_extracts/', mass_def = 'mean'):


		if also_90:
			#fname = '%s/tSZ_CIB_radio_extracts_090_M200min1.3_zmin0.0_boxsize50.0am_dx0.5am.pkl.gz' %(ipfolder)
			fname = '%s/tSZ_CIB_radio_extracts_090_M200min1.3_zmin0.1_boxsize50.0am_dx0.5am.pkl.gz' %(ipfolder)
			sehgal_dic_90 = pickle.load(gzip.open(fname, 'rb'))

		#fname = '%s/tSZ_CIB_radio_extracts_148_M200min1.3_zmin0.0_boxsize50.0am_dx0.5am.pkl.gz' %(ipfolder)
		fname = '%s/tSZ_CIB_radio_extracts_148_M200min1.3_zmin0.1_boxsize50.0am_dx0.5am.pkl.gz' %(ipfolder)
		sehgal_dic = pickle.load(gzip.open(fname, 'rb'))

		if also_90: assert sehgal_dic.keys() == sehgal_dic_90.keys()

		keynames = sehgal_dic.keys()
		M200m_sehgal_sims = np.asarray(keynames)[:,3]
		redshift_sehgal_sims = np.asarray(keynames)[:,2]

		### from IPython import embed; embed()

		if also_90:
			sel_cutouts_90 = []
		sel_cutouts = []

		for mmcnt, mm in enumerate(halo_mass):
			mm = mm/1e14
			mtoldel = mm * percent_mass_tolerance
			closestinds = np.where( (M200m_sehgal_sims>=mm-mtoldel) & (M200m_sehgal_sims<=mm+mtoldel) & (redshift_sehgal_sims>=zmin))[0]

			### print mm, len(closestinds)

			if len(closestinds) == 0: 
				from IPython import embed; embed()
				print 'No sehgal sims in the given mass range: %s, %s' %(mm - mtoldel, mm + mtoldel)
				sys.exit()
			'''
			while len(closestinds) == 0:
				loopcnt = 2 #increase tolerance
				mtoldel = mm * percent_mass_tolerance * loopcnt
				closestinds = np.where( (M200m_sehgal_sims>=mm-mtoldel) & (M200m_sehgal_sims>=mm+mtoldel))[0]
				loopcnt += 1
			'''
			### from IPython import embed; embed()
			#randomly pick one of the closest inds
			selkey = keynames[closestinds[np.random.randint(len(closestinds))]]

			### print mm, selkey

			iterrange = 1
			if also_90: 
				iterrange = 2

			### from IPython import embed; embed()

			for nn in range(2):
				#print selkey
				if nn == 0:
					tmp_sel_cutout = np.asarray( sehgal_dic[selkey] )* self.K_uK
				else:
					tmp_sel_cutout = np.asarray( sehgal_dic_90[selkey] )* self.K_uK

				if which_comps == 'all': s,e = 0,4
				elif which_comps == 'tSZ': s,e = 0,1
				elif which_comps == 'ir': s,e = 1,2
				elif which_comps == 'radio': s,e = 2,3
				elif which_comps == 'ksz': s,e = 3,4

				tmp_sel_cutout = np.sum(tmp_sel_cutout[s:e], axis = 0)
				tmp_sel_cutout -= np.mean(tmp_sel_cutout)

				if nn == 0:
					sel_cutouts.append(tmp_sel_cutout)
				else:
					sel_cutouts_90.append(tmp_sel_cutout)

		if also_90:
			return np.asarray( sel_cutouts ), np.asarray( sel_cutouts_90 )
		else:
			return np.asarray( sel_cutouts )

	def fn_pick_cutout(self, nu, halo_mass, mapparams, return_as_dic = 0):

		#if self.sims.is_seq(self.tSZ_cutouts) == None:
		if self.tSZ_cutouts == None:
			tSZ = self.fn_read_cutouts()
			self.tSZ_M200 = tSZ['M200m']
			self.tSZ_cutouts = tSZ['y_cutouts']
			self.tSZ_reso_arcmin = tSZ['reso_arcmin']
			self.tSZ_npix = tSZ['npix']

		y_to_T_conv = self.compton_y_to_delta_Tcmb(nu)
		nx, ny, dx, dy = mapparams

		deg_fac = dx/self.tSZ_reso_arcmin
		self.tSZ_npix /= deg_fac

		tSZ_nx = self.tSZ_cutouts.shape[0]
		if nx > tSZ_nx:
			nx = ny = tSZ_nx

		ex1, ex2 = int(self.tSZ_npix/2) - int(nx/2), int(self.tSZ_npix/2) + int(nx/2)
		ey1, ey2 = int(self.tSZ_npix/2) - int(ny/2), int(self.tSZ_npix/2) + int(ny/2)

		sel_cutouts = []
		sel_cutouts_dic = {}

		np.random.seed(self.cmbrandomseedval)
		for mmcnt, mm in enumerate(halo_mass):
			
			closestinds = self.fn_get_closest_mass_halos(mm)

			#randomly pick one of the closest inds
			selind = closestinds[np.random.randint(len(closestinds))]

			sel_cutout = self.tSZ_cutouts[:,:,selind]
			sel_cutout = self.sims.downsample_map(sel_cutout,deg_fac)
			sel_cutout = sel_cutout[ex1:ex2, ey1:ey2]

			sel_cutout =  sel_cutout * y_to_T_conv * self.K_uK

			#imshow(sel_cutout);colorbar();show();quit()
			sel_cutouts.append(sel_cutout)

			sel_cutouts_dic[(mmcnt,mm)] = sel_cutout

			#imshow(sel_cutout);colorbar();grid(1);show();quit()

		#print sel_cutouts_dic.keys()
		#quit()

		if not return_as_dic:
			return np.asarray( sel_cutouts )
		else:
			return sel_cutouts_dic


	#def fn_beta_model_for_tSZ(self, MAP):

		
		

class SPTpol_map_using_SPTpol_pipeline():

	def __init__(self):

		#T->P monopole deprojection values
		self.project_out_t = {'qScale':0.0139216, 'uScale':0.0078030}
		self.clustercatfile = 'data/500sqdeg_bleem/sptpol_4sigma_matched_sptsz.txt' #provided by Linsey Bleem
		self.map_centre = [0.0, -59.033333]
		self.map_resol = 0.5 #arcmins
		self.mapshape = np.asarray( [2640, 5040] )
		self.map_proj = 5 #
		self.mapfile = '/data57/jhenning/ra0hdec-57p5_V3/20160518/bundles/150/no_cuts/coadd_all.h5'
		self.degrade_fac = None
		self.cutout_size = 60 #arcmins #should be the same as COV matrix
		self.clus_snrlevel = 4.
		self.ang_dist_tol_from_source = 0.3#(self.cutout_size / 60. ) * 3. #degrees
		self.ang_dist_tol_from_other_cutouts = 0.1 #if point sources are already masked

	def _ini_single_params(self,param_name, param_val):
		cmd = '%s = %s' %(param_name, param_val)
		exec(cmd)

        def fn_ang_dist(self,ip_x1, ip_y1, ip_x2, ip_y2): #all numpy arrays

                # CONVERT INTO RADIANS
                x1 = np.radians(ip_x1)
                y1 = np.radians(ip_y1)
                x2 = np.radians(ip_x2)
                y2 = np.radians(ip_y2)

                ang_dist_rad = np.arccos(np.sin(y1)*np.sin(y2)+(np.cos(y1)*np.cos(y2)*np.cos(x1 - x2)))
                ang_dist_deg = np.degrees(ang_dist_rad) # Converting radians to degrees

                return ang_dist_deg

		
	def fn_read_process_SPTpol_maps(self, stored_map_file, mapfile = None, maskfile = None, project_out_t = None, degrade_fac = None):
		"""
		adopted from email from Zhen Hou on Nov. 4, 2015
		"""

		import sptpol_software.util.files as files
		import sptpol_software.analysis.maps as maps

		smap = files.read(mapfile)

		#Remove weight and flatten pol angles
		smap.removeWeight().flattenPol()

		#Converting units to uK_CMB.
		smap*=1e6

		#Actually deproject T from P.
		if project_out_t == None:
			project_out_t = self.project_out_t
		smap = maps.projectTFromP(smap, qScale=project_out_t['qScale'], uScale=project_out_t['uScale'])

		if not maskfile == None:
			mask = files.read(maskfile)
			smap *= mask

		if degrade_fac <> None:
			self.degrade_fac = degrade_fac

		if self.degrade_fac <> None:
			smap = smap.degradedMap(self.degrade_fac)
			stored_map_file = stored_map_file.replace('.pkl.gz','_degraded_%s.pkl.gz' %(self.degrade_fac))

		MAPS = np.asarray( [smap['T'].map,smap['Q'].map,smap['U'].map] )

		pickle.dump(MAPS, gzip.open(stored_map_file,'wb'), protocol = 2 )

		#return MAPS

	def fn_read_stored_map(self, mapfile):

		if self.degrade_fac <> None:
			mapfile = mapfile.replace('.pkl.gz','_degraded_%s.pkl.gz' %(self.degrade_fac))
		MAPS = pickle.load( gzip.open(mapfile,'r') )

		return MAPS, mapfile

	def fn_extract_map_cutouts_for_noise_foreground(self, mapfile, no_cutouts, MAPS = None, clustercatfile = None, maskfile = None, point_source_file = None, cutout_dic_name = None, ret_cutout_dic = 0, field_size = 500., degrade_cutout = None, weight_map = None, minweight = None, maxweight = None):

		srcra, srcdec = [], []
		if clustercatfile <> None: #do not pick regions near clusters
			cra,cdec,cz,cconf = self.fn_read_cluster_catalogue(clustercatfile, get_pixel_info = 0)
			srcra, srcdec = np.copy(cra), np.copy(cdec)

		if point_source_file <> None: #do not pick regions near point sources
			pdata = np.loadtxt(point_source_file,usecols = [1,2])
			pra,pdec = np.asarray( zip(*pdata) )
			srcra = np.concatenate( (srcra,pra) )
			srcdec = np.concatenate( (srcdec,pdec) )


		#read map
		if MAPS == None:
			MAPS,mapfile = self.fn_read_stored_map(mapfile)
		MAPSHAPE = MAPS.shape[1:]


		if not maskfile == None:
			total_mask = pickle.load(open(maskfile,'r'))
			#assert MAPS[0].shape == total_mask.shape
			try:
				MAPS *= total_mask
			except ValueError as errvals:
				print '\n Mask not applied. Becaause %s' %(errvals)

		map_resol = self.map_resol
		if self.degrade_fac<>None:
			 map_resol *= self.degrade_fac

		cutout_dic = {}
		cutout_dic['map_source'] = mapfile
		cutout_dic['cluster_catalogue'] = clustercatfile
		cutout_dic['cutout_size_arcmins'] = self.cutout_size
		cutout_dic['resol_arcmins'] = map_resol
		cutout_dic['degraded_from_ori'] = self.degrade_fac
		cutout_dic['ang_dist_tol_from_source_arcmins'] = self.ang_dist_tol_from_source * 60.
		cutout_dic['ang_dist_tol_from_other_cutouts_arcmins'] = self.ang_dist_tol_from_other_cutouts * 60.

		cutout_dic['cutouts'] = {}

		if cutout_dic_name == None:
			cutout_dic_name = '%s_noise_foreground_%s_cutouts' %(mapfile,no_cutouts)

		'''
		if field_size == 500.:
			delta_ra_field, delta_dec_field = 18., 18. #500 sq-deg field = 500 **0.5 = 22.36 degrees; being conservative and adopting 18 deg.
		else:
			delta_ra_field, delta_dec_field = 5., 5. #500 sq-deg field = 500 **0.5 = 22.36 degrees; being conservative and adopting 18 deg.
		'''

		delta_ra_field = delta_dec_field = field_size**0.5 - 5.

		minra, maxra = self.map_centre[0] - delta_ra_field, self.map_centre[0] + delta_ra_field
		mindec, maxdec = self.map_centre[1] - delta_dec_field, self.map_centre[1] + delta_dec_field

		import sky_local

		map_resol = self.map_resol
		mapshape = self.mapshape
		if self.degrade_fac <> None:
			map_resol *= self.degrade_fac
			mapshape /= self.degrade_fac
		
		start = time.time()
		picked_ra,picked_dec = [], []

		#MAPFORCOV = []

		log_file = 'tmp/noise_foreground_cutout_extraction_%s.txt' %(no_cutouts)
		n = 0
		while len( cutout_dic['cutouts'].keys() ) < no_cutouts:

			randra, randdec = np.random.uniform(minra, maxra), np.random.uniform(mindec, maxdec)

			ang_dist_deg = self.fn_ang_dist(srcra, srcdec, np.tile(randra,len(srcra)), np.tile(randdec,len(srcra)) )
			if min(ang_dist_deg) < self.ang_dist_tol_from_source:
				#print 'Ingored as the pixel is too close to some cluster/ pointsource'
				continue

			
			if len(picked_ra)>0:
				ang_dist_deg = self.fn_ang_dist(np.asarray(picked_ra), np.asarray(picked_dec), np.tile(randra,len(picked_ra)), np.tile(randdec,len(picked_ra)) )
				if min(ang_dist_deg) < self.ang_dist_tol_from_other_cutouts:
					continue

			pixels = sky_local.ang2Pix([randra, randdec], self.map_centre, map_resol, mapshape, proj = self.map_proj, round=True, bin_center_zero=True, return_validity=False, use_c_code=False)

			y, x = pixels

			if weight_map <> None:
				#print weight_map
				W_CUTOUT = self.fn_pick_cutout(x,y,np.asarray( [weight_map] ),self.cutout_size,map_resol,MAPSHAPE, perform_mean_sub = 0, perform_checks = 0)[0]				
				avg_weight = np.mean(W_CUTOUT)
				if avg_weight<minweight or avg_weight>maxweight:
					continue

			"""
			ignedges = 700
			minval, maxval = ignedges, MAPSHAPE[0] - ignedges
			x, y = np.random.randint(minval, maxval), np.random.randint(minval, maxval)
			randra, randdec = x,y
			#CUTOUTS = MAPS[:,x-nx/2:x+nx/2, y-ny/2:y+ny/2]
			"""

			CUTOUTS = self.fn_pick_cutout(x,y,MAPS,self.cutout_size,map_resol,MAPSHAPE)
			if not hasattr(CUTOUTS, '__len__'):
				#print randra, randdec, 'not picked'
				continue

			if degrade_cutout <> None:
				sims = simulations()
				TMP = []
				for cutout in CUTOUTS: 
					cutout = sims.downsample_map(cutout,degrade_cutout)
					TMP.append(cutout)

				CUTOUTS = np.asarray(TMP)

			if len(CUTOUTS) == 1: #only T then just make Q, U to be zero
				Q = U = np.zeros(CUTOUTS[0].shape)
				CUTOUTS = np.asarray( [CUTOUTS[0], Q, U] )

			keyname = (round(randra,3),round(randdec,3))
			cutout_dic['cutouts'][keyname] = CUTOUTS

			#MAPFORCOV.append(CUTOUTS[0])

			#print randra, randdec, 'picked'
			logline = 'Picked %s of %s; Time spent so far = %s seconds' %(len( cutout_dic['cutouts'].keys() ), no_cutouts, time.time()-start)
			logfile = open(log_file,'a')
			logfile.writelines('%s\n' %(logline))
			logfile.close()
			print logline

			'''
			#add this guy to srcra, srcdec to not repeat the cutout
			srcra = srcra.tolist(); srcra.append(randra); srcra = np.asarray(srcra)
			srcdec = srcdec.tolist(); srcdec.append(randdec); srcdec = np.asarray(srcdec)
			'''

			#commenting this - as the process becomes extremely slow
			#picked_ra.append(randra)
			#picked_dec.append(randdec)

			'''
			#if np.random.randint(no_cutouts)%100 == 0: #randomly show some of these
			if n<=10:
				subplot(1,len(MAPS),1);imshow(CUTOUTS[0]);cbar=colorbar();cbar.ax.tick_params(labelsize=8);title('T - %s,%s, x = %s, y = %s' %(n,str(keyname),x,y), fontsize = 8)
				if len(MAPS)>1:
					subplot(1,len(MAPS),2);imshow(CUTOUTS[1]);cbar=colorbar();cbar.ax.tick_params(labelsize=8);title('Q - %s,%s' %(n,str(keyname)), fontsize = 8)
					subplot(1,len(MAPS),3);imshow(CUTOUTS[2]);cbar=colorbar();cbar.ax.tick_params(labelsize=8);title('U - %s,%s' %(n,str(keyname)), fontsize = 8);
				show()#;quit()
			'''

			n += 1


		'''
		### check the cov now
		MAPFORCOV = np.asarray(MAPFORCOV)
		nx, ny = MAPFORCOV.shape[1], MAPFORCOV.shape[2]
		npixels = nx * ny
		sims = simulations()
		COV = sims.calcCov(MAPFORCOV, len(MAPFORCOV), npixels)
		imshow(COV, vmin=-200, vmax=1200);colorbar();show();quit()
		'''

		#plot(picked_ra,'bo');plot(picked_dec,'go');show();quit()
		end = time.time()
		print 'All cutouts obtained. Total time taken: %s minutes' %( (end-start)/60.)

		pickle.dump(cutout_dic, gzip.open(cutout_dic_name, 'wb') , protocol = 2)

		if ret_cutout_dic:
			return cutout_dic			

	def fn_check_for_scan_edges(self, CUTOUT):

		#clusters are not really near edge of scanned region - still good to check
		#make sure cluster is not outside scanned area
		tmp = np.asarray( CUTOUT.ravel() )
		tol_fac = 0.05 #atleat (100 - tol_fac) = (100 - 5) = 95 percent of data should not be zero

		#if len(np.where(tmp==0.)[0]) > tol_fac * len(tmp):
		if len(np.where(abs(tmp)<1.)[0]) > tol_fac * len(tmp): #1. is such a small number as we are dealing with uK maps here
			return 1

		return 0	

	def fn_check_for_map_edge(self, x1, x2, y1, y2, MAPSHAPE):

		if x1<0. or x2>MAPSHAPE[1] or y1 <0. or y2>MAPSHAPE[0]:
			return 1

		return 0

	def fn_pick_cutout(self,x,y,MAPS,cutout_size,map_resol,MAPSHAPE,perform_mean_sub = 0, perform_checks = 1):

		nx = ny = cutout_size/map_resol
		#CUTOUT = np.zeros( (nx,ny) )
		x1,x2 = x-nx/2,x+nx/2
		y1,y2 = y-ny/2,y+ny/2

		if not perform_checks:
			tqulen = len(MAPS)
			CUTOUT = np.asarray( [MAPS[tqucnt][y1:y2,x1:x2] for tqucnt in range(tqulen)] )

			return CUTOUT


		#clusters are not really near map edge - still good to check
		map_edge_check = self.fn_check_for_map_edge(x1,x2,y1,y2,MAPSHAPE)

		#clusters are not really near map edge - still good to check
		map_edge_check = self.fn_check_for_map_edge(x1,x2,y1,y2,MAPSHAPE)
		if map_edge_check:
			#print 'Ignoring clusters at the map edges';n+=1
			return None

		#clusters are not really near edge of scanned region - still good to check
		#make sure cluster is not outside scanned area
		T_CUTOUT = MAPS[0][y1:y2,x1:x2]
		scan_edge_check = self.fn_check_for_scan_edges(T_CUTOUT)
		if scan_edge_check:
			#print 'Ignoring clusters near the edges of scanned region';n+=1
			return None
		tqulen = len(MAPS)
		CUTOUT = np.asarray( [MAPS[tqucnt][y1:y2,x1:x2] for tqucnt in range(tqulen)] )
		if perform_mean_sub:
			CUTOUT = np.asarray( [TMP - np.mean(TMP.ravel()) for TMP in CUTOUT] )

		return CUTOUT


	########################################################################################################################
	########################################################################################################################
	def fn_tf_deconv(self, CUTOUT, ra, dec):

		sims = simulations()
		try:
			TWODTF = self.TWODTF
		except:
			dx = dy = self.map_resol
			nx = ny = int(self.cutout_size / dx)
			mapparams = [nx, ny, dx, dy]
			TWODTF = sims.fn_get_SPTpol_TWODTF(mapparams)[0]
			self.TWODTF = TWODTF

		#clf();subplot(121);imshow(np.fft.fftshift(TWODTF));colorbar()
		TWODTF_rot, rot_angle = sims.rotate_tf(TWODTF, ra, dec, in_elspace = 1, return_angle = 1)
		#subplot(122);imshow(np.fft.fftshift(TWODTF_rot));colorbar();title(rot_angle); show();quit()
		DECONV = 1./TWODTF_rot
		DECONV[DECONV<0.] = 0.
		DECONV[DECONV == inf] = 0.

		#clf();imshow(DECONV);colorbar();show();quit()

		for tqu in range(len(CUTOUT)):

			DUMMY = np.fft.fft2(CUTOUT[tqu], s=(TWODTF_rot.shape))
			TF_DECONV = DUMMY * DECONV
			DUMMY_IFFT = np.fft.ifft2( TF_DECONV * TWODTF ).real
			CUTOUT[tqu] = DUMMY_IFFT[0:CUTOUT[tqu].shape[1], 0:CUTOUT[tqu].shape[0]]

		return CUTOUT

	########################################################################################################################

	def fn_extract_map_cutouts(self, mapfile, MAPS = None, pixels = None, clustercatfile = None, maskfile = None, cutout_dic_name = None, ret_cutout_dic = 0, degrade_cutout = None, weight_map = None):

		if pixels == None:
			assert clustercatfile <> None

		if pixels == None:
			ra,dec,z,conf,x_inds,y_inds = self.fn_read_cluster_catalogue(clustercatfile, get_pixel_info = 1)
		else:
			x_inds,y_inds = np.asarray( zip(*pixels) )


		if MAPS == None:
			MAPS,mapfile = self.fn_read_stored_map(mapfile)

		'''
		clf();imshow(MAPS[0],vmin=-150., vmax = 150.);
		plot(x_inds,y_inds,'ko');colorbar()
		xlim(0,MAPS[0].shape[1])
		ylim(0,MAPS[0].shape[0])
		show();quit()
		'''

		if not maskfile == None:
			total_mask = pickle.load(open(maskfile,'r'))
			#assert MAPS[0].shape == total_mask.shape
			try:
				MAPS *= total_mask
			except ValueError as errvals:
				print '\n Mask not applied. Because %s' %(errvals)

		#imshow(MAPS[0],vmin=-150, vmax=150.);colorbar();show();quit()

		map_resol = self.map_resol
		if self.degrade_fac<>None:
			 map_resol *= self.degrade_fac
		
		cutout_dic = {}
		cutout_dic['map_source'] = mapfile
		cutout_dic['cluster_catalogue'] = clustercatfile
		cutout_dic['cutout_size_arcmins'] = self.cutout_size
		cutout_dic['resol_arcmins'] = map_resol
		cutout_dic['degraded_from_ori'] = self.degrade_fac

		cutout_dic['cutouts'] = {}
		MAPSHAPE = MAPS.shape[1:]

		'''
		subplot(111);imshow(MAPS[0]);colorbar();title('T');show();quit()
		subplot(132);imshow(MAPS[1]);colorbar();title('Q')
		subplot(133);imshow(MAPS[2]);colorbar();title('U');show();quit()
		'''

		if cutout_dic_name == None:
			cutout_dic_name = '%s_cutouts' %(mapfile)

		n = 0
		for rr,dd,zz,cc,x,y in zip(ra,dec,z,conf,x_inds,y_inds): #ra,dec,z,conf,x_inds,y_inds loop

			'''
			if clustercatfile <> None:
				keyname = (round(rr,3),round(dd,3),round(zz,3),round(cc,3))
			else:
				keyname = n
			'''

			CUTOUTS = self.fn_pick_cutout(x,y,MAPS,self.cutout_size,map_resol,MAPSHAPE)
			if not hasattr(CUTOUTS, '__len__'):
				n+=1
				continue

			'''
			## deconvolve TF
			subplot(121);imshow(CUTOUTS[0]);colorbar()
			CUTOUTS = self.fn_tf_deconv(CUTOUTS, rr, dd)
			subplot(122);imshow(CUTOUTS[0]);colorbar();show();quit()
			'''

			avg_weight = 0.
			if weight_map <> None:
				#print weight_map
				W_CUTOUT = self.fn_pick_cutout(x,y,np.asarray( [weight_map] ),self.cutout_size,map_resol,MAPSHAPE, perform_mean_sub = 0, perform_checks = 0)[0]
				avg_weight = np.mean(W_CUTOUT)

			keyname = (round(rr,3),round(dd,3),round(zz,3),round(cc,3), round(avg_weight,3))

			debug = 0
			if debug:
				sims = simulations()
				T = CUTOUTS[0]
				PSD = np.fft.fftshift( abs(np.fft.fft2(T)) )
				nx, ny = T.shape
				mapparams = [nx, ny, 0.25, 0.25]
				lx, ly = sims.get_lxly(mapparams)
				clf();subplot(131);imshow(PSD,extent = [np.min(lx),np.max(lx),np.min(ly),np.max(ly)]);colorbar()

			if degrade_cutout <> None:
				sims = simulations()
				TMP = []
				for cutout in CUTOUTS: 
					cutout = sims.downsample_map(cutout,degrade_cutout)
					TMP.append(cutout)

				CUTOUTS = np.asarray(TMP)

			if debug:
				T = CUTOUTS[0]
				PSD = np.fft.fftshift( abs(np.fft.fft2(T)) )
				nx, ny = T.shape
				mapparams = [nx, ny, 0.5, 0.5]
				lx, ly = sims.get_lxly(mapparams)
				subplot(132);imshow(PSD,extent = [np.min(lx),np.max(lx),np.min(ly),np.max(ly)]);colorbar();

				subplot(133);imshow(T);colorbar();show();quit()


			if len(CUTOUTS) == 1: #only T then just make Q, U to be zero
				Q = U = np.zeros(CUTOUTS[0].shape)
				CUTOUTS = np.asarray( [CUTOUTS[0], Q, U] )


			print n, keyname
			if np.random.randint(len(x_inds))%20 == 0: #randomly show some of these
			#if n == 0 or n == 10 or n == 32:# some examples in 100 sq. deg. field fot tSZ_free maps
				#print keyname
				subplot(1,len(MAPS),1);imshow(CUTOUTS[0]);cbar=colorbar();cbar.ax.tick_params(labelsize=8);title('T - %s,%s, x = %s, y = %s' %(n,str(keyname),x,y), fontsize = 8);grid(1,ls='solid')
				if len(MAPS)>1:
					subplot(1,len(MAPS),2);imshow(CUTOUTS[1]);cbar=colorbar();cbar.ax.tick_params(labelsize=8);title('Q - %s,%s' %(n,str(keyname)), fontsize = 8)
					subplot(1,len(MAPS),3);imshow(CUTOUTS[2]);cbar=colorbar();cbar.ax.tick_params(labelsize=8);title('U - %s,%s' %(n,str(keyname)), fontsize = 8);
				show()#;quit()
			cutout_dic['cutouts'][keyname] = CUTOUTS
			n+=1

		print 'Total cutouts picked = %s' %(len(cutout_dic['cutouts'].keys()))
		pickle.dump(cutout_dic, gzip.open(cutout_dic_name, 'w') )

		if ret_cutout_dic:
			return cutout_dic	

	def fn_ang2pix(self, angles, ref_angles = None, mapshape = None, map_resol = None): #pass angles = [ra,dec] in degrees

		"""
		Based on: http://mathworld.wolfram.com/LambertAzimuthalEqual-AreaProjection.html
		Should be very similar to the one on SPTpol pipeline
		"""
		
		sin,cos = np.sin, np.cos

		if ref_angles == None:
			ref_angles = self.map_centre

		if mapshape == None:
			mapshape = self.mapshape

		if  map_resol == None:
			 map_resol = self. map_resol

		#convert to radians
		angles = np.radians(angles)
		ref_angles = np.radians(ref_angles)

		#get dec = phi, and ra = lamb
		lamb, phi = angles
		lamb_0, phi_1 = ref_angles


		k_prime_nr = 2.
		k_prime_dr_1 = 1.
		k_prime_dr_2 = sin(phi_1) * sin(phi) + cos(phi_1) * cos(phi) * cos(lamb - lamb_0)
		k_prime_dr = k_prime_dr_1 + k_prime_dr_2

		k_prime = np.sqrt( k_prime_nr/k_prime_dr )

		x = k_prime * cos(phi) * sin(lamb - lamb_0)
		y = k_prime * ( cos(phi_1) * sin(phi) - sin(phi_1) * cos(phi) * cos(lamb-lamb_0) )

		'''
		#do the 90 deg. adjustments for dec
		phi = np.radians(90.) - phi
		phi_1 = np.radians(90.) - phi_1

		k_prime_nr = 2.
		k_prime_dr_1 = 1.
		#k_prime_dr_2 = sin(phi_1) * sin(phi) + cos(phi_1) * cos(phi) * cos(lamb - lamb_0)
		k_prime_dr_2 = cos(phi_1) * cos(phi) + sin(phi_1) * sin(phi) * cos(lamb - lamb_0)
		k_prime_dr = k_prime_dr_1 + k_prime_dr_2

		x = k_prime * sin(phi) * sin(lamb - lamb_0)
		y = k_prime * (sin(phi_1) * cos(phi) - cos(phi_1) * sin(phi) * cos(lamb-lamb_0) )
		'''

		#picked from SPTpol
		reso_rad = np.radians(map_resol/60.)
		min_pos = -0.5*mapshape*reso_rad  # minimum position (left or bottom edge) in radians
		min_pos_y, min_pos_x = min_pos  # gets overwritten in BICEP projection since x and y resolutions are different

		# Then shift the radians position and divide by the resolution
		round = 1
		if round:
			# without the floor function, int conversion just chops off fractional part.
			# slightly out-of-range negative values will get incorrectly put in pixel 0,
			#  and the validity checks below will never find them.
			x = np.asarray(np.floor((x - min_pos_x)/reso_rad), dtype=int)
			y = np.asarray(np.floor((y - min_pos_y)/reso_rad), dtype=int)
		else:
			x = ne.evaluate("((x - min_pos_x)/reso_rad) - bin_center_zero*0.5")
			y = ne.evaluate("((y - min_pos_y)/reso_rad) - bin_center_zero*0.5")


		return np.asarray( [x, y] )

	def fn_pix2ang(self, pixels, ref_angles = None):

		"""
		Based on: http://mathworld.wolfram.com/LambertAzimuthalEqual-AreaProjection.html
		Should be very similar to the one on SPTpol pipeline
		"""
		
		sin,cos = np.sin, np.cos

		if ref_angles == None:
			ref_angles = self.map_centre

		#convert to radians
		ref_angles = np.radians(ref_angles)


	def fn_read_cluster_catalogue(self,clustercatfile = None, SNRlevel = None, get_pixel_info = None):

		"""
		reads the cluster catalogue file and returns RA, DEC information
		"""

		if clustercatfile == None:
			clustercatfile = self.clustercatfile

		if clustercatfile.split('/')[-1] == self.clustercatfile.split('/')[-1]:
			cluscat = np.loadtxt(clustercatfile, usecols = [2,3,4,6,7])
			ra,dec,conf,photz,specz = np.asarray( zip(*cluscat) )
		elif clustercatfile.find('combined_cluster_catalogue.pkl.gz')>-1:
			clusdic = pickle.load(gzip.open(clustercatfile,'rb'))
			cluscat = np.asarray( clusdic['final_clus_cat'].values() )
			ra,dec,conf,photz,specz = cluscat[:,0], cluscat[:,1], cluscat[:,2], cluscat[:,4], cluscat[:,5]
			cluscat = np.asarray( [ra,dec,conf,photz,specz] ).T
		elif clustercatfile.find('2500d_cluster_sample_fiducial_cosmology.fits')>-1:
			sptszfits = fits.open(clustercatfile)[1]
			sptsz = sptszfits.data
			ra = np.asarray( map(lambda x: x[1], sptsz) )
			dec = np.asarray( map(lambda x: x[2], sptsz) )
			conf = np.asarray( map(lambda x: x[3], sptsz) )
			specz = np.asarray( map(lambda x: x[8], sptsz) )
			photz = specz
			cluscat = np.asarray( [ra,dec,conf,photz,specz] ).T

		if SNRlevel == None:
			SNRlevel = self.clus_snrlevel

		cluscat = cluscat[conf>=SNRlevel]
		ra,dec,conf,photz,specz = np.asarray( zip(*cluscat) )

		#look for specz, then photoz
		z = np.copy(specz)
		no_spec_z = np.where(z==0)
		z[no_spec_z] = photz[no_spec_z]

		if not get_pixel_info:
			return ra,dec,z,conf

		"""
		tries to find the pixels corresponding to cluster centers - useful for cluster cutout extraction
		More info here:
		http://mathworld.wolfram.com/LambertAzimuthalEqual-AreaProjection.html
		"""

		#pixels = self.fn_ang2pix([ra,dec]) #yet to be complete

		import sky_local

		if self.degrade_fac <> None:
			map_resol *= self.degrade_fac
			mapshape /= self.degrade_fac
			
		pixels = sky_local.ang2Pix([ra,dec], self.map_centre, self.map_resol, self.mapshape, proj = self.map_proj, round=True, bin_center_zero=True, return_validity=False, use_c_code=False)

		y, x = pixels

		return ra, dec, z, conf, x, y

##########################################################################################
##########################################################################################
##########################################################################################

class FB():

	def __init__(self):
		ombh2 = 0.022
		omb0 = ombh2/param_dict['h']**2.
		#param_dict['omega_m'] = 0.3089
		from astropy.cosmology import FlatLambdaCDM
		cosmo = FlatLambdaCDM(H0 = param_dict['h']*100., Om0 = param_dict['omega_m'], Ob0 = omb0)

		self.rho_bar = cosmo.critical_density(0.).to('M_sun/Mpc3').value * cosmo.Om0
		self.kmin = 1e-4
		self.kmax = 40.
		self.h = cosmo.h
		self.nk = 1000
		self.ombh2 = 0.022#cosmo.Ob0 * cosmo.h**2.
		self.omch2 = cosmo.Om0 * cosmo.h**2.
		self.ns = 0.965

		pars = camb.CAMBparams()
		pars.set_cosmology(H0=self.h * 100, ombh2=self.ombh2, omch2=self.omch2)
		pars.set_dark_energy() #re-set defaults
		pars.InitPower.set_params(ns=self.ns)
		self.cambmatterpower = camb.get_matter_power_interpolator(pars, hubble_units = 0, kmax = self.kmax, k_hunit = 0, nonlinear = 0)

	#def pkz(self, z, k):
	#	return self.cambmatterpower.P(z, k)
	
	def sigma_Rz(self, R, z=0.):
		""" 
		Computes the square root of the variance of matter fluctuations within a sphere of radius R [Mpc]

		\sigma^2(R)= \frac{1}{2 \pi^2} \int_0^\infty \frac{dk}{k} k^3 P(k,z) W^2(kR)
	
		where
	
		W(kR) = \\frac{3j_1(kR)}{kR}
		"""
		import scipy.integrate as integrate
		if np.isscalar(R) or (np.size(R) == 1):

			def W_k_tophat(k):
				""" 
				Returns the Fourier Transform of a tophat window function. 
				"""
				return 3./k**3*(np.sin(k) - k*np.cos(k))

			def int_sigma(logk):


				k  = np.exp(logk)
				kR = k * R
				W  = W_k_tophat(kR)#3.0 * (np.sin(kR) - kR * np.cos(kR)) / kR**3
				#pk = self.pkz.P(z,k)
				
				#20171010	
				pk = self.cambmatterpower.P(z, k)
				return k**3 * W**2 * pk


			lnks   = np.linspace(np.log(self.kmin), np.log(self.kmax), self.nk)
			#lnks   = np.logspace(-4.,np.log10(40.), 10000)
			#pks = fb.cambmatterpower.P(0, lnks)
			sigma2 = integrate.simps( int_sigma(lnks), x=lnks)
			sigma2 /= (2.*np.pi**2)

			return np.sqrt(sigma2)
			# return np.sqrt(1.0/(2.0*np.pi**2.0) * integrate.romberg(int_sigma, np.log(self.kmin), np.log(self.kmax)))
		else:
			return np.asarray([ self.sigma_Rz(Rs) for Rs in R ])

	def sigma_Mz(self, M, z=0.):
		"""         
		Computes the square root of the variance of matter fluctuations within a sphere of mass M [M_sun]
		"""
		R = self.M2R(M)#, z=z)

		return self.sigma_Rz(R, z=z)

	def M2R(self, M):#, z=None): # [Mpc]
		"""
		Lagrangian scale R [Mpc] of density of mass M [M_sun]
		FIXME: check h factor
		"""
		#return ((3.* M)/(4.*np.pi * self.rho_bar(0.)))**(1./3.)
		return ((3.* M)/(4.*np.pi * self.rho_bar))**(1./3.)

	def R2M(self, R):#, z=None): # [M_sun]
		"""
		FIXME: check h factor
		"""
		#return 4./3.*np.pi * self.rho_bar(0.) * R**3
		return 4./3.*np.pi * self.rho_bar * R**3

	def delta_c(self):
		# FIXME: valid for flat-only cosmology, see NFW97
		#return 0.15*(12.0*np.pi)**(2.0/3.0)
		return 1.686
	
	def nu_M(self, M, z=0.):
		"""
		Returns the normalized mass overdensity as function of mass [M_sun] at a given redshift (z=0 by default)
		"""
		return self.delta_c() / self.sigma_Mz(M, z=z)


class FB():

	def __init__(self, param_dict):
		import camb
		ombh2 = 0.022
		omb0 = ombh2/param_dict['h']**2.
		#param_dict['omega_m'] = 0.3089
		from astropy.cosmology import FlatLambdaCDM
		cosmo = FlatLambdaCDM(H0 = param_dict['h']*100., Om0 = param_dict['omega_m'], Ob0 = omb0)

		self.rho_bar = cosmo.critical_density(0.).to('M_sun/Mpc3').value * cosmo.Om0
		self.kmin = 1e-4
		self.kmax = 40.
		self.h = cosmo.h
		self.nk = 1000
		self.ombh2 = 0.022#cosmo.Ob0 * cosmo.h**2.
		self.omch2 = cosmo.Om0 * cosmo.h**2.
		self.ns = 0.965

		pars = camb.CAMBparams()
		pars.set_cosmology(H0=self.h * 100, ombh2=self.ombh2, omch2=self.omch2)
		pars.set_dark_energy() #re-set defaults
		pars.InitPower.set_params(ns=self.ns)
		self.cambmatterpower = camb.get_matter_power_interpolator(pars, hubble_units = 0, kmax = self.kmax, k_hunit = 0, nonlinear = 0)

	#def pkz(self, z, k):
	#	return self.cambmatterpower.P(z, k)
	
	def sigma_Rz(self, R, z=0.):
		""" 
		Computes the square root of the variance of matter fluctuations within a sphere of radius R [Mpc]

		\sigma^2(R)= \frac{1}{2 \pi^2} \int_0^\infty \frac{dk}{k} k^3 P(k,z) W^2(kR)
	
		where
	
		W(kR) = \\frac{3j_1(kR)}{kR}
		"""
		import scipy.integrate as integrate
		if np.isscalar(R) or (np.size(R) == 1):

			def W_k_tophat(k):
				""" 
				Returns the Fourier Transform of a tophat window function. 
				"""
				return 3./k**3*(np.sin(k) - k*np.cos(k))

			def int_sigma(logk):


				k  = np.exp(logk)
				kR = k * R
				W  = W_k_tophat(kR)#3.0 * (np.sin(kR) - kR * np.cos(kR)) / kR**3
				#pk = self.pkz.P(z,k)
				
				#20171010	
				pk = self.cambmatterpower.P(z, k)
				return k**3 * W**2 * pk


			lnks   = np.linspace(np.log(self.kmin), np.log(self.kmax), self.nk)
			#lnks   = np.logspace(-4.,np.log10(40.), 10000)
			#pks = fb.cambmatterpower.P(0, lnks)
			sigma2 = integrate.simps( int_sigma(lnks), x=lnks )
			sigma2 /= (2.*np.pi**2)

			return np.sqrt(sigma2)
			# return np.sqrt(1.0/(2.0*np.pi**2.0) * integrate.romberg(int_sigma, np.log(self.kmin), np.log(self.kmax)))
		else:
			return np.asarray([ self.sigma_Rz(Rs) for Rs in R ])

	def sigma_Mz(self, M, z=0.):
		"""         
		Computes the square root of the variance of matter fluctuations within a sphere of mass M [M_sun]
		"""
		R = self.M2R(M)#, z=z)

		return self.sigma_Rz(R, z=z)

	def M2R(self, M):#, z=None): # [Mpc]
		"""
		Lagrangian scale R [Mpc] of density of mass M [M_sun]
		FIXME: check h factor
		"""
		#return ((3.* M)/(4.*np.pi * self.rho_bar(0.)))**(1./3.)
		return ((3.* M)/(4.*np.pi * self.rho_bar))**(1./3.)

	def R2M(self, R):#, z=None): # [M_sun]
		"""
		FIXME: check h factor
		"""
		#return 4./3.*np.pi * self.rho_bar(0.) * R**3
		return 4./3.*np.pi * self.rho_bar * R**3

	def delta_c(self):
		# FIXME: valid for flat-only cosmology, see NFW97
		#return 0.15*(12.0*np.pi)**(2.0/3.0)
		return 1.686
	
	def nu_M(self, M, z=0.):
		"""
		Returns the normalized mass overdensity as function of mass [M_sun] at a given redshift (z=0 by default)
		"""
		return self.delta_c() / self.sigma_Mz(M, z=z)

############################################################################################################
############################################################################################################
############################################################################################################

from colossus.cosmology import cosmology
from colossus.halo import concentration, mass_defs
cosmology.setCosmology('planck15')
import Nikhel_Cosmology1 as nc
import Nikhel_Cosmology1_sp as nc1

from scipy import integrate
from scipy.io import readsav
import numpy as np, pickle, sys, gzip, os, glob, time
import scipy.optimize as optimize
import astropy.io.fits as fits
if str(os.getcwd()).find('sraghunathan') > -1:
	from pylab import *
	from IPython import embed;
	import py_ini; cmap_planck = py_ini.get_planck_cmap()

