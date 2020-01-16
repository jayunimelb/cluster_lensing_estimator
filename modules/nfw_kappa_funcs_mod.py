import numpy as np
from scipy import integrate

class kappa_stuff():

	def __init__(self):
		#self.h = 0.7
		#self.omega_lambda = 0.7
		#self.omega_m = 0.3
		#self.omega_r = 0.
		#self.omega_k = 0.
		#self.z_cmb = 1100.0
		#self.z_lss = 1089.0
                self.Mpc_to_metres = 3.086e22
		self.solar_mass_to_kgs = 1.9891e+30 #kg
		#self.vellight = 299797.0 #km/s
		self.vellight = 299792.458 #km/s
		self.G = 6.674e-11 #Nm^2 / kg^2

	def fn_update_cosmo_params(self,param_dict):
		
		"""
		supply a dic here to update cosmology
		example: param_dic = {'h':0.74,'omega_lambda':0.73}
		"""

		for keys in sorted( param_dict.keys() ):
			cmd = 'self.%s = param_dict[keys]' %(keys)
			exec(cmd)

        def fn_cosmo_dist(self,z,return_in_metres = 0):

		"""
		routine to calculate cosmological distances in Mpc. 
		Can also return in standard units. Just turn return_in_metres = 1
		based on: Hogg 2000 (arXiv:astro-ph/9905116)
		"""

	        vellight = self.vellight
                omega_lambda = self.omega_lambda #D.E
                omega_m = self. omega_m #matter
		try: #radiation
	                self.omega_r
		except AttributeError:
	                self.omega_r = 0.

                omega_r = self.omega_r

                omega = omega_lambda + omega_m + omega_r #total density
                omega_k = 1 - omega #curvature

		#hubble constant stuff
                H_0 = 100 * self.h
                H_0_std = (H_0/(3.08568025 * 10**19)) #s-1
		self.H_0_std = H_0_std
		self.H_0 = H_0

	        # HUBBLE DISTANCE IN Mpc
	        d_h = vellight/H_0#_std

	        # INTEGRATION FUNCTION E(z)
                # SCIPY RETURNS THE INTEGRATION RESULT AS A TWO DIMENSIONAL ARRAY. FIRST ELEMENT IS THE RESULT. SECOND ELEMENT IS THE ERROR ASSOCIATED
	        def e_z(z):
		        return (1.0/np.sqrt(omega_m*((1+z)**3)+ omega_k*((1+z)**2) + omega_lambda))

	        e_z_int, e_z_int_err = integrate.quad(e_z,0.,z)
	
	        #ALL DISTANCES ARE IN Mpc (SINCE I'VE CALCULATED HUBBLE DISTANCE ABOVE IN Mpc)
	        # TOTAL LINE-OF-SIGHT COMOVING DISTANCE
		
	        d_c = d_h * e_z_int

	        # TRANSVERSE COMOVING DISTANCE
	        if (omega_k==0.0): #FLAT
		        d_t = d_c
	        elif (omega_k>0.0): #POSITIVELY CURVED
		        d_t = d_h/np.sqrt(omega_k) * np.sinh(np.sqrt(omega_k)*d_c/d_h)
	        else: #NEGATIVELY CURVED
		        d_t = d_h/np.sqrt(abs(omega_k)) * np.sinh(np.sqrt(abs(omega_k))*d_c/d_h)

	        if (omega_lambda==0.0): #UNIVERSE WITHOUT COS. CONSTANT / DARK ENERGY
		        d_t = d_h * 2 *(2 - (omega_m *(1-z)) - ((2-omega_m) * (np.sqrt(1+(omega_m*z))))) / (omega_m**2 * (1+z))

	        # ANGULAR DIAMETER DISTANCE
	        d_a = d_t / (1+z)

	        # LUMINOSITY DISTANCE
	        d_l = (1+z) * d_t

	        # DISTANCE MODULUS
	        dm = 5.0 * np.log10(d_l*10**6/10) # 1 Mpc = 10**6 pc
	
                distances = np.asarray([d_h,d_c,d_t,d_a,d_l,dm])
                if return_in_metres:
			distances = distances * self.Mpc_to_metres

	        return distances

	def fn_d_a12(self,z1,z2,dist_reqd = 'comoving',return_in_metres = 0):

		"""
		Comoving / Angular diameter distance between 2 objects zt z1, z2 - important for gravitational lensing
		"""

		try:
			self.omega_k
		except AttributeError:
			self.omega_k = 0.

		omega_k = self.omega_k

		distances_z1 = self.fn_cosmo_dist(z1)
		if dist_reqd == 'comoving':
			d_m1 = distances_z1[1]
		else:
			d_m1 = distances_z1[3]
                d_h1 = distances_z1[0]

                distances_z2 = self.fn_cosmo_dist(z2)
		if dist_reqd == 'comoving':
			d_m2 = distances_z2[1]
		else:
			d_m2 = distances_z2[3]
		d_h2 = distances_z2[0]

		assert d_h1==d_h2
		d_h = d_h1

		assert dist_reqd == 'comoving' or dist_reqd == 'ang_dia'

		if dist_reqd == 'comoving':
			t1 = 1.
		else: #angular diamater distance
			t1 = 1. / (1+z2)
		t2 = d_m2 * np.sqrt( 1 +  (omega_k * d_m1**2 / d_h**2.)  )
		t3 = d_m1 * np.sqrt( 1 +  (omega_k * d_m2**2 / d_h**2.)  )

		d_a12 = t1 * (t2 - t3)
                if return_in_metres:
			distances = d_a12 * self.Mpc_to_metres

		return d_a12

	def fn_rho_cric(self,z):
		"""
		rho_cric = 3 * H**2. / ( 8 * pi * G )
		"""
                import numpy as np

		G = self.G
		H = self.fn_H(z)
		rho_cric = 3 * H**2. / ( 8 * np.pi * G )

		return rho_cric	

	def fn_r200(self,M_200,z,M_def = 200,dens = 'mean'): #M_200 in solar mass
		"""
		r_200 = M_200 / ( 4*pi/3 * 200 * rho_cric(z) )
		"""
		import numpy as np

		rho_cric = self.fn_rho_cric(z) #kg/m^3
		if dens == 'mean':
			r_200_cubed = M_200 / ( (4*np.pi/3) * M_def * rho_cric * self.omega_m)
		else:
	                r_200_cubed = M_200 / ( (4*np.pi/3) * M_def * rho_cric)

                r_200 = (r_200_cubed)**(1/3.) #metres

		return r_200

	def fn_g(self,x):
		
		t1 = 1/(x**2. - 1.)
		if x>1:
			t3 = 2./np.sqrt(x**2. - 1.) * np.arctan( np.sqrt( (x-1 ) / (x+1) ) )
			t2 = 1 - t3

			fval = t1 * t2

		elif x<1:
			t3 = 2./np.sqrt(1. - x**2.) * np.arctanh( np.sqrt( (1-x) / (1+x) ) )
			t2 = 1 - t3

			fval = t1 * t2

		elif x == 1:
			fval = 1./3.

		return fval

	def fn_f(self,c):

		return np.log(1+c) - (c / (1+c))

	def fn_kappa(self,theta,M_200,c,z_L,z_S):

		M_200_kg = M_200 * self.solar_mass_to_kgs
                r_200 = self.fn_r200(M_200_kg,0) / self.Mpc_to_metres #Mpc
		r_s = r_200 / c

                distances_z_S = self.fn_cosmo_dist(z_S)
                distances_z_L = self.fn_cosmo_dist(z_L)

		d_S = distances_z_S[1] #comoving LOS / transeverse - same for flat LCDM #14.345 Gpc - correct # same as Whu
		d_L = distances_z_L[1] #comoving LOS / transeverse - same for flat LCDM
		d_LS = self.fn_d_a12(z_L,z_S)

		'''
		from astropy.cosmology import FlatLambdaCDM
		cosmo = FlatLambdaCDM(H0 = self.h*100., Om0 = self.omega_m)
		d_L = cosmo.comoving_distance(z_L).value
		d_S = cosmo.comoving_distance(z_S).value
		d_LS = (cosmo.comoving_distance(z_S)-cosmo.comoving_distance(z_L)).value
		'''

		t1 = 3/(4 * np.pi * self.vellight**2.) 
		t2 = self.H_0 * d_L * d_LS  * (1+z_L) / d_S

		rho_cric_kg_Mp3 = self.fn_rho_cric(0) / (1./self.Mpc_to_metres)**3.

		t3 = M_200_kg * self.H_0 / (rho_cric_kg_Mp3 * r_s**2.)

		theta_s = r_s/d_L #~0.94 arcmins for 2.7e14 @z = 0.7 -> get the same answer at Hu et al. 2007

		if len(theta.shape)>1:
			nx,ny = theta.shape
			t4 = np.asarray( map(lambda x: self.fn_g(x/theta_s) / self.fn_f(c), theta.ravel()) ).reshape(nx,ny)
		else:
			t4 = np.asarray( map(lambda x: self.fn_g(x/theta_s) / self.fn_f(c), theta.ravel()) )

		return t1 * t2 *t3 * t4

	####################################################################################################################################
	####################################################################################################################################
	####################################################################################################################################
	## below functions are not being used

	def fn_A(self,M_200,c):
		"""
		A = (M_200 * c^2) / ( 4 * pi * [ln(1+c) - (c/1+c)] )
		"""
		nr = M_200 * (c**2.)
		dr = 4 * np.pi * ( np.log(1+c) - (c / (1+c)) )

		return nr / dr

	def fn_H(self,z):
		"""
		H(a) = H_0 * np.sqrt ( omega_lambda + omega_m / a^3 + omega_r / a^4 - (omega - 1) / a^2)
		"""

		omega_lambda = self.omega_lambda
		omega_m = self.omega_m
		try:
	                self.omega_r
		except AttributeError:
	                self.omega_r = 0.

		omega_r  = self.omega_r 

		omega = omega_lambda + omega_m + omega_r
		omega_k = 1 - omega

		H_0 = 100 * self.h
		H_0_std = (H_0/(3.08568025 * 10**19)) #s-1
		a  = 1./ (1+z)
	
		H = H_0_std * np.sqrt ( omega_lambda + (omega_m /a**3.) + (omega_r / a**4.) + (omega_k/a**2.) )

		return H


