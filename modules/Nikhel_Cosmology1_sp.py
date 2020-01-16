import numpy as np
##import matplotlib.pyplot as pl
##import pyfits as pyf
from math import *
##import cosmolopy.distance as cd
from scipy.integrate import quad

"""
This module contains following tools
1) M500_to_R500
2) R500_to_M500
3) Angular_diamter_distance using cosmolopy
4) Generalized NFW profile with Arnaud et al. 2010 parameters
5) Integrated SZ Y parameter (M500_to_Yx500_steradian)
6) Integrated SZ Y parameter (M500_to_Yx500_arcm2)
7) SZE flux in mJy for a given survey frequency neu
8) Radiu_to_theta_arcm
9) myrange
"""


#  -------1
# M500 to R500
def M500_to_R500 (Mass,  redshift,  Om_M,   Om_L,  H0):
		a_3              =       (1.+redshift)**3.
		Hz               = H0*np.sqrt((Om_M*a_3)+Om_L)
		G_c              = 4.302e-9                                                     # Gravitational constant
		cri_den  = 3.*Hz*Hz/8./pi/G_c                   # Critical density (redshift)
		R500     = (3.*Mass/(500.*cri_den*4.*pi))**(1./3.)
		return R500


#  -------2     
# R500 to M500
def R500_to_M500 (R500,  redshift,  Om_M,   Om_L,  H0):
		a_3              =       (1.+redshift)**3.
		Hz               = H0*np.sqrt((Om_M*a_3)+Om_L)
		G_c              = 4.302e-9                                                     # Gravitational constant
		cri_den  = 3.*Hz*Hz/8./pi/G_c                   # Critical density (redshift)
		Mass     = 4./3. * pi * R500**3. *500.*cri_den
		return Mass

#  -------3     
# Angular diameter distance for flat cosmlogy   
def Angular_diameter_distance (redshift,  Om_M,   Om_L,  H0):
		'''
		cosmo = {'omega_M_0' : Om_M, 'omega_lambda_0' : Om_L, 'h' : H0/100.}
		cosmo = cd.set_omega_k_0(cosmo)
		d_M = cd.angular_diameter_distance(redshift, **cosmo)
		'''

		#SR edits
		def e_z(z):
			omega_k = 0.
			return (1/np.sqrt(Om_M*((1+z)**3)+ omega_k*((1+z)**2) + Om_L))

		e_z_int, e_z_int_err = quad(e_z, 0.,redshift)
		h  =H0/100.
		d_c, d_t, d_a, d_l, dm, v_c, dv_c = fn_distances_volume(e_z_int, redshift, h, Om_M, 0., Om_L)
		d_M = d_a

		return d_M


#  -------4     
# Generalized Pressure profile [P(u)] Arnaud et al. 2010 
def Arnaud_GNFW (u, redshift, Om_M, Om_L, H0):
		h_70     = H0/70.
		a_3              =       (1.+redshift)**3.
		Hz               = H0*np.sqrt((Om_M*a_3)+Om_L)
		P0               = 8.403* (h_70)**(-3./2.)
		c_500    = 1.177
		gamma = 0.3081
		alph     = 1.0510
		beta     = 5.4905
		h_z        = Hz/H0
		P_u      = P0/( (c_500*u)**gamma * (1. + (c_500*u)**alph)**( (beta-gamma)/alph) )
		return P_u



#  -------5
# SZE Y(x R500) in Steradian from Arnaud et al. 2010 eqtn (25)
# Define the integrand I_u = 3 P(u) u^2 
def I_u (u,     redshift, Om_M, Om_L, H0):
		return 3.* Arnaud_GNFW (u, redshift, Om_M, Om_L, H0)*u**2
# Define the integrand I_k = 3 P(u) sqrt(u^2 - x^2) u   
def I_k (u, x, redshift, Om_M, Om_L, H0):
		return 3.* Arnaud_GNFW (u, redshift, Om_M, Om_L, H0)*np.sqrt(u**2. - x**2.)*u

### SR edits
def fn_distances_volume(e_z_int, z, h, omega_m, omega_k, omega_lambda, dz = 0.01):

	#H0 = 100.0    # HUBBLE CONSTANT = 100h Km/s/Mpc

	# HUBBLE CONSTANT IN STANDARD UNITS - DIMENSIONS sec-1
	H0 = h * 100.
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
### SR edits

def M500_to_Yx500_steradian (x,  Mass,  redshift,  Om_M,   Om_L,  H0):


		h_70    =       H0/70.
		a_3             =       (1.+redshift)**3.
		I_5             =       quad (I_u, 0., 5., args = (   redshift,Om_M, Om_L, H0) )
		I_x             =  quad (I_k,  x, 5., args = (x,redshift,Om_M, Om_L, H0) )
		J_x             =       I_5[0] - I_x[0]
		B_x             =       2.925e-5 * J_x / h_70
		alph            =       1.78
		Dist            =       Angular_diameter_distance (redshift, Om_M, Om_L, H0)

		Ez                      =       np.sqrt((Om_M*a_3)+(Om_L))
		Yx500   =       B_x * (Mass/(3.e14/h_70))**alph / Dist**2. / Ez**(-2./3.)
		return Yx500


#  -------6     
# SZE Y(x R500) in arcmin^2 from Arnaud et al. 2010 eqtn (25)
def M500_to_Yx500_arcm2 (x,  Mass):
		redshift  = 0.7
		Om_M = 0.307320
		Om_L = 0.692680
		H0 = 0.6774* 100.
		
		return M500_to_Yx500_steradian (x,  Mass,  redshift,  Om_M,   Om_L,  H0) * 3282.8 * 3600.



#  -------7
#  SZE flux in mJy for a given survey frequency neu     
def SZE_flux_frequency_mJy (neu,  x,  Mass,  redshift,  Om_M,   Om_L,  H0):
		#I_CMB          = 2. * (1.38e-16*2.725)**3. / (1.e-27*4.136*1.602*299792*100.)**2.
		I_CMB           = 2.7033e8                      # Jy/str.
		f_neu           = (neu/56.78)**4. * np.exp(neu/56.78) / (np.exp(neu/56.78)-1.)**2.
		g_neu           = (neu/56.78)*(1. / np.tanh( (neu/56.78 ) / 2.) ) - 4.
		Yx500           = M500_to_Yx500_steradian (x, Mass, redshift, Om_M, Om_L, H0)
		SZflux          = I_CMB * f_neu * g_neu * Yx500 * 1.e3                          # mJy
		return SZflux



#  -------8
# Radius to theta
def radius_to_theta_arcm (radius,  redshift,  Om_M,   Om_L,  H0):
		Dist                    =       Angular_diameter_distance (redshift, Om_M, Om_L, H0)
		theta           =       radius / Dist
		return theta * 180. / pi * 60.



#  -------9
# my range
def my_range(start, end, step):
	while start <= end:
		yield start
		start += step

