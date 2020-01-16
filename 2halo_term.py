def fn_get_halo_bias_terms(delta):
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



def fn_get_halo_bias(M, z, delta = 200.):

	"""
	http://iopscience.iop.org/article/10.1088/0004-637X/724/2/878/pdf
	Eq. 5
	"""

	nu = fb.nu_M(M, z)
	delta_c =fb.delta_c()

	A, B, C, alpha, beta, gamma = fn_get_halo_bias_terms(delta)
	t1 = 1.
	t2 = A * ( (nu**alpha) / ( nu**alpha + delta_c**alpha) ) 
	t3 = B * nu**beta
	t4 = C * nu**gamma

	b_nu = t1 - t2 + t3 + t4

	return b_nu, nu

def fn_kappa_two_halo_term(theta, M, z, lmin = 0, lmax = 20000):
	import astropy.constants as const
	import scipy.special as special
	import scipy.integrate as integrate

	b_M_z, nu = fn_get_halo_bias(M, z)
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

	kappa_two_halo = np.zeros(theta.shape)
	t2_nr = rho_matter_z * b_M_z
	t2_dr = (1. + z)**3. * sigma_c.value * D_L.value**2.
	t2 = t2_nr / t2_dr

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

		def fn_for_int_for_simps(el, bessel_order = 0):
			t1 = el * special.jn(bessel_order, el * theta_val) /2. /np.pi
			kl = el / (1.+z) / D_L.value
			t3 = fb.cambmatterpower.P(z, kl)
			return t1 * t2 * t3

		els = np.arange(lmin, lmax)
		kappa_two_halo[tcnt] = integrate.simps( fn_for_int_for_simps(els), x=els )
		#print kappa_two_halo[tcnt];sys.exit()

	return kappa_two_halo
