import numpy as np, sys
from astropy import constants as const
from astropy import units as u
from astropy import coordinates as coord
import pdb
from scipy import integrate
from pylab import *

#Put multiple clusters into kappa map
def get_NFW_kappa_fullmap(param_dict, ra_map, dec_map, ra_list, dec_list, M_200_list, c_200_list, z_L_list, z_LSS, mass_def, rho_def, theta_max = -1., Tcmb = 2.725, profile = 'NFW', ret_nfw_dic = 0, just_surface_mass_density = 0, des_offsets_correction = 0, richval = None):

    from astropy.cosmology import FlatLambdaCDM
    cosmo = FlatLambdaCDM(H0 = param_dict['h']*100., Tcmb0 = Tcmb, Om0 = param_dict['omega_m'])
    N_clusters = len(ra_list)
    kappa_map = np.zeros(ra_map.shape)
    map_coords = coord.SkyCoord(ra = ra_map*u.degree, dec = dec_map*u.degree)
    for ci in xrange(0,N_clusters):

        #print "ci = ", ci
        cluster_coords = coord.SkyCoord(ra = ra_list[ci]*u.degree, dec = dec_list[ci]*u.degree)
        theta_map = map_coords.separation(cluster_coords).value*(np.pi/180.)

        ##kappa_map += get_NFW_kappa(cosmo, theta_map, M_200_list[ci], c_200_list[ci], z_L_list[ci], z_LSS, mass_def, rho_def, theta_max)

	if profile == 'NFW' and just_surface_mass_density == 0: #original - uses analytic expressions for nFW kappa
		if not ret_nfw_dic:
		        kappa_map += get_NFW_kappa(cosmo, theta_map, M_200_list[ci], c_200_list[ci], z_L_list[ci], z_LSS, mass_def, rho_def, theta_max)
		else:
		        RET = get_NFW_kappa(cosmo, theta_map, M_200_list[ci], c_200_list[ci], z_L_list[ci], z_LSS, mass_def, rho_def, theta_max, ret_nfw_dic = ret_nfw_dic, des_offsets_correction = des_offsets_correction, richval = richval, param_dict = param_dict)
			kappa_map += RET[0]
			nfw_dic = RET[1]
	else: #20161231 - kappa from numeric integration methods
		class_kappa = class_convergence()
		retvals = class_kappa.fn_get_kappa_general(cosmo, theta_map, M_200_list[ci], c_200_list[ci], z_L_list[ci], z_LSS, mass_def, rho_def, theta_max = theta_max, profile = profile, just_surface_mass_density = just_surface_mass_density)
		if not just_surface_mass_density:
			kappa_map_tmp, nfw_dic = retvals
		else:
			just_surface_mass_density = retvals
			return just_surface_mass_density
		kappa_map += kappa_map_tmp


	'''
	#from pylab import *
	#imshow(kappa_map);colorbar();show();quit()

	print get_NFW_kappa(cosmo, [np.radians(0.94/60.)], M_200_list[ci], c_200_list[ci], z_L_list[ci], z_LSS, mass_def, rho_def)

	from pylab import *
	subplot(121);imshow(kappa_map);colorbar();title('$k_{EBAX}$')

	import nfw_kappa_funcs_mod as nfw_mod
	kappa_stuff = nfw_mod.kappa_stuff()
	kappa_stuff.fn_update_cosmo_params(param_dict)
	whu_kappa = kappa_stuff.fn_kappa(np.asarray( [np.radians(0.94/60.)] ), M_200_list[ci], c_200_list[ci], z_L_list[ci], z_LSS)
	print whu_kappa
	kappa_map_2 = kappa_stuff.fn_kappa(theta_map, M_200_list[ci], c_200_list[ci], z_L_list[ci], z_LSS)
	subplot(122);imshow(kappa_map_2);colorbar();title('$k_{SR}$')

	show();quit()
	'''
    if ret_nfw_dic == 0:
		return kappa_map
    else:
		return kappa_map, nfw_dic

def get_NFW_kappa(cosmo, theta, mass, concentration, z_L, z_S, mass_def, rho_def, theta_max = -1, ret_nfw_dic = 0, des_offsets_correction = 0, richval = None, param_dict = None, totalclus = None):
    nfw_dic = {}
    mass = mass*u.Msun
    if (rho_def == 'crit'):
        rho_c_z = cosmo.critical_density(z_L)
    elif (rho_def == 'mean'):
        rho_c_z = cosmo.Om(z_L)*cosmo.critical_density(z_L)
    else:
        print "rho definition not specified correctly in cluster profile"
        assert(0)

    #NFW profile properties
    delta_c = (mass_def/3.)*(concentration**3.)/( np.log(1.+concentration) - concentration/(1.+concentration) )
    r_v = (((mass/(mass_def*4.*np.pi/3.))/rho_c_z)**(1./3.)).to('Mpc')
    #r_v = (((mass/(mass_def*4.*np.pi/3.))/ (.3*rho_c_z) )**(1./3.)).to('Mpc') #change of 20160531

    r_s = r_v/concentration
    #Angular diameter distances
    D_L = cosmo.comoving_distance(z_L)/(1.+z_L)
    D_S = cosmo.comoving_distance(z_S)/(1.+z_S)
    D_LS = (cosmo.comoving_distance(z_S)-cosmo.comoving_distance(z_L))/(1.+z_S)
    #Normalization of kappa
    Sigma_c = (((const.c.cgs**2.)/(4.*np.pi*const.G.cgs))*(D_S/(D_L*D_LS))).to('M_sun/Mpc2')
    '''
    D_L = cosmo.comoving_distance(z_L)
    D_S = cosmo.comoving_distance(z_S)
    D_LS = (cosmo.comoving_distance(z_S)-cosmo.comoving_distance(z_L))
    #Normalization of kappa
    Sigma_c = (((const.c.cgs**2.)/(4.*np.pi*const.G.cgs))*(D_S/(D_L*D_LS * (1+z_L)))).to('M_sun/Mpc2')
    '''
    #Useful variables
    def fn_get_g_theta(x):
        g_theta = np.zeros(x.shape)
        gt_one = np.where(x > 1.0)
        lt_one = np.where(x < 1.0)
        #eq_one = np.where(np.abs(x - 1.0) < 1.0e-5)
        eq_one = np.where(x == 1.)
        ##eq_zeros = np.where(x == 0.)
        g_theta[gt_one] = (1./(x[gt_one]**2. - 1))*(1. - (2./np.sqrt(x[gt_one]**2. - 1.))*np.arctan(np.sqrt((x[gt_one]-1.)/(x[gt_one]+1.))) )#.value)
        g_theta[lt_one] = (1./(x[lt_one]**2. - 1))*(1. - (2./np.sqrt(1. - x[lt_one]**2.))*np.arctanh(np.sqrt((1. - x[lt_one])/(x[lt_one]+1.))) )#.value)
        ##g_theta[eq_zeros] = (1./(x[eq_zeros]**2. - 1))*(1. - (2./np.sqrt(1. - x[eq_zeros]**2.))*np.arctanh(np.sqrt((1. - x[eq_zeros])/(x[eq_zeros]+1.))) )#.value)
        g_theta[eq_one] = 1./3.

        return g_theta


    R = D_L*theta
    x = (R/r_s).value
    ##from IPython import embed; embed()
    g_theta = fn_get_g_theta(x)
    ###g_theta[g_theta == np.inf] = 0.

    ###from IPython import embed; embed()

    #Functional form of profile
    #Projected mass
    Sigma = ((2.*r_s*delta_c*rho_c_z)*g_theta).to('M_sun/Mpc2')

    #Final kappa
    kappa = Sigma/Sigma_c

    if des_offsets_correction: ##DES offset correction
    	print '\n\t\t\tKappa offsets correction'
    	assert richval <> None and param_dict<>None and totalclus<>None

        def fn_sigma_offsets_correction_int(phi, Rval, Roffval):
            #print Rval, Roffval
            Rmis = np.sqrt(Rval**2. + Roffval**2. + 2 * Rval * Roffval * np.cos(phi))
            x = np.asarray([(Rmis/r_s).value])
            g_theta_mis = fn_get_g_theta(x)[0]
            Sigma_mis_val = (1./2./np.pi) * ( ((2.*r_s*delta_c*rho_c_z)*g_theta_mis).to('M_sun/Mpc2') ).value
            return Sigma_mis_val

        param_dict['fmis'] = 0.22
        param_dict['lncmis'] = -1.13
        fmis = param_dict['fmis'] ##0.22 #+/- 0.11
        lncmis = param_dict['lncmis'] ##-1.13 #+/-0.22
        cmis = np.exp(lncmis)
        R_lambda = (richval/100.)**0.2 / param_dict['h']
        offset_scatter = cmis * R_lambda
        #from IPython import embed; embed()
        Roff_dist = np.random.rayleigh(offset_scatter, totalclus) ## * theta
        Roff = Roff_dist[np.random.randint(len(Roff_dist))]

        R = (D_L*theta).value.flatten()
        if len(Roff)>1:
        	Roff = Roff.flatten()

        Sigma_mis_int = np.zeros(R.shape)
        phis = np.linspace(0, 2*np.pi, 1000)
        for ii in range(len(R)):
            ##if ii % 500 == 0: print ii
            ## Sigma_mis_int[ii] = integrate.quad(fn_sigma_offsets_correction_int, 0.,2 * np.pi, args = (R[ii], Roff[ii]))[0]
            #Sigma_mis_int[ii] = integrate.simps( fn_sigma_offsets_correction_int(phis, R[ii], Roff[ii]), x=phis )
            Sigma_mis_int[ii] = integrate.simps( fn_sigma_offsets_correction_int(phis, R[ii], Roff), x=phis )

        Sigma_mis = Sigma_mis_int.reshape(theta.shape)

        #Final kappa
        kappa_mis = Sigma_mis/Sigma_c.value

        kappa_full = (1-fmis) * kappa + ( fmis * kappa_mis )

        kappa = np.copy(kappa_full)

	#from IPython import embed; embed()
	#sys.exit()
    ##imshow(kappa, interpolation = 'bicubic');colorbar();show();sys.exit()

    '''
    from pylab import *
    subplot(221);plot(np.degrees(theta[len(theta)/2,:])*60., np.cumsum(kappa[len(theta)/2,:]),'r');
    subplot(222);plot(np.degrees(theta[len(theta)/2,:])*60., kappa[len(theta)/2,:],'g')
    subplot(223);pcolor(kappa);colorbar()
    show();quit()
    '''

    if (theta_max > 0):
        beyond_theta_max = np.where(theta > theta_max)
        kappa[beyond_theta_max] = 0.0
    nfw_dic['r_s'] = r_s
    nfw_dic['D_L'] = D_L

    if not ret_nfw_dic:
	    return kappa.value
    else:
	    return kappa.value, nfw_dic

################################################################################################################################################
################################################################################################################################################
################################################################################################################################################

class class_convergence():

    """
    added on 20161231 - SR
    Perform numeric integration of surface mass density to obtain kappa
    """

    def __init__(self):

	self.epsrel = 1.49e-6
	self.epsabs = 0.

    def fn_get_surface_mass_density(self, theta, D_L, profile = 'NFW', alpha_ein = 0.18, theta_max = -1.):

	assert profile == 'NFW' or profile == 'Einasto' or profile == 'gNFW_Moore' or profile == 'DandK'

	R = theta * D_L
	R = R.value.flatten()
        Sigma_int = np.zeros(len(R)) + 1e-10 #some small number
	mean_Sigma_num_int = np.zeros(len(R)) #not yet implemented

	concentration = self.concentration
	mass_def = self.mass_def
	rho_c_z = self.rho_c_z
	mass = self.mass

	#scale radius and R_200
	r_v = (((mass/(mass_def*4.*np.pi/3.))/rho_c_z)**(1./3.)).to('Mpc')
	r_s = r_v/concentration

	if profile == 'NFW':

	    """
	    #working with R instead of x = R/r_s. The profile is modified accordingly
	    rho_R_nr = delta_c * rho_c(z) * r_s**3
	    rho_R_dr = R * (r_s + R)**2.
	    rho_R = rho_R_nr / rho_R_dr

	    \begin{eqnarray}
	       \Sigma(x) & = &  2 \delta_c \rho_{crit}^z R_s^3 \int_0^{\infty} \rho(x, z)\ dz
	    \end{eqnarray}

	     Here $r = \sqrt{R^2 + z^2}$. Using change of variables we get $z = \sqrt{r^2 - R^2 }$ and $dz = \frac{rdr}{\sqrt{r^2 - R^2}}$

	    \begin{eqnarray}
	      \Sigma(x) =  2 \delta_c \rho_{crit}^z R_s^3 \int_R^{\infty} \frac{1}{r (R_s + r)^2} \frac{r\ dr}{\sqrt{r^2 - R^2}}
	    \end{eqnarray}

	    """

	    #NFW profile properties
	    delta_c = (mass_def/3.)*(concentration**3.)/(np.log(1.+concentration)-concentration/(1.+concentration))

	    #make sure that we get same delta_c using num. integration as well
            t1 = (mass_def/3.) * r_v.value**3.
	    t2 = 1./ r_s.value**3. / integrate.quad(lambda r : r / (r + r_s.value)**2., 0., r_v.value, epsabs=self.epsabs, epsrel=self.epsrel)[0]
	    delta_c_2 = t1 * t2 #should be same as delta_c and it is same!
	    assert round(delta_c,1) == round(delta_c_2,1)

	    norm_for_int = 2. * delta_c * rho_c_z.value * r_s.value**3.

	    for ii in range(len(R)):
		Sigma_int[ii] =  integrate.quad(lambda r : norm_for_int / ( r * (r_s.value + r )**2. ) * ( r / np.sqrt(r**2. - R[ii]**2.) ), R[ii], np.inf, epsabs=self.epsabs, epsrel=self.epsrel)[0]
	     #try to using map(lambda ) function. Minimal gain in processing time.
	     #Sigma_int = map(lambda z: integrate.quad(lambda r : norm_for_int / ( r * (r_s.value + r )**2. ) * ( r / np.sqrt(r**2. - z**2.) ), z, np.inf, epsabs=self.epsabs, epsrel=self.epsrel)[0], R)

	elif profile == 'gNFW_Moore': #arXiv:1005.0411 (Eq. 27); https://arxiv.org/abs/astro-ph/9903164

	    concentration = 2.5
	    r_s = r_v/concentration

	    #alpha, beta, gamma = 1., 3., 1. #will be similar to NFW
	    alpha, beta, gamma = 1.5, 3., 1.5 #arXiv:1005.0411 (Eq. 27); https://arxiv.org/abs/astro-ph/9903164

	    t1 = (mass_def/3.) * r_v.value**3.
	    t2 = 1. / integrate.quad(lambda r : r**2. * (2. / ( (r/r_s.value)**gamma * ( 1 + (r/r_s.value)**(3-gamma) ) ) ), 0., r_v.value, epsabs=self.epsabs, epsrel=self.epsrel)[0]
	    delta_c = t1 * t2

	    norm_for_int = 2. * delta_c * rho_c_z.value

	    for ii in range(len(R)):
		Sigma_int[ii] =  integrate.quad(lambda r : norm_for_int * (2. / ( (r/r_s.value)**gamma * ( 1 + (r/r_s.value)**(3-gamma) ) ) ) * ( r / np.sqrt(r**2. - R[ii]**2.) ), R[ii], np.inf, epsabs=self.epsabs, epsrel=self.epsrel)[0]

	elif profile == 'Einasto':

		#NFW profile properties
		r_minus2 = r_s.value
		#alpha = 0.18 #gives something similar to NFW
		#alpha = .3 #galaxy cluster #https://arxiv.org/pdf/1302.0288v2.pdf
		alpha = alpha_ein

		t1 = (mass_def/3.) * r_v.value**3.
		t2 = 1. / integrate.quad(lambda r : r**2 * np.exp ( - (2/alpha) * ( (r/r_minus2)**alpha - 1)  ), 0., r_v.value, epsabs=self.epsabs, epsrel=self.epsrel)[0]

		delta_c = t1 * t2
		norm_for_int = 2. * delta_c * rho_c_z.value
		### import time; start = time.time(); print '\n\n', theta_max, start
		### from IPython import embed; embed()
		if theta_max<>-1.:
			### print '\nEntering faster execution\n'
			theta_arcmins_flatten = np.degrees(theta.flatten())*60.

			for ii in range(len(R)):
				if abs( theta_arcmins_flatten[ii] )< theta_max:
					#print ii, len(R), theta_arcmins_flatten[ii], abs( theta_arcmins_flatten[ii] )< theta_max, "True"
					Sigma_int[ii] =  integrate.quad(lambda r : (norm_for_int * np.exp ( - (2/alpha) * ( (r/r_minus2)**alpha - 1)  ) ) * ( r / np.sqrt(r**2. - R[ii]**2.) ), R[ii], np.inf, epsabs=self.epsabs, epsrel=self.epsrel)[0]

				#else:
				#	print ii, len(R), theta_arcmins_flatten[ii], abs( theta_arcmins_flatten[ii] )< theta_max, "False"
		else:
			### print '\nEntering normal execution\n'
			for ii in range(len(R)):
				Sigma_int[ii] =  integrate.quad(lambda r : (norm_for_int * np.exp ( - (2/alpha) * ( (r/r_minus2)**alpha - 1)  ) ) * ( r / np.sqrt(r**2. - R[ii]**2.) ), R[ii], np.inf, epsabs=self.epsabs, epsrel=self.epsrel)[0]
		### end = time.time()
		### print end, 'time taken = %s' %(end-start), '\n\n'#;sys.exit()

	return Sigma_int.reshape(theta.shape), mean_Sigma_num_int.reshape(theta.shape)

	"""
	elif profile == 'DandK': #Diemer & Kravtsov


	    #NFW profile properties
	    r_minus2 = r_s.value
	    alpha = 0.18 #gives something similar to NFW
	    #alpha = .3 #galaxy cluster #https://arxiv.org/pdf/1302.0288v2.pdf

	    nu = ((alpha - 0.155)/0.0095)**0.5 #Eq. 5 of arXiv:1401.1216
	    r_t = (1.9 - 0.18 * nu) * r_s.value #Eq. 6 of arXiv:1401.1216
	    beta, gamma = 4., 8.
	    b_e, s_e = 1., 1.5 #Fig. 18 for nu ~ 4.0

	    t1 = (mass_def/3.) * r_v.value**3.
	    r_t = np.inf
	    rho_m = 0.
	    t2 = 1. / integrate.quad(lambda r : r**2 * ( np.exp ( - (2/alpha) * ( (r/r_minus2)**alpha - 1)  ) * (1 + (r/r_t)**beta)**-(gamma/beta) + ( rho_m * (b_e * (r/5/r_s.value)**-s_e + 1 ) ) ), 0., r_v.value, epsabs=self.epsabs, epsrel=self.epsrel)[0]

	    delta_c = t1 * t2
	    norm_for_int = 2. * delta_c * rho_c_z.value

	    for ii in range(len(R)):
	        Sigma_int[ii] =  integrate.quad(lambda r : (norm_for_int * np.exp ( - (2/alpha) * ( (r/r_minus2)**alpha - 1)  ) * (1 + (r/r_t)**beta)**-(gamma/beta) + ( rho_m * (b_e * (r/5/r_s.value)**-s_e + 1 ) )  ) * ( r / np.sqrt(r**2. - R[ii]**2.) ), R[ii], np.inf, epsabs=self.epsabs, epsrel=self.epsrel)[0]

	return Sigma_int.reshape(theta.shape), mean_Sigma_int.reshape(theta.shape)
	"""

    def fn_get_kappa_symmetric(self, sigma, sigma_norm):

	return sigma / sigma_norm

    def fn_get_kappa_general(self,cosmo, theta, mass, concentration, z_L, z_S, mass_def, rho_def, theta_max = -1, profile = 'NFW', just_surface_mass_density = 0):

	mass = mass*u.Msun
	if (rho_def == 'crit'):
	    rho_c_z = cosmo.critical_density(z_L)
	elif (rho_def == 'mean'):
	    rho_c_z = cosmo.Om(z_L)*cosmo.critical_density(z_L)
	else:
	    print "rho definition not specified correctly in cluster profile"
	    assert(0)

        rho_c_z = rho_c_z.to('M_sun/Mpc3') #convert critical density into M_sun/MPc3

	self.nfw_dic = {}
	self.mass = mass
	self.rho_c_z = rho_c_z
	self.concentration = concentration
	self.mass_def = mass_def
	self.rho_m = cosmo.Om(z_L)*cosmo.critical_density(z_L)

	#Angular diameter distances
	D_L = cosmo.comoving_distance(z_L)/(1.+z_L)
	D_S = cosmo.comoving_distance(z_S)/(1.+z_S)
	D_LS = (cosmo.comoving_distance(z_S)-cosmo.comoving_distance(z_L))/(1.+z_S)


	self.nfw_dic['D_L'] = D_L

	#surface mass density of the cluster
	Sigma_num_int, mean_Sigma_num_int = self.fn_get_surface_mass_density(theta, D_L, profile = profile, theta_max = theta_max)
	if just_surface_mass_density: return Sigma_num_int

	#Normalization of kappa
	Sigma_c = (((const.c.cgs**2.)/(4.*np.pi*const.G.cgs))*(D_S/(D_L*D_LS))).to('M_sun/Mpc2').value

	kappa_num_int = self.fn_get_kappa_symmetric(Sigma_num_int, Sigma_c)

	#def_angle = theta * mean_Sigma_num_int / Sigma_c
	#from pylab import *
	#imshow(kappa_num_int);colorbar();show();quit()

	return kappa_num_int, self.nfw_dic

################################################################################################################################################
################################################################################################################################################
################################################################################################################################################
