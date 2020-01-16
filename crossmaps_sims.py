#!/usr/bin/env python
########################################################################################################
########################################################################################################
########################################################################################################
import scipy.optimize as optimize

def fn_powerlaw_fitting_func(p, x, y = None, return_fit = 0):

	def fn_law(p, xp):
		return p[0]*(xp**p[1])

	if not return_fit:
		return np.ravel(fn_law(p, x) - y)
	else:
		return fn_law(p, x)

def fn_get_extra_variance_dic(richness_arr_for_sz_variance, use_lgmca = 0):
##def fn_get_extra_variance_dic(weights_folder = 'data/sims/sptpol_des/crossmaps/sptpol_map_like_sims/DES_RM_like_sample/year3_lgt_5_with_FG/CR_maps/sehgal_sims/crossmaps_nomask_noFG_extraszvariance_weights_for_year3_lgt5_vl/SNR_20.0_100.0/0.5_tSZ_sptpol/v1'):

	print '\n\tgetting SZ extra variance now from sims'
	if (1):#not use_lgmca:
		weights_folder = 'data/sims/sptpol_des/crossmaps/sptpol_map_like_sims/DES_RM_like_sample/year3_lgt_5_with_FG/CR_maps/sehgal_sims/crossmaps_nomask_noFG_extraszvariance_weights_for_year3_lgt5_vl/SNR_20.0_100.0/0.5_tSZ_sptpol/'
		#weights_folder = 'data/sims/sptpol_des/crossmaps/sptpol_map_like_sims/DES_RM_like_sample/year3_lgt_5_with_FG/CR_maps/sehgal_sims/crossmaps_nomask_noFG_extraszvariance_weights_for_year3_lgt5_vl/SNR_20.0_1000.0/0.5_tSZ_white_0.1/'
	else:
		weights_folder = 'data/sims/sptpol_des/crossmaps/sptpol_map_like_sims/DES_RM_like_sample/year3_lgt_5_with_FG/CR_maps/sehgal_sims/crossmaps_nomask_noFG_extraszvariance_weights_for_year3_lgt5_vl/lgmca/SNR_20.0_100.0/0.5_tSZ_sptpol/'
	weights_files = sorted( glob.glob('%s/stacked_kappa*' %(weights_folder)) )##[0:5]

	extra_variance_weights_dic = {}
	for wfcnt, wfname in enumerate( weights_files ):
		###print wfcnt
		dic = pickle.load(gzip.open(wfname, 'rb'))['SZ_signal_extravar']
		for kkk in dic.keys():
			if kkk in extra_variance_weights_dic:
				extra_variance_weights_dic[kkk].append(dic[kkk])
			else:
				extra_variance_weights_dic[kkk] = [ dic[kkk] ]

	extra_variance_weights_keys = np.asarray(extra_variance_weights_dic.keys())
	extra_variance_weights_noise_values = np.asarray(extra_variance_weights_dic.values())
	extra_variance_weights_noise_mean = np.mean(extra_variance_weights_noise_values, axis = 1)

	##from IPython import embed; embed()

	richarr = extra_variance_weights_keys[:,3]
	x, y = richarr, extra_variance_weights_noise_mean
	p0 = [y[0]/x[0]**2., 2.]
	p1, pcov, infodict, errmsg, success = optimize.leastsq(fn_powerlaw_fitting_func, p0[:], args=(x, y), full_output=1)
	extra_variance_weights_noise_fit = fn_powerlaw_fitting_func(p1, richness_arr_for_sz_variance, return_fit = 1)

	if (0):
		plot(x, y, color = 'orange', marker ='.', ls ='None', label = r'\textbf{$\sigma_{SZ}$}')
		plot(richness_arr_for_sz_variance, extra_variance_weights_noise_fit, 'k.', ls = 'None')
		legend(loc = 2);show()

	return extra_variance_weights_noise_fit, p1

	'''
	#set default weight using max. richness values
	#this is because Sehgal sims did not have clusters more than lambda>~100
	defaultnoiseind = np.argmax(extra_variance_weights_keys[:,3])
	defaultnoisevalue = extra_variance_weights_noise_mean[defaultnoiseind]

	return extra_variance_weights_keys, extra_variance_weights_noise_mean, defaultnoisevalue
	'''

def fn_get_weights(RA, DEC, CLUS_IDENTIFIER_full, kappa_qe_arr_full, kappa_weights = 1, extra_variance_weights = 1, clus_radius_am = 10., testing = 0, use_lgmca = 0):

	print '\n\tgetting weights now'

	RADIUS = np.sqrt( RA**2. + DEC**2. ) * 60.
	inds_outer = np.where( (RADIUS>=clus_radius_am) & (RADIUS<=30.) )

	if extra_variance_weights:
		##extra_variance_weights_keys, extra_variance_weights_noise_mean, defaultnoisevalue = fn_get_extra_variance_dic()
		richness_arr_for_sz_variance = CLUS_IDENTIFIER_full[:,3]
		extra_variance_weights_noise_fit, extra_var_params_plaw_fit = fn_get_extra_variance_dic(richness_arr_for_sz_variance, use_lgmca = use_lgmca)
		##from IPython import embed; embed()

	WEIGHTS = []
	WEIGHTS_IND = []
	WEIGHTS_DIC = {}
	kappa_qe_arr_full_corrected = []
	for cntr in range( len(CLUS_IDENTIFIER_full) ):

		keyname = CLUS_IDENTIFIER_full[cntr] #get cluster info.
		kappa_qe = kappa_qe_arr_full[cntr] #get curr kappa qe for weight estimation

		std1, std2 = 0., 0. #default weights			
		#convergence weight ignoring cluster region
		if kappa_weights:
			std1 = np.std(kappa_qe[inds_outer])

		if extra_variance_weights:
			'''
			r,d,z,l = keyname[0:4]
			wind = np.where( (extra_variance_weights_keys[:,0] == r) & (extra_variance_weights_keys[:,1] == d) \
					& (extra_variance_weights_keys[:,2] == z) & (extra_variance_weights_keys[:,3] == l))[0]

			if len(wind) == 0: 
				std2 = defaultnoisevalue #clusters for which we don't have a value from sims (roughly 5 clusters in vl sample)
			else:
				std2 = extra_variance_weights_noise_mean[wind[0]]
			'''
			std2 = extra_variance_weights_noise_fit[cntr] #using a power law fit from the vl sample + sehgal SZ sims

		###print cntr, std1, std2

		#total error
		totstd = np.sqrt( std1**2 + std2**2.)

		#calc. weights now
		if totstd == 0.: totweight = 1.
		else: totweight = 1./totstd**2.

		if testing: 
			richval = keyname[3]

			##subplot(111);
			if cntr == 0:
				plot(richval, std1, '.', color = 'k', label = '$\sigma_{\kappa}$')
				plot(richval, std2, '.', color = cm.gnuplot(200), label = '$\sigma_{SZ}$');
			else:
				plot(richval, std1, '.', color = 'k');plot(richval, std2, '.', color = cm.gnuplot(200));
			##subplot(122);plot(richval, totweight, '.', color = 'k')

		WEIGHTS.append(totweight)
		if extra_variance_weights:
			WEIGHTS_IND.append( [1/std1**2., 1/std2**2.] )

		WEIGHTS_DIC[tuple(keyname)] = totweight

	###from IPython import embed; embed()
	WEIGHTS = np.asarray( WEIGHTS )
	if extra_variance_weights:
		WEIGHTS_IND = np.asarray( WEIGHTS_IND )

	if testing: 
		legend(loc =  1, fontsize = 10)
		xlabel(r'\textbf{Richness}', fontsize = 12); ylabel(r'\textbf{Error}', fontsize = 12); title(r'\texbf{Noise in kappa maps}')
		show();sys.exit()

	if not extra_variance_weights:
		return WEIGHTS, WEIGHTS_DIC##, WEIGHTS_IND, extra_var_params_plaw_fit
	else:
		return WEIGHTS, WEIGHTS_DIC, WEIGHTS_IND, extra_var_params_plaw_fit
################################################################################################################
################################################################################################################
################################################################################################################
def gaussian_fitting_func(p, p0, X, Y, MAP, lbounds = None, ubounds = None, fixed = None, return_fit = 0):
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

def fitting_func(p, p0, X, cluster_mass, DATA = None, lbounds = None, ubounds = None, fixed = None, return_fit = 0, commit_20161221_change = 0):
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

def fn_get_QE_sims_for_COV_MF(noofsims, totalclus, mapparams, param_dict, randomseedval = None, pos_dep_MASK = None, return_all_kappa_sims = 0, store_sims_folder = None):

	nx, ny, dx, dy = mapparams
	STACKED_KAPPA_QE_SIMS = np.zeros( (noofsims, nx, ny) )
	for simno in range(noofsims):

		## print randomseedval
		timeval = '%.6f' %(time.time())
		randomseedval = int(timeval.split('.')[-1])

		#get lensed CMB for specific cluster M, c, z + add beam; imshow(CMB_SIMS[0,0]);colorbar();show()
		CMB_SIMS = sims.fn_perform_lensing_nfw(mapparams, totalclus, param_dict, param_dict['use_beam'], degrade_resol_fac = param_dict['degrade_resol_fac'], random_seed = randomseedval, no_lensing = 1)
		#from IPython import embed; embed()

		'''
		#from IPython import embed; embed()
		#print CMB_SIMS
		#20171114 - FG additions
		nu = 150.
		FG_uncorr = np.asarray( [sims.fn_get_foreground_power('all_Gaussian', mapparams, 1, perform_lensing = 0, nu1=nu, nu2=nu)[0] for nn in range(totalclus)] )
		FG_uncorr *= 1e6		
		CMB_SIMS += FG_uncorr
		'''

		### print '\n\n90 to be performed now\n\n'

		if param_dict['cross_maps']:
			### GRAD_SIMMAPS = sims.fn_perform_lensing_nfw(mapparams, totalclus, param_dict, param_dict['use_beam'], degrade_resol_fac = param_dict['degrade_resol_fac'], random_seed = randomseedval, no_lensing = 1, Bl = Bl_gradient_map)
			GRAD_SIMMAPS = sims.fn_perform_lensing_nfw(mapparams, totalclus, param_dict, param_dict['use_beam'], degrade_resol_fac = param_dict['degrade_resol_fac'], random_seed = randomseedval, no_lensing = 1, Bl = Bl_gradient_map, nu=90)
			'''
			nu = 90.
			FG_uncorr_90 = np.asarray( [sims.fn_get_foreground_power('all_Gaussian', mapparams, 1, perform_lensing = 0, nu1=nu, nu2=nu)[0] for nn in range(totalclus)] )
			FG_uncorr_90 *= 1e6		
			GRAD_SIMMAPS += FG_uncorr_90
			'''
		#add noise, TF
		#CMB_SIMS = sims.fn_tf_noise(CMB_SIMS, mapparams, noofsims, reqdbox = None, add_TF = param_dict['add_TF'], add_noise = add_noise, expnoiselevel = expnoiselevel)
		CMB_SIMS = np.asarray( map(lambda x, y: sims.fn_tf_noise(x, mapparams, 1, reqdbox = None, add_TF = param_dict['add_TF'], add_noise = add_noise, expnoiselevel = expnoiselevel, cluscnt = y, noise_random_seed = 0), CMB_SIMS, np.arange(totalclus)) )
		CMB_SIMS = CMB_SIMS[:,0]/1e6

		if param_dict['cross_maps']:
			if expnoiselevel<>None:
				expnoiselevel_grad = np.tile(param_dict['noise_level_grad_map'], len(expnoiselevel))
			else:
				expnoiselevel_grad = None

			#GRAD_SIMMAPS = sims.fn_tf_noise(GRAD_SIMMAPS, mapparams, noofsims, reqdbox = None, add_TF = param_dict['add_TF'], add_noise = add_noise, expnoiselevel = expnoiselevel_grad, special_map = 1)
			GRAD_SIMMAPS = np.asarray( map(lambda x, y: sims.fn_tf_noise(x, mapparams, 1, reqdbox = None, add_TF = param_dict['add_TF'], add_noise = add_noise, expnoiselevel = expnoiselevel_grad, special_map = 1, cluscnt = y, nu = 90., noise_random_seed = 0), GRAD_SIMMAPS, np.arange(totalclus)) )
			GRAD_SIMMAPS = GRAD_SIMMAPS[:,0]/1e6

		'''
		subplot(121);imshow(CMB_SIMS[0,0]);colorbar()
		if 	param_dict['cross_maps']:
			subplot(122);imshow(GRAD_SIMMAPS[0,0]);colorbar();
		show();sys.exit()
		'''

		logline = '\t\t cluster sims complete for sim = %s. next kappa estimation' %(simno)
		logfile = open(log_file,'a');logfile.writelines('%s\n' %(logline));logfile.close()
		print logline

		if param_dict['use_mask']:
			try:
				MASK = apodMASK
			except:
				apodMASK = fn_apodize(CMB_SIMS[0], mapparams, mask = 'circle', just_return_mask = 1)
				MASK = apodMASK

			#logline = '\t\t apodmask obtained'
			#if position dependent mask is supplied the multiply it to apdomask here
			if pos_dep_MASK <> None: #do not use this for stacked sample
				#from IPython import embed; embed()
				MASK = np.copy(MASK) * pos_dep_MASK

			CMB_SIMS_MASKED = np.copy(CMB_SIMS) * MASK
			if param_dict['cross_maps']:
				GRAD_SIMMAPS_MASKED = np.copy(GRAD_SIMMAPS) * MASK
		else:
			CMB_SIMS_MASKED = np.copy(CMB_SIMS)
			if param_dict['cross_maps']:
				GRAD_SIMMAPS_MASKED = np.copy(GRAD_SIMMAPS)

		#GRAD_SIMMAPS_MASKED, CMB_SIMS_MASKED = np.copy(GRAD_SIMMAPS), np.copy(CMB_SIMS)
		#20171207 - remove mean of the cutouts here before kappa_qe
		CMB_SIMS_MASKED = CMB_SIMS_MASKED - np.mean(CMB_SIMS_MASKED)
		if param_dict['cross_maps']:
			GRAD_SIMMAPS_MASKED = GRAD_SIMMAPS_MASKED - np.mean(GRAD_SIMMAPS_MASKED)

		"""
		print CMB_SIMS_MASKED.shape
		subplot(121);imshow(CMB_SIMS_MASKED[0, 0]);colorbar()
		subplot(122);imshow(GRAD_SIMMAPS_MASKED[0, 0]);colorbar();show();
		sys.exit()
		"""

		if param_dict['cross_maps']: 
			#KAPPA_QE = sims.fn_get_kappa_QE(OBSMAP, simmapparams, Dls_len, Dls_unlen, tszfree = 0, OBSMAP2 = OBSMAP2)
			KAPPA_QE_SIMS = np.asarray( map(lambda x, y: sims.fn_get_kappa_QE([x], mapparams, Dls_len, Dls_unlen, tszfree = 0, OBSMAP2 = [y]).real[0], CMB_SIMS_MASKED, GRAD_SIMMAPS_MASKED) )
		else:
			KAPPA_QE_SIMS = np.asarray( map(lambda x: sims.fn_get_kappa_QE([x], mapparams, Dls_len, Dls_unlen, tszfree = 0).real[0], CMB_SIMS_MASKED) )

		KAPPA_QE_SIMS = np.asarray( [KQES - np.mean(KQES) for KQES in KAPPA_QE_SIMS] )

		#from IPython import embed; embed()

		#logline = '\t\t KAPPA_QE_SIMS otained'
		#logfile = open(log_file,'a');logfile.writelines('%s\n' %(logline));logfile.close()
		#print logline

		if 0 == 1:#mf_sub:
			#20171123 - remove mean feild from sims
			MEAN_FIELD = np.mean(KAPPA_QE_SIMS, axis = 0)
			### kappa_qe = np.copy(kappa_qe) - MEAN_FIELD[0] #no mean-field sub here
			KAPPA_QE_SIMS = KAPPA_QE_SIMS - MEAN_FIELD[0]
		else:
			KAPPA_QE_SIMS = np.asarray( [KQES - np.mean(KQES) for KQES in KAPPA_QE_SIMS] )
			### kappa_qe = kappa_qe - np.mean(kappa_qe) #20171129 - removing mean from estimated kappa

		STACKED_KAPPA_MAP = np.mean(KAPPA_QE_SIMS, axis = 0)

		#get mean-field also from sims here
		try:
			MEAN_FIELD_FROM_STACKED_KAPPA_QE_SIMS += np.sum(KAPPA_QE_SIMS, axis = 0)
			totstacked_for_mean_field += len(KAPPA_QE_SIMS)
		except:
			MEAN_FIELD_FROM_STACKED_KAPPA_QE_SIMS = np.sum(KAPPA_QE_SIMS, axis = 0)
			totstacked_for_mean_field = len(KAPPA_QE_SIMS)

		if store_sims_folder<>None: #then store each sim
			sims_store_file_for_cov = '%s/SIM_%sclusters_%s.pkl.gz' %(store_sims_folder, totalclus, randomseedval)
			pickle.dump(STACKED_KAPPA_MAP, gzip.open(sims_store_file_for_cov, 'wb'), protocol = 2)

		STACKED_KAPPA_QE_SIMS[simno] = STACKED_KAPPA_MAP

	MEAN_FIELD_FROM_STACKED_KAPPA_QE_SIMS /= totstacked_for_mean_field

	return STACKED_KAPPA_QE_SIMS, MEAN_FIELD_FROM_STACKED_KAPPA_QE_SIMS

def fitting_func(p, p0, X, cluster_mass, DATA = None, lbounds = None, ubounds = None, fixed = None, return_fit = 0, commit_20161221_change = 0):
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


def fn_calc_delta_chisq_one_width(Marr, del_chi_sq_arr, cluster_mass, perform_ip = 0, fit_parabola = 1):

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
			p1, success = optimize.leastsq(fitting_func, p0, args=(p0, Marr,  cluster_mass, del_chi_sq_arr, lbounds, None, fixed))#, ubounds))
		else:
			p1, success = optimize.leastsq(fitting_func, p0, args=(p0, Marr,  cluster_mass, del_chi_sq_arr, lbounds))#, ubounds))

		#print p0
		#print p1, success;quit()

		'''
		if success == 1:
			width = p1[0]# * 2.
		else:
			width = 0.
#		if p1[0] > 2 * p0:
#			width = 0.
		'''

		#print lbounds, ubounds, p1
		#x_ip = np.arange(min(Marr),max(Marr),0.01)
		x_ip = np.arange(min(Marr),max(Marr),cluster_mass/100.)
		#p1 = p1/2.
		y_ip = fitting_func(p1,p1,x_ip, cluster_mass,return_fit = 1)

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

		if 0==1:
			clf()
			y_ip_ori = fitting_func(p0,p0,x_ip, cluster_mass,return_fit = 1)

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


def fn_get_det_significance(Marr, res_arr, massval_for_L_peak = None):
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

def fn_apodize(MAP, mapparams, mask = 'square', just_return_mask = 1):

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

	return apodMASKf

def fn_calc_likelihood(MAP, MODEL, C, modelmass, onedfit = 0, kappa_two_halo_term = None, plot_onedfit = 0):

	#imshow(C);colorbar();show();quit()
	Cinv = sc.linalg.pinv2(C)
	C = np.mat(C)

	sign, logdetval = np.linalg.slogdet(C)
	logdetval = logdetval * sign

	#t1 = -.5 * npixels * np.log( 2 * np.pi )
	#t2 = -0.5 * logdetval

	if onedfit:
		#raprfmap = sims.fn_radial_profile(MAP, RADEC)
		raprfmap = sims.fn_radial_profile(MAP, RADEC, bin_size=binsize, minbin=0.0, maxbin=maxbin)
		d = raprfmap[:,1].flatten()###  - np.mean(raprfmap[:,1].flatten())

		#raprfmodel = sims.fn_radial_profile(MODEL, RADEC)
		raprfmodel = sims.fn_radial_profile(MODEL, RADEC, bin_size=binsize, minbin=0.0, maxbin=maxbin)
		dp = raprfmodel[:,1].flatten()### - np.mean(raprfmodel[:,1].flatten())

	else:

		d= MAP.flatten()
		dp = MODEL.flatten()

	if 0==1:#sims.is_seq(kappa_two_halo_term):

		dp += kappa_two_halo_term

	#subplot(221);imshow(MAP);colorbar();title(np.mean(MAP))
	#subplot(222);imshow(MODEL);colorbar();title(np.mean(MODEL))
	#subplot(223)
	if plot_onedfit and local:
		##plot(raprfmap[:,0], d, 'ro')
		if modelmass == Aarr[-1]:
			errors = np.diag(C)**0.5
			errorbar(raprfmap[:,0], d, yerr = errors, marker = 'o', color='k', capsize = 2., elinewidth = 1., label=r'\textbf{Data}')
		try:
			Marr
		except:
			Marr = np.copy(Aarr)
		colorinds = np.linspace(0,255,len(Marr)); colorarr = np.asarray( [cm.jet( int(c) ) for c in colorinds] )
		ccurrcolind = np.where(Marr == modelmass)[0][0];colorval = colorarr[ccurrcolind]
		plot(raprfmap[:,0], dp, color = colorval,alpha=0.5)#;title('SIM = %s; Model mass = %s' %(simcnt, modelmass))
		#show()

	d = d-dp

	#Hartlap correction factor
	correction_factor = (float(noofsims) - len(d) - 1.)/noofsims
	#print correction_factor, float(noofsims), len(d)
	Cinv = Cinv * correction_factor

	t3 = -0.5 * np.asarray( np.dot(d.T, np.dot( Cinv, d ))).squeeze()

	logLval = t3# + t2

	return logLval
	
def fn_get_CLUS_IDENIFIER(totalclus, preserve_radec = 0, minrich = 20, maxrich = 100):
	#deprecated!
	return None
	mapfolder = 'data/sanjay_maps_201705xx/final_maps/20171027_downsample_x2_tSZfree_plus_150_v2_no_masking'
	data_file = '%s/map_150.pkl.gz_cutouts_150ghz_no_downsampling_y3_v6.4.21_lgt5_redmapper_clusters_full_TQUforMLE' %(mapfolder)
	CLUSKEYNAMES = pickle.load(gzip.open(data_file, 'rb'))['cutouts'].keys()
	if minrich == None:
		selectedinds = np.arange(totalclus) #sorted( np.random.choice(np.arange(len(CLUSKEYNAMES)), size = totalclus, replace = 0) )
	else:
		richarr_RM = np.asarray( map(lambda x: x[3], CLUSKEYNAMES) )
		selectedinds = np.where( (richarr_RM>= minrich) & (richarr_RM <= maxrich))[0][0:totalclus]

	CLUSKEYNAMES = map(lambda x: CLUSKEYNAMES[x], selectedinds)
	CLUS_IDENTIFIER = []

	for k in range(len(CLUSKEYNAMES)):
		currval = CLUSKEYNAMES[k]
		clra, cldec, z_val, richval = currval[0], currval[1], currval[2], currval[3]
		if not preserve_radec:
			clra, cldec = 0., 0.
		CLUS_IDENTIFIER.append([clra, cldec, z_val, richval, 1.])

	CLUS_IDENTIFIER = np.asarray( CLUS_IDENTIFIER )

	return CLUS_IDENTIFIER

def fn_select_clusters(kappa_qe_arr, CLUS_IDENTIFIER, rich1, rich2, z1 = None, z2 = None):

	passed_inds = []
	for kcnt, cluskey in enumerate( CLUS_IDENTIFIER ):
		ra, dec, z_val, rich, weight = cluskey

		passed = 0
		if rich >= rich1 and rich<rich2:
			passed = 1

		if z1<>None:
			if z_val>=z1: 
				passed = 1
			else: 
				passed = 0

		if z2>None:
			if z_val<z1: 
				passed = 1
			else: 
				passed = 0

		if passed: passed_inds.append(kcnt)

	if len(passed_inds) == 0: return None, None

	passed_inds = np.asarray( passed_inds )

	return kappa_qe_arr[passed_inds], CLUS_IDENTIFIER[passed_inds]


def fn_add_tf_beam_to_model(kappa_model):

	lx, ly = sims.get_lxly(mapparams)
	L = (lx**2. + ly**2.)**0.5
	ngaus = int( np.sqrt(8. * np.log(2.)) / np.radians(sims.exp_beam/60.) )
	above_beam_scale = np.where(L>=ngaus)

	#20180407 - removing the next 3 steps: beam and TF convolution and removing modes similar to data in QE. 
	#20180407 - Adding beam after this step
	'''
	kappa_model_fft = np.fft.fft2(kappa_model)
	kappa_model_fft[above_beam_scale] = 0.
	kappa_model_old = np.fft.ifft2( kappa_model_fft * TWODTF ).real #twodtf similar to data here
	'''
	
	beam_tf_decov = 1./np.copy(sims.Bl * TWODTF)
	bad = np.where( (beam_tf_decov==np.inf) | (beam_tf_decov>=1e10) | (beam_tf_decov<0.) )
	kappa_model_fft = np.fft.fft2(kappa_model)
	kappa_model_fft[above_beam_scale] = 0.
	kappa_model_fft[bad] = 0.
	kappa_model = np.fft.ifft2( kappa_model_fft ).real

	return kappa_model

########################################################################################################
########################################################################################################
########################################################################################################
def fn_get_kappa_model_A_dic_from_stored_models(AA, alpha_fit_arr, minrich, maxrich, inv_var_weights = 0):
	#20180411 - read all kappadics and store. pick the model for specified richness range from the stored dics

	### from IPython import embed; embed()
	searchstr = '%s/stacked_kappa_model_forA%s*.pkl.gz' %(opfolder, AA)
	kappa_model_A_file_arr = glob.glob(searchstr)
	if len(kappa_model_A_file_arr) == 0: print '\n\nNo model files found. aborting here\n\n'

	tmprich_arr_for_models = []
	pick_this_file_arr = []
	for mmm in kappa_model_A_file_arr:
		tmp = mmm.split('_')[-2:]
		r1, r2 = float( tmp[0].replace('minrich','') ), float( tmp[1].replace('maxrich','').replace('.pkl.gz','') )
		if r1 >= minrich and r2 <= maxrich:
			pick_this_file_arr.append(mmm)
		tmprich_arr_for_models.append(r1)
		tmprich_arr_for_models.append(r2)
	tmprich_arr_for_models = np.unique(tmprich_arr_for_models)
	assert minrich in tmprich_arr_for_models and maxrich in tmprich_arr_for_models


	#which files to combine to produce the model is stored in "pick_this_file_arr"
	tmp_dic = {}
	for tmpfname in pick_this_file_arr:
		dic = pickle.load(gzip.open(tmpfname, 'rb'))
		### print tmpfname
		for alpha_fit in alpha_fit_arr:
			keyname = round(alpha_fit, 3)
			'''
			howmanyclus = dic[keyname][2]
			if keyname in tmp_dic:
				tmp_dic[keyname][0] += ( modelsum * howmanyclus) # * howmanyclus --> becuase I have already calculated the mean instead of just a sum
				tmp_dic[keyname][1] += howmanyclus
			else:

				tmp_dic[keyname] = [dic[keyname][0] * howmanyclus, dic[keyname][2]]
			'''
			### from IPython import embed; embed()
			if not inv_var_weights:
				howmanyclus = dic[keyname][1]
				modelsum = dic[keyname][0]
			else:
				howmanyclus = dic[keyname][3]
				modelsum = dic[keyname][2]

			if keyname in tmp_dic:
				tmp_dic[keyname][0] += modelsum # * howmanyclus --> becuase I have already calculated the mean instead of just a sum
				tmp_dic[keyname][1] += howmanyclus
			else:

				tmp_dic[keyname] = [modelsum, howmanyclus]

	kappa_model_A_dic = {}
	for keyname in tmp_dic:
		model_for_this_alpha = tmp_dic[keyname][0]/tmp_dic[keyname][1]
		kappa_model_A_dic[keyname] = model_for_this_alpha

	return kappa_model_A_dic

########################################################################################################
########################################################################################################
########################################################################################################

################################################################################################
################################################################################################
################################################################################################

import numpy as np, pickle, gzip, sys, glob, time, os, scipy as sc, pdb
sys.path.append('modules')
import modules
import sky_local
from scipy.interpolate import interp1d
#import modules.radialProfile as radialProfile
#import healpy as H

local = 1
if str(os.getcwd()).find('sri')>-1:
	local = 0
if local:
	#import matplotlib;matplotlib.use('Agg')
	from pylab import *

sims = modules.scl_cmb.simulations()
args = sys.argv[1:]

onedfit = 1

#ipfolder, clustype, noofsims, paramfile, totalclus, mf_sub, start, end = args
#start, end = int(start), int(end)

ipfolder, noofsims, minrich, maxrich, mf_sub = args
noofsims = int(noofsims)
mf_sub = int(mf_sub)
minrich = float(minrich)
maxrich = float(maxrich)

timedatestr='%s' %(time.strftime('%Y%m%d_%H%M%S', time.localtime()))
log_file = 'tmp/logs_%s_fit_kappa.txt' %(timedatestr)
sims.log_file = log_file

sims.quiet = 1

minbin = 0.
binsize, maxbin = 1.0, 10.
sims.CR_maps = 1

########################################################################################
########################################################################################
#what are we fitting
#1. Just masss
if 0==1:
	minM, maxM, delM = 0., 6.0, 0.2
	#minM, maxM, delM = 2., 6., 0.4
	Marr = np.arange(minM,maxM+delM,delM) * 1e14

#2. M-lambda relation
#model: M = A * (lambda / lambda_pivot)**alpha * ( (1+z) / (1+z_pivot) )**beta
#now we are just fitting A and fixing alpha = 1.225 (1.12 M17 and 1.33 S16); beta = 0.18 from M17 (see page 7 of 1708.01360)

alpha_fit = 1.12 #from M17 ###1.225
lambda_pivot = 30.
beta_fit = 0.18
z_pivot = 0.5

minA, maxA, delA = 0.0, 10.0, 0.5
#minA, maxA, delA = 0.2, 6., 0.4
Aarr = np.arange(minA,maxA+delA,delA) * 1e14
########################################################################################
########################################################################################


if ipfolder.find(',')>-1:
	ipfolder, extra = ipfolder.split(',')
else:
	extra = ''
searchstr = '%s/*stacked*%s*.pkl.gz*' %(ipfolder,extra)
### fnamearr = glob.glob(searchstr)#[0:2]
fnamearr = sorted(glob.glob(searchstr), key=os.path.getsize)
### fnamearr = fnamearr[0:5]

if len(fnamearr) == 0:
	print 'No files found like %s' %(searchstr)
	sys.exit()

if ipfolder.find('null_test')>-1:
	null_test = 1
else:
	null_test = 0

################################################################################################
################################################################################################
#20180325 - first pick the file with kappa_qe_arr and get the kappa_COV
#from IPython import embed;embed()
#sys.exit()
if (1):
	cov_file = '%s/kappa_COV_%sJK_1d_%sam_%sdx_richnessbins_%s_%s.pkl.gz' %(ipfolder, noofsims, maxbin, binsize, minrich, maxrich)
	kappa_qe_arr = None
	if not os.path.exists(cov_file):
		print '\n\n getting kappa_COV using a JK approch\n'
		for fname in fnamearr:
			### print fname
			kappadic = pickle.load(gzip.open(fname,'rb'))
			if not 'kappa_qe_arr' in kappadic:
				continue
			kappa_qe_arr = np.asarray( kappadic['kappa_qe_arr'] )
			CLUS_IDENTIFIER = np.asarray( kappadic['CLUS_IDENTIFIER'] )
			break
		if sims.is_seq(kappa_qe_arr):
			totalclus = len(kappa_qe_arr)
			### CLUS_IDENTIFIER = fn_get_CLUS_IDENIFIER(totalclus, preserve_radec = 1)
			kappa_qe_arr, CLUS_IDENTIFIER = fn_select_clusters(kappa_qe_arr, CLUS_IDENTIFIER, minrich, maxrich)
			kappa_COV = sims.fn_get_COV_from_JK(kappa_qe_arr, CLUS_IDENTIFIER, noofsims, [ (minrich, maxrich) ] )
			pickle.dump(kappa_COV, gzip.open(cov_file,'wb'), protocol = 2)
		else:
			print 'No kappa_qe_arr for JK cov'
			sys.exit()
		sys.exit()

################################################################################################
################################################################################################

BEST_FIT_MODELS = {}
if local: colorarr = [cm.jet(int(c)) for c in np.linspace(10, 250, len(fnamearr))]
#run s2_2_fit_kappa_crossmaps.py data/sims/sptpol_des/crossmaps/sptpol_map_like_sims/data/rich_20_26/0.5_no_tSZ_sptpol/ sptpol_des 10 data/sims/sptpol_des/crossmaps/sptpol_map_like_sims/data/rich_20_26/0.5_no_tSZ_sptpol/params_planck_r_0.0_2015_cosmo_lensed_LSS_crossmap_sptpol_used.txt 1
det_arr = np.zeros(len(fnamearr))

for fcnt,fname in enumerate(fnamearr):

	if fname.find('for_cov_file')>-1 or fname.find('randoms')>-1: continue

	#if fname.find('81840l1')==-1:continue

	opfname = fname.replace('data/sims','pkls').replace('.pkl.gz','_results.pkl.gz')
	if not mf_sub:
		tmp = opfname.split('/')
		tmp[-1] = 'no_mf_sub/%s' %(tmp[-1])
		opfname = '/'.join(tmp)
	if not onedfit:
		tmp = opfname.split('/')
		tmp[-1] = '2dfitting/%s' %(tmp[-1])
		opfname = '/'.join(tmp)
	opfname = '%s_' %(opfname)

	logline = '%s, %s' %(fname, fcnt)
	print logline

	#### read kappa that needs to be fit
	kappadic = pickle.load(gzip.open(fname,'rb'))
	

	try:
		kappa_qe_arr_full = np.asarray( kappadic['kappa_qe_arr'] )
	except:
		kappa_qe_arr_full = None
		pass
	CLUS_IDENTIFIER_full = np.asarray( kappadic['CLUS_IDENTIFIER'] )
	z_L_list = CLUS_IDENTIFIER_full[:,2]
	totalclus = len(CLUS_IDENTIFIER_full)

	print '\n\t\ttotalclus = %s' %(totalclus)

	kappa_qe_arr, CLUS_IDENTIFIER = kappa_qe_arr_full, CLUS_IDENTIFIER_full
	CLUS_IDENTIFIER = np.asarray( CLUS_IDENTIFIER )

	param_dict = kappadic['param_dict']
	nx,ny = kappadic['stacked_kappa_qe'].shape
	mapparams = nx,ny,0.5,0.5
	for AA in Aarr:
		cc = self.c_Duffy(MM, z_val, param_dict['h'], profile_name = 'NFW')
		kappa_model = sims.fn_lensing_ini(mapparams,param_dict, RA, DEC, [clra], [cldec], [MM], [cc], [z_val], param_dict['z_lss'], param_dict['mass_def'], param_dict['rho_def'], truncate_kappa_at_radius = param_dict['truncate_kappa_at_radius'], profile_name = profile_name, theta_max = theta_max)
		lx, ly = self.get_lxly(mapparams)
		L = (lx**2. + ly**2.)**0.5
		ngaus = int( np.sqrt(8. * np.log(2.)) / np.radians(self.exp_beam/60.) )
		above_beam_scale = np.where(L>=ngaus)
	beam_tf_decov = 1./np.copy(Bl * TWODTF)
	bad = np.where( (beam_tf_decov==np.inf) | (beam_tf_decov>=1e10) | (beam_tf_decov<0.) )
	kappa_model_fft = np.fft.fft2(kappa_model)
	kappa_model_fft[above_beam_scale] = 0.


	if sims.is_seq(kappa_qe_arr): 
		kappa_qe_arr, CLUS_IDENTIFIER = fn_select_clusters(kappa_qe_arr, CLUS_IDENTIFIER, minrich, maxrich)
		if not sims.is_seq(kappa_qe_arr):
			print '\n\n\t\tNothing passed here\n\n'
			continue
		stacked_kappa_qe = np.mean( kappa_qe_arr, axis = 0)
	else:
		stacked_kappa_qe = kappadic['stacked_kappa_qe']
	raarr, decarr = CLUS_IDENTIFIER[:,0], CLUS_IDENTIFIER[:,1]

	########################################################################################################
	#remove mean field now based on random locations
	### randpatch_kappa_files = glob.glob('%s/*randoms' %(ipfolder) )
	###tmpfolder = '/'.join(ipfolder.split('/')[:7])
	###tmpfolder = '%s/SNR_%s_%s/*/' %(tmpfolder, minrich, maxrich)
	tmpfolder = '/'.join(ipfolder.split('/')[:9])

	#picking from data randoms
	if 0==1:
		tmpfolder = 'data/sims/sptpol_des/crossmaps/sptpol_map_like_sims/data/year3_lgt_5/SNR_20.0_1000.0/0.5_no_tSZ_sptpol/'
		randpatch_kappa_files = glob.glob('%s/*randoms' %(tmpfolder) )
		if len(randpatch_kappa_files)>0:
			randpatch_kappa = randpatch_kappa_files[0]
			MEAN_FIELD = pickle.load(gzip.open(randpatch_kappa, 'rb'))['stacked_kappa_qe']
		else:
			MEAN_FIELD = None


	randpatch_kappa_files = glob.glob('%s/MEANFIELD_%ssims*.pkl.gz' %(ipfolder, noofsims))
	if len(randpatch_kappa_files)>0:
		randpatch_kappa = randpatch_kappa_files[0]
		MEAN_FIELD = pickle.load(gzip.open(randpatch_kappa, 'rb'))
	else:
		MEAN_FIELD = None


	#20180408 - picking mean field from no_lensing sim runs
	splname = ipfolder.strip('/').split('/')[-3]
	if splname.find('ngrad')>-1:
		tmpfolder = '/'.join(ipfolder.strip('/').split('/')[:-3])
	splname = '%s_no_lensing/' %(splname)

	searchstr = '%s/%sSNR_20.0_40.0/0.5_no_tSZ_sptpol/*_05000_clusters*' %(tmpfolder, splname)
	### searchstr = '%s/%sSNR_20.0_40.0/0.5_no_tSZ_sptpol/*_25000_clusters*' %(tmpfolder, splname)
	randpatch_kappa_files = glob.glob( searchstr )
	param_dict = kappadic['param_dict']

	if param_dict['use_mask']:
		if len(randpatch_kappa_files) == 0: print '\n\n\tno mean-field file. aborting\n\n'; sys.exit()

		for mfcnt, mffile in enumerate( randpatch_kappa_files ):
			mfdic = pickle.load(gzip.open(mffile, 'rb')) 
			print mffile
			CURR_STACK = mfdic['stacked_kappa_qe'] * mfdic['totalclus']
			if mfcnt == 0: 
				totformf = mfdic['totalclus']
				MEAN_FIELD = CURR_STACK
			else:
				totformf += mfdic['totalclus']
				MEAN_FIELD += CURR_STACK

		MEAN_FIELD /= totformf


	### subplot(121);imshow(stacked_kappa_qe[80:120, 80:120]);colorbar()
	#if sims.is_seq(MEAN_FIELD): stacked_kappa_qe = stacked_kappa_qe - MEAN_FIELD
	### subplot(122);imshow(stacked_kappa_qe[80:120, 80:120]);colorbar();show()
	#sys.exit()
	########################################################################################################

	########################################################################################################
	### 20171214_2026: replace param_dict from the kappadic - this has all params stored
	param_dict = kappadic['param_dict']
	param_dict['rho_def'] = 'crit' #change it to mean density according to DES definitions
	#20180109 - modify other cosmo to match DES definition - Melchoir (arXiv: 1610.06890)
	#but note we must get kappa again as Cls will change
	"""
	param_dict['omega_lambda'] = 0.7
	param_dict['omega_m'] = 0.3
	param_dict['h'] = 0.7
	"""
	#20180318 - change l1 in param_dict if l1_for_qe is present
	try:
		param_dict['l1'] = param_dict['l1_for_qe']
	except:
		pass
	sims.inidic = param_dict
	########################################################################################################

	########################################################################################################
	#### Cls
	Dlfile_len = param_dict['Dlfile_len']
	#Dlfile_len = 'data/output_spt_r_0.07_lensedCls.dat'
	Dls_len = np.loadtxt(Dlfile_len,usecols=[0,1])#,2,3])

	Dlfile_unlen = param_dict['Dlfile_unlen']
	#Dlfile_unlen = 'data/output_spt_r_0.07_scalCls.dat'
	Dls_unlen = np.loadtxt(Dlfile_unlen,usecols=[0,1])#,2,3])
	########################################################################################################

	########################################################################################################
	dx = dy = kappadic['reso_arcmin']
	add_noise = kappadic['add_noise']
	try:
		expnoiselevel = [kappadic['expnoiselevel'][0]]
	except:
		expnoiselevel = None
	### boxsize = kappadic['boxsize']
	nx, ny = stacked_kappa_qe.shape	
	boxsize = nx * dx
	randomseedval = param_dict['random_seed_val_for_cov']### * int(time.time()/1e5)
	mapparams = [nx, ny, dx, dy]
	########################################################################################################

	########################################################################################################
	#### cluster stuff
	clra, cldec = 0., 0.

	minval, maxval = clra-boxsize/2/60.,  clra+boxsize/2/60.
	ra = dec = np.linspace(minval, maxval, nx)
	RA, DEC = np.meshgrid(ra,dec)
	RADEC = [RA, DEC]
	RADPRF = sims.fn_radial_profile(stacked_kappa_qe, RADEC, bin_size=binsize, minbin=0.0, maxbin=maxbin)
	### errorbar(RADPRF[:,0], RADPRF[:,1], yerr = RADPRF[:,2], marker = 'o', color = 'r');show()#;sys.exit()

	if 0==1:
		radprf = sims.fn_radial_profile(stacked_kappa_qe, RADEC, bin_size=binsize, minbin=0.0, maxbin=maxbin)
		plot(radprf[:,0], radprf[:,1]);
		if sims.is_seq(MEAN_FIELD):
			radprf_mf = sims.fn_radial_profile(MEAN_FIELD, RADEC, bin_size=binsize, minbin=0.0, maxbin=maxbin)
			plot(radprf_mf[:,0], radprf_mf[:,1]);

		show();sys.exit()

	########################################################################################################

	########################################################################################################
	#### simulation details			
	if add_noise == 2:
		whichmap = param_dict['whichmap']
		sims.whichmap = whichmap
	########################################################################################################

	########################################################################################################
	#### initialize simulations
	if sims.inidic['add_TF'] == 3:
		TWODTF = sims.fn_get_EBAX_2016_anal(mapparams, l1=sims.inidic['l1'], l2=sims.inidic['l2'], l3=sims.inidic['l3'])[0]
	elif sims.inidic['add_TF'] == 0:
		TWODTF = sims.fn_get_HPF(mapparams, ideal = 1)[0]
	sims.exp_beam = param_dict['exp_beam']
	sims.Dls2map(Dls_len, mapparams, CMB_outputscale = param_dict['CMB_outputscale'])
	if param_dict['perform_tSZ_removal']:
		sims.fn_beam_stuff(mapparams, use_beam = param_dict['use_beam'], nu = 90) #90 GHz beam fot tSZfree maps
	else:
		sims.fn_beam_stuff(mapparams)

	try:
		param_dict['cross_maps']
	except:
		param_dict['cross_maps'] = 0
	if param_dict['cross_maps'] == 1: #then two maps will be used: gradient and lensing map for QE
		assert param_dict['exp_beam_gradient_map']
		assert param_dict['noise_level_grad_map']
		#Bl_gradient_map = sims.fn_beam_stuff(mapparams, use_beam = param_dict['use_beam'], exp_beam = param_dict['exp_beam_gradient_map'], return_beam = 1)
		Bl_90 = sims.fn_beam_stuff(mapparams, use_beam = param_dict['use_beam'], nu = 90, return_beam = 1)
		Bl_150 = sims.fn_beam_stuff(mapparams, use_beam = param_dict['use_beam'], nu = 150, return_beam = 1)
		Bl_gradient_map = Bl_90 ### / Bl_150
		'''
		Bl_gradient_map[Bl_gradient_map == np.inf] = 0.
		Bl_gradient_map[np.where(np.isnan(Bl_gradient_map))] = 0.

		if param_dict['use_beam'] == 2: #remove this later
			Bl_gradient_map = Bl_90

		sims.Bl_gradient_map = Bl_gradient_map
		'''

		"""
		subplot(121);imshow(sims.Bl);colorbar()
		subplot(122);imshow(sims.Bl_gradient_map);colorbar()
		show();sys.exit()
		"""
	########################################################################################################
	#sys.exit()
	### kappa_model_file = '%s/stacked_kappa_model_minrich%s_maxrich%s.pkl.gz' %(ipfolder, minrich, maxrich)
	"""
	totalclus = kappadic['totalclus']
	use_real_Mz_dist = param_dict['use_real_Mz_dist']
	mapfolder = 'data/sanjay_maps_201705xx/final_maps/20171027_downsample_x2_tSZfree_plus_150_v1_no_masking'
	data_file = '%s/map_150.pkl.gz_cutouts_150ghz_no_downsampling_year3_lgt5_redmapper_clusters_TQUforMLE' %(mapfolder)
	CLUSKEYNAMES = pickle.load(gzip.open(data_file, 'rb'))['cutouts'].keys()
	selectedinds = np.arange(totalclus) #sorted( np.random.choice(np.arange(len(CLUSKEYNAMES)), size = totalclus, replace = 0) )
	richarr_RM = np.asarray( map(lambda x: x[3], CLUSKEYNAMES) )[selectedinds]
	z_L_list = np.asarray( map(lambda x: x[2], CLUSKEYNAMES) )[selectedinds]
	M_200_list = sims.fn_rich_mass_M17(richarr_RM, z_L_list)
	"""
	#get true mass here
	richarr_RM = np.asarray( map(lambda x: x[3], CLUS_IDENTIFIER) )
	z_L_list = np.asarray( map(lambda x: x[2], CLUS_IDENTIFIER) )
	M_200_list = sims.fn_rich_mass_M17(richarr_RM, z_L_list)
	truemass = np.mean(M_200_list)
	totalclus = len(CLUS_IDENTIFIER)

	opfolder = '%s/M_lambda_modelv1.0' %(ipfolder)
	kappa_model_file = '%s/stacked_kappa_model_minrich%s_maxrich%s.pkl.gz' %(opfolder, minrich, maxrich)

	if 1==1:
		#opfolder = '%s/M_lambda_model' %(ipfolder)
		## opfolder = '%s/M_lambda_modelv1.0' %(ipfolder)
		if not os.path.exists(opfolder): os.system('mkdir %s' %(opfolder))
		if os.path.exists(kappa_model_file):
			kappa_model_dic = pickle.load(gzip.open(kappa_model_file, 'rb'))
		else:
			kappa_model_dic = {}

		for acnt, AA in enumerate( Aarr ):

			keyname = AA/1e14
			logline = '\tGetting stacked kappa model for mass = %s e14 solar mass' %(keyname)
			logfile = open(log_file,'a');logfile.writelines('%s\n' %(logline));logfile.close()
			print logline

			if keyname not in kappa_model_dic: 
				### STACKED_KAPPA_MODEL = fn_get_STACKED_KAPPA_MODEL_FN_MASS(MM, CLUS_IDENTIFIER)
				### STACKED_KAPPA_MODEL = fn_get_STACKED_KAPPA_MODEL_FN_MASS(AA, CLUS_IDENTIFIER)
				### kappa_model_dic[keyname] = STACKED_KAPPA_MODEL

				'''
				returns = sims.fn_get_STACKED_KAPPA_MODEL_FN_MASS(AA, alpha_fit, RA, DEC, CLUS_IDENTIFIER, param_dict, mapparams, TWODTF)
				if len(returns) == 3:
					STACKED_KAPPA_MODEL, STACKED_KAPPA_MODEL_SUM, totalsummed = returns
				else:
					STACKED_KAPPA_MODEL = returns
					STACKED_KAPPA_MODEL_SUM, totalsummed = None, None
				kappa_model_dic[keyname] = [STACKED_KAPPA_MODEL, STACKED_KAPPA_MODEL_SUM, totalsummed]
				'''
				if param_dict['add_TF'] == 1:
					TWODTF = sims.fn_get_HPF(mapparams, minel = param_dict['min_el'], maxel = param_dict['max_el'])[0]
				rad_profile_info = [binsize, minbin, maxbin]
				#from IPython import embed;embed()
				KAPPA_MODELS_rad_profiles = sims.fn_get_STACKED_KAPPA_MODEL_FN_MASS(AA, alpha_fit, RA, DEC, CLUS_IDENTIFIER \
				, param_dict, mapparams, Bl = sims.Bl, TWODTF = TWODTF, kappa_model_fixed = 0 \
				, rad_profile_info = rad_profile_info, use_TF_beam = 1, perform_offsets_correction = 0)

				kappa_model_dic[keyname] = KAPPA_MODELS_rad_profiles

				pickle.dump(kappa_model_dic, gzip.open(kappa_model_file, 'wb'), protocol = 2)

		#from IPython import embed;embed()

		raprfmap = sims.fn_radial_profile(stacked_kappa_qe, RADEC, bin_size=binsize, minbin=0.0, maxbin=maxbin)
		d = raprfmap[:,1].flatten() ### - np.mean(raprfmap[:,1].flatten())
		DATA_VECTOR = np.copy(d)
		cov_file = '%s/kappa_COV_%sJK_1d_%sam_%sdx_richnessbins_%s_%s.pkl.gz' %(ipfolder, noofsims, maxbin, binsize, minrich, maxrich)
		kappa_COV = pickle.load(gzip.open(cov_file, 'rb'))
		kappa_std = np.diag(kappa_COV)**0.5
		kappa_COV = np.mat(kappa_COV)
		sign, logdetval = np.linalg.slogdet(kappa_COV)
		logdetval = logdetval * sign

		Cinv = sc.linalg.pinv2(kappa_COV)
		correction_factor = (float(noofsims) - len(DATA_VECTOR) - 1.)/noofsims
	#print correction_factor, float(noofsims), len(d)
		Cinv = Cinv * correction_factor
		logLarr = []
		for acnt, AA in enumerate( Aarr ):
		
			keyname = AA/1e14
			model_arr = kappa_model_dic[keyname]
			model_arr_1d= sims.fn_radial_profile(model_arr, RADEC, bin_size=binsize, minbin=0.0, maxbin=maxbin)
			model_arr_1d = model_arr_1d[:,1].flatten() 
			diff_arr = DATA_VECTOR - model_arr_1d
			##from IPython import embed; embed()

			logLval = -0.5 * np.asarray( np.dot(diff_arr.T, np.dot( Cinv, diff_arr ))).squeeze()

			#print AA, logLval

			logLarr.append( logLval )
		tmp = logLarr - max(logLarr); L = np.exp(tmp); L = L/max(L)
		plot(Aarr, L)
		recov_mass = Aarr[np.argmax(logLarr)]

		res_arr = np.copy(logLarr)
		res_arr *= 2.

		res_arr_at_peak_mass = max(res_arr)
		delta_chisq = res_arr_at_peak_mass - res_arr
		width = fn_calc_delta_chisq_one_width(Aarr, delta_chisq, recov_mass)
	
		det_sig = fn_get_det_significance(Aarr, logLarr)		
		det_arr[fcnt] = det_sig
		print det_sig
		if fcnt == len(fnamearr)-1:
			from IPython import embed;embed()	
		continue
		#sys.exit()
	#from IPython import embed;embed()
	'''
	#get model from data as these are the same clusters
	opfolder = '%s/M_lambda_model_with_slopev2.0' %(ipfolder)

	try:
		kappa_model_dic
	except:
		kappa_model_dic = {}
		for acnt, AA in enumerate( Aarr ):
			keyname = AA/1e14
			print keyname
			kappa_model_dic[keyname] = [fn_get_kappa_model_A_dic_from_stored_models(AA, [alpha_fit], minrich, maxrich, inv_var_weights = 0)[alpha_fit]]
	'''

	########################################################################################################
	########################################################################################################
	########################################################################################################
	#now get kappa_COV from sims
	if 0==1:#param_dict['use_mask']:
		totthreads = 4
		os.putenv('OMP_NUM_THREADS',str(totthreads))
		totalclus_for_cov = 100 #kappadic['totalclus']
		sims_cov_file = '%s/kappa_COV_%ssims_%sam_%sdx_%scluters.pkl.gz' %(ipfolder, noofsims, maxbin, binsize, totalclus_for_cov)
		MF_file = '%s/MEANFIELD_%ssims_%sam_%sdx_%scluters.pkl.gz' %(ipfolder, noofsims, maxbin, binsize, totalclus_for_cov)
		if param_dict['use_mask'] and not os.path.exists(MF_file):
			STACKED_KAPPA_QE_SIMS, MEAN_FIELD_FROM_STACKED_KAPPA_QE_SIMS = fn_get_QE_sims_for_COV_MF(noofsims, totalclus_for_cov, mapparams, param_dict, randomseedval, return_all_kappa_sims = 1)

			#also get mean field from these sims
			#imshow(MEAN_FIELD_FROM_STACKED_KAPPA_QE_SIMS);colorbar();show()
			#sys.exit()

			RADPROFILES = np.asarray( map(lambda x: sims.fn_radial_profile(x, RADEC, bin_size=binsize, minbin=0.0, maxbin=maxbin), STACKED_KAPPA_QE_SIMS) )
			RADPRF = RADPROFILES[:,:,1]
			totbins = np.shape(RADPRF)[1]
			noofsims = len(STACKED_KAPPA_QE_SIMS)
			sims_kappa_COV = sims.calcCov(RADPRF, noofsims, npixels = totbins)
			#imshow(kappa_COV);colorbar();show();sys.exit()
			pickle.dump(kappa_COV, gzip.open(sims_cov_file, 'wb'), protocol = 2)
			pickle.dump(MEAN_FIELD_FROM_STACKED_KAPPA_QE_SIMS, gzip.open(MF_file, 'wb'), protocol = 2)
			sys.exit()
		else:
			sims_kappa_COV = pickle.load(gzip.open(sims_kappa_COV, 'rb'))

		scalefactor = kappadic['totalclus']/totalclus_for_cov
		sims_kappa_COV = sims_kappa_COV/scalefactor
		###imshow(sims_kappa_COV);colorbar();show();sys.exit()
		MEAN_FIELD = pickle.load(gzip.open(MF_file, 'rb'))

	else:
		cov_file = '%s/kappa_COV_%sJK_1d_%sam_%sdx_richnessbins_%s_%s.pkl.gz' %(ipfolder, noofsims, maxbin, binsize, minrich, maxrich)
		kappa_COV = pickle.load(gzip.open(cov_file, 'rb'))

	if sims.is_seq(MEAN_FIELD):
		stacked_kappa_qe = stacked_kappa_qe - MEAN_FIELD
	BEST_FIT_MODELS['MEAN_FIELD'] = MEAN_FIELD
	BEST_FIT_MODELS['kappa_COV'] = kappa_COV

	#kappa_COV = np.eye(kappa_COV.shape[0]) * np.diag(kappa_COV)
	#clf();imshow(kappa_COV);colorbar();show()#;sys.exit()

	########################################################################################################
	########################################################################################################
	########################################################################################################
	#get weights here

	inv_var_weights = 1
	if inv_var_weights:
		extra_variance_weights = 1
		INV_VAR_WEIGHTS_full, INV_VAR_WEIGHTS_full_dic, INV_VAR_WEIGHTS_SPLIT_full, extra_var_params_plaw_fit = fn_get_weights(RA, DEC, CLUS_IDENTIFIER_full, kappa_qe_arr_full, kappa_weights = inv_var_weights, extra_variance_weights = extra_variance_weights, testing = 0)
	else:
		extra_variance_weights = 0
		INV_VAR_WEIGHTS_full, INV_VAR_WEIGHTS_full_dic = fn_get_weights(RA, DEC, CLUS_IDENTIFIER_full, kappa_qe_arr_full, kappa_weights = inv_var_weights, extra_variance_weights = extra_variance_weights, use_lgmca = use_lgmca)##, testing = 1)

	kappa_qe_arr, CLUS_IDENTIFIER, INV_VAR_WEIGHTS = sims.fn_select_clusters(kappa_qe_arr_full, CLUS_IDENTIFIER_full, minrich, maxrich, z1 = 0., weights_full = INV_VAR_WEIGHTS_full)
	if inv_var_weights:
		stacked_kappa_qe = np.sum( np.asarray( map(lambda x, y: x * y, kappa_qe_arr, INV_VAR_WEIGHTS) ), axis = 0) / np.sum( INV_VAR_WEIGHTS )
	else:
		stacked_kappa_qe = np.mean( kappa_qe_arr, axis = 0)

	truemass = np.sum( np.asarray( map(lambda x, y: x * y, M_200_list, INV_VAR_WEIGHTS) ), axis = 0) / np.sum( INV_VAR_WEIGHTS )

	########################################################################################################
	########################################################################################################
	########################################################################################################

	#now use the above models to get the likelihood curve
	#sys.exit()


	fix_centroid = 0
	if fix_centroid:
		
		dummy_apodMASK = fn_apodize(stacked_kappa_qe, mapparams, mask = 'circle', just_return_mask = 1)
		tmpdata = stacked_kappa_qe * dummy_apodMASK
		height, amp = 0, np.max(tmpdata)
		ypini, xpini = np.unravel_index(np.argmax(tmpdata), tmpdata.shape)
		xpini, ypini = (xpini-nx/2), (ypini-ny/2)
		x_cen, y_cen = xpini*dx, ypini*dx
		wx, wy = 1., 1. #arcmins
		rot = 0. #degrees
		p0 = [height, amp, x_cen, y_cen, wx, wy, rot]
		p1, pcov, infodict, errmsg, success = optimize.leastsq(gaussian_fitting_func, p0[:], args=(p0, RA*60., DEC*60., stacked_kappa_qe), full_output=1)
		xcen, ycen = p1[2], p1[3]
		### imshow(stacked_kappa_qe, interpolation = 'bicubic',extent = [min(ra)*60., max(ra)*60., min(dec)*60., max(dec)*60.], origin = 'lower')
		### plot(xcen, ycen, 'k+');show()
		xp, yp = xcen / dx, ycen / dx
		xp = int(round(xp, 0))
		yp = int(round(yp, 0))
		if xp>3 or yp>3:
			xp,yp = xpini, ypini
		'''
		from IPython import embed; embed()
		stacked_kappa_qe_2 = np.roll(np.roll(stacked_kappa_qe, -yp, axis = 0), -xp, axis = 1)
		subplot(121);imshow(stacked_kappa_qe, origin = 'lower', interpolation = 'bicubic');colorbar(); grid(1)
		subplot(122);imshow(stacked_kappa_qe_2, origin = 'lower', interpolation = 'bicubic');colorbar(); grid(1)
		show();sys.exit()
		'''
		stacked_kappa_qe = np.roll(np.roll(stacked_kappa_qe, -yp, axis = 0), -xp, axis = 1)

	if not onedfit:
		sn = 10.
		boxsize = nx * dx
		e1, e2 = int(nx/2 - sn/dx/2), int(nx/2 + sn/dx/2)
		stacked_kappa_qe = stacked_kappa_qe[e1:e2, e1:e2]

	###############################################################
	#data vector here
	raprfmap = sims.fn_radial_profile(stacked_kappa_qe, RADEC, bin_size=binsize, minbin=0.0, maxbin=maxbin)
	d = raprfmap[:,1].flatten() ### - np.mean(raprfmap[:,1].flatten())
	DATA_VECTOR = np.copy(d)

	###############################################################
	#inverse cov. matrix and corrections
	kappa_std = np.diag(kappa_COV)**0.5
	kappa_COV = np.mat(kappa_COV)
	sign, logdetval = np.linalg.slogdet(kappa_COV)
	logdetval = logdetval * sign

	Cinv = sc.linalg.pinv2(kappa_COV)
	correction_factor = (float(noofsims) - len(DATA_VECTOR) - 1.)/noofsims
	#print correction_factor, float(noofsims), len(d)
	Cinv = Cinv * correction_factor
	###############################################################

	def fn_get_kappa_model_A_dic_from_stored_models(kappa_model_file, weights = None):

		dic = pickle.load(gzip.open(kappa_model_file, 'rb'))
		kappa_model_dic = {}
		for keyname in dic.keys():
			modelrafprf_arr = dic[keyname]

			if weights <> None:
				modelrafprfsum = np.sum( np.asarray( map(lambda x, y: x * y, modelrafprf_arr, weights) ), axis = 0) 
				weightsum = np.sum( weights )
				kappa_model_dic[keyname] = modelrafprfsum/weightsum

			else:
				kappa_model_dic[keyname] = np.mean(modelrafprf_arr, axis = 0)

		return kappa_model_dic

	kappa_model_dic = fn_get_kappa_model_A_dic_from_stored_models(kappa_model_file, weights = INV_VAR_WEIGHTS)

	plot_onedfit = 0
	if plot_onedfit and local: clf()
	logLarr = []
	### for mcnt, MM in enumerate( Marr ):
	for acnt, AA in enumerate( Aarr ):
		### keyname = MM/1e14
		keyname = AA/1e14

		'''
		STACKED_KAPPA_MODEL = kappa_model_dic[keyname][0]
		
		#now incorporate beam and TF into the model - 20180430
		STACKED_KAPPA_MODEL = fn_add_tf_beam_to_model(STACKED_KAPPA_MODEL)

		#if len(STACKED_KAPPA_MODEL) == 3: STACKED_KAPPA_MODEL = STACKED_KAPPA_MODEL[0]

		#subplot(121);imshow(STACKED_KAPPA_MODEL[80:120, 80:120]);colorbar();title(AA)
		#subplot(122);imshow(stacked_kappa_qe[80:120, 80:120]);colorbar();show()
		if not onedfit:
			STACKED_KAPPA_MODEL = STACKED_KAPPA_MODEL[e1:e2, e1:e2]

		if 0==1:#AA>0:
			subplot(121);imshow(stacked_kappa_qe);colorbar();#title(sims.inidic['l1_for_qe'])
			subplot(122);imshow(STACKED_KAPPA_MODEL);colorbar();title(AA);show()#;sys.exit()
		'''

		#logL = fn_calc_likelihood(stacked_kappa_qe, STACKED_KAPPA_MODEL, kappa_COV, AA, onedfit = onedfit, plot_onedfit = plot_onedfit)#, kappa_two_halo_term = kappa_two_halo_term)
		model_arr = kappa_model_dic[keyname]
		diff_arr = DATA_VECTOR - model_arr
		##from IPython import embed; embed()

		logLval = -0.5 * np.asarray( np.dot(diff_arr.T, np.dot( Cinv, diff_arr ))).squeeze()

		print AA, logLval

		logLarr.append( logLval )


	#plot(theta*60., kappa_two_halo_term, 'm', label = r'\textbf{$\kappa_{2h}$}')
	recov_mass = Aarr[np.argmax(logLarr)]

	res_arr = np.copy(logLarr)
	res_arr *= 2.

	res_arr_at_peak_mass = max(res_arr)
	delta_chisq = res_arr_at_peak_mass - res_arr
	width = fn_calc_delta_chisq_one_width(Aarr, delta_chisq, recov_mass)
	### det_sig = fn_get_det_significance(Marr, logLarr)
	det_sig = fn_get_det_significance(Aarr, logLarr)

	if 0 == 1:
		RES_DIC = {}
		RES_DIC['recov_mass'] = recov_mass
		RES_DIC['Aarr'] = Aarr
		RES_DIC['logLarr'] = logLarr
		RES_DIC['width'] = width
		RES_DIC['det_sig'] = det_sig
		RES_DIC['delta_chisq'] = delta_chisq
		RES_DIC['stacked_kappa_qe'] = stacked_kappa_qe
		RES_DIC['kappa_COV'] = kappa_COV
		RES_DIC['noofsims'] = noofsims
		RES_DIC['minrich'] = minrich
		RES_DIC['maxrich'] = maxrich
		totalclus = len(kappa_qe_arr)

		RES_DIC['totalclus'] = totalclus
		fd= '%s/for_20171229_reports' %(ipfolder)
		if not os.path.exists(fd): os.system('mkdir %s' %(fd))
		res_dic_file = '%s/%sclusters_%sminrich_%smaxrich_%ssims.pkl.gz' %(fd, totalclus, minrich, maxrich, noofsims )
		pickle.dump(RES_DIC, gzip.open( res_dic_file, 'wb'), protocol = 2)

	if plot_onedfit and local:
		xlim(0.,10.); legend(loc=1)
		#xlabel(r'\textbf{$\theta$}', fontsize = 12)
		xlabel(r'\textbf{arcmin}', fontsize = 12)
		ylabel(r'\textbf{Convergence $\kappa$}', fontsize = 12)
		tit = r'\textbf{Mass = %s $\pm$ %.3f; Det. significance = %.3f; Min. rich = %s; Max. rich = %s}' %(recov_mass/1e14, width/1e14, det_sig, minrich, maxrich)
		title(tit, fontsize = 10)
		show()

	#tit = 'Mass = %s + %.3f; Det. significance = %.3f; Min. rich = %s; Max. rich = %s' %(recov_mass/1e14, width/1e14, det_sig, minrich, maxrich)
	tmp = logLarr - max(logLarr); L = np.exp(tmp); L = L/max(L)
	if plot_onedfit and local: 
		clf();
		plot(Aarr, L, color = cm.gnuplot(200)); #title(fname, fontsize = 10);
		axvline(2.35e14)
		show();sys.exit()

	#store the best fits for plots
	recov_mass_keyname = recov_mass/1e14
	'''
	BF_STACKED_KAPPA_MODEL = kappa_model_dic[recov_mass_keyname][0]
	#if len(BF_STACKED_KAPPA_MODEL) == 3: BF_STACKED_KAPPA_MODEL = BF_STACKED_KAPPA_MODEL[0]
	BF_raprfmodel = sims.fn_radial_profile(BF_STACKED_KAPPA_MODEL, RADEC, bin_size=binsize, minbin=0.0, maxbin=maxbin)
	BEST_FIT_MODELS[fname] = [recov_mass, BF_STACKED_KAPPA_MODEL, BF_raprfmodel, width, det_sig, Aarr, L, totalclus, stacked_kappa_qe, logLarr]
	'''
	BF_raprfmodel = kappa_model_dic[recov_mass_keyname][0]
	BEST_FIT_MODELS[fname] = [recov_mass, None, BF_raprfmodel, width, det_sig, Aarr, L, totalclus, stacked_kappa_qe, logLarr, INV_VAR_WEIGHTS]

	lab = r'\textbf{%s + %.3f; %.3f$\sigma$}' %(recov_mass/1e14, width/1e14, det_sig)
	if local:
		plot(Aarr, L, color = colorarr[fcnt], label = lab, lw = 1.)#; title(tit, fontsize = 10);
		### truemass = 2.35e14
		if fcnt == 0: vlines(truemass, 0., 1., lw = 2., linestyle = '--', color = 'r')
		xlabel(r'\textbf{Normalisation - $M-\lambda$ relation}')
		tit = 'SPTpol-like simulations: DES RM year 3 sample. Clusters = %s' %(kappadic['totalclus'])
		if  param_dict['use_beam'] == 1:
			plname = 'sptpol_desrm_sims_gaussian_beam.png'
			tit = '%s; Beam = %s\' Gaussian' %(tit, param_dict['exp_beam'])
		else:
			tit = '%s; Beam = SPTpol' %(tit)
			plname = 'sptpol_desrm_sims_sptpol_beam.png'
		if ipfolder.find('with_FG')>-1: 
			plname = '%s_with_FG.png' %(plname.replace('.png',''))
			tit = '%s; FG = 1' %(tit)
		else:
			tit = '%s; FG = 0' %(tit)
		tit = r'\textbf{%s}' %(tit)
		title(tit, fontsize = 10);

if local:
	legend(loc = 2, fontsize = 6.5)
	#show()
	plfolder = 'reports/20180125_beamfixed'
	plname = '%s/%s' %(plfolder, plname)
	savefig(plname, bbox_to_inches = 'tight', pad_inches = 0.1)
	close()


## now plot all best-fits

'''
if truemass in kappa_model_dic.keys():
	trueprofile = kappa_model_dic[truemass/1e14]
else:
	truemass_rounded_1 = Aarr[np.argmin(abs(round(truemass, 1) - Aarr))]
	truemass_rounded_2 = Aarr[np.argmin(abs(round(truemass, 1) + (delA * 1e14) - Aarr))]

	trueprofile_1 = kappa_model_dic[truemass_rounded_1/1e14][0]#[0]
	trueprofile_2 = kappa_model_dic[truemass_rounded_2/1e14][0]#[0]
	trueprofile = (trueprofile_1 + trueprofile_2)/2.
trueprofile_radial = sims.fn_radial_profile(trueprofile, RADEC, bin_size=binsize, minbin=0.0, maxbin=maxbin)
'''
trueprofile = None
trueprofile_radial = None
trueA = 2.35e14

BEST_FIT_MODELS['true'] = [truemass, trueA, trueprofile, trueprofile_radial]
### BEST_FIT_MODELS['CLUS_IDENTIFIER'] = CLUS_IDENTIFIER

op_file_sims = '%s/sim_results_for_plotting_JK_%s_%s.pkl.gz' %(ipfolder, float(minrich), float(maxrich))
pickle.dump(BEST_FIT_MODELS, gzip.open(op_file_sims, 'wb'), protocol = 2)

if not local: sys.exit()
####################################################################################################
####################################################################################################
####################################################################################################
### op_file_sims = '%s/sim_results_for_plotting.pkl.gz' %(ipfolder)
BEST_FIT_MODELS = pickle.load(gzip.open(op_file_sims, 'rb'))

kappa_COV = BEST_FIT_MODELS['kappa_COV']
MEAN_FIELD = BEST_FIT_MODELS['MEAN_FIELD']
true = BEST_FIT_MODELS['true']
true_m, trueA, true_stack, true_radprf = true

from scipy.interpolate import interp1d
true_x, true_y = true_radprf[:,0], true_radprf[:,1]
ip_x = np.arange(0.5, maxbin - 0.5, 0.1)
interp_type = 'cubic'
fninterp = interp1d(true_x, true_y, kind = interp_type, bounds_error = 0, fill_value = 0.)
### ip_y = np.interp(ip_x, x, y)
true_ip_y = fninterp(ip_x)

clf()
subplots_adjust(wspace=0.0)
labfs = 11
cov_std = np.diag(kappa_COV)**0.5
### colorval = cm.gnuplot(200)
colorval = 'darkorchid'
alphaval, lwval = 1., 0.2
###plot(ip_x, true_ip_y, 'k', lw = 2.)
colorarr = [cm.jet(int(d)) for d in np.linspace(0,255,len(BEST_FIT_MODELS))]

for kcnt, keyname in enumerate( BEST_FIT_MODELS.keys() ):
	if keyname == 'kappa_COV' or keyname == 'MEAN_FIELD' or keyname == 'true': continue

	M, stack, radprf, width, det_sig, Aarr, L, totalclus, stacked_kappa_qe  = 	BEST_FIT_MODELS[keyname]

	ax  = subplot(121)
	x, y = radprf[:,0], radprf[:,1]
	fninterp = interp1d(x, y, kind = interp_type, bounds_error = 0, fill_value = 0.)
	ip_y = fninterp(ip_x)
	plot(ip_x, ip_y, color = colorval, alpha = alphaval, lw = 0.5)

	ax  = subplot(122)

	plot(Aarr/1e14, L, lw = 0.5, color = colorval, alpha = alphaval)#colorarr[kcnt])
	try:
		totalL.append(L)
	except:
		totalL = [L]

totalL = np.asarray(totalL)
combL = np.prod(totalL, axis = 0); combL/= max(combL)
plot(Aarr/1e14, combL, 'k')#colorarr[kcnt])

ax  = subplot(121)
labval = r'\textbf{Recovered}'
plot(1e10, 0., lw = 6., alpha = alphaval, label = labval, color = colorval)
axhline(color = 'k', lw = 1., linestyle = 'solid')
#labval_true = r'\textbf{Input: M$_{200,m}$ = %.2f $\times$ 10 $^{14}$ M$_{\odot}$}' %(true_m/1e14)
#labval_true = r'\textbf{Input: %.2f $\times$ 10 $^{14}$ M$_{\odot}$}' %(true_m/1e14)
labval_true = r'\textbf{M$_{\rm true}$: %.2f $\times$ 10 $^{14}$ M$_{\odot}$}' %(true_m/1e14)
errorbar(true_radprf[:,0], true_radprf[:,1], yerr = cov_std, color = 'k', ecolor = 'k', ms = 6., capsize = 1., ls = 'None', marker = 'o', label = labval_true, alpha = 1., zorder=100)
legend(loc = 1, fontsize = labfs - 3, framealpha = 1., numpoints = 2)
xlim(0., 10.)
xlabel(r'\textbf{Radial distance [arcmin]}', fontsize = labfs)
ylabel(r'\textbf{Convergence $\kappa$}', fontsize = labfs)
ax.xaxis.set_major_locator(MaxNLocator(nbins=5, prune='upper'))
#title(r'\textbf{Simulations - DES redMaPPer Year 3-like sample}', fontsize = 12)

ax  = subplot(122)
### axvline(true_m/1e14, linestyle = '-.', lw = 2., color = 'k', label = labval_true)#, zorder = 10)
axvline(trueA/1e14, linestyle = '-.', lw = 2., color = 'k', label = labval_true)#, zorder = 10)
legend(loc = 2, fontsize = labfs - 3, framealpha = 1., numpoints = 2)

figtext(0.18,0.9, r'\textbf{Simulations - DES redMaPPer Year 3-like sample}', fontsize = 14)
##suptitle(r'\textbf{Simulations - DES redMaPPer Year 3-like sample}', fontsize = 14, verticalalignment = 'center')
xlabel(r'\textbf{Normalisation A [10$^{14}$ M$_{\odot}$]}', fontsize = labfs)
ylabel(r'\textbf{Normalised $\mathcal{L}$}', fontsize = labfs)
ax.yaxis.tick_right(); ax.yaxis.set_label_position('right')
ylim(-0.01,None)
xlim(0.,4.)
show();sys.exit()

plname = 'reports/2018028_fittingslope_polarisation/sims_kappa_model_MF_null_y3_lgt5_v6_4_21_full.pdf'
savefig(plname, bbox_to_inches = 'tight', pad_inches = 0.1)
sys.exit()



########################################################################################################
########################################################################################################
########################################################################################################

