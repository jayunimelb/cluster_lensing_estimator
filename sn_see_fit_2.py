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

def fn_see_goodness_of_fit(KAPPA_QE, recov_mass, foldername, colorval, maxbin = 10., binsize = 1.):

	recov_mass = recov_mass*1e14

	foldername = foldername.replace('pkls/', 'data/sims/').replace('no_mf_sub','')
	paramfile = glob.glob('%s/params_planck_r_0.0_2015*' %(foldername))[0]
	param_dict = sims.fn_get_param_dict_from_paramfile(paramfile)

	nx, ny = KAPPA_QE.shape
	dx = 0.5
	boxsize = nx * dx
	clra, cldec = 0., 0.
	minval, maxval = clra-boxsize/2/60.,  clra+boxsize/2/60.
	ra = dec = np.linspace(minval, maxval, nx)
	RA, DEC = np.meshgrid(ra,dec)
	cc = 3. #change this
	zz = 0.6 #change this
	mapparams = [nx, ny, dx, dx]

	sims.fn_lensing_ini(mapparams,param_dict, RA, DEC, [clra], [cldec], [recov_mass], [cc], [zz], param_dict['z_lss'], param_dict['mass_def'], param_dict['rho_def'], truncate_kappa_at_radius = param_dict['truncate_kappa_at_radius'])
	kappa_true = sims.KAPPA

	TWODTF = sims.fn_get_EBAX_2016_anal(mapparams, l1=param_dict['l1'], l2=param_dict['l2'], l3=param_dict['l3'])[0]

	lx, ly = sims.get_lxly(mapparams)
	L = (lx**2. + ly**2.)**0.5
	sims.exp_beam = param_dict['exp_beam']
	try:
		ngaus = int( np.sqrt(8. * np.log(2.)) / np.radians(sims.exp_beam/60.) )
	except:
		ngaus = int( np.sqrt(8. * np.log(2.)) / np.radians(sims.inidic['exp_beam']/60.) )
	above_beam_scale = np.where(L>=ngaus)
	kappa_true_fft = np.fft.fft2(kappa_true)
	kappa_true_fft[above_beam_scale] = 0.
	kappa_true = (kappa_true_fft * TWODTF) #twodtf similar to data here
	kappa_true = np.fft.ifft2(kappa_true_fft).real

	subplot(121);imshow(KAPPA_QE);colorbar()
	subplot(122);imshow(kappa_true);colorbar();show()

	RADEC = (RA,DEC)
	KAPPA_QE -= np.mean(KAPPA_QE)
	#kappa_true -= np.mean(kappa_true)
	raprfmap = sims.fn_radial_profile(KAPPA_QE, RADEC, bin_size=binsize, minbin=0.0, maxbin=maxbin)
	d = raprfmap[:,1]# - np.mean(raprfmap[:,1])
	raprfmodel = sims.fn_radial_profile(kappa_true, RADEC, bin_size=binsize, minbin=0.0, maxbin=maxbin)
	dp = raprfmodel[:,1]# - np.mean(raprfmodel[:,1])

	#clf()
	errorbar(raprfmap[:,0], d, yerr = raprfmap[:,2], color = colorval, marker = 'o', ls ='None')
	plot(raprfmap[:,0], dp, color = colorval)
	show();sys.exit()


################################################################################################################################
################################################################################################################################
################################################################################################################################

import pickle, gzip, numpy as np, sys, glob
from pylab import *
from matplotlib import rc;rc('text', usetex=True);rc('font', weight='bold');matplotlib.rcParams['text.latex.preamble'] = [r'\boldmath']
from scipy.interpolate import interp1d
import modules
import modules.scl_cmb as scl_cmb

sims = scl_cmb.simulations()

#folders = sys.argv[1:2]

redmapper = int(sys.argv[1])
if redmapper:
	#folders = ['pkls/sptpol_des/crossmaps/sptpol_map_like_sims/data/0.5_no_tSZ_sptpol']
	#folders = ['pkls/sptpol_des/crossmaps/sptpol_map_like_sims/data/year3/0.5_no_tSZ_sptpol']
	#folders = ['pkls/sptpol_des/crossmaps/sptpol_map_like_sims/data/year3/0.5_no_tSZ_sptpol/no_mf_sub/']
	#folders = ['pkls/sptpol_des/crossmaps/sptpol_map_like_sims/data/v6.4.17_full/SNR_20.0_1000.0/0.5_no_tSZ_sptpol/no_mf_sub/']
	#folders = ['pkls/sptpol_des/crossmaps/sptpol_map_like_sims/data/v6.4.17_full/SNR_0.0_1000.0/0.5_no_tSZ_sptpol/no_mf_sub/']
	folders = ['pkls/sptpol_des/crossmaps/sptpol_map_like_sims/data/year3/SNR_20.0_1000.0/0.5_no_tSZ_sptpol/no_mf_sub/']
	### folders = ['pkls/sptpol_des/crossmaps/sptpol_map_like_sims/data/v6.4.17_full/SNR_30.0_40.0/0.5_no_tSZ_sptpol/no_mf_sub/']
else:
	folders = ['pkls/sptpol_des/crossmaps/sptpol_map_like_sims/data/sptpol_clusters/0.5_no_tSZ_sptpol']

for fd in folders:
	fnames = glob.glob('%s/final_*' %(fd))

	colorinds = np.linspace(0,255,len(fnames)); colorarr = np.asarray( [cm.jet( int(c) ) for c in colorinds] )
	#richarr = np.asarray( map(lambda x: float( x.rsplit('_')[-2] ), fnames) )
	richarr = map(lambda x: x.rsplit('_')[-2:], fnames)
	richarr = np.asarray( sorted( map(lambda x: ( float(x[0]) + float(x[1]) )/2., richarr) ) ) 

	for f in fnames:#[1:2]:
		dic = pickle.load(gzip.open(f, 'rb'))
		M = dic['Marr']/1e14
		logl = dic['logl'].T[0]
		L = dic['L']
		recov_mass = M[np.argmax(L)]
		stacked_kappa_qe=dic['stacked_kappa_qe']
		tit = f.rsplit('/')[-1].replace('_','\_')
		#imshow(stacked_kappa_qe, interpolation = 'bicubic');colorbar();title(tit);show()

		richval = (float( f.rsplit('_')[-2] ) + float( f.rsplit('_')[-1] ) )/2.
		colorval = colorarr[ np.where(richarr == richval)[0] ][0]
		

		#subplot(121);fn_see_goodness_of_fit(stacked_kappa_qe, recov_mass, fd, colorval)

		res_arr = np.copy(logl)
		res_arr *= 2.

		res_arr_at_peak_mass = max(res_arr)
		delta_chisq = res_arr_at_peak_mass - res_arr
		width = fn_calc_delta_chisq_one_width(M, delta_chisq, recov_mass)
		det_sig = fn_get_det_significance(M, logl)

		lab = '%s: ' %(dic['label'])
		lab = '%s$%.3f \pm %.2f$: SNR: %.2f' %(lab, recov_mass, width, det_sig)
		lab = lab.replace('redmapper:', '').strip()
		lab = lab.replace('sptpol:', '').strip()

		#subplot(122);
		plot(M,L, label =r'\textbf{%s}' %(lab), lw = 1., color = colorval)

		print recov_mass, width, det_sig

leg = legend(loc=4,fontsize = 12, framealpha = 1.)
leg.set_alpha(1.)
if redmapper:
	title(r'\textbf{DES redMaPPer year 1 masses}', fontsize = 14)
else:
	title(r'\textbf{SPTpol SZ-clusters}', fontsize = 14)
xlabel(r'\textbf{Mass} [$10^{14}\ M_{\odot}$]', fontsize = 14)
ylabel(r'\textbf{Likelihood} $\mathcal{L}$', fontsize = 14)
if redmapper:
	pl_name = '/Users/sraghunathan/Research/SPTPol/analysis/2016_11/QE_SPTpol/reports/20171109_crossmap_results/RM_richness_bins.png'
else:
	pl_name = '/Users/sraghunathan/Research/SPTPol/analysis/2016_11/QE_SPTpol/reports/20171109_crossmap_results/SPTpol_clusters.png'
#savefig(pl_name, bbox_inches = 'tight', pad_inches=0.1);close()
#sys.exit()
show()