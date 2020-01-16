"""
Takes in folder with stacked kappa as input

Output:
Likelihood fitted to the observed kappa


*Note : it has to go through changes for different richness bins
"""

import numpy as np
import pickle,gzip,glob
import sys
import modules.scl_cmb as scl
import numpy as np, pickle, gzip, sys, glob, time, os, scipy as sc, pdb

def fitting_func_gaussian(p, p0, X, DATA = None, lbounds = None, ubounds = None, fixed = None, return_fit = 0, gauss = 1, cluster_mass = None, commit_20161221_change = 0):
    if hasattr(fixed, '__len__'):
        p[fixed] = p0[fixed]

    if hasattr(lbounds, '__len__'):
        linds = abs(p)<abs(lbounds)
        p[linds] = lbounds[linds]

    if hasattr(ubounds, '__len__'):
        uinds = abs(p)>abs(ubounds)
        p[uinds] = ubounds[uinds]

    if not gauss:
        if commit_20161221_change:
            fitfunc = lambda p, X: p[2] + ((X - p[1])) ** 2. / p[0] ### added on 20161221 for better fitting
        else:

            fitfunc = lambda p, X: ((X - cluster_mass)) ** 2. / p[0]
        if not return_fit:
            return fitfunc(p, X) - DATA
        else:

            return fitfunc(p, X)
    else:

        fitfunc = lambda p, x: p[1]*(np.exp(-(x-p[2])**2/(2*p[3]**2)))
        if not return_fit:
            return fitfunc(p, X) - DATA
        else:
            return fitfunc(p, X)


def fn_gaussian_1d(xarr,yarr, p0 = None, return_just_fit = 0):

    import scipy.optimize as optimize
    import scipy.special

    #fitfunc = lambda p, x: p[1]*(np.exp(-(x-p[2])**2/(2*p[3]**2)))
    fitfunc = lambda p, x: p[1]*(np.exp(-(x-p[2])**2/(2*p[3]**2)))
    #fitfunc = lambda p, x: p[1] * x**(p[0]-1) * np.exp(-x**2/2) / (2**(p[0]/2-1) * scipy.special.gamma(p[0]/2))
    errfunc = lambda p, x, y: fitfunc(p, x) - y #minimization

    if return_just_fit: #then yarr is param array
        p1 = yarr
        return fitfunc(p1,xarr)
    else:

        #if p0 == None:
        #    p0 = [np.std(yarr),np.max(yarr),xarr[np.argmax(yarr)],np.std(xarr)*2]
        p1, success = optimize.leastsq(errfunc, p0[:], args=(xarr, yarr))
        return p1

def fn_likelihood_finer_resol(M, L, cluster_mass = None,delM=0.001, max_M = 5.):

    import scipy.optimize as optimize

    if max(M)>100.:
        M/=1e14

    if cluster_mass == None: cluster_mass = M[np.argmax(L)]

    #print min(M),max(M)
    M_ip = np.arange(min(M),max(M),delM) #5e11
    if max(M)<max_M:
        M_ip = np.arange(min(M),max_M,delM) #5e11

    #first guess a good parameter
    gau_width = abs(cluster_mass - M[np.argmin(abs(L))])/2.35 * 2.
    p0 = [0.,np.max(L),cluster_mass,gau_width]
    lbounds = np.asarray([0., 1., cluster_mass*0.03, gau_width * 2.])
    ubounds = np.asarray([1e-2, 1., cluster_mass*1.03, gau_width * 0.1])

    #fixed = np.asarray([0.,1.,0.,0.])
    p1, success = optimize.leastsq(fitting_func_gaussian, p0, args=(p0, M, L, lbounds, ubounds))
    gauss_params=fn_gaussian_1d(M, L, p1)
    L_ip = fn_gaussian_1d(M_ip, gauss_params, return_just_fit = 1)

    return M_ip, L_ip, gauss_params




def fn_get_STACKED_KAPPA_MODEL_FN_MASS(MM,RA, DEC, CLUS_IDENTIFIER, param_dict, mapparams, Bl = None, TWODTF = None, clra = 0., cldec = 0., profile_name = 'NFW', kappa_model_fixed = 0, null_test = 0, theta_max = -1., weights = None, use_TF_beam = 0, rad_profile_info = None, perform_offsets_correction = 1, reqd_shape = None):

    zvals = np.asarray( map(lambda x: x[2], CLUS_IDENTIFIER) )
    if len(np.unique(zvals)) == 1:
        kappa_model_fixed = 1
    totalclus = len(zvals)

    for kcnt, cluskey in enumerate(CLUS_IDENTIFIER):
        if not null_test:
            try:
                ra, dec, z_val, rich, weight, weights_norm = cluskey
            except:
                ra, dec, z_val, rich, weight = cluskey
        else:
            ra, dec, z_val, rich, weight, weights_norm = cluskey

        if kappa_model_fixed:
            z_val = np.mean(zvals) # set to 0.7 for fixed kappa
            print "Model is fixed to redshift of %s"%(z_val)
            cc = sims.c_Duffy(MM, z_val, param_dict['h'], profile_name = profile_name)
            sims.fn_lensing_ini(mapparams,param_dict, RA, DEC, [clra], [cldec], [MM], [cc], [z_val], param_dict['z_lss'], param_dict['mass_def'], param_dict['rho_def'], truncate_kappa_at_radius = param_dict['truncate_kappa_at_radius'], profile_name = profile_name, theta_max = theta_max)
            kappa_model = sims.KAPPA
        else:
            cc = sims.c_Duffy(MM, z_val, param_dict['h'], profile_name = profile_name)
            sims.fn_lensing_ini(mapparams,param_dict, RA, DEC, [clra], [cldec], [MM], [cc], [z_val], param_dict['z_lss'], param_dict['mass_def'], param_dict['rho_def'], truncate_kappa_at_radius = param_dict['truncate_kappa_at_radius'], profile_name = profile_name, theta_max = theta_max)
            kappa_model = sims.KAPPA

        # cut beam and TF modes
        if use_TF_beam:
            try:
                above_beam_scale = np.where(L>=ngaus)
            except:
                lx, ly = sims.get_lxly(mapparams)
                L = (lx**2. + ly**2.)**0.5
                ngaus = int( np.sqrt(8. * np.log(2.)) / np.radians(param_dict['exp_beam']/60.) )
                above_beam_scale = np.where(L>=ngaus)

            beamtf_lmat = Bl * TWODTF
            to_deconvolve_lmat = np.copy(beamtf_lmat)
            bad = np.where(to_deconvolve_lmat==0.)
            to_deconvolve_lmat[bad] = 1.
            deconvolve_lmat = 1./to_deconvolve_lmat
            highly_filtered = np.where((deconvolve_lmat > 1.0e10) | (deconvolve_lmat < 0) )

            kappa_model_fft = np.fft.fft2(kappa_model)
            kappa_model_fft[above_beam_scale] = 0.
            kappa_model_fft[highly_filtered] = 0.
            kappa_model_fft[bad] = 0
            kappa_model = np.fft.ifft2( kappa_model_fft ).real
            # hard coded change to stacked kappas shape
            kappa_model = kappa_model[nx/2 - int(reqd_shape[0]/2.): nx/2+int(reqd_shape[0]/2.),ny/2-int(reqd_shape[1]/2.):ny/2 + int(reqd_shape[1]/2.)]
        if kappa_model_fixed:
            return kappa_model
        else:
            if kcnt ==0:
                stacked_kappa = kappa_model
            else:
                stacked_kappa =stacked_kappa+kappa_model
    return stacked_kappa/float(totalclus)

def fitting_func(p, p0, X, cluster_mass, DATA = None, lbounds = None, ubounds = None, fixed = None, return_fit = 0):
        if hasattr(fixed, '__len__'):
                p[fixed] = p0[fixed]

        if hasattr(lbounds, '__len__'):
                linds = abs(p)<abs(lbounds)
                p[linds] = lbounds[linds]

        if hasattr(ubounds, '__len__'):
                uinds = abs(p)>abs(ubounds)
                p[uinds] = ubounds[uinds]

        fitfunc = lambda p, X: ((X - cluster_mass)) ** 2. / p[0]
        if not return_fit:
                return fitfunc(p, X) - DATA
        else:
                return fitfunc(p, X)

def fn_calc_delta_chisq_one_width(Marr, del_chi_sq_arr, cluster_mass, perform_ip = 0, fit_parabola = 1):
    """
    Fits a parabola to the data and calculates one chisquare

    """
    from scipy.interpolate import interp1d


    default_width = 1.
    if 1==1:
        reqdinds = np.where( (del_chi_sq_arr<4*25.) )
        Marr = Marr[reqdinds]
        del_chi_sq_arr = del_chi_sq_arr[reqdinds]
    value_for_error = 1.
    interp_type = 'linear'
    linds, uinds = np.where(Marr<cluster_mass)[0], np.where(Marr>cluster_mass)[0]
    minreqdindsforfitting = 2
    if len(linds) < minreqdindsforfitting  or len(uinds) < minreqdindsforfitting:
        width = default_width
        reqdinds = np.where( (del_chi_sq_arr<4*25.) )

        Marr = Marr[reqdinds]
        del_chi_sq_arr = del_chi_sq_arr[reqdinds]
        value_for_error = 1.
        interp_type = 'linear'
        linds, uinds = np.where(Marr<cluster_mass)[0], np.where(Marr>cluster_mass)[0]
        minreqdindsforfitting = 2

    fninterp = interp1d(del_chi_sq_arr[linds], Marr[linds], kind = interp_type, bounds_error = 0, fill_value = 0.)
    l_err = fninterp(value_for_error)
    fninterp = interp1d(del_chi_sq_arr[uinds], Marr[uinds], kind = interp_type, bounds_error = 0, fill_value = 0.)
    u_err = fninterp(value_for_error)
    width = 0.001
    default_width = 0.001
    if width == 0.:
        width = .05

    if fit_parabola:
        import scipy.optimize as optimize
        p0 = np.asarray([width,Marr[np.argmin(del_chi_sq_arr)], min(del_chi_sq_arr)])

        lbounds = np.asarray([p0[0] * 1e-4, 0., min(del_chi_sq_arr)])
        ubounds = np.asarray([p0[0] * 10., 10., min(del_chi_sq_arr)*1.1])
        fixed = np.asarray([False, False,True])
        p1, success = optimize.leastsq(fitting_func, p0, args=(p0, Marr,  cluster_mass, del_chi_sq_arr, lbounds))
        x_ip = np.arange(min(Marr),max(Marr),cluster_mass/100.)
        y_ip = fitting_func(p1,p1,x_ip, cluster_mass,return_fit = 1)
        linds, uinds = np.where(x_ip<=cluster_mass)[0], np.where(x_ip>=cluster_mass)[0]
        value_for_error = 1.
        interp_type = 'linear'
        fninterp = interp1d(y_ip[linds], x_ip[linds], kind = interp_type, bounds_error = 0, fill_value = 0.)
        l_err = fninterp(value_for_error)
        fninterp = interp1d(y_ip[uinds], x_ip[uinds], kind = interp_type, bounds_error = 0, fill_value = 0.)
        u_err = fninterp(value_for_error)
        width = (u_err - l_err)/2.

    if width == 0:
        width = default_width
    return width


sims = scl.simulations()
sims.tqulen = 1
ipfolder = sys.argv[1]
noofsims = int(sys.argv[2])
minrich = float(sys.argv[3])
maxrich = float(sys.argv[4])



# Covariance matrix calculation
files = sorted(glob.glob('%s/st*'%(ipfolder)))
cov_file = '%s/kappa_COV_%s_simsJK.pkl.gz' %(ipfolder,noofsims)
data_for_cov = pickle.load(gzip.open(files[0]))
boxsize = data_for_cov['param_dict']['boxsize']
dx,dy = data_for_cov['param_dict']['reso_arcmin'],data_for_cov['param_dict']['reso_arcmin']
nx,ny = int(boxsize/dx), int(boxsize/dy)
minval, maxval = -boxsize/2/60.,  boxsize/2/60.
ra = np.linspace(minval, maxval, nx)
dec = np.linspace(minval, maxval, ny)
RA, DEC = np.meshgrid(ra,dec)
simmaparams = mapparams = [nx,ny,dx,dy]
param_dict = data_for_cov['param_dict']
kappa_qe_arr  = np.asarray(data_for_cov['kappa_qe_arr'])
CLUS_IDENTIFIER = np.asarray(data_for_cov['CLUS_IDENTIFIER'])
#kappa_qe_arr, CLUS_IDENTIFIER = fn_select_clusters(kappa_qe_arr, CLUS_IDENTIFIER, minrich, maxrich)
if not os.path.exists(cov_file):
    kappa_COV = sims.fn_get_COV_from_JK(kappa_qe_arr, CLUS_IDENTIFIER, noofsims, [ (minrich, maxrich) ] )
    pickle.dump(kappa_COV, gzip.open(cov_file,'wb'), protocol = 2)
else:
    kappa_COV = pickle.load(gzip.open(cov_file, 'rb'))

stacked_kappa_shape = data_for_cov['stacked_kappa_qe'].shape

kappa_COV = np.mat(kappa_COV)
sign, logdetval = np.linalg.slogdet(kappa_COV)
logdetval = logdetval * sign
Cinv = sc.linalg.pinv2(kappa_COV)

minM, maxM, delM = 0, 12.0, 0.1
Marr = np.arange(minM,maxM+delM,delM) * 1e14
#load transfer and beam function

if param_dict['add_TF'] == 0:
    TWODTF = sims.fn_get_HPF(mapparams, ideal = 1)[0]
elif param_dict['add_TF'] == 1:
    TWODTF = sims.fn_get_HPF(mapparams, minel = param_dict['min_el'], maxel = param_dict['max_el'])[0]
elif param_dict['add_TF'] == 3:
    TOWDTF = sims.fn_get_EBAX_2016_anal(mapparams, l1=param_dict['l1'], l2=param_dict['l2'], l3=param_dict['l3'])[0]
beam_error = 0
simmapparams = mapparams
Bl = sims.fn_beam_stuff(simmapparams, beam_error = beam_error,use_beam = param_dict['use_beam'],exp_beam = param_dict['exp_beam'], return_beam = 1)

# Generate models

model_file = '%s/Model_for_all_z.pkl.gz'%(ipfolder)
if os.path.exists(model_file):
    kappa_dic = pickle.load(gzip.open(model_file))
else:
    kappa_dic = {}


for MM in Marr:
    keyname = MM/1e14
    if not keyname in kappa_dic.keys():
        # fix kappa_model_fixed = 1 to generate models at a particular redshift
        model = fn_get_STACKED_KAPPA_MODEL_FN_MASS(MM,RA, DEC, CLUS_IDENTIFIER, param_dict, mapparams, Bl = Bl, TWODTF = TWODTF, clra = 0., cldec = 0., profile_name = 'NFW', kappa_model_fixed = 0, null_test = 0, theta_max = -1., weights = None, use_TF_beam = 1,  perform_offsets_correction = 0, reqd_shape = stacked_kappa_shape)
        kappa_dic[keyname] = model#/5145.

pickle.dump(kappa_dic,gzip.open(model_file,'w'))

#fit observed data to model
total_logLarr = np.ones(len(Marr))
width_arr = []
percentage_arr = []
res = {}

# Go over all the stacked simulated CMB cluster lensing maps
for ff in files:
    # Get one dimensional data vector of the stacked kappa profiles
    data = pickle.load(gzip.open(ff))
    stkp = data['stacked_kappa_qe']
    boxsize_kappa = stkp.shape[0]*data['param_dict']['reso_arcmin']
    minval, maxval = -boxsize_kappa/2/60.,  boxsize_kappa/2/60.
    nx,ny = data['stacked_kappa_qe'].shape[0], data['stacked_kappa_qe'].shape[1]
    ra = np.linspace(minval, maxval, nx)
    dec = np.linspace(minval, maxval, ny)
    RA, DEC = np.meshgrid(ra,dec)
    RADEC = [RA,DEC]
    binsize = 1
    maxbin = 10
    DATA_VECTOR = sims.fn_radial_profile(stkp, RADEC, bin_size = binsize, minbin = 0., maxbin = maxbin)[:,1]
    correction_factor = (float(noofsims) - len(DATA_VECTOR) - 1.)/noofsims
    Cinv = Cinv * correction_factor
    logLarr = []
    for acnt, MM in enumerate( Marr ):
        keyname = MM/1e14
        model_arr = kappa_dic[keyname]
        model_arr_1d= sims.fn_radial_profile(model_arr, RADEC, bin_size=binsize, minbin=0.0, maxbin=maxbin)
        model_arr_1d = model_arr_1d[:,1].flatten()
        diff_arr = DATA_VECTOR - model_arr_1d
        logLval = -0.5 * np.asarray( np.dot(diff_arr.T, np.dot( Cinv, diff_arr ))).squeeze()
        logLarr.append( logLval )

    logLarr_tmp = np.copy(logLarr) - max(logLarr)
    L_for_mass = np.exp(logLarr_tmp)
    L_for_mass/=np.max(L_for_mass)
    L = L_for_mass
    Aarr = Marr/1e14

    Aip, Lip, p = fn_likelihood_finer_resol(np.copy(Aarr), L, max_M = 10., delM = 0.0005)

    #Aip, Lip, p = fn_likelihood_finer_resol(np.copy(Aarr), L)
    recov_values = np.asarray( sims.fn_get_width_from_sampling(Aip, Lip) )
    recov_values = np.asarray( sims.fn_get_width_from_sampling(Marr/1e14, L_for_mass) )
    recov_mass, recov_mass_low_err, recov_mass_high_err = recov_values
    print recov_mass
    width = (recov_mass_low_err +recov_mass_high_err)*0.5
    width_arr.append(width)
    total_logLarr = total_logLarr + logLarr
    percentage_uncen = (width/(recov_mass/1e14))*100.
    percentage_arr.append(percentage_uncen)
final_mass = Marr[np.argmax(total_logLarr)]

res['final_mass'] = final_mass
res['width_arr'] = width_arr
res['percentage_uncen'] = percentage_arr
res['total_logLarr'] = total_logLarr
res['Marr'] = Marr
quit()
pickle.dump(res,gzip.open('%s/new_results.pkl.gz'%(ipfolder),'w'))
print width_arr
