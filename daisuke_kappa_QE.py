"""
Takes the simulated CMB lensing maps and extracts lensing convergence profiles



"""
import sys, pickle, gzip, os, time, glob, numpy as np, argparse
import astropy.io.fits as fits
import modules,scipy.ndimage as ndimage, pickle,gzip
from scipy import interpolate as intrp
if str(os.getcwd()).find('sri')>-1: import matplotlib; matplotlib.use('Agg')
from colossus.cosmology import cosmology
from colossus.halo import concentration, mass_defs
cosmology.setCosmology('planck15')
import time
timedatestr='%s' %(time.strftime('%Y%m%d_%H%M%S', time.localtime()))

# Load classes where most the functions are written
sims = modules.scl_cmb.simulations()
noisesims = modules.scl_cmb.noise_sims()
tSZ_fitting = modules.scl_cmb.tSZ_fitting()
cluster_stuff = modules.scl_cmb.cluster_stuff()
def fn_apodize(MAP, mapparams, mask = 'square', just_return_mask = 1):

	import numpy as np, scipy.ndimage as ndimage, scipy.signal as signal

	nx, ny, dx, dy = mapparams

	start =  time.time()
	pix = dx
	radius = (nx * pix)/10. #(nx * pix)/3. #changed on 20171206
	npix_cos = int(radius/pix)
	ker=np.hanning(npix_cos)
	ker2d=np.asarray( np.sqrt(np.outer(ker,ker)) )

	##from IPython import embed; embed()
	MASKf=np.zeros((nx,ny))
	minval, maxval = -(nx*pix)/2,  (nx*pix)/2
	x = y = np.linspace(minval, maxval, nx)
	X, Y = np.meshgrid(x,y)
	xc, yc = 0., 0.
	if mask == 'circle':
		if padding_zeros>1:
			radius = (nx * dx/2) - (nx/3) ###(nx/50) #changed on 20171206
		else:
			radius = (nx * dx/2) - (nx/50) ###(nx/50) #changed on 20171206
		inds=np.where((X-xc)**2. + (Y-yc)**2. <= radius**2.) #all in arcmins
	elif mask == 'square':
		radius = (nx * dx/2) - 1.
		inds=np.where((abs(X)<=radius) & (abs(Y)<=radius)) #all in arcmins

	MASKf[inds]=1.

	apodMASKf=ndimage.convolve(MASKf, ker2d)#, mode='wrap')
	apodMASKf/=apodMASKf.max()

	###imshow(apodMASKf);colorbar();show()

	return apodMASKf


# allot the arguments
parser = argparse.ArgumentParser(description='')
parser.add_argument('-cmbrandomseedval', dest='cmbrandomseedval', action='store', help='cmbrandomseedval', type=int, default=8128) #some random number
parser.add_argument('-clustype', dest='clustype', action='store', help='clustype', type=str, default='sptpol_des')
parser.add_argument('-paramfile', dest='paramfile', action='store', help='paramfile', type=str, required=True)
parser.add_argument('-prefix', dest='prefix', action='store', help='prefix', type=str, required=True)
parser.add_argument('-noiselevels', dest='noiselevels', action='store', help='noiselevels', type=str, required = True)

parser.add_argument('-extra_name', dest='extra_name', action='store', help='extra_name', type=str, default = None)# = '_y3_v6.4.21_lgt5_redmapper_clusters_vl')
parser.add_argument('-totalclus', dest='totalclus', action='store', help='totalclus', type=int)#, required=True)
parser.add_argument('-redmapper', dest='redmapper', action='store', help='redmapper', type=int, default = 1)

parser.add_argument('-just_SZfreemap', dest='just_SZfreemap', action='store', help='just_SZfreemap', type=int, default = 0)

parser.add_argument('-use_data', dest='use_data', action='store', help='use_data', type=int, default=0)
parser.add_argument('-CR_maps', dest='CR_maps', action='store', help='CR_maps', type=int, default = 1)
parser.add_argument('-no_lensing', dest='no_lensing', action='store', help='no_lensing', type=int, default=0)
parser.add_argument('-dx', dest='dx', action='store', help='dx', type=float, default=0.5)
parser.add_argument('-add_tSZ', dest='add_tSZ', action='store', help='add_tSZ', type=int, default=0)
parser.add_argument('-minSNR', dest='minSNR', action='store', help='minSNR', type=str, default = None)
parser.add_argument('-maxSNR', dest='maxSNR', action='store', help='maxSNR', type=str, default = None)
parser.add_argument('-null_test', dest='null_test', action='store', help='null_test', type=int, default = 0)
parser.add_argument('-curl_test', dest='curl_test', action='store', help='curl_test', type=int, default = 0)
parser.add_argument('-lx_for_qe', dest='lx_for_qe', action='store', help='lx_for_qe', type=int, default = 0)
parser.add_argument('-beam_error', dest='beam_error', action='store', help='beam_error', type=float, default = 0)
parser.add_argument('-mean_field', dest='mean_field', action='store', help='mean_field', type=int, default = 0)
parser.add_argument('-which_rand', dest='which_rand', action='store', help='which_rand', type=str, default = None)
parser.add_argument('-howmanyrands', dest='howmanyrands', action='store', help='howmanyrands', type=int, default = 5500)
parser.add_argument('-store_kappa_array', dest='store_kappa_array', action='store', help='store_kappa_array', type=int, default = 0)
parser.add_argument('-min_sep_from_psource', dest='min_sep_from_psource', action='store', help='min_sep_from_psource', type=float, default = 10.)
parser.add_argument('-use_lgmca', dest='use_lgmca', action='store', help='use_lgmca', type=int, default = 0)
parser.add_argument('-padding_zeros', dest='padding_zeros', action='store', help='padding_zeros', type=int, default = 1)
parser.add_argument('-add_psource', dest='add_psource', action='store', help='add_psource', type=int, default = 0)
parser.add_argument('-matched_filter', dest='matched_filter', action='store', help='matched_filter', type=int, default = 0)
parser.add_argument('-szextravariance', dest='szextravariance', action='store', help='szextravariance', type=int, default = 0)
parser.add_argument('-boxsize_extract', dest='boxsize_extract', action='store', help='boxsize_extract in am', type=float, default = 30.)
parser.add_argument('-nx_large', dest='nx_large', action='store', help='nx_large in am', type=int, default = None)
parser.add_argument('-ny_large', dest='ny_large', action='store', help='ny_large in am', type=int, default = None)
parser.add_argument('-make_sptpol_mask_homo', dest='make_sptpol_mask_homo', action='store', help='make_sptpol_mask_homo', type=int, default = 1)
parser.add_argument('-cluster_mass', dest='cluster_mass', action='store', help='cluster_mass', type=float, default = 2.8)
parser.add_argument('-fit_arnaud', dest='fit_arnaud', action='store', help='fit_arnaud', type=float, default = 0)
parser.add_argument('-fit_gaussian', dest='fit_gaussian', action='store', help='fit_gaussian', type=float, default = 0)
parser.add_argument('-sptsz_sim', dest='sptsz_sim', action='store', help='sptsz_sim', type=int, default = 0)
parser.add_argument('-sehgal_sims', dest='sehgal_sims', action='store', help='sehgal_sims', type=int, default = 0)
parser.add_argument('-takahashi_sims', dest='takahashi_sims', action='store', help='takahashi_sims', type=int, default = 0)
parser.add_argument('-arnaud_profile', dest='arnaud_profile', action='store', help='arnaud_profile', type=int, default = 0)
parser.add_argument('-arnaud_scatter', dest='arnaud_scatter', action='store', help='arnaud_scatter', type=int, default = 0)
args = parser.parse_args()
args_keys = args.__dict__
for kargs in args_keys:
        param_value = args_keys[kargs]
	if isinstance(param_value, str):
		cmd = '%s = "%s"' %(kargs, param_value)
	else:
		cmd = '%s = %s' %(kargs, param_value)
	exec(cmd)

# pass the random seed value to Simulations class
sims.cmbrandomseedval = cmbrandomseedval
if add_tSZ:
    if sehgal_sims==0 and arnaud_profile == 0 and takahashi_sims ==0:
        print "Please specify which SZ profile to add"
        sys.exit()

pref_z_for_lnlike = None
params = np.recfromtxt(paramfile,usecols=[0],delimiter = '=')
paramvals = np.recfromtxt(paramfile,usecols=[1],delimiter = '=')

# pass all the paramsfile values in param_dict (dictionary)
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

# pass the param dictionary file to Simulations class
sims.inidic = param_dict
sims.cmb_noise_randomseedval = cmbrandomseedval
# create a log file to log all the steps
timedatestr='%s' %(time.strftime('%Y%m%d_%H%M%S', time.localtime()))
log_file = 'tmp/logs_%s_%s.txt' %(timedatestr,cmbrandomseedval)
if os.path.exists(log_file):
	os.system('rm %s' %(log_file))
sims._ini_log_file(log_file)


# Use generated CAMB power spectrum
Dlfile_len = param_dict['Dlfile_len']
Dlfile_unlen = param_dict['Dlfile_unlen']
# Load Dls
if param_dict['do_pol']:
	Dls_len = np.loadtxt(Dlfile_len,usecols=[0,1,2,3,4])
	Dls_unlen = np.loadtxt(Dlfile_unlen,usecols=[0,1,2,3,4])
else:
	Dls_len = np.loadtxt(Dlfile_len,usecols=[0,1])
	Dls_unlen = np.loadtxt(Dlfile_unlen,usecols=[0,1])
tqulen = Dls_len.shape[1] - 1
sims.tqulen = tqulen


redshift = 0.7
M_200_list = np.ones(totalclus)*float(cluster_mass)*1e14
z_L_list = np.ones(totalclus)*redshift



# get beams
nx,ny,dx,dy = 200,200,0.5,0.5
simmapparams = [nx,ny, dx,dy]
Bl_150 = sims.fn_beam_stuff(simmapparams, use_beam = param_dict['use_beam'], exp_beam = param_dict['exp_beam'],nu = 150, return_beam = 1)
Bl_90 = sims.fn_beam_stuff(simmapparams, use_beam = param_dict['use_beam'], exp_beam = param_dict['exp_beam_90'],nu = 90, return_beam = 1)
Bl_gradient_map = Bl_90
# pass the beams to simulations class
sims.Bl_gradient_map = Bl_gradient_map
sims.Bl_90 = Bl_90
sims.Bl_150 = Bl_150




# add_noise variable will determine whether added noise will be a white noise or a realisation of a spectrum
if noiselevels.lower() == 'sptpol':
	add_noise = 2
	noise_folder = 'sptpol'
elif noiselevels.lower() == 'ilc_res_file':
	add_noise = 0
	noise_folder = 'ilc_res_file'
else:
	add_noise = 1
sims.noise_present = add_noise

# noiselevel variable determines the white noise level
if noiselevels.split(',')[0] == 'None' or float(noiselevels.split(',')[0]) == 0.:
    logline = '\n\t#######################################################\n\tNoise cannot be zero / None. Setting it to 0.1uK-arcmin to COV regularisation'
    logfile = open(log_file,'a')
    logfile.writelines('%s\n' %(logline))
    logfile.close()
    print logline
    noiselevels = noiselevels.replace('None','0.1')
expnoiselevel = np.asarray(noiselevels.split(',')).astype('float')
noise_folder = 'white_%s' %(expnoiselevel[0])

# extract the CMB lensing cutouts from the Takahashi SIMS
cmb_map_arr = sims.fn_extract_takahashi_cmb_maps(M_200_list, z_L_list)

# convolve CMB with beam function and add experimental noise
kappa_qe_arr,CLUS_IDENTIFIER =[],[]
res_dic = {}
for ii in range(totalclus):
    CMB = cmb_map_arr[ii]
    # subtract the mean off
    CMB = CMB - np.mean(CMB)
    CMB = np.fft.ifft2( np.fft.fft2( CMB ) * Bl_150 ).real
    SIMMAPS_L_BEAM_TF = sims.fn_tf_noise(np.asarray([CMB]), simmapparams, 1, reqdbox = None, add_TF = param_dict['add_TF'], add_noise = add_noise, expnoiselevel = expnoiselevel, cluscnt = ii, return_full_T_map = 0)[0]
    TMAP = (SIMMAPS_L_BEAM_TF[0]/1e6)
    OBSMAP = np.asarray([TMAP])

    # get the apodisation mask
    apodMASK = fn_apodize(OBSMAP[0], simmapparams, mask = 'circle', just_return_mask = 1)
    OBSMAP = OBSMAP * apodMASK
    OBSMAP = OBSMAP - np.mean(OBSMAP)
    if param_dict['cross_maps']:
        TMAP2 = np.asarray([(GRAD_SIMMAPS_L_BEAM_TF[0]/1e6)])
        OBSMAP2 = TMAP2 - np.mean(TMAP2)
    else:
        OBSMAP2 = None
    #hardcoded
    richval = None; noise_weight = 1;sims.just_SZfreemap = 0;sims.nx_large = None; sims.ny_large = None
    sims.matched_filter = 0; sims.sehgal_fg_added = 0
    KAPPA_QE = sims.fn_get_kappa_QE(OBSMAP, simmapparams, Dls_len, Dls_unlen, tszfree = 1, OBSMAP2 = OBSMAP2, richval = richval, use_data = use_data, curl_test = curl_test, noise_weight = noise_weight)
    # extract only boxsize_extract arcminutes to save space
    nx_small, nx_large  = int(nx/2 -  boxsize_extract), int(nx/2 + boxsize_extract)
    ny_small, ny_large  = int(ny/2 -  boxsize_extract), int(ny/2 + boxsize_extract)
    KAPPA_QE = KAPPA_QE.real - np.mean(KAPPA_QE.real)
    kappa_qe_arr.append(KAPPA_QE[nx_small:nx_large,ny_small:ny_large])
    CLUS_IDENTIFIER.append([0.0, 0.0, 0.69999999999999996, 100.0, 1.0])

opfolder = 'data/Takahashi_sims/cluster_mass_%s'%(cluster_mass)
opfolder_split = opfolder.split('/')
tmp = ''
for ff in opfolder_split:
	tmp = '%s/%s' %(tmp,ff);tmp = tmp.strip('/')
	if not os.path.exists(tmp):
		os.system('mkdir %s' %(tmp))
opfile = '%s/stacked_kappa_%s_%05d_clusters_%s.pkl.gz'%(opfolder,timedatestr,totalclus, cmbrandomseedval)
res_dic['CLUS_IDENTIFIER'] =CLUS_IDENTIFIER
res_dic['param_dict'] = param_dict
res_dic['kappa_qe_arr'] = kappa_qe_arr
stacked_kappa_qe = np.mean(kappa_qe_arr, axis = 0)
res_dic['stacked_kappa_qe'] = stacked_kappa_qe

pickle.dump(res_dic,gzip.open(opfile,'w'))
print opfile
