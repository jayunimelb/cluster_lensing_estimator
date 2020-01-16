def fn_ang_dist_fast(ip_x1, ip_y1, ip_x2, ip_y2): #all numpy arrays

        import numpy as np

        if isinstance(ip_x1,str):
                ip_x1=float(ip_x1)
                ip_y1=float(ip_y1)
                ip_x2=float(ip_x2)
                ip_y2=float(ip_y2)


        # CONVERT INTO RADIANS
        x1 = np.radians(ip_x1)
        y1 = np.radians(ip_y1)
        x2 = np.radians(ip_x2)
        y2 = np.radians(ip_y2)

        ang_dist_rad = np.arccos(np.sin(y1)*np.sin(y2)+(np.cos(y1)*np.cos(y2)*np.cos(x1 - x2)))
        ang_dist_deg = np.degrees(ang_dist_rad) # Converting radians to degrees

        return ang_dist_deg

def fn_get_radec_M_z_cl(clusatcen = 1, noofclus = 1):
	clus_ra, clus_dec = [], []
	M_200_cl, z_L_cl = [], []
	c_200_cl = []
	for ii in range(noofclus):
		M_val = (np.random.random() * 10)
		while M_val>5. or M_val<.5:
			M_val = (np.random.random() * 10)
		M_val *= 1e14

		z_val =  z_L_list[np.random.randint(len(z_L_list))]


		if noofclus == 1 and clusatcen == 1:
			clus_ra_val, clus_dec_val = 0., 0.
		else:
			clus_ra_val = np.float(raarr[np.random.randint(edge_protection,len(raarr)-edge_protection)])
			clus_dec_val = np.float(decarr[np.random.randint(edge_protection,len(decarr)-edge_protection)])

		#if testing:
		if 1==1:
			M_val, z_val = 2e14, 0.6

		z_L_cl.append(z_val)
		M_200_cl.append(M_val)
		clus_ra.append(clus_ra_val)
		clus_dec.append(clus_dec_val)

		c_200_cl.append(3.)

	return clus_ra, clus_dec, M_200_cl, z_L_cl, c_200_cl

################################################################################################
################################################################################################
################################################################################################

import numpy as np
import sys, pickle, gzip, os, time
####
import modules
from modules import scl_cmb
from pylab import *

sims = modules.scl_cmb.simulations()

timedatestr='%s' %(time.strftime('%Y%m%d_%H%M%S', time.localtime()))
log_file = 'tmp/logs_%s.txt' %(timedatestr)
if not os.path.exists(log_file):
        os.system('rm %s' %(log_file))
sims._ini_log_file(log_file)

testing = 0

####
#noofsims = int(sys.argv[1])
totalclus = int(sys.argv[1])

####
beamfwhmarcmins = float(sys.argv[2])
noiselevels_T_arr = [float(sys.argv[3])]

####

isolated_clus = int(sys.argv[4])
max_noniso_clus = 3 #maximum of 3 clusters in a 100 arcmin box

####
#dx = float(sys.argv[5])
degraded_dx = float(sys.argv[5])
### dx = 0.05
dx = degraded_dx

####
Dlfile_len = 'data/fd.cllens'
Dlfile_unlen = 'data/fd.cl'
CMB_outputscale = 1
Dls_len = np.loadtxt(Dlfile_len,usecols=[0,1])
Dls_unlen = np.loadtxt(Dlfile_unlen,usecols=[0,1])
#Dls_unlen = Dls_len


###
param_dict = {}
param_dict['h'] = 0.7
param_dict['omega_m'] = 0.3
param_dict['z_lss'] = 1089.
param_dict['T_cmb'] = 2.73

c_200_list = 3.0
mass_def = 200.0
rho_def = 'crit'
theta_max = -1.

## create unlesned sims
boxsize = 40
#dx, dy = .25, .25 #arcmins
dy = dx
####

nx, ny = int(boxsize/dx), int(boxsize/dy)

boxsize = int(nx * dx)

simmapparams = [nx, ny, dx, dy]

minval, maxval = -boxsize/2/60.,  boxsize/2/60.
raarr = np.linspace(minval, maxval, nx)
decarr = np.linspace(minval, maxval, ny)
RA, DEC = np.meshgrid(raarr,decarr)

reqd_box = 6. #arcmins
reqd_box_pix = reqd_box/dx

edge_protection = reqd_box_pix

####

opfolder = 'detections/boxsize_%s_%samresol' %(boxsize,dx)


####

z_del = 0.01
minz, maxz = 0.3, 1.
z_L_list = np.arange(minz,maxz,z_del)

"""
SIMMAPS = sims.Dls2map(Dls_len, simmapparams, nosims = 1, CMB_outputscale = CMB_outputscale)
#SIMMAPS *= 1e6
SIMMAPS_BEAM_TF = sims.fn_beam_tf(SIMMAPS, simmapparams, 1, reqdbox = None, use_beam = usebeam, exp_tmat=None, exp_noise_level=None)

radecdic = {}
clus_ra, clus_dec, M_200_cl, z_L_cl, c_200_cl = [], [], [], [], []
for bb in range(totalclus):

	start = time.time()
	clus_ra_val = np.float(raarr[np.random.randint(edge_protection,len(raarr)-edge_protection)])
	clus_dec_val = np.float(decarr[np.random.randint(edge_protection,len(decarr)-edge_protection)])

	#clus_ra_val,clus_dec_val = 10.,10.

	clus_ra.append(clus_ra_val)
	clus_dec.append(clus_dec_val)

	M_val = (np.random.random() * 10.) * 1e14
	M_200_cl.append(M_val)

	z_val =  z_L_list[np.random.randint(len(z_L_list))]
	z_L_cl.append(z_val)

	print bb,clus_ra_val,clus_dec_val,M_val,z_val

	c_val = 3.0
	c_200_cl.append(c_val)

	radecdic[(bb,clus_ra_val,clus_dec_val)] = [M_val,z_val,c_val]


del_noise = 1.
noiselevels_T_arr = np.arange(1.5,5,del_noise)

#noiselevels_T_arr = [1., 1., 3., 5.]
#noiselevels_T_arr = [1e-5, 1., 3., 5.]
#noiselevels_T_arr = [1e-10, 1., 3., 5.]
#noiselevels_T_arr = [1e-3,1.,2.5,4.]
"""
OPDIC = {}
for nn, noiselevel_T in enumerate(noiselevels_T_arr):

	print '\nNoise = %s uk-arcmin' %(noiselevel_T)

	STACKED = np.zeros((nx,ny))#, dtype='complex')

	for bb in range(totalclus):

		'''
		clus_ra=[.0]
		clus_dec=[.0]

		#M_val = (np.random.randint(1,5e3)/1e3) * 1e14 #limiting mass to 5e13 - 5e14 solar mass
		M_val = (np.random.random() * 10)
		while M_val>5. or M_val<.5:
			M_val = (np.random.random() * 10)
					
		M_val *= 1e14
		M_200_cl=[M_val]

		z_val =  z_L_list[np.random.randint(len(z_L_list))]
		z_L_cl=[z_val]
		'''
		clus_ra, clus_dec, M_200_cl, z_L_cl, c_200_cl = fn_get_radec_M_z_cl(clusatcen = 1, noofclus = 1) #get one cluster at the centre of the field

		logline = '\n\t\t########################################\n\t\tCluster number: %s; Mass = %s, z = %s' %(bb,M_200_cl[0],z_L_cl[0])
		logfile = open(log_file,'a')
		logfile.writelines('%s\n' %(logline))
		logfile.close()
		print logline

		#CMB sims
		SIMMAPS = sims.Dls2map(Dls_len, simmapparams, nosims = 1, CMB_outputscale = CMB_outputscale)#, random_seed = 123)

		SIMMAPS[0,0] = SIMMAPS[0,0] - np.mean(SIMMAPS[0,0])
		imshow(SIMMAPS[0,0] * 1e6);colorbar();show();quit()
		print SIMMAPS.shape;quit()
		if testing:
			SIMMAPS_BEAM_TF = sims.fn_beam_tf(SIMMAPS, simmapparams, beamfwhmarcmins, 1, reqdbox = None, use_beam = 1, add_TF=1, exp_noise_level=noiselevel_T)

		#perform lensing
		SIMMAPS_L, KAPPA = sims.fn_perform_lensing_nfw(np.copy(SIMMAPS), simmapparams, param_dict, RA, DEC, clus_ra, clus_dec, M_200_cl, c_200_cl, z_L_cl, param_dict['z_lss'], mass_def, rho_def, return_kappa = 1)#, smooth_kappa_resol = degraded_dx)
		#SIMMAPS_L, KAPPA = sims.fn_perform_lensing_nfw_finer(np.copy(SIMMAPS), simmapparams, param_dict, RA, DEC, clus_ra, clus_dec, M_200_cl, c_200_cl, z_L_cl, param_dict['z_lss'], mass_def, rho_def, return_kappa = 1, truncate_kappa_at_radius=None)

		if bb == 0 and not testing:
			lx, ly = sims.get_lxly(simmapparams)
			l2d = np.sqrt(lx**2. + ly**2.)

			fwhm = np.radians(1./60.)
			filter2 = sims.gauss_beam(fwhm, l2d)[0]
			KAPPA_FFT = np.fft.fft2(KAPPA)
			KAPPA_FFT = KAPPA_FFT * filter2
			KAPPA = np.fft.ifft2(KAPPA_FFT).real
			op_file = '%s/KAPPA_m_%s_z_%s.pkl.gz' %(opfolder,M_200_cl[0],z_L_cl[0])
			pickle.dump(KAPPA,gzip.open(op_file,'w'))
			#quit()

#		if 0 == 1:
		if not isolated_clus: #then add more clusters in the field
			moreclus = np.random.randint(max_noniso_clus)
			logline = '\t\tIncluding %s more clusters in the box' %(moreclus)
			logfile = open(log_file,'a')
			logfile.writelines('%s\n' %(logline))
			logfile.close()
			print logline

			clus_ra, clus_dec, M_200_cl, z_L_cl, c_200_cl = fn_get_radec_M_z_cl(clusatcen = 0, noofclus = moreclus)
			SIMMAPS_L, KAPPA_2 = sims.fn_perform_lensing_nfw(SIMMAPS_L, simmapparams, param_dict, RA, DEC, clus_ra, clus_dec, M_200_cl, c_200_cl, z_L_cl, param_dict['z_lss'], mass_def, rho_def, return_kappa = 1)
			KAPPA = KAPPA + KAPPA_2


		#convolve beam and exp. noise level
		SIMMAPS_L_BEAM_TF = sims.fn_beam_tf(SIMMAPS_L, simmapparams, beamfwhmarcmins, 1, reqdbox = None, use_beam = 1, add_TF=1, exp_noise_level=noiselevel_T)

		'''
		if testing:
			#display plots
			subplot(221);imshow(SIMMAPS[0,0],origin='lower');colorbar();grid(True,ls='solid')
			subplot(222);imshow(SIMMAPS_L[0,0],origin='lower');colorbar();grid(True,ls='solid')			
			subplot(223);css=imshow(np.copy(KAPPA)/max(KAPPA.ravel()),origin='lower');css.set_clim(-.1, .3);colorbar();grid(True,ls='solid')
			#subplot(224);imshow(SIMMAPS[0,0] - SIMMAPS_L[0,0],origin='lower');colorbar();grid(True,ls='solid')
			subplot(224);imshow(SIMMAPS_BEAM_TF[0,0] - SIMMAPS_L_BEAM_TF[0,0],origin='lower');colorbar();grid(True,ls='solid')
			show()#;quit()
		'''

		#get KAPPA_QE
		KAPPA_QE = sims.fn_get_kappa_QE(SIMMAPS_L[0,0], simmapparams, Dls_len, Dls_unlen, noiselevel_T)

		STACKED += np.copy(KAPPA_QE.real)

	STACKED/=totalclus

	keyname = '(%s,%s,%s)' %(beamfwhmarcmins,noiselevel_T,totalclus)
	OPDIC[keyname] = STACKED

	if testing and totalclus<100: #otherwise python crashes for some reason
		reqd_box = 10. #arcmins
		reqd_box_pix = reqd_box/dx
		x1,x2 = int(len(STACKED)/2-reqd_box_pix/2), int(len(STACKED)/2+reqd_box_pix/2)
		y1,y2 = int(len(STACKED)/2-reqd_box_pix/2), int(len(STACKED)/2+reqd_box_pix/2)
		imshow(STACKED[y1:y2,x1:x2]);colorbar();title(keyname);show();quit()
	

	"""
	subplot(111);
	css=imshow(STACKED,origin='lower')#;css.set_clim(0.,.1)
	grid(True,ls='solid',lw=0.3)
	#plot(STACKED[30,:],'r')
	colorbar()
	title('Noise = %s uk-arcmin; Clusters = %s' %(noiselevel_T, totalclus), fontsize = 12)
	show();quit()
	"""
	
	if not isolated_clus:
		op_file = '%s/non_iso_stacked_%s_arcmins_%s_noise_%s_clusters.pkl.gz_%s' %(opfolder,beamfwhmarcmins,noiselevel_T,totalclus,timedatestr)
	else:
		op_file = '%s/stacked_%s_arcmins_%s_noise_%s_clusters.pkl.gz_%s' %(opfolder,beamfwhmarcmins,noiselevel_T,totalclus,timedatestr)
	pickle.dump(OPDIC,gzip.open(op_file,'w'))

print '\nDone ..\n'


