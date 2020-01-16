################################################################################################
#Comments: This is an update to helper_funcs/qe_funcs.py to include the polarisation estimators
#Original Author: Eric Baxter (EB)
#Modified by: Srini Raghunathan (SR)
#Date start/finish: 20161107 / 201611xx

#new functions:
#1. fn_get_filters to get the gradient and inv. variance filters based on the estimator
#2,3. fn_get_EB_betas, fn_get_EB_gammas for the kappa normalization. Check Eric Baxter's notes on implementing QE.

#testing:
#tested TT estimator and the results match the original qe_funcs.py module

#status:
#work in progress

#pending:
#Extra multiplication factors for G_l^XY, L_l^Y fields
################################################################################################
import numpy as np
import pdb, sys
from pylab import *

#import helper_funcs.fourier_funcs as ff

#Get angle for specified vectors
def get_angle(x,y):
#return angle in range [-pi,pi)
	tempx = np.atleast_1d(x)
	tempy = np.atleast_1d(y)
	angle = np.zeros(tempx.shape)
	zero_x_zeroy = np.where((tempx == 0.) & (tempy == 0.))
	zero_x_posy  = np.where((tempx == 0.) & (tempy > 0.))
	zero_x_negy  = np.where((tempx == 0.) & (tempy < 0.))
	neg_x_zeroy = np.where((tempx < 0.) & (tempy == 0.))
	
	nonzero_x = np.where((tempx != 0.))
	negx_negy = np.where((tempx < 0) & (tempy < 0))
	negx_posy = np.where((tempx < 0) & (tempy > 0))
	
	angle[nonzero_x] = np.arctan(tempy[nonzero_x]/tempx[nonzero_x])
	angle[negx_posy] = np.pi+angle[negx_posy]
	angle[negx_negy] = -np.pi+angle[negx_negy]

	angle[zero_x_zeroy] = 0.
	angle[zero_x_posy] = np.pi/2.
	angle[zero_x_negy] = -np.pi/2.
	angle[neg_x_zeroy] = -np.pi
	return angle

#get QE of kappa
def get_kappa(lx, ly, beamtf_lmat, weight_gradient_lmat, weight_inv_var_lmat, T_mat, norm_lmat, dxdy, ngaus, estimator, T_mat2 = None):

	# Following steps in 1411.7999

	if T_mat2 == None: T_mat2 = np.copy(T_mat)

	'''for i in range(3):
		subplot(1,3,i+1);imshow(T_mat[i]);colorbar()
	show();sys.exit()'''

	#Useful later
	l2d = np.sqrt(lx**2. + ly**2.)

	# Subtract mean
	obs_T_mat_meansub = np.asarray([T_mat[i] - np.mean(T_mat[i]) for i in range(len(T_mat))])
	obs_T_mat_meansub2 = np.asarray([T_mat2[i] - np.mean(T_mat2[i]) for i in range(len(T_mat2))])
	
	#from pylab import *
	#imshow(obs_T_mat_meansub);colorbar();show();quit()

	#For doing deconvolution
	to_deconvolve_lmat = np.copy(beamtf_lmat)
	bad = np.where(to_deconvolve_lmat==0.)
	to_deconvolve_lmat[bad] = 1.
	deconvolve_lmat = 1./to_deconvolve_lmat
	deconvolve_lmat[bad] = 0.

	highly_filtered = np.where((deconvolve_lmat > 1.0e10) | (deconvolve_lmat < 0) )
	deconvolve_lmat[highly_filtered] = 0.0

	# Do filtering and deconvolution
	T_lmat = dxdy*np.fft.fft2(obs_T_mat_meansub)*deconvolve_lmat
	T_lmat2 = dxdy*np.fft.fft2(obs_T_mat_meansub2)*deconvolve_lmat

	loc_dic = {'T':0, 'E': 1, 'B': 2}
	if 1==1:#len(T_mat.shape)>2: #check 20171117
		T_lmat2 = T_lmat2[loc_dic[estimator[0]]] #gradient map
		T_lmat = T_lmat[loc_dic[estimator[1]]]
	else:
		T_lmat2 = T_lmat

	'''
	#Eq. 14 of Hu 2007
	if estimator[1] == 'E' or estimator[1] == 'B':
		psi = 2. * get_angle(lx,ly)
		X_lmat = X_lmat * ( np.cos(psi) + 1j * np.sin(psi) )
		Y_lmat = Y_lmat * ( np.cos(psi) + 1j * np.sin(psi) )
	'''
	if estimator[1] == 'E' or estimator[1] == 'B':
		#psi = 2. * get_angle(lx,ly)
		psi = 2*np.arctan2(-lx, ly)
		T_lmat2 = T_lmat2 * ( np.cos(psi) + 1j * np.sin(psi) )
		T_lmat = T_lmat * ( np.cos(psi) + 1j * np.sin(psi) )

	#from pylab import *
	#imshow(np.fft.ifftshift(np.fft.ifft2(T_lmat2).real));colorbar();show();quit()
	
	#Filtered gradient field
	filtered_gradient_x_lmat = 1.j*lx*weight_gradient_lmat*T_lmat2
	filtered_gradient_y_lmat = 1.j*ly*weight_gradient_lmat*T_lmat2


	#W^T = weight_inv_var_l
	T_inv_var_lmat = weight_inv_var_lmat*T_lmat

	#imshow(T_inv_var_lmat.real);colorbar();show();sys.exit()

	#Convert to configuration space
	filtered_gradient_x_mat = (1./dxdy)*np.fft.ifft2(filtered_gradient_x_lmat)
	filtered_gradient_y_mat = (1./dxdy)*np.fft.ifft2(filtered_gradient_y_lmat)
	T_inv_var_mat = (1./dxdy)*np.fft.ifft2(T_inv_var_lmat)

	#imshow(T_inv_var_mat.real);colorbar();show();sys.exit()

	#Multiply filtered fields in configuration space and fourier transform
	GL_x_lmat = dxdy*np.fft.fft2(filtered_gradient_x_mat*T_inv_var_mat)
	GL_y_lmat = dxdy*np.fft.fft2(filtered_gradient_y_mat*T_inv_var_mat)

	#Filtered divergence, and applying normalization
	temp_norm_lmat = np.copy(norm_lmat)
	bad = np.where(temp_norm_lmat == 0.0)
	temp_norm_lmat[bad] = 1.0
	ampval = 1./temp_norm_lmat

	#WHAT DOES ngaus DO???77777777
	#filter2 = np.exp(-0.5*(l2d/float(ngaus))**2)
	filter2 = np.ones(l2d.shape) 
	above_beam_scale = np.where(l2d>=ngaus)
	filter2[above_beam_scale] = 0.

	teye = -ampval*l2d*1.j*(lx*GL_x_lmat + ly*GL_y_lmat)
	kappa_lmat = filter2*teye*l2d/2.
	kappa_lmat[bad] = 0.0
	kappa_mat = (1./dxdy)*np.fft.ifft2(kappa_lmat)

	#imshow(np.fft.fftshift(norm_lmat.real), extent = [np.min(lx), np.max(lx), np.min(ly), np.max(ly)]);colorbar();show();quit()

	if (0):
		import matplotlib.pyplot as pl
		pdb.set_trace()
	
		# 777777777 testing
		threshold = 0.5
		filtered_by_beam_and_tf = np.where(((to_deconvolve_lmat) < threshold*np.max(to_deconvolve_lmat)))
		plot_filter = 1.0 + np.zeros(l2d.shape)
		plot_filter[filtered_by_beam_and_tf] = 0.0
		kappa_lmat_filtered = kappa_lmat*plot_filter
		kappa_filtered = (1./dxdy)*np.fft.ifft2(kappa_lmat_filtered)
		

		fig, ax = pl.subplots(1,1)
		ax.imshow(np.real(kappa_filtered))

		fig, ax = pl.subplots(1,1)
		ax.imshow(np.real(filtered_gradient_x_mat*T_inv_var_mat))

		fig, ax = pl.subplots(1,1)
		ax.imshow(np.real(teye))
		
		# 777777 for testing
		T_mat_deconvolve = (1./dxdy)*np.fft.ifft2(T_lmat)


		figt, axt = pl.subplots(1,1)
		axt.imshow(np.real(T_mat_deconvolve))

		fig_all, ax_all = pl.subplots(1,3)
		ax_all[0].imshow(np.log(np.abs(dxdy*np.fft.fft2(obs_T_mat_meansub))))
		ax_all[1].imshow(np.log(np.abs(dxdy*np.fft.fft2(obs_T_mat_meansub*window_func))))
		ax_all[2].imshow(np.log(np.abs(dxdy*np.fft.fft2(obs_T_mat_meansub*window_func)*deconvolve_lmat)))

		fig_real, ax_real = pl.subplots(1,3)
		ax_real[0].imshow(obs_T_mat_meansub)
		ax_real[1].imshow(obs_T_mat_apod)
		ax_real[2].imshow(np.real(T_mat_deconvolve))
		
		fig_test, ax_test = pl.subplots(1,3)
		ax_test[0].imshow(np.real(filtered_gradient_x_mat))
		ax_test[1].imshow(np.log(np.abs(T_inv_var_mat)))
		ax_test[2].imshow(np.log(np.real(filtered_gradient_x_mat*T_inv_var_mat)))
		
		fig_beam, ax_beam = pl.subplots(1,4)
		ax_beam[0].imshow(beam_lmat)
		ax_beam[1].imshow(transfer_lmat)
		ax_beam[2].imshow(transfer_lmat*beam_lmat)
		ax_beam[3].imshow(to_deconvolve_lmat)
		pdb.set_trace()
		
	return kappa_mat


def fn_get_filters(l2d, C_l_len_lmat_dic, C_l_unl_lmat_dic, N_l_lmat_dic, estimator, l_G = 1800, lmax_inv_var = 15000):

	######################################
	# weight_gradient_lmat (gradient filter)
	######################################
	#unlensed Cls (Eq. 21 of Hu)
	if estimator[1] == 'T' or estimator[1] == 'E':
		if estimator in C_l_unl_lmat_dic.keys(): #TT, EE, TE
			field_for_unlensedcls = estimator
		else: #ET		
			field_for_unlensedcls = estimator[::-1] #just reverse the estimator

	else: #any of the "B" lensing estimator
		field_for_unlensedcls = '%sE' %(estimator[0])

	C_l_unl_lmat = C_l_unl_lmat_dic[field_for_unlensedcls]

	#lensed Cls / Nls (Eq. 21/22 of Hu)
	field_for_lensedcls = '%s%s' %(estimator[0],estimator[0])
	C_l_len_lmat = C_l_len_lmat_dic[field_for_lensedcls]
	N_l_lmat = N_l_lmat_dic[field_for_lensedcls]

	weight_gradient_lmat = C_l_unl_lmat / (C_l_len_lmat + N_l_lmat)

	above_lg = np.where(l2d > l_G)
	weight_gradient_lmat[above_lg] = 0.0

	######################################
	# weight_inv_var_lmat (lensing filter)
	######################################
	field_for_lensedcls = '%s%s' %(estimator[1],estimator[1])
	C_l_len_lmat = C_l_len_lmat_dic[field_for_lensedcls]
	N_l_lmat = N_l_lmat_dic[field_for_lensedcls]
	
	weight_inv_var_lmat = 1. / (C_l_len_lmat + N_l_lmat)

	above_inv_var_lmax = np.where(l2d > lmax_inv_var)
	weight_inv_var_lmat[above_inv_var_lmax] = 0.0

	#subplot(121);imshow(weight_gradient_lmat);colorbar()
	#subplot(122);imshow(weight_inv_var_lmat);colorbar();show()
	#sys.exit()

	return weight_gradient_lmat, weight_inv_var_lmat

#Calculate QE normalization

def get_lxly(mapparams):

	nx, ny, dx, dy = mapparams
	dx = np.radians(dx/60.)
	dy = np.radians(dy/60.)

	lx, ly = np.meshgrid( np.fft.fftfreq( nx, dx ), np.fft.fftfreq( ny, dy ) )
	lx *= 2* np.pi
	ly *= 2* np.pi

	return lx, ly

def fn_get_EB_gammas(l2d, theta_lmat, C_l_unl_lmat_dic, estimator):

	nlx, nly = l2d.shape
	psi = 2*theta_lmat

	"""
	dx = dy = 0.5 #arcmins
	mapparams = [nlx, nly, dx, dy]
	lx, ly = get_lxly(mapparams)
	psi = 2.*np.arctan2(lx, -ly)
	"""
	

	#if estimator == 'EB' or 'TB':
	if estimator[1] == 'B': #EB, TB
		if estimator == 'EB':
			C_l_unl_lmat = C_l_unl_lmat_dic['EE']
		elif estimator == 'TB':
			C_l_unl_lmat = C_l_unl_lmat_dic['TE']

		R_gamma = np.zeros((2,nlx,nly), dtype = np.cfloat)
		P_gamma = np.zeros((2,nlx,nly), dtype = np.cfloat)
		Q_gamma = np.zeros((2,nlx,nly), dtype = np.cfloat)

		R_gamma[0,:,:] = 0.5*l2d*np.exp(1.j*theta_lmat) 
		R_gamma[1,:,:] = 0.5*l2d*np.exp(-1.j*theta_lmat)
		P_gamma[0,:,:] = l2d*np.exp(-1.j*theta_lmat)
		P_gamma[1,:,:] = l2d*np.exp(1.j*theta_lmat)
		Q_gamma[0,:,:] = C_l_unl_lmat* np.sin(psi)
		Q_gamma[1,:,:] = C_l_unl_lmat* np.sin(psi)

	#if estimator == 'TT' or estimator == 'EE':
	else: #TT, TE, EE

		if estimator in C_l_unl_lmat_dic.keys(): #TT, EE, TE
			field_for_unlensedcls = estimator
		else: #ET		
			field_for_unlensedcls = estimator[::-1] #just reverse the estimator

		C_l_unl_lmat = C_l_unl_lmat_dic[field_for_unlensedcls] #TT or EE

		if estimator == 'EE':
			extra_gamma_term = np.cos(psi)
		else: #TT, TE, ET
			extra_gamma_term = 1.

		#extra_gamma_term = 1.

		R_gamma = np.zeros((5,nlx,nly), dtype = np.cfloat)
		P_gamma = np.zeros((5,nlx,nly), dtype = np.cfloat)
		Q_gamma = np.zeros((5,nlx,nly), dtype = np.cfloat)

		R_gamma[0,:,:] =  0.5*l2d*np.exp(1.j*theta_lmat)
		R_gamma[1,:,:] =  0.5*l2d*np.exp(-1.j*theta_lmat)
		R_gamma[2,:,:] = -0.5*l2d*np.exp(1.j*theta_lmat)
		R_gamma[3,:,:] = -0.5*l2d*np.exp(-1.j*theta_lmat)
		R_gamma[4,:,:] = l2d**2.

		P_gamma[0,:,:] = l2d*np.exp(-1.j*theta_lmat)*C_l_unl_lmat * extra_gamma_term
		P_gamma[1,:,:] = l2d*np.exp(1.j*theta_lmat)*C_l_unl_lmat * extra_gamma_term
		P_gamma[2,:,:] = l2d*np.exp(-1.j*theta_lmat) * extra_gamma_term
		P_gamma[3,:,:] = l2d*np.exp(1.j*theta_lmat) * extra_gamma_term
		P_gamma[4,:,:] = 1.0* extra_gamma_term

		Q_gamma[0,:,:] = 1.0 * extra_gamma_term
		Q_gamma[1,:,:] = 1.0 * extra_gamma_term
		#if estimator == 'TE'  or estimator == 'ET': #replace extra_Q_gamma_term now by 1 as the last three terms are not multiplied by extra_Q_gamma_term for TE
		#	extra_Q_gamma_term = 1. 
		Q_gamma[2,:,:] = C_l_unl_lmat * extra_gamma_term
		Q_gamma[3,:,:] = C_l_unl_lmat * extra_gamma_term
		Q_gamma[4,:,:] = C_l_unl_lmat * extra_gamma_term


	return R_gamma, P_gamma, Q_gamma

def fn_get_EB_betas(l2d, theta_lmat, weight_gradient_lmat, weight_inv_var_lmat, estimator):

	nlx, nly = l2d.shape
	psi = 2.*theta_lmat #np.arctan2(-lx, ly)


	R_beta = np.zeros((1,nlx,nly), dtype = np.cfloat)
	P_beta = np.zeros((1,nlx,nly), dtype = np.cfloat)
	Q_beta = np.zeros((1,nlx,nly), dtype = np.cfloat)

	if estimator[1] == 'T':
		extra_beta_term = 1.
	elif estimator[1] == 'E':
		extra_beta_term = np.cos(psi)
	elif estimator[1] == 'B':
		extra_beta_term = np.sin(psi)

	#extra_beta_term = 1.

	R_beta[:,:] = 1.0
	P_beta[:,:] = weight_gradient_lmat#* extra_beta_term
	Q_beta[:,:] = weight_inv_var_lmat#* extra_beta_term

	return R_beta, P_beta, Q_beta

def get_norm(lx, ly, C_l_unl_lmat_dic, weight_gradient_lmat, weight_inv_var_lmat, dxdy, estimator, mapparams):

	l2d = np.sqrt(lx**2. + ly**2.)
	theta_lmat = np.arctan2(-lx, ly)


	if estimator in C_l_unl_lmat_dic.keys(): #TT, EE, TE
		field_for_unlensedcls = estimator
	else: #ET		
		field_for_unlensedcls = estimator[::-1] #just reverse the estimator

	C_l_unl_lmat = C_l_unl_lmat_dic[field_for_unlensedcls] #TT or EE

	#lx, ly = np.unique(lx.flatten()), np.unique(ly.flatten())
	L = np.max(l2d) #maxl to go
	norm = np.zeros( (l2d.shape) )


	#for y1c, l1y in 


	"""
	
	for l1cnt, l1 in np.ndenumerate(l2d):
		l2 = L - l1
		l2cnt = np.where(l2d == l2)
		psl1 = theta_lmat[l1cnt]
		psl2 = theta_lmat[l2cnt]

		print l1cnt, l2cnt

		'''
		t2 = weight_gradient_lmat[lcnt]
		t3 = np.interp(l2, np.unique( l2d.flatten() ), np.unique( C_l_unl_lmat.flatten() ))
		dummy = np.interp(l1, np.unique( l2d.flatten() ), np.unique( C_l_unl_lmat.flatten() ))
		#plot(l2d.flatten(), theta_lmat.flatten())
		#show();sys.exit()
		print dummy, C_l_unl_lmat[lcnt]
		'''
	"""

	sys.exit()

	
def get_norm_v1(lx, ly, C_l_unl_lmat_dic, weight_gradient_lmat, weight_inv_var_lmat, dxdy, estimator, mapparams):

	"""
	updating Eric Baxter's code to include polarisation estimators
	additional parameters:
	C_l_unl_lmat replaced by C_l_unl_lmat_dic that contains the C_l_unl_lmat for all possible xcorr (TT,EE,TE)
	estimator = TT, EE, ET, EB, etc.
	"""

	l2d = np.sqrt(lx**2. + ly**2.)

	#compute normalization
	'''#20171122 - changing this to below definition for angles
	#theta_lmat = get_angle(lx,ly)
	'''
	nlx = lx.shape[0]
	nly = lx.shape[1]
	theta_lmat = np.arctan2(-lx, ly)

	R_alpha = np.zeros((2,nlx,nly), dtype = np.cfloat)
	P_alpha = np.zeros((2,nlx,nly), dtype = np.cfloat)
	Q_alpha = np.zeros((2,nlx,nly), dtype = np.cfloat)

	R_alpha[0,:,:] = 0.5*l2d*np.exp(1.j*theta_lmat)
	R_alpha[1,:,:] = 0.5*l2d*np.exp(-1.j*theta_lmat)
	P_alpha[0,:,:] = l2d*np.exp(-1.j*theta_lmat)
	P_alpha[1,:,:] = l2d*np.exp(1.j*theta_lmat)
	Q_alpha[0,:,:] = 1.0
	Q_alpha[1,:,:] = 1.0
	
	### SR: 20161107 - modifying R, P, and Q betas/gammas to include the extra fields for polarisation estimators
	R_beta, P_beta, Q_beta = fn_get_EB_betas(l2d, theta_lmat, weight_gradient_lmat, weight_inv_var_lmat, estimator)
	R_gamma, P_gamma, Q_gamma = fn_get_EB_gammas(l2d, theta_lmat, C_l_unl_lmat_dic, estimator)
	
	alpha_len, beta_len, gamma_len = len(R_alpha), len(R_beta), len(R_gamma)
	#print alpha_len, beta_len, gamma_len;sys.exit()
	rpq_sum = np.zeros((nlx, nly), dtype = np.cfloat)
	for ii in xrange(0,alpha_len):
		for kk in xrange(0,beta_len):
			for jj in xrange(0,gamma_len):
				R_tot = R_alpha[ii,:,:]*R_beta[kk,:,:]*R_gamma[jj,:,:]
				P_tot = P_alpha[ii,:,:]*P_beta[kk,:,:]*P_gamma[jj,:,:]
				Q_tot = Q_alpha[ii,:,:]*Q_beta[kk,:,:]*Q_gamma[jj,:,:]
				Px = (1./dxdy)*np.fft.ifft2(P_tot)
				Qx = (1./dxdy)*np.fft.ifft2(Q_tot)
				rpq_sum += (R_tot)*dxdy*np.fft.fft2(Px*Qx)

	norm = rpq_sum

	#imshow(norm.real);colorbar();title(estimator);show();sys.exit()

	return norm

