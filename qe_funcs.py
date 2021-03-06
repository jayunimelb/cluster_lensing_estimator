import numpy as np, modules
import pdb
from pylab import *

#import helper_funcs.fourier_funcs as ff
sims = modules.scl_cmb.simulations()

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
#def get_kappa(lx, ly, beamtf_lmat, weight_gradient_lmat, weight_inv_var_lmat, T_mat, norm_lmat, dxdy, ngaus):
def get_kappa(lx, ly, beamtf_lmat, weight_gradient_lmat, weight_inv_var_lmat, T_mat, norm_lmat, dxdy, ngaus, T_mat2 = None, curl_test = 0, beam_lmat_transfer_lmat_grad = None, l_HPF = 100, unlensed_CMB = None):

	if T_mat2 is None: T_mat2 = np.copy(T_mat)

	#print T_mat, T_mat2, 'hi', curl_test; sys.exit()

	# Following steps in 1411.7999

	#Useful later
	l2d = np.sqrt(lx**2. + ly**2.)

	if (0):###l_HPF<>None: #remove long wavelength modes

		# Subtract mean
		T_mat -= np.mean(T_mat)
		T_mat2 -= np.mean(T_mat2)


		#HPF gradient to remove low-ell modes that can affect mean subtraction
		lowelHPF = np.ones(l2d.shape)
		lowelHPF[l2d<=l_HPF] = 0.
		lowelHPF[abs(lx)<=l_HPF] = 0.
		#from IPython import embed; embed()
		weight_inv_var_lmat[abs(lx)<=l_HPF] = 0.
		weight_gradient_lmat[abs(lx)<=l_HPF] = 0.

		T_mat = np.fft.ifft2( np.fft.fft2(T_mat) * lowelHPF).real
		T_mat2 = np.fft.ifft2( np.fft.fft2(T_mat2) * lowelHPF).real

	# Subtract mean
	obs_T_mat_meansub = T_mat - np.mean(T_mat)
	obs_T_mat_meansub2 = T_mat2 - np.mean(T_mat2)

	###print np.mean(obs_T_mat_meansub), np.mean(obs_T_mat_meansub2)

	### from IPython import embed; embed()

	if np.allclose(obs_T_mat_meansub, obs_T_mat_meansub2) == False:
		print '\n\t\tMaps for gradient and lensing are different\n'
	else:
		print '\n\t\tMaps for gradient and lensing are the same\n'
		pass

	##subplot(121);imshow(obs_T_mat_meansub);colorbar();subplot(122);imshow(obs_T_mat_meansub2);colorbar();title('in QE funcs');show();quit()

	#For doing deconvolution
	to_deconvolve_lmat = np.copy(beamtf_lmat)
	bad = np.where(to_deconvolve_lmat==0.)

	to_deconvolve_lmat[bad] = 1.
	deconvolve_lmat = 1./to_deconvolve_lmat
	#from IPython import embed;embed()
	#sys.exit()
	deconvolve_lmat[bad] = 0.

	highly_filtered = np.where((deconvolve_lmat > 1.0e10) | (deconvolve_lmat < 0) )
	deconvolve_lmat[highly_filtered] = 0.0

	if hasattr(beam_lmat_transfer_lmat_grad, '__len__'):
		deconvolve_lmat_grad = np.copy(beam_lmat_transfer_lmat_grad)
		bad = np.where(deconvolve_lmat_grad==0.)
		deconvolve_lmat_grad[bad] = 1.
		deconvolve_lmat_grad = 1./deconvolve_lmat_grad
		deconvolve_lmat_grad[bad] = 0.

		highly_filtered = np.where((deconvolve_lmat_grad > 1.0e10) | (deconvolve_lmat_grad < 0) )
		deconvolve_lmat_grad[highly_filtered] = 0.0
	else:
		deconvolve_lmat_grad = deconvolve_lmat

	#imshow(deconvolve_lmat);colorbar();show();quit()

	# Do filtering and deconvolution
	T_lmat = dxdy*np.fft.fft2(obs_T_mat_meansub)*deconvolve_lmat

	T_lmat2 = dxdy*np.fft.fft2(obs_T_mat_meansub2)*deconvolve_lmat_grad
	filtered_gradient_x_lmat = 1.j*lx*weight_gradient_lmat*T_lmat2 #20171012 - different map for gradient
	filtered_gradient_y_lmat = 1.j*ly*weight_gradient_lmat*T_lmat2 #20171012 - different map for gradient


	#W^T = weight_inv_var_l
	T_inv_var_lmat = weight_inv_var_lmat*T_lmat
	#### T_inv_var_lmat = T_inv_var_lmat.conjugate()
	#T_inv_var_lmat = weight_inv_var_lmat*T_lmat2

	#Convert to configuration space
	filtered_gradient_x_mat = (1./dxdy)*np.fft.ifft2(filtered_gradient_x_lmat)
	filtered_gradient_y_mat = (1./dxdy)*np.fft.ifft2(filtered_gradient_y_lmat)
	T_inv_var_mat = (1./dxdy)*np.fft.ifft2(T_inv_var_lmat)

	#for checking
	filtered_gradient_mat = (1./dxdy)*np.fft.ifft2(weight_gradient_lmat * T_lmat2)


	#subplot(121);imshow((filtered_gradient_mat.real));colorbar()
	#subplot(122);imshow((T_mat));colorbar();show()
	##subplot(121);imshow(np.fft.fftshift(weight_gradient_lmat));colorbar()
	##subplot(122);imshow(np.fft.fftshift(weight_inv_var_lmat));colorbar();show();

	'''
	#plot and check maps
	map_1 = np.fft.ifft2(T_lmat * weight_gradient_lmat).real
	map_2 = np.fft.ifft2(T_lmat * weight_inv_var_lmat).real
	subplot(121);imshow(map_1);colorbar()
	subplot(122);imshow(map_2);colorbar();show();quit()
	#imshow(T_inv_var_mat.real);colorbar();show();quit()
	'''
	#imshow(T_mat);colorbar();show();sys.exit()
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

	###filter2 = np.exp(-0.5*(l2d/float(ngaus))**2) #should it be a Gaussian smoothing or LPF at ngaus #20161117
	filter2 = np.ones(l2d.shape)
	above_beam_scale = np.where(l2d>=ngaus)
	filter2[above_beam_scale] = 0.
	
	#imshow(np.fft.fftshift(filter2), extent = [np.min(lx), np.max(lx), np.min(ly), np.max(ly)]);colorbar();show();quit()

	if not curl_test:
		teye = -ampval*l2d*1.j*(lx*GL_x_lmat + ly*GL_y_lmat)
	else:
		teye = -ampval*l2d*1.j*(ly*GL_x_lmat - lx*GL_y_lmat)

	kappa_lmat = filter2*teye*l2d/2.
	kappa_lmat[bad] = 0.0
	kappa_mat = (1./dxdy)*np.fft.ifft2(kappa_lmat)

	#print kappa_mat.real

	'''
	imshow(kappa_mat.real, interpolation = 'bicubic');colorbar();
	nx, ny = kappa_mat.shape
	axvline(nx/2., color = 'k');axhline(nx/2., color = 'k');show();sys.exit()
	'''
	if (0):
		kappa_mat = kappa_mat.real
		u, t, v = np.linalg.svd(kappa_mat)
		t[t>20.] = 0
		kappa_mat = np.dot(u * t, v)

	return kappa_mat

	if (0):

		figtext(0.35,0.95,r'\textbf{10 $\times$ 10 arcmin cutouts}', fontsize = 14)
		def get_planck_cmap():
			from matplotlib.colors import ListedColormap
			colombi1_cmap = ListedColormap(np.loadtxt("Planck_Parchment_RGB.txt")/255.)
			colombi1_cmap.set_bad("gray") # color of missing pixels
			colombi1_cmap.set_under("white") # color of background, necessary if you want to use
			# this colormap directly with hp.mollview(m, cmap=colombi1_cmap)
			return colombi1_cmap
		mycmap = get_planck_cmap()

		nx, ny = kappa_mat.shape
		snx = 10
		s, e = (nx/2)-snx, (nx/2)+snx
		intrp_val = 'bicubic'
		fsval = 10.

		"""
		from IPython import embed; embed()
		##gradient_dir = np.mean( np.degrees( np.arctan2(filtered_gradient_x_mat[s:e,s:e].real, filtered_gradient_y_mat[s:e,s:e].real) ) )

		'''
		dot_1 = np.dot(filtered_gradient_x_mat[s:e,s:e].real,filtered_gradient_x_mat[s:e,s:e].real)
		dot_2 = np.dot(filtered_gradient_y_mat[s:e,s:e].real,filtered_gradient_y_mat[s:e,s:e].real)
		grad_mag = np.sqrt( dot_1 * dot_2)
		dot_product = np.dot( filtered_gradient_x_mat[s:e,s:e].real, filtered_gradient_y_mat[s:e,s:e].real )
		angle = np.degrees( np.arccos ( dot_product / grad_mag ) )
		angle = np.nan_to_num(angle)
		'''
		xgrad, ygrad = filtered_gradient_x_mat[s:e,s:e].flatten().real, filtered_gradient_y_mat[s:e,s:e].flatten().real
		cross_product = map( lambda x, y: np.cross( [x], [y]), xgrad, ygrad )

		gradient_dir = np.mean( angle[angle<>0.])

		###gradient_dir = np.mean( np.degrees( np.arctan2(filtered_gradient_y_mat[s:e,s:e].real, filtered_gradient_x_mat[s:e,s:e].real) ) )
		##if gradient_dir<0: gradient_dir += 360.

		'''
		xgrad, ygrad = filtered_gradient_x_mat.real[s:e, s:e], filtered_gradient_y_mat.real[s:e, s:e]
		subplot(221);imshow(xgrad, origin = 'lower', cmap  =mycmap);colorbar()
		subplot(222);imshow(ygrad, origin = 'lower', cmap  =mycmap);colorbar();
		subplot(223);imshow( filtered_gradient_mat.real[s:e, s:e], origin = 'lower' , cmap  =mycmap);colorbar();title(gradient_dir)
		x, y = np.meshgrid(np.arange(e-s), np.arange(e-s))
		quiver( x, y, ygrad, xgrad)
		show(); sys.exit()
		'''
		"""


		subplot(331);imshow(filtered_gradient_mat.real[s:e,s:e], interpolation = intrp_val, cmap = mycmap);colorbar();
		title(r'\textbf{gradient}', fontsize = fsval);axis('off')
		axvline((e-s)/2., color = 'k');axhline((e-s)/2., color = 'k')
		subplot(332);imshow((T_inv_var_mat).real[s:e,s:e], interpolation = intrp_val, cmap = mycmap);colorbar();
		title(r'\textbf{fil. Map2}', fontsize = fsval);axis('off')
		axvline((e-s)/2., color = 'k');axhline((e-s)/2., color = 'k')

		subplot(333);imshow(kappa_mat.real[s:e,s:e], interpolation = intrp_val, cmap = mycmap);colorbar();
		title(r'\textbf{kappa QE}', fontsize = fsval);axis('off')
		axvline((e-s)/2., color = 'k');axhline((e-s)/2., color = 'k')
		##show()

		'''
		#rotate it in grad. direction
		import scipy.ndimage as ndimage
		##clf()
		angle_arr = np.arange(45., 360., 90.)
		angle_arr = np.concatenate( (angle_arr,[gradient_dir]) )
		for acnt, a in enumerate(angle_arr):
			kappa_mat_rot = ndimage.interpolation.rotate(kappa_mat.real, -a, reshape = False, mode = 'reflect')
			subplot(3,3,acnt+4);imshow(kappa_mat_rot.real[s:e,s:e], interpolation = intrp_val, cmap = mycmap);colorbar();
			title(r'\textbf{kappa QE}', fontsize = fsval);axis('off')
			axvline((e-s)/2., color = 'k');axhline((e-s)/2., color = 'k')
			title(r'\textbf{%.2f}' %a, fontsize = 10)
		show();sys.exit()

		subplot(336);imshow(kappa_mat_rot.real[s:e,s:e], interpolation = intrp_val, cmap = mycmap);colorbar();
		title(r'\textbf{kappa QE}', fontsize = fsval);axis('off')
		axvline((e-s)/2., color = 'k');axhline((e-s)/2., color = 'k')
		title(r'\textbf{%.2f}' %gradient_dir, fontsize = 10)
		'''

		subplot(334);imshow((filtered_gradient_x_mat).real[s:e,s:e], interpolation = intrp_val, cmap = mycmap);colorbar();
		title(r'\textbf{filtered gradX}', fontsize = fsval);axis('off')
		axvline((e-s)/2., color = 'k');axhline((e-s)/2., color = 'k')
		subplot(335);imshow((filtered_gradient_y_mat).real[s:e,s:e], interpolation = intrp_val, cmap = mycmap);colorbar();
		title(r'\textbf{filtered gradY}', fontsize = fsval);axis('off')
		axvline((e-s)/2., color = 'k');axhline((e-s)/2., color = 'k')
		subplot(337);imshow((filtered_gradient_x_mat*T_inv_var_mat).real[s:e,s:e], interpolation = intrp_val, cmap = mycmap);colorbar();
		title(r'\textbf{gradX x map2}', fontsize = fsval);axis('off')
		axvline((e-s)/2., color = 'k');axhline((e-s)/2., color = 'k')
		subplot(338);imshow((filtered_gradient_y_mat*T_inv_var_mat).real[s:e,s:e], interpolation = intrp_val, cmap = mycmap);colorbar();
		title(r'\textbf{gradY x map2}', fontsize = fsval);axis('off')
		axvline((e-s)/2., color = 'k');axhline((e-s)/2., color = 'k')

		'''
		kappa_psd = np.fft.fftshift( abs(np.fft.fft2(kappa_mat))**2. );
		subplot(336);imshow(kappa_psd, interpolation = 'bicubic', extent = [np.min(lx), np.max(lx), np.min(ly), np.max(ly)], cmap = mycmap);colorbar();title('kappa PSD')
		simmapparams = [nx, ny, 0.5, 0.5]
		cls_lensing = sims.fn_plot_pow_spec(simmapparams, [kappa_mat.real])[0][0]
		dls_fac = 1.#cls_lensing[:,0] * (cls_lensing[:,0] + 1) / 2 / np.pi
		subplot(339, yscale = 'log');plot(cls_lensing[:,0], dls_fac * cls_lensing[:,1]);title('kappa 1d power')
		xlim(0,ngaus); ylim(1e-10, 1e-5)#ylim(1e-5, 1.)
		'''

		show();sys.exit()

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


#Calculate QE normalization
def get_norm(lx, ly, C_l_unl_lmat, weight_gradient_lmat, weight_inv_var_lmat, dxdy):
	theta_lmat = get_angle(lx, ly)
	l2d = np.sqrt(lx**2. + ly**2.)

	#compute normalization
	nlx = lx.shape[0]
	nly = lx.shape[1]
	theta_lmat = get_angle(lx,ly)
	R_alpha = np.zeros((2,nlx,nly), dtype = np.cfloat)
	P_alpha = np.zeros((2,nlx,nly), dtype = np.cfloat)
	Q_alpha = np.zeros((2,nlx,nly), dtype = np.cfloat)

	R_alpha[0,:,:] = 0.5*l2d*np.exp(1.j*theta_lmat)
	R_alpha[1,:,:] = 0.5*l2d*np.exp(-1.j*theta_lmat)
	P_alpha[0,:,:] = l2d*np.exp(-1.j*theta_lmat)
	P_alpha[1,:,:] = l2d*np.exp(1.j*theta_lmat)
	Q_alpha[0,:,:] = 1.0
	Q_alpha[1,:,:] = 1.0

	R_beta = np.zeros((nlx,nly), dtype = np.cfloat)
	P_beta = np.zeros((nlx,nly), dtype = np.cfloat)
	Q_beta = np.zeros((nlx,nly), dtype = np.cfloat)

	R_beta[:,:] = 1.0
	P_beta[:,:] = weight_gradient_lmat
	Q_beta[:,:] = weight_inv_var_lmat

	R_gamma = np.zeros((5,nlx,nly), dtype = np.cfloat)
	P_gamma = np.zeros((5,nlx,nly), dtype = np.cfloat)
	Q_gamma = np.zeros((5,nlx,nly), dtype = np.cfloat)

	R_gamma[0,:,:] =  0.5*l2d*np.exp(1.j*theta_lmat)
	R_gamma[1,:,:] =  0.5*l2d*np.exp(-1.j*theta_lmat)
	R_gamma[2,:,:] = -0.5*l2d*np.exp(1.j*theta_lmat)
	R_gamma[3,:,:] = -0.5*l2d*np.exp(-1.j*theta_lmat)
	R_gamma[4,:,:] = l2d**2.

	P_gamma[0,:,:] = l2d*np.exp(-1.j*theta_lmat)*C_l_unl_lmat
	P_gamma[1,:,:] = l2d*np.exp(1.j*theta_lmat)*C_l_unl_lmat
	P_gamma[2,:,:] = l2d*np.exp(-1.j*theta_lmat)
	P_gamma[3,:,:] = l2d*np.exp(1.j*theta_lmat)
	P_gamma[4,:,:] = 1.0

	Q_gamma[0,:,:] = 1.0
	Q_gamma[1,:,:] = 1.0
	Q_gamma[2,:,:] = C_l_unl_lmat
	Q_gamma[3,:,:] = C_l_unl_lmat
	Q_gamma[4,:,:] = C_l_unl_lmat

	rpq_sum = np.zeros((nlx, nly), dtype = np.cfloat)
	for ii in xrange(0,2):
		for jj in xrange(0,5):
			R_tot = R_alpha[ii,:,:]*R_beta[:,:]*R_gamma[jj,:,:]
			P_tot = P_alpha[ii,:,:]*P_beta[:,:]*P_gamma[jj,:,:]
			Q_tot = Q_alpha[ii,:,:]*Q_beta[:,:]*Q_gamma[jj,:,:]
			Px = (1./dxdy)*np.fft.ifft2(P_tot)
			Qx = (1./dxdy)*np.fft.ifft2(Q_tot)
			rpq_sum += (R_tot)*dxdy*np.fft.fft2(Px*Qx)

	norm = rpq_sum

	return norm
