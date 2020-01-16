import numpy as np
from pylab import *
import pickle; import gzip
from scipy import ndimage
import scipy as sc
import scipy.signal
import scl_cmb
import scipy.fftpack as fft

#import pdb

#Get angle for specified vectors
def get_angle(x,y):
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
def get_kappa_old(lx, ly, weight_gradient_lmat, weight_inv_var_lmat, T_mat, norm_lmat, dxdy, ngaus, **kwargs):

   
    sims = scl_cmb.simulations() #babbloo

    theta_lmat = get_angle(lx, ly)
    l2d = np.sqrt(lx**2. + ly**2.)

    #Convert observed map to Fourier space
    #check this -- apodization necessary?
    #x, y = np.meshgrid( np.arange(0,ny), np.arange(0,nx) )
    #apod = 1.0#np.sin( np.pi/nx*x )*np.sin( np.pi/nx*y )

    '''
    #### T_lmat = fft.fft2(T_mat)
    ### babbloo padding zeros
    pnx = pny = T_mat.shape[0]*2
    nx, ny = T_mat.shape
    PADDING = np.zeros((pnx,pny))
    p1,p2 = int(pnx/2-nx/2), int(pny/2+ny/2)
    #PADDING[p1:p2,p1:p2] = T_mat
    PADDING[0:nx,0:ny] = T_mat

    T_lmat = fft.fft2(PADDING)#[p1:p2,p1:p2]
    subplot(131);css = imshow(fft.fftshift(T_lmat).real);css.set_clim(np.median(T_lmat.ravel())/2., np.median(T_lmat.ravel())*2.);colorbar()
    subplot(132);css = imshow(T_lmat_1.real);css.set_clim(np.median(T_lmat.ravel())/2., np.median(T_lmat.ravel())*2.);colorbar()
    T_lmat_2 = scfns.pad_ft(T_lmat_1)
    subplot(133);css = imshow(T_lmat_2.real);css.set_clim(np.median(T_lmat.ravel())/2., np.median(T_lmat.ravel())*2.);colorbar();show();quit()
#   T_lmat =  np.roll(np.roll(T_lmat,-nx/2,axis=0),-ny/2,axis=1)
    #T_lmat = scfns.pad_ft(T_lmat)
    #T_lmat = scfns.unpad_ft(T_lmat) #unpad now

    #T_lmat =  np.roll(np.roll(T_lmat,-nx,axis=0),-ny,axis=1)

    css=imshow(fft.ifft2(T_lmat).real);css.set_clim(min(T_mat.ravel()),max(T_mat.ravel()));colorbar()
    #subplot(122);css=imshow(fft.ifft2(PADDING2).real - fft.ifft2(T_lmat).real);colorbar()
    show()
    quit()
    '''

    T_lmat = fft.fft2(T_mat)
    T_lmat = dxdy*T_lmat

    padd_please = 0
    if padd_please:
	    T_lmat = scfns.pad_ft(T_lmat)
	    weight_gradient_lmat = scfns.pad_ft(weight_gradient_lmat)
	    weight_inv_var_lmat = scfns.pad_ft(weight_inv_var_lmat)

    diff_grad = 0
    if ('T_map_for_gradient' in kwargs):
        diff_grad = 1
        T_map_for_gradient = kwargs['T_map_for_gradient']
        T_lmat_forgrad = dxdy*fft.fft2(T_map_for_gradient)
    #W^TT = weight_gradient_l
    if (diff_grad == 0):
	GRAD_WEI_MAP_FFT = weight_gradient_lmat*T_lmat
	#### babbloo changing this
	if padd_please:
		#GRAD_WEI_MAP_FFT = sc.signal.fftconvolve(fft.ifft2(weight_gradient_lmat),T_mat, mode = 'same')
		GRAD_WEI_MAP_FFT = sims.fn_convolve_using_QL(weight_gradient_lmat,T_lmat)
		GRAD_WEI_MAP_FFT = scfns.unpad_ft(GRAD_WEI_MAP_FFT) #unpad again

        filtered_gradient_x_lmat = 1.j*lx*GRAD_WEI_MAP_FFT
        filtered_gradient_y_lmat = 1.j*ly*GRAD_WEI_MAP_FFT
    if (diff_grad == 1):
        filtered_gradient_x_lmat = 1.j*lx*weight_gradient_lmat*T_lmat_forgrad
        filtered_gradient_y_lmat = 1.j*ly*weight_gradient_lmat*T_lmat_forgrad

    #W^T = weight_inv_var_l    
    T_inv_var_lmat = weight_inv_var_lmat*T_lmat

    '''
    MAP_ORI = np.fft.ifft2(T_lmat).real/dxdy * 1e6
    MAP_ORI_INV_VAR = np.fft.ifft2(T_inv_var_lmat).real/dxdy * 1e6
    MAP_ORI_GRAD_WEI = np.fft.ifft2(GRAD_WEI_MAP_FFT).real/dxdy * 1e6
    subplot(131);imshow( MAP_ORI );colorbar()
    subplot(132);imshow(MAP_ORI_INV_VAR);colorbar()
    subplot(133);imshow(MAP_ORI - MAP_ORI_GRAD_WEI);colorbar()

    show();quit()
    '''
    ####babbloo changing above to pad zeros before convolution
    if padd_please:
		#T_inv_var_lmat = sc.signal.fftconvolve(fft.ifft2(weight_inv_var_lmat),T_mat, mode = 'same')
		T_inv_var_lmat = sims.fn_convolve_using_QL(weight_inv_var_lmat,T_lmat)
		T_inv_var_lmat = scfns.unpad_ft(T_inv_var_lmat)
    
    #subplot(221);imshow( fft.ifftshift( fft.ifft2(weight_gradient_lmat).real ) ); colorbar()
    #subplot(222);imshow( fft.ifftshift( fft.ifft2(weight_inv_var_lmat).real ) ); colorbar()
    #subplot(223);imshow( fft.ifft2(GRAD_WEI_MAP_FFT).real); colorbar()
    #subplot(224);imshow( fft.ifft2(T_inv_var_lmat).real); colorbar()
    #subplot(131);imshow(DUMMY1);colorbar()
    #subplot(132);imshow(DUMMY2);colorbar()
    #subplot(133);imshow(DUMMY1 - DUMMY2);colorbar()
    #show();quit()

    if 1 == 1:
	    #Convert to configuration space
	    filtered_gradient_x_mat = (1./dxdy)*fft.ifft2(filtered_gradient_x_lmat)
	    filtered_gradient_y_mat = (1./dxdy)*fft.ifft2(filtered_gradient_y_lmat)
	    T_inv_var_mat = (1./dxdy)*fft.ifft2(T_inv_var_lmat)

	    #Multiply filtered fields in configuration space and fourier transform
	    GL_x_lmat = dxdy*fft.fft2(filtered_gradient_x_mat*T_inv_var_mat)
	    GL_y_lmat = dxdy*fft.fft2(filtered_gradient_y_mat*T_inv_var_mat)
    else:

	    GL_x_lmat = fft.fft2(fft.ifft2(filtered_gradient_x_lmat)*fft.ifft2(T_inv_var_lmat))
	    GL_y_lmat = fft.fft2(fft.ifft2(filtered_gradient_y_lmat)*fft.ifft2(T_inv_var_lmat))

    #subplot(221);imshow((GL_x_lmat).real);colorbar()
    #subplot(222);imshow((GL_y_lmat).real);colorbar()
    #subplot(223);imshow((filtered_gradient_x_mat*T_inv_var_mat).real);colorbar()
    #subplot(224);imshow((filtered_gradient_y_mat*T_inv_var_mat).real);colorbar()
    #show();quit()

    #Filtered divergence, and applying normalization
    #WHAT DOES ngaus DO???77777777
    ### ngaus = 20000
    ### filter2 = np.exp(-0.5*(l2d/float(ngaus))**2)

    #subplot(121);imshow(np.fft.fftshift(filter2));colorbar()
    #### get the healpy type beam
    fwhm = np.radians(1./60.)
    filter2 = sims.gauss_beam(fwhm, l2d)[0]
    #subplot(122);imshow(np.fft.fftshift(filter2));colorbar()
    #show();quit()
    
    temp_norm_lmat = np.copy(norm_lmat)
    bad = np.where(temp_norm_lmat == 0.0)
    temp_norm_lmat[bad] = 1.0

    ampval = 1./temp_norm_lmat
    #ampval = ampval * filter2

    #css=imshow(fft.ifftshift(fft.ifft2(temp_norm_lmat).real));colorbar();show();quit()

    teye = -ampval*l2d*1.j*(lx*GL_x_lmat + ly*GL_y_lmat)
    kappa_lmat = filter2*teye*l2d/2.
    ### babbloo removing filter2
    #kappa_lmat = teye*l2d/2.

    #imshow(fft.ifft2(kappa_lmat).real);colorbar();grid(True,ls='solid');show();quit()

    kappa_lmat[bad] = 0.0
    kappa_mat = (1./dxdy)*fft.ifft2(kappa_lmat)

    return kappa_mat


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
            Px = (1./dxdy)*fft.ifft2(P_tot)
            Qx = (1./dxdy)*fft.ifft2(Q_tot)
            rpq_sum += (R_tot)*dxdy*fft.fft2(Px*Qx)

    norm = rpq_sum
    return norm


