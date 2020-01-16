import numpy as np
import pdb
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
def get_kappa(lx, ly, beamtf_lmat, weight_gradient_lmat, weight_inv_var_lmat, T_mat, norm_lmat, dxdy, ngaus):

    # Following steps in 1411.7999

    #Useful later
    l2d = np.sqrt(lx**2. + ly**2.)

    # Subtract mean
    obs_T_mat_meansub = T_mat - np.mean(T_mat)

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

    #imshow(deconvolve_lmat);colorbar();show();quit()

    # Do filtering and deconvolution
    T_lmat = dxdy*np.fft.fft2(obs_T_mat_meansub)*deconvolve_lmat
    
    #Filtered gradient field
    filtered_gradient_x_lmat = 1.j*lx*weight_gradient_lmat*T_lmat
    filtered_gradient_y_lmat = 1.j*ly*weight_gradient_lmat*T_lmat

    #W^T = weight_inv_var_l
    T_inv_var_lmat = weight_inv_var_lmat*T_lmat

    #Convert to configuration space
    filtered_gradient_x_mat = (1./dxdy)*np.fft.ifft2(filtered_gradient_x_lmat)
    filtered_gradient_y_mat = (1./dxdy)*np.fft.ifft2(filtered_gradient_y_lmat)
    T_inv_var_mat = (1./dxdy)*np.fft.ifft2(T_inv_var_lmat)

    '''
    #plot and check maps
    map_1 = np.fft.ifft2(T_lmat * weight_gradient_lmat).real
    map_2 = np.fft.ifft2(T_lmat * weight_inv_var_lmat).real
    subplot(121);imshow(map_1);colorbar()
    subplot(122);imshow(map_2);colorbar();show();quit()
    #imshow(T_inv_var_mat.real);colorbar();show();quit()
    '''

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

    teye = -ampval*l2d*1.j*(lx*GL_x_lmat + ly*GL_y_lmat)
    kappa_lmat = filter2*teye*l2d/2.
    kappa_lmat[bad] = 0.0
    kappa_mat = (1./dxdy)*np.fft.ifft2(kappa_lmat)

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

