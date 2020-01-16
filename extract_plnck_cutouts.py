'''
Extracts galaxy cluster cutouts from planck LGMCA maps
Which are eventually be used as gradient maps for SPT-SZ clusters


'''

def rotate_tf(tf_map_center, ra0, dec0, ra, dec, tf_rot_angle = None):
        '''                                                                                                                                                                                               
        Based on code by Srini Raghunathan                                                                                                                                                                
        tf_map_center is transfer funciton evaluated at the center of the map                                                                                                                             
        ra0, dec0 specify the map center                                                                                                                                                                  
        ra, dec are the coordinates at which we want to evaluate the TF                                                                                                                                   
        '''
        import scipy.ndimage as ndimage
        if tf_rot_angle == None:
                tf_rot_angle = get_tf_rot_angle(ra,dec,ra0,dec0) #radians                                                                                                                                 
        tf_rotated = ndimage.interpolation.rotate(tf_map_center, tf_rot_angle*(180./np.pi), reshape = False, mode = 'nearest')

        return tf_rotated

def get_tf_rot_angle(ra_orig, dec_orig, ra0_orig, dec0_orig):
        """                                                                                                                                                                                               
        Author: Srini Raghunathan                                                                                                                                                                         
        Modified: Eric Baxter                                                                                                                                                                             
                                                                                                                                                                                                          
        this is only valid for ZEA (proj 5) projection                                                                                                                                                    
        http://arxiv.org/pdf/1111.7245v1.pdf - section 5.4 for reference                                                                                                                                  
        (although note that formulae therein seem to have an error)                                                                                                                                       
        this is roughly the angle required to convert the fourier space filtering to real space filtering                                                                                                 
        Should match output of get_proj5_angle.pro in the SPT analysis repo                                                                                                                               
        """

        ra_val = np.radians(ra_orig)
        dec_val = np.radians(dec_orig)

        #do the 90deg. subtraction for declinations                                                                                                                                                       
        dec_val = np.radians(90.) - dec_val
        dec0 = np.radians(90.) - np.radians(dec0_orig)
        ra0 = np.radians(ra0_orig)

        #nr - numerator; dr - denominator; _1 - first term; _2 - second term, etc. in the equation.                                                                                                       
        gamma_dr = 1 + ( np.cos(dec0) * np.cos(dec_val) ) + ( np.sin(dec0) * np.sin(dec_val) * np.cos(ra0 - ra_val) )
        gamma = 0.5 / gamma_dr

        A_1 = gamma * np.sin(dec0) * np.sin(ra_val - ra0) * ( ( np.sin(dec0)*np.cos(dec_val) ) - ( np.sin(dec_val)*np.cos(dec0)*np.cos(ra_val - ra0) ) )
        A_2 = np.cos(dec0) * np.sin(ra_val-ra0)
        A = A_1 + A_2

        B_1 = gamma * np.sin(dec0) * (np.sin(ra_val-ra0))**2. * np.sin(dec_val)
        B_2 = np.cos(ra_val-ra0)

        B = B_1 + B_2

        alpha = np.arctan(-1.*A/B)
        return alpha


import numpy as np
import healpy as hp, pickle,gzip
from pylab import *
import scipy.ndimage as ndimage
catalog_file = 'data_sp/sptsz/data/spt/catalog.pkl.gz' #'data_sp/sptsz/catalog.pkl.gz'
lgmca_map_file = 'data_sp/sptsz/data/LGMCA/WPR2_CMB_muK.fits'#'data_sp/sptsz/WPR1_CMB_muK.fits'
lgmca_map_file = 'compton_sptsz_nside2048.fits'
mask_file = '/Users/sanjaykumarp/transfer/spt/MASK.fits'
op_fname = 'data_sp/sptsz/lgmca2.pkl.gz'
op_fname = 'ysz_cts.pkl.gz'
catalog = pickle.load(gzip.open(catalog_file))
sz_ra, sz_dec, conf, zvals,field_name = catalog['RA'],catalog['DEC'],catalog['XI'],catalog['REDSHIFT'],catalog['field_name']
lgmca_map = hp.read_map(lgmca_map_file)
#mask_map = hp.read_map(mask_file)
#lgmca_map *= mask_map
cutouts = {}
dx = 0.5 # arcminutes

boxsize = 100/dx# total of 100 arcminutes
ra0 ,dec0 = 180.,0
plnk_cutouts = {}
for i,(ra,dec) in enumerate(zip(sz_ra,sz_dec)):
    keyname = ra,dec,zvals[i],conf[i],1,field_name[i]
    if zvals[i] == 0:
	    continue
    ##sys.exit()
    lgmca_gnom = hp.gnomview(lgmca_map, rot = [ra, dec,-90.],flip = 'geo', xsize = boxsize, coord = ['G','C'], reso = dx, return_projected_map = 1)
    
    close()
    ##imshow(lgmca_gnom);colorbar();show();sys.exit()
    tf_rot_angle = get_tf_rot_angle(ra,dec,ra0,dec0)
    tf_rot_angle *= -1
    lgmca_gnom = rotate_tf(lgmca_gnom, ra0, dec0, ra, dec, tf_rot_angle = tf_rot_angle)
    lgmca_gnom = lgmca_gnom/1e6 # to convert it into kelvin
    cutouts[keyname]= np.asarray([lgmca_gnom])
    print i
plnk_cutouts['cutouts'] = cutouts
plnk_cutouts['resol_arcmins'] = dx
plnk_cutouts['cutout_size_arcmins'] = 100
plnk_cutouts['cluster_catalogue'] = 'data_sp/sptsz/list_for_sanjay.txt'
plnk_cutouts['map_source'] = lgmca_map_file

pickle.dump(plnk_cutouts,gzip.open(op_fname,'w'))
