def wrapAzAround(az, az_center=180., modify_in_place=True):
    """
    Like wrapAz, but instead of always returning an array with values in [0,360),
    wrapAzAround can change the "branch cut" so that all values returned are
    within 180. of az_center.
    
    INPUTS
       az: A numpy array of azimuth values, in degrees.

       az_center [180.]: wrap around this value.
          The default of 180 makes this equivalent to the simple wrapAz.
          In general, the values returned are in [az_center-180., az_center+180.)
          e.g. if az_center is zero, values returned are in [-180.,+180.).
          e.g. if we are making a map centered around az0, set az_center=az0
       
       modify_in_place [True]: If True, the input array is 
          overwritten with the result.

    OUTPUT
       fall in the range [0,360).
       
    EXAMPLE OF USE

    In [1]: from sptpol_software.util.tools import wrapAzAround

    In [2]: arr = np.array([-1.,0.,1.,179.,180.,181.,359.,360.,361.])
    
    In [3]: wrapAzAround(arr, az_center=180., modify_in_place=False)
    Out[3]: array([ 359.,    0.,    1.,  179.,  180.,  181.,  359.,    0.,    1.])

    In [4]: wrapAzAround(arr, az_center=0., modify_in_place=False)
    Out[4]: array([  -1.,    0.,    1.,  179., -180., -179.,   -1.,    0.,    1.])
       An array of azimuth values, adjusted so that all values
    """
    az_shift = 180.-az_center
    if modify_in_place:
        az += az_shift
        wrapAz(az, modify_in_place=True)
        az -= az_shift
        return az
    else:
        return wrapAz(az+az_shift, modify_in_place=False) - az_shift
    

def wrapAz(az, modify_in_place=True):
    """
    Looks through an input array for values < 0 or > 360,
    and wraps them into the [0,360) range.
    
    INPUTS
       az: A numpy array of azimuth values, in degrees.

       modify_in_place [True]: If True, the input array is 
          overwritten with the result.

    OUTPUT
       An array of azimuth values, adjusted so that all values
       fall in the range [0,360).
       
    EXAMPLE OF USE

    In [1]: from sptpol_software.util.tools import wrapAz

    In [2]: arr = np.array([-722, -4, -3, -2, -1, 0,   1,   2,   3, 359, 360, 361, 1132])
    
    In [3]: wrapAz(arr, modify_in_place=False)
    Out[3]: array([358, 356, 357, 358, 359,   0,   1,   2,   3, 359,   0,   1,  52])

    In [4]: arr
    Out[4]: 
    array([-722,   -4,   -3,   -2,   -1,    0,    1,    2,    3,  359,  360,
            361, 1132])

    In [5]: wrapAz(arr)
    Out[5]: array([358, 356, 357, 358, 359,   0,   1,   2,   3, 359,   0,   1,  52])

    In [6]: arr
    Out[6]: array([358, 356, 357, 358, 359,   0,   1,   2,   3, 359,   0,   1,  52])
    """
    if np.isscalar(az):
        scalar=True
        az = np.asarray([az])
    else:
        scalar=False
    
    
    if modify_in_place: az2use = np.asarray(az)
    else: az2use = np.asarray(az).copy()

    # Shift all values outside of [0,360) into the range. 
    # Use a while loop to allow for the possibility that some
    # values could be further away then 360 degrees from the range.
    # i.e., correctly treat cases like [-743, 1132].
    while az2use.min() < 0:
        az2use[az2use < 0] += 360.
    while az2use.max() >= 360:
        az2use[az2use >= 360] -= 360.

    if scalar:
        return az2use[0]
    else:
        return az2use

########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################

"""
Classes and functions relating directly to the sky. Contains, for example, Map objects
(which hold data organized by location on the sky) and functions to convert
between different sky projections. 

Classes
    Map
    PolarizedMap
    TrueSky
    ObservedSky
    PolarizedTrueSky
    PolarizedObservedSky
    SkyStream
    PolarizedSkyStream

Functions
    pix2Ang
    ang2Pix

Non-Core Dependencies
   NumPy
   SciPy
"""

__metaclass__ = type  #Use "new-style" classes
__author__    = "Stephen Hoover"
__email__     = "hoover@kicp.uchicago.edu"
__version__   = "1.08"
__date__      = "2013-09-30" #Date of last modification

import warnings, pdb, unittest
import os
import numexpr as ne
import numpy as np
import scipy as sp
from scipy import interpolate#, weave, ndimage
from numpy import sin, cos, sqrt
from glob import glob
import matplotlib.pyplot as plt
import matplotlib
"""
import telescope
import receiver
from .. import constants
from ..constants import DTOR, RTOD
from .. import float_type
from ..util import tools, math, fits, hdf, time
from ..util.tools import struct
from ..data import c_interface
from copy import deepcopy
from mpl_toolkits.axes_grid1 import make_axes_locatable
"""

DTOR = np.pi/180.
RTOD = 180./np.pi

ne.set_num_threads(int(os.environ.get('OMP_NUM_THREADS',4))) # Don't use all the processors!

tqu=['T','Q','U']
teb=['T','E','B']
#pol_index = {'T':0, 'Q':1, 'U':2, 'I':0, 'V':3, 'T':0, 'E':1, 'B':2, 'I':0} # For translating letters into array indices

proj_name_to_index = {'Sanson-Flamsteed':0, 'CAR':1, 'SIN':2, 'Healpix':3, 'Sterographic':4, 'Lambert':5, 'CAR00':7, 'rot00':8, 'BICEP':9}
__proj_name_to_index = dict( [(key.lower(), value) for key, value in proj_name_to_index.iteritems()]) # Lower-case version for comparisons.

# Projections tested in testAllProj:
proj_index_in_idl = [0,1,2,4,5]
proj_to_test = [0,1,2,4,5,7]

##########################################################################################
def read(filename, file_type=None, verbose=True):
    """
    Read out a Map or PolarizedMap which has been written to an output file.
    
    INPUTS
        filename : (string or list of strings) The filename(s) to read. May be given 
            as a fully-specified filename, a string parseable by glob, or a list which
            may contain both kinds of strings.
        
        file_type [None]: (string) If the "type" input is "FITS" or "HDF5" 
            (case insensitive), then we read the file as being in that format, 
            regardless of extension.
            
        verbose [True]: (bool) Extra screen output.
        
    OUTPUT
        A Map or PolarizedMap object.
    """
    # If the filename is a string, proceed. Otherwise, treat it as a list of files, and read each separately.
    try: 
        filename+''
        # In case the file name is a string with wildcards, feed it through glob.
        globized_name = glob(filename)
        if len(globized_name)==1:
            # The Map class's read function will return a PolarizedMap if that's what the file contains.
            return Map.read(globized_name[0], file_type=file_type)
        else:
            return read(globized_name, file_type=file_type, verbose=verbose)
    except TypeError: 
        if verbose: print "Reading %d files." % len(filename)
        return map(lambda x: read(x, file_type=file_type), filename)
    
##########################################################################################
def projNameToIndex(proj):
    """
    Converts a proj name to the corresponding index. Lets us use more descriptive function calls.
    
    Supported projections are:

     0, "Sanson[-Flamsteed]":  Sanson-Flamsteed projection (x = ra*cos(dec), y = dec)
     1, "CAR":  CAR projection (x = ra, y = dec)
     2, "SIN":  SIN projection
     3, "Healpix":  Healpix (not a projection at all, but pixels on the sphere) [NOT IMPLEMENTED]
     4, "Stereo[graphic]":  stereographic projection [SAME AS proj5??]
     5, "Lambert":  Lambert azimuthal equal-area projection
     7, "CAR00": CAR projection (x = ra, y = dec), if the field is first 
         rotated so that the RA, dec center is 0,0.
     8, "rot00"
     9, "BICEP":  ra,dec but pixel size in ra is set to be square at the mean dec: roughly width reso/cos(mean_dec), 
                  but rounded to make integer pixels for ra -55->+55 and dec -70->-45
                  Basically http://en.wikipedia.org/wiki/Equirectangular_projection, with phi1~middle_dec (but rounded)
    INPUTS
        proj: (string, int): Returns ints directly, or converts strings to a corresponding proj index.
        
    OUTPUT
        An integer proj index.
    """
    # In case 'proj' is a complete name or index, return right away.
    try: proj=__proj_name_to_index.get(proj.lower(), proj)
    except AttributeError: pass
     
    if proj in proj_name_to_index.values(): return proj

    # If 'proj' isn't a complete name, look for acceptable partial names.
    proj = str(proj).lower()
    if proj.startswith('sanson') or proj=='sf': proj=0
    elif proj=='car00': proj=7
    elif proj.startswith('car'): proj=1
    elif proj.startswith('sin'): proj=2
    elif proj.startswith('healpix'): proj=3
    elif proj.startswith('stereo'): proj=4
    elif proj.startswith('lambert'): proj=5
    elif proj.startswith('bicep'): proj=9
    else:
        raise ValueError('I don\'t know what projection "'+proj+'" is!') 

    return proj

##########################################################################################
def pix2Ang(pixel_coords, ra_dec_center, reso_arcmin, map_pixel_shape, proj=0, wrap=True):
    """
    Supported projections are:

     0:  Sanson-Flamsteed projection (x = ra*cos(dec), y = dec)
     1:  CAR projection (x = ra, y = dec)
     2:  SIN projection
     3:  Healpix (not a projection at all, but pixels on the sphere) [NOT IMPLEMENTED]
     4:  stereographic projection
     5:  Oblique Lambert azimuthal equal-area projection  (ref p. 185, Snyder, J. P. 1987, Map Projections-A Working Manual (Washington, DC: U.S. Geological Survey))
     7:  CAR projection, with map rotated so that the center is at RA, dec = (0, 0)
     8:  "rot00"
     9:  BICEP projection
    INPUTS
        pixel_coords : (2-element tuple of arrays) A tuple or list of arrays, [y_coord, x_coord].
            Note the order! The first element is the "y", the mostly-dec coordinate, and the 
            second element is the "x", the mostly-RA coordinate.
        
        ra_dec_center : (2-element array) The [RA, declination] of the center of the map.
        
        reso_arcmin : (float) The width of each pixel in arcminutes, assumed to be the same for
            both x and y directions.
            
        map_pixel_shape : (2-element array) The height and width of the map, in pixels.
        
        proj [0]: (int or string) Which map projection should I use to turn pixel coordinates
            on a flat map into angles on the curved sky? May be an integer index, or string
            name of the projection.
            
    OUTPUT
        (ra, dec): A 2-tuple of arrays. 'ra' and 'dec' are each arrays with the same number of elements as
        the arrays in the 'pixel_coords' input. 'ra' is the right ascension in degrees 
        (wrapped to the range [0, 360) if wrap==True), and 'dec' is the declination in degrees.
        returns the ra,dec of the CENTER of the pixels
    """
    
    # Convert the "proj" input to an index.
    proj = projNameToIndex(proj)
    
    # Break out the x and y coordinates, and cast them as floats.
    float_type = np.float32
    pixel_coords = (y_coord, x_coord) = float(pixel_coords[0]), float(pixel_coords[1])
    n_pixels = map_pixel_shape#.astype(float_type)
    #pixel_coords = (y_coord, x_coord) = float(pixel_coords[0]), float(pixel_coords[1])
    #n_pixels = map_pixel_shape#.astype(float_type)

    # shift to the center of the pixel, subtract off npix/2 to center around 0, then convert to degrees
    y_coord = (y_coord + 0.5 - 0.5*n_pixels[0])*reso_arcmin/60
    x_coord = (x_coord + 0.5 - 0.5*n_pixels[1])*reso_arcmin/60

    ra_dec_center_rad = ra_dec_center * DTOR # Convert to radians

    if proj==0:
        dec = ra_dec_center[1] - y_coord
        ra = x_coord / np.cos(dec*DTOR) + ra_dec_center[0]
    elif proj==1:
        dec = ra_dec_center[1] - y_coord
        ra = x_coord + ra_dec_center[0]
    elif proj==2:
        rho = sqrt( x_coord**2 + y_coord**2) * DTOR
        c = np.arcsin(rho)
        phi_temp = np.arcsin(cos(c)*sin(ra_dec_center_rad[1]) - DTOR*y_coord*cos(ra_dec_center_rad[1]))
        lambda_temp = RTOD * np.arctan2(x_coord*DTOR*sin(c),  
                                        rho*cos(ra_dec_center_rad[1])*cos(c) + y_coord*DTOR*sin(ra_dec_center_rad[1])*sin(c) )
        bad_rho = rho < 1e-8
        phi_temp[bad_rho] = ra_dec_center_rad[1]
        lambda_temp[bad_rho] = 0.
        
        dec = phi_temp * RTOD
        ra = ra_dec_center[0] + lambda_temp
    elif proj==5 or proj==4:
        rho = sqrt( x_coord**2 + y_coord**2) * DTOR
        if proj==5:
          c = 2*np.arcsin(rho/2)
        if proj==4:
          c = 2*np.arctan(rho/2)
        phi_temp = np.arcsin(cos(c)*sin(ra_dec_center_rad[1]) - DTOR*y_coord*sin(c)/rho*cos(ra_dec_center_rad[1])) 
        lambda_temp = RTOD * np.arctan2(x_coord*DTOR*sin(c),
                                          rho*cos(ra_dec_center_rad[1])*cos(c) + y_coord*DTOR*sin(ra_dec_center_rad[1])*sin(c))
        bad_rho = rho < 1e-8

        ####phi_temp[bad_rho] = ra_dec_center_rad[1] #what is this?
        ####lambda_temp[bad_rho] = 0.
        
        dec = phi_temp * RTOD
        ra = ra_dec_center[0] + lambda_temp
    elif proj==7:
        ra, dec = applyPointingOffset(ra_dec_center, [x_coord, y_coord], offset_units='degrees', as_azel=False)
    elif proj==9:
        # BICEP  ra_dec_center=(0.0, -57.5)
        ny, nx = n_pixels
        #nx, ny = (236, 100) # for BICEP
        pixsize = reso_arcmin/60.  # pixsize = 0.25 for BICEP  ->  reso_arcmin=15.
        sy  = pixsize  # dec spacing in degrees
        asx = pixsize/cos(ra_dec_center_rad[1])  # approximate spacing in ra
        sx  = np.round(asx*nx)/nx  # round so an integer number of degrees fit in the given number of pixels.
        dec = ra_dec_center[1] - y_coord
        ra  = x_coord*sx/sy + ra_dec_center[0]  # scale x_coord to make it square in center (undo ang2Pix)
    else:
        raise ValueError("I don't know what to do with proj "+str(proj)+".")
    
    if wrap:
        # Wrap RA values to the [0, 360) range.
        #tools.wrapAz(ra)
	wrapAz(ra)
    else:
        # make sure branch cut for ra is opposite map center
        too_low = np.where(ra-ra_dec_center[0] < -180.)
        ra[too_low] += 360.
        too_high = np.where(ra-ra_dec_center[0] >= 180.)
        ra[too_high] -= 360.
        
    return ra, dec

    if return_validity:
        n_pix_y, n_pix_x = n_pixels
        pixel_y_is_good = ne.evaluate("(y_pixel_coord>=0) & (y_pixel_coord<n_pix_y)")
        pixel_x_is_good = ne.evaluate("(x_pixel_coord>=0) & (x_pixel_coord<n_pix_x)")
    
        return telescope.PointingSequence(pixel_coords, map_shape=map_pixel_shape, proj=proj,
                                          pixel_resolution_arcmin=reso_arcmin), (pixel_y_is_good, pixel_x_is_good)
    else:
        return telescope.PointingSequence(pixel_coords, map_shape=map_pixel_shape, proj=proj,
                                          pixel_resolution_arcmin=reso_arcmin)

##########################################################################################

def ang2Pix(ra_dec, ra_dec_center, reso_arcmin, map_pixel_shape, proj=0, 
            round=True, bin_center_zero=True, return_validity=False, use_c_code=False):
    """
    Supported projections are:

     0:  Sanson-Flamsteed projection (x = ra*cos(dec), y = dec)
     1:  CAR projection (x = ra, y = dec)
     2:  SIN projection
     3:  Healpix (not a projection at all, but pixels on the sphere) [NOT IMPLEMENTED]
     4:  stereographic projection
     5:  Oblique Lambert azimuthal equal-area projection (ref p. 185, Snyder, J. P. 1987, Map Projections-A Working Manual (Washington, DC: U.S. Geological Survey))
     7:  CAR projection, with map rotated so that the center is at RA, dec = (0, 0)
     8:  "rot00"
     9:  BICEP projection
     
    INPUTS
        ra_dec : (2-element tuple of arrays) A tuple or list of arrays, [ra, dec]. Can also 
            be a PointingSequence object. In degrees.
        
        ra_dec_center : (2-element array) The [RA, declination] of the center of the map. In degrees.
        
        reso_arcmin : (float) The width of each pixel in arcminutes, assumed to be the same for
            both x and y directions.
            
        map_pixel_shape : (2-element array) The height and width of the map, in pixels.  (n_pix_y, n_pix_x)
        
        proj [0]: (int or string) Which map projection should I use to turn angles on the
            curved sky into pixels on a flat map? May be an integer index, or string
            name of the projection.
        
        round [True]: (bool)
        
        bin_center_zero [True]: (bool)  only applies if round==False.  If bin_center_zero==True, shift output by 0.5
        
        return_validity [True]: (bool) If True, return a 2-tuple stating which pixel coordinates are good.
            If False, return only the PointingSequence object.
        
        use_c_code [False]: (bool) If True, use compiled C code to find the pointing. ~10-20% faster than
            the Python implementation.
        
    OUTPUT
        A PointingSequence object and a 2-tuple. The PointingSequence has pixel coordinates, in 
        the form (y_coord, x_coord), where each element is an array of pixel indices.
        The second output tuple designates which pixel coordinates are "good", i.e., within
        bounds. It is of the form (good_ys, good_xs), where each element is
        an array of booleans.
        
        Note:  y_coord is more like elevation than dec: index zero corresponds to the least negative dec 
               (closest to horizon at south pole)
               
     EXAMPLE to test 360 wrap:
        from sptpol_software.observation import sky
        ra =  np.array([-1.,0.,1., 89.,90.,91., 179.,180.,181., 269.,270.,271., 359.,360.,361.])
        dec = np.array([0]*len(ra))  # equator
        sky.ang2Pix([ra, dec], ra_dec_center=[0,    0], reso_arcmin=60., map_pixel_shape=(360, 360), return_validity = True, use_c_code = False, proj = 1)
        sky.ang2Pix([ra, dec], ra_dec_center=[180., 0], reso_arcmin=60., map_pixel_shape=(360, 360), return_validity = True, use_c_code = False, proj = 1)
        sky.ang2Pix([ra, dec], ra_dec_center=[180., 0], reso_arcmin=60., map_pixel_shape=(360, 360), return_validity = True, use_c_code = False, proj = 0)
        dec = np.array([-60]*len(ra))  # -60 (so cos = 0.5)
        sky.ang2Pix([ra, dec], ra_dec_center=[180., -60], reso_arcmin=60., map_pixel_shape=(360, 360), return_validity = True, use_c_code = False, proj = 1)
        sky.ang2Pix([ra, dec], ra_dec_center=[180., -60], reso_arcmin=60., map_pixel_shape=(360, 360), return_validity = True, use_c_code = False, proj = 0)  # y_coord's get squished
        sky.ang2Pix([ra, dec], ra_dec_center=[0., -60], reso_arcmin=60., map_pixel_shape=(360, 360), return_validity = True, use_c_code = False, proj = 5)  # proj5 not great for all-sky maps
        sky.ang2Pix([360-0.01, -40], ra_dec_center=[15.02715, -35.03219], reso_arcmin=60., map_pixel_shape=(360, 360), return_validity = True, use_c_code = False, proj = 0)
    """
    # If the input is just a single RA, dec, convert it to an array of arrays, as this function expects.
    if np.size(np.asarray(ra_dec))==2:
        ra_dec = np.asarray([ [ra_dec[0]], [ra_dec[1]] ])
    
    # Convert the "proj" input to an index.
    proj = projNameToIndex(proj)

    assert type(bin_center_zero)==bool

    # Check if we can access the C code before launching into the calculation.
    if use_c_code:
        try:
            clibs = c_interface.initializeCLibs()
        except OSError,err:
            warnings.warn("Failed to load sptpol_cdevel python library. Reverting to Python ang2Pix.", RuntimeWarning)
            use_c_code=False

    ## babbloo
    float_type = np.float
    n_pixels = np.asarray(map_pixel_shape, dtype=float_type)  # used below if return_validity for both C and python
    
    if use_c_code:
        # Set up some of the variables needed for the C code.
        C = c_interface.C
        data_type = c_interface.typeindex if round else np.float64
        x_pixel_coord = np.empty(len(ra_dec[0]), dtype=data_type)
        y_pixel_coord = np.empty(len(ra_dec[1]), dtype=data_type)
        ra = np.ascontiguousarray(ra_dec[0], dtype=np.float64)
        dec = np.ascontiguousarray(ra_dec[1], dtype=np.float64)
        map_shape_deg = np.ascontiguousarray(map_pixel_shape) * reso_arcmin/60.
        mapinfo = c_interface.MAPINFO(proj=proj, map_center=ra_dec_center, map_shape=map_shape_deg, reso_arcmin=reso_arcmin)
        
        # Call the appropriate function (it differs by data type returned).
        if round:
            clibs.ang2Pix(C.byref(mapinfo), 
                          ra.ctypes.data_as(C.POINTER(C.c_double)),
                          dec.ctypes.data_as(C.POINTER(C.c_double)),
                          len(x_pixel_coord), 
                          x_pixel_coord.ctypes.data_as(C.POINTER(c_interface.typeindex)),
                          y_pixel_coord.ctypes.data_as(C.POINTER(c_interface.typeindex)))
        else:
            clibs.ang2Pix_nonrounded(C.byref(mapinfo), 
                                     ra.ctypes.data_as(C.POINTER(C.c_double)),
                                     dec.ctypes.data_as(C.POINTER(C.c_double)),
                                     len(x_pixel_coord), 
                                     x_pixel_coord.ctypes.data_as(C.POINTER(C.c_double)),
                                     y_pixel_coord.ctypes.data_as(C.POINTER(C.c_double)),
                                     bin_center_zero)
    else:
        # Unpack some variables into forms that we'll use in this function.
        ra = ra_dec[0].copy()  # copy because we unwrap in place below
        elev = -ra_dec[1]  # assume geographic south pole
        ra0, dec0 = ra_dec_center
        reso_rad = reso_arcmin/60. * DTOR

        # Wrap ra around ra0.
        # After this, all values of ra are within 180 degrees of ra0.
        # and below, (phi-phi0) will be in [-pi,+pi).
        # phi and phi0 always get subtracted -- 
        # they appear together as (phi-phi0) because of rotational invariance.
        # And we want to handle arbitrarily big maps centered around any ra0.
        ###tools.wrapAzAround(ra, ra0, modify_in_place=True)
	wrapAzAround(ra, ra0, modify_in_place=True)

        min_pos = -0.5*n_pixels*reso_rad  # minimum position (left or bottom edge) in radians
        min_pos_y, min_pos_x = min_pos  # gets overwritten in BICEP projection since x and y resolutions are different
        phi, phi0 = ra*DTOR, ra0*DTOR  # ra and ra_center in radians
        theta = (90 - elev)*DTOR  # theta polar angle in radians
        theta0 = (90 + dec0)*DTOR
        # First, get an x and y position in radians centered around zero
        # This is position relative to the map center
        if proj==0:
            x_pos = ne.evaluate("(phi-phi0)*cos(elev*DTOR)")
            y_pos = ne.evaluate("(elev + dec0)*DTOR")
        elif proj==1:
            x_pos = ne.evaluate("phi-phi0")
            y_pos = ne.evaluate("(elev + dec0)*DTOR")
        elif proj==2:
            x_pos = ne.evaluate("sin(phi-phi0)*sin(theta)")
            y_pos = ne.evaluate("-sin(theta)*cos(theta0)*cos(phi-phi0) + sin(theta0)*cos(theta)")
        elif proj==4:
            k = ne.evaluate("2/(1+cos(theta0)*cos(theta) + sin(theta0)*sin(theta)*cos(phi-phi0))")
            x_pos = ne.evaluate("k * sin(theta)*sin(phi-phi0)")
            y_pos = ne.evaluate("k * (-sin(theta)*cos(theta0)*cos(phi-phi0)+sin(theta0)*cos(theta))")
        elif proj==5:
            k = ne.evaluate("sqrt(2/(1+cos(theta0)*cos(theta)+sin(theta0)*sin(theta)*cos(phi-phi0)))")
            x_pos = ne.evaluate("k * sin(theta)*sin(phi-phi0)")
            y_pos = ne.evaluate("k * (-sin(theta)*cos(theta0)*cos(phi-phi0)+sin(theta0)*cos(theta))")
        elif proj==7:
            y_pos, x_pos = calculatePointingOffset([ra0, dec0], [ra, ra_dec[1]], offset_units='radians', as_azel=False)
        elif proj==8:
            x_pos, y_pos = rotateToZeroZero(np.array([ra, ra_dec[1]])*constants.RTOD, np.array([ra0,dec0])*constants.RTOD)
            x_pos*=constants.DTOR
            y_pos*=constants.DTOR
        elif proj==9:
            # BICEP default of 0.25 deg pixels, ra -55->+55, and dec -70->-45, gives a 236 by 100 pixel map.
            ny, nx = n_pixels
            #nx, ny = (236, 100) # for BICEP
            pixsize = reso_arcmin/60.  # pixsize = 0.25 for BICEP  ->  reso_arcmin=15.
            sy  = pixsize  # dec spacing in degrees
            asx = pixsize/cos(dec0*DTOR)  # approximate spacing in ra
            sx  = np.round(asx*nx)/nx  # ra spacing in degrees.  round so an integer number of degrees fit in the given number of pixels.
            #sx  = 0.4661016949152542  # for BICEP
            x_pos = ne.evaluate("(phi-phi0)*sy/sx")  # multiplying by the ratio sy/sx makes the ra pixels square-on-sky at dec0
            y_pos = ne.evaluate("(elev + dec0)*DTOR")
            min_pos_x = -0.5*nx*sy*DTOR  # OVERWRITE with new, scaled position.
                                         # Note this is scaled by sy! It is no longer the 'minimum ra in radians', but done this way
                                         # so we can subtract it from x_pos and divide by reso_rad to get pixel number.
        else:
            raise ValueError("I don't know what to do with proj "+str(proj))
        

        # Then shift the radians position and divide by the resolution
        if round:
            # without the floor function, int conversion just chops off fractional part.
            # slightly out-of-range negative values will get incorrectly put in pixel 0,
            #  and the validity checks below will never find them.
            x_pixel_coord = np.asarray(np.floor((x_pos - min_pos_x)/reso_rad), dtype=int)
            y_pixel_coord = np.asarray(np.floor((y_pos - min_pos_y)/reso_rad), dtype=int)
        else:
            x_pixel_coord = ne.evaluate("((x_pos - min_pos_x)/reso_rad) - bin_center_zero*0.5")
            y_pixel_coord = ne.evaluate("((y_pos - min_pos_y)/reso_rad) - bin_center_zero*0.5")
         
    pixel_coords = (y_pixel_coord, x_pixel_coord)
    return pixel_coords

    if return_validity:
        n_pix_y, n_pix_x = n_pixels
        pixel_y_is_good = ne.evaluate("(y_pixel_coord>=0) & (y_pixel_coord<n_pix_y)")
        pixel_x_is_good = ne.evaluate("(x_pixel_coord>=0) & (x_pixel_coord<n_pix_x)")
    
        return telescope.PointingSequence(pixel_coords, map_shape=map_pixel_shape, proj=proj,
                                          pixel_resolution_arcmin=reso_arcmin), (pixel_y_is_good, pixel_x_is_good)
    else:
        return telescope.PointingSequence(pixel_coords, map_shape=map_pixel_shape, proj=proj,
                                          pixel_resolution_arcmin=reso_arcmin)


