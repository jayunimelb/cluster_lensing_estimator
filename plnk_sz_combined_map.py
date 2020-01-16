import healpy as hp
import numpy as np

lmax = 10000
nside = 8192
fwhm_arcmin = 1.85 # final resolution of the maps

band = 150 # can be 95, 150 or 220

# I/O:
# Read combined data alms
# These alms have already had the 2d filter and 1.85 arcmin FWHM Gaussian applied
alm = hp.read_alm('alm_combined_data_%iGHz_ptsrcmasked_50mJy.fits'%(band,))
from IPython import embed;embed()
# Read 2d filters for each frequency
Mlm = np.fromfile('Mlm_lmax10000_%iGHz.bin'%(band,))

# Operations on the alms:
# Undo 2d filters from alms
alm[np.nonzero(Mlm)] /= Mlm[np.nonzero(Mlm)]

# Make a new 2d filter
# e.g. zero out a box of modes from ell = ell0 to ell1 and m = m0 to m1
Mlm = np.ones(hp.Alm.getsize(lmax))
ells, ms = hp.Alm.getlm(lmax)
ell0 = 1000
ell1 = lmax
m0 = 0
m1 = 300
for ell in xrange( ell0, ell1+1 ):
    Mlm[hp.Alm.getidx(lmax, ell, np.arange(m0, m1+1))] = 0.

# Apply 2d filter
alm *= Mlm

# Make a 1d filter (m-independent)
# e.g. a Gaussian beam with FWHM = fwhm_arcmin
fwhm_arcmin = 3.0
bl = hp.gauss_beam( fwhm_arcmin * np.pi/180./60., lmax=lmax )

# Apply 1d filter
ells, ms = hp.Alm.getlm(lmax)
alm *= bl[ells]

# lower lmax of alms
def lower_lmax(alm_in, lmax_out):
    lmax_in = hp.Alm.getlmax(alm_in.size)
    assert (lmax_in > lmax_out), "New lmax must be less than lmax of input alms"
    alm_out = np.zeros(hp.Alm.getsize(lmax_out), dtype=np.complex)
    for il in xrange(0, lmax_out+1):
        alm_out[hp.Alm.getidx(lmax_out, il, np.arange(0,il+1))] = alm_in[hp.Alm.getidx(lmax_in, il, np.arange(0,il+1))]
    return alm_out

# Usage:
new_lmax = 8000
alm = lower_lmax(alm, new_lmax)

# Make map from alms
tmap = hp.alm2map( alm, nside )

# Look at map interactively
hp.mollzoom(tmap)


# SPT-only products:
# Read filter transfer functions
tlm = hp.read_alm( 'tlm_150GHz_lmax10000.fits' )

# Read beam file
bl = np.fromfile( 'bl_gauss_fwhm1p75am_lmax10000.bin' )

# Create new beam:
fwhm_arcmin_new = 1.85 # FWHM of new beam
bl_new = hp.gauss_beam( fwhm_arcmin_new / 60. * np.pi/180., lmax=10000)
