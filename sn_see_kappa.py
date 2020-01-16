def gaussian(height, center_x, center_y, width_x, width_y):
    """Returns a gaussian function with the given parameters"""
    width_x = float(width_x)
    width_y = float(width_y)
    return lambda x,y: height*exp(
                -(((center_x-x)/width_x)**2+((center_y-y)/width_y)**2)/2)

def moments(data):
    """Returns (height, x, y, width_x, width_y)
    the gaussian parameters of a 2D distribution by calculating its
    moments """
    total = data.sum()
    X, Y = indices(data.shape)
    x = (X*data).sum()/total
    y = (Y*data).sum()/total
    col = data[:, int(y)]
    width_x = sqrt(abs((arange(col.size)-y)**2*col).sum()/col.sum())
    row = data[int(x), :]
    width_y = sqrt(abs((arange(row.size)-x)**2*row).sum()/row.sum())
    height = data.max()
    return height, x, y, width_x, width_y

def fitgaussian(data):
    """Returns (height, x, y, width_x, width_y)
    the gaussian parameters of a 2D distribution found by a fit"""
    params = moments(data)
    errorfunction = lambda p: ravel(gaussian(*p)(*indices(data.shape)) -
                                 data)
    p, success = optimize.leastsq(errorfunction, params)
    return p

def fitting_func(p, p0, X, Y, MAP, lbounds = None, ubounds = None, fixed = None, return_fit = 0):
	'''
	if hasattr(fixed, '__len__'):
		p[fixed] = p0[fixed]

	if hasattr(lbounds, '__len__'):
		linds = abs(p)<abs(lbounds)
		p[linds] = lbounds[linds]

	if hasattr(ubounds, '__len__'):
		uinds = abs(p)>abs(ubounds)
		p[uinds] = ubounds[uinds]
	'''

	def fn_gaussian(p, xp, yp):
		#rota = np.radians(p[6])
		#xp = xp * np.cos(rota) - yp * np.sin(rota)
		#yp = xp * np.sin(rota) + yp * np.cos(rota)
		g = p[0]+p[1]*np.exp( -(((p[2]-xp)/p[4])**2+ ((p[3]-yp)/p[5])**2)/2.)

		return g

	if not return_fit:
		return np.ravel(fn_gaussian(p, X, Y) - MAP)
	else:
		return fn_gaussian(p, X, Y)



import numpy as np, pickle, gzip, glob, os
import scipy.optimize as optimize
from pylab import *
import modules
from matplotlib import rc;rc('text', usetex=True);rc('font', weight='bold');matplotlib.rcParams['text.latex.preamble'] = [r'\boldmath']

sims = modules.scl_cmb.simulations()

folderpath = sys.argv[1]
plfolder = folderpath.replace('data/sims/','plots/')

plfoldersplit = plfolder.split('/')

for plcnt, pl in enumerate(plfoldersplit):
	if plcnt == 0:
		plfoldername = pl
	else:
		plfoldername = '%s/%s' %(plfoldername, pl)
	if not os.path.exists(plfoldername): os.system('mkdir %s' %(plfoldername))

files = glob.glob('%s/*stacked*' %(folderpath))
#files = files[0:5]

clf()
sn = 10.
for fcnt, fname in enumerate(files):
	kappadic = pickle.load(gzip.open(fname,'rb'))
	PSRC_DIST =kappadic['PSRC_DIST']
	#kappa_qe_arr = kappadic['kappa_qe_arr']
	kappa_qe = kappadic['stacked_kappa_qe']
	dx = dy = kappadic['reso_arcmin']
	add_noise = kappadic['add_noise']
	#expnoiselevel = [kappadic['expnoiselevel'][0]]
	#### boxsize = kappadic['boxsize']
	nx, ny = kappa_qe.shape
	boxsize = nx * dx

	clra, cldec = 0., 0.

	minval, maxval = clra-boxsize/2/60.,  clra+boxsize/2/60.
	#minval, maxval = 0,  clra+boxsize/60.
	ra = dec = np.linspace(minval, maxval, nx)

	RA, DEC = np.meshgrid(ra,dec)

	e1, e2 = int(nx/2 - sn/dx/2), int(nx/2 + sn/dx/2)
	kappa_qe = kappa_qe[e1:e2, e1:e2]
	RA = RA[e1:e2, e1: e2]# * 60.
	DEC = DEC[e1:e2, e1: e2]# * 60.

	RADEC = [RA, DEC]
	radprofile = sims.fn_radial_profile(kappa_qe, RADEC, bin_size = 0.5, minbin = 0., maxbin = sn)
	#radprofile = sims.fn_radial_profile(kappa_qe, RADEC)
	

	#fit Gaussian to the map to get the centroid
	height, amp = 0,np.max(kappa_qe)
	x_cen, y_cen = 0.,0.
	wx, wy = 1., 1. #arcmins
	rot = 0. #degrees
	p0 = [height, amp, x_cen, y_cen, wx, wy, rot]
	p1, pcov, infodict, errmsg, success = optimize.leastsq(fitting_func, p0[:], args=(p0, RA*60., DEC*60., kappa_qe), full_output=1)
	xcen, ycen = p1[2], p1[3]

	MAP_FIT = fitting_func(p1,p1,RA*60.,DEC*60.,kappa_qe,return_fit = 1)
	#sys.exit()



	'''
	smooth_beam_arcmins = 1.
	fwhm = np.radians(smooth_beam_arcmins/60.)
	mapparams = [nx, ny, dx, dy]
	lx, ly = sims.get_lxly(mapparams)
	EL = (lx ** 2. + ly ** 2.)**0.5	
	LPF = sims.gauss_beam(fwhm, EL)[0]
	kappa_qe = np.fft.ifft2( np.fft.fft2(kappa_qe) * LPF).real
	'''
	fname = fname.split('/')[-1]
	"""
	new_nx, new_ny = kappa_qe.shape
	minval, maxval = None, None
	clf()
	subplot(111); imshow(kappa_qe, vmin = minval, vmax = maxval, extent = [np.min(RA)*60., np.max(RA)*60., np.min(DEC)*60., np.max(DEC)*60.], cmap = cm.jet, interpolation = 'bicubic');colorbar();grid(1, ls = 'solid')#;plot(0.,0.,'k+', ms = 10.)#;xlim(0, new_nx); ylim(0, new_ny);
	title(fname, fontsize = 8)
	#subplot(122)#, xscale = 'log'); #plot(radprofile[:,0], radprofile[:,1], 'ro-');
	#errorbar(radprofile[:,0],radprofile[:,1], yerr = [radprofile[:,2],radprofile[:,2]], color='r', marker = 'o')#, lw = 2., ls = 'solid')
	#xlim(0, 10.)
	show();#sys.exit()
	"""
	clf()
	#minval, maxval = None, None#-0.25,0.25#None, None#-0.4,0.1		
	minval, maxval = None, None#-0.4,0.1
	imshow(kappa_qe, origin='lower', extent = [np.min(RA)*60., np.max(RA)*60., np.min(DEC)*60., np.max(DEC)*60.], vmin = minval, vmax = maxval, interpolation = 'bicubic', cmap = cm.jet);colorbar()
	#contour(MAP_FIT, extent = [np.min(RA)*60., np.max(RA)*60., np.min(DEC)*60., np.max(DEC)*60.], cmap = cm.gray)
	plot(xcen, ycen, 'k+', ms = 8.); text(xcen*2., ycen*2., r'\textbf{(%.3f, %.3f)}' %(xcen, ycen), fontsize = 14, color = 'k')
	axhline(0., color = 'k');axvline(0., color = 'k')
	totalclus = int(fname.split('_')[3])
	extra= '_'.join(fname.split('_')[5:7])
	#print plfoldername;sys.exit()
	pl_name = '%s/%05d_clusters_%s.png' %(plfoldername, totalclus, extra)
	title(r'\textbf{Total clusters = %s}' %(totalclus))
	xlabel(r'\textbf{ARCMIN}', fontsize =14)
	ylabel(r'\textbf{ARCMIN}', fontsize =14)
	#show();sys.exit()
	savefig(pl_name, bbox_inches = 'tight', pad_inches=0.1);close()
	



sys.exit()

