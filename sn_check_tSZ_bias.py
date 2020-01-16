import numpy as np, pickle, gzip, glob, os
from pylab import *

opfolder = 'outputs/cmass/'

noisearr = [0.1,10.]
tszarr = [0,1]
dx = dy = 0.5
totalclus = 4000

sbpl = 1
for ncnt, noise in enumerate(noisearr):
	for tszcnt, tsz in enumerate(tszarr):
		if tsz == 0:
			extra = 'no_tSZ'
		else:
			extra = 'tSZ'
		foldername = '%s/%s_%s_white_%.1f' %(opfolder, dx, extra, noise)
		fname = glob.glob('%s/*_%s_clusters*.pkl.gz' %(foldername,totalclus))[0]
		kappadic = pickle.load(gzip.open(fname,'rb'))
		stacked_kappa_qe = kappadic['stacked_kappa_qe']

		nx, ny = stacked_kappa_qe.shape
		snx = sny = 20
		ex1, ex2 = (nx/2) - int(snx/dx/2), (nx/2) + int(snx/dx/2)
		ey1, ey2 = (ny/2) - int(sny/dy/2), (ny/2) + int(sny/dy/2)
		stacked_kappa_qe = stacked_kappa_qe[ex1:ex2, ey1: ey2]
		

		ax = subplot(2,2,sbpl)
		imshow(stacked_kappa_qe);colorbar();title('Noise = %s uK\'; Clusters = %s; tSZ = %s' %(noise, totalclus, tsz))
		sbpl += 1
		#show();quit()

show()
				
