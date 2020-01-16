import pickle, gzip, numpy as np, sys
from pylab import *

fname = sys.argv[1]

RES_DIC = pickle.load(gzip.open(fname,'rb'))
keys = np.asarray( sorted(RES_DIC.keys()) )
Marr = keys[:,0]
carr = keys[:,1]

M,C = np.meshgrid(Marr,carr)
logLarr = np.zeros( M.shape )
for mcnt, MM in enumerate(Marr):
	for ccnt, cc in enumerate(carr):
		keyname = (round(MM,3),round(cc,2))
		logLarr[ccnt,mcnt] = RES_DIC[keyname]

#imshow(logLarr);colorbar();show();quit()
imshow(logLarr, extent=[np.min(Marr),np.max(Marr),np.min(carr),np.max(carr)], aspect = 'auto');colorbar()
#pcolor(M, C, logLarr);colorbar();xlim( min(Marr), max(Marr) );ylim( min(carr), max(carr) )
show()

