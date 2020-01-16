import pickle, gzip, numpy as np, sys, glob, os
make_plot = 1
if str(os.getcwd()).find('sri')>-1:
	make_plot = 0

if make_plot: from pylab import *

folders = sys.argv[1:2]#:]
#minsnr, maxsnr = None, None#20.,40.#1000.#None, None
#minsnr, maxsnr = 26.,1000.#None, None
minsnr, maxsnr = int(sys.argv[2]), int(sys.argv[3])
use_weights = int(sys.argv[4])

def fn_normalize(data, arrmin = None, arrmax = None, Nmin = 0, Nmax = 1.):

        #data -= max(data*1.1) #to make all negative
        if arrmin == None:
                arrmin = min(data)
        if arrmax == None:
                arrmax = max(data)
        normed_data=(data-arrmin)/(arrmax-arrmin)
        normed_data=normed_data*Nmax+Nmin

        return normed_data

totweight = 0.
for fd in folders:
	fnames = sorted( glob.glob('%s/data*' %(fd)) )
	if fd.find('null_test')>-1:
		null_test = 1
	else:
		null_test = 0
	totclus = 0
	for fcnt, f in enumerate( fnames ):
		dic = pickle.load(gzip.open(f, 'rb'))
		print f, len(dic.keys())
		if len(dic.keys()) == 0: continue

		if use_weights: weightarr = np.asarray(dic.keys())[:,4]

		for kcnt, keyname in enumerate(dic.keys()):
			if not null_test:
				ra, dec, z, snr, weight = keyname
			else:
				ra, dec = keyname
				minsnr, maxsnr = None, None
			if minsnr <> None:
				if snr<minsnr: continue

			if maxsnr <> None:
				if snr>maxsnr: continue


			#kappa_qe = np.zeros((60, 60))
			#M, c, logL, kappa_qe, MF = dic[keyname]
			if len(dic[keyname]) == 5:
				M, c, logL, kappa_qe, MF = dic[keyname]
				#kappa_qe -= MF
			elif len(dic[keyname]) == 4:
				M, c, logL, kappa_qe = dic[keyname]
			else:
				M, c, logL = dic[keyname]
			#subplot(121);imshow(kappa_qe);colorbar();subplot(122);imshow(MF);colorbar();show();sys.exit()
			try:
				if totclus == 0: stacked_kappa_qe = np.copy(kappa_qe)
				else: stacked_kappa_qe+= np.copy(kappa_qe)
			except:
				stacked_kappa_qe = None
				print 'kappa_qe not stored'

			currlogL = logL

			if use_weights:
				weight = np.log( snr )#np.random.random()
				#weight = fn_normalize(weight, min(weightarr), max(weightarr))
				#weight = 1.
				currlogL = currlogL + weight
				totweight = totweight * weight

			#print currlogL;sys.exit()

			if totclus == 0:
				logl_arr = np.copy(currlogL)
			else:
				logl_arr += currlogL				

			totclus += 1

	if use_weights:
		logl_arr-=totweight

	if stacked_kappa_qe<>None:
		stacked_kappa_qe /= totclus

		nx, ny = stacked_kappa_qe.shape
		sn, dx = 10., 0.5
		boxsize = nx * dx
		clra, cldec = 0., 0.
		minval, maxval = clra-boxsize/2/60.,  clra+boxsize/2/60.
		ra = dec = np.linspace(minval, maxval, nx)
		e1, e2 = int(nx/2 - sn/dx/2), int(nx/2 + sn/dx/2)
		stacked_kappa_qe = stacked_kappa_qe[e1:e2, e1:e2]

		if make_plot: imshow(stacked_kappa_qe, interpolation = 'bicubic');colorbar();show()#;sys.exit()
		
	print 'Total clusters stacked = %s; SNR range = (%s,%s)' %(totclus, minsnr, maxsnr)
	#sys.exit()

	#print logl_arr;sys.exit()
	logl_arr_tmp = np.copy(logl_arr) - max(logl_arr)
	L = np.exp(logl_arr_tmp); L = L/max(L)

	tmp = f.strip('/').split('/')
	f1, f2 = '/'.join(tmp[0:-1]), tmp[-1]
	f2 = '_'.join(f2.split('_')[0:-2])
	res_file_name='%s/final_%s_%s_%s' %(f1, f2, minsnr, maxsnr)
	#sys.exit()
	if minsnr <> None and maxsnr <> None:
		lab = '%s: %s to %s' %(totclus, int(minsnr), int(maxsnr))
	else:
		lab = '%s: Full' %(totclus)
	opdic = {}
	opdic['label'] = lab
	opdic['Marr'] = M
	opdic['logl'] = logl_arr
	opdic['L'] = L
	opdic['stacked_kappa_qe'] = stacked_kappa_qe
	pickle.dump(opdic, gzip.open(res_file_name,'wb'), protocol = 2)
	if make_plot: plot(M,L)
if make_plot: show()
sys.exit()


######
"""
import pickle, gzip, numpy as np, sys, glob
from pylab import *

folders = sys.argv[1:]
minsnr, maxsnr = None, None#20.,40.#1000.#None, None
#minsnr, maxsnr = 26.,1000.#None, None

for fd in folders:
	fnames = glob.glob('%s/data*' %(fd))
	print fnames;sys.exit()
	if fd.find('null_test')>-1:
		null_test = 1
	else:
		null_test = 0
	for f in fnames:
		dic = pickle.load(gzip.open(f, 'rb'))
		for kcnt, keyname in enumerate(dic.keys()):
			if not null_test:
				ra, dec, z, snr, weight = keyname
			else:
				ra, dec = keyname
			if minsnr <> None:
				if snr<minsnr: continue

			if maxsnr <> None:
				if snr>maxsnr: continue

			M, c, logL = dic[keyname]
			
			try:
				logl_arr += logL
			except:
				logl_arr = np.copy(logL)

		logl_arr_tmp = np.copy(logl_arr) - max(logl_arr)
		L = np.exp(logl_arr_tmp); L = L/max(L)

		tmp = f.strip('/').split('/')
		f1, f2 = '/'.join(tmp[0:-1]), tmp[-1]
		res_file_name='%s/final_%s_%s_%s' %(f1, f2, minsnr, maxsnr)
		if minsnr <> None and maxsnr <> None:
			lab = '%s: %s to %s' %(f.split('_')[-2], int(minsnr), int(maxsnr))
		else:
			lab = '%s: Full' %(f.split('_')[-2])
		opdic = {}
		opdic['label'] = lab
		opdic['Marr'] = M
		opdic['logl'] = logl_arr
		opdic['L'] = L
		pickle.dump(opdic, gzip.open(res_file_name,'wb'), protocol = 2)
		plot(M,L)
show()
#sys.exit()
"""