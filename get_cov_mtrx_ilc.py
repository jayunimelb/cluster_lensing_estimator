"""


"""

import numpy as np
import pickle, gzip

cutouts_90 = pickle.load(gzip.open('data_sp/sptsz/min_variance/lensed_cmb_cutouts_500_clusters_90GHz_Arnaud_tSZ_fg_noise.pkl.gz','r'))
cutouts_150 = pickle.load(gzip.open('data_sp/sptsz/min_variance/lensed_cmb_cutouts_500_clusters_150GHz_Arnaud_tSZ_fg_noise.pkl.gz','r'))
cutouts_220 = pickle.load(gzip.open('data_sp/sptsz/min_variance/lensed_cmb_cutouts_500_clusters_220GHz_Arnaud_tSZ_fg_noise.pkl.gz','r'))

boxsize_extract = 10.
nx,ny = cutouts_90.shape[2:]
dx,dy = 0.5,0.5
sn = boxsize_extract
e1, e2 = int(nx/2. - sn/dx/2.), int(nx/2. + sn/dx/2.)

totalclus = len(cutouts_90)
data_90,data_150, data_220 = [],[],[]

for i in range(totalclus):
	extract_90 = (cutouts_90[i,0][e1:e2, e1:e2])/(1e6)
	extract_150 = (cutouts_150[i,0][e1:e2, e1:e2])/(1e6)
	extract_220 = (cutouts_220[i,0][e1:e2, e1:e2])/(1e6)
	if i >0:
		data_90 = np.append(data_90,extract_90.flatten())
		data_150 = np.append(data_150,extract_150.flatten())
		data_220 = np.append(data_220,extract_220.flatten())
	else:
		data_90 = extract_90.flatten()
		data_150 = extract_150.flatten()
		data_220 = extract_220.flatten()


from IPython import embed;embed()
