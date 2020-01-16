import numpy as np
import pickle,gzip
import sys
from pylab import *
data_file = sys.argv[1]

redmapper_data = pickle.load(gzip.open(data_file))

z_cut = 0.5
lmda_cut = 30
kappa_1_1 =[];kappa_2_1=[];kappa_1_2=[];kappa_2_2=[]
CLUS_IDENTIFIER = redmapper_data['CLUS_IDENTIFIER']
mean_field = pickle.load(gzip.open('/Users/sanjaykumarp/transfer/spt/data_stacked_kappa_21829_clusters.pkl.gz_20180123_212904_redmapper_randoms'))['stacked_kappa_qe']

for j,clus_iden in enumerate(redmapper_data['CLUS_IDENTIFIER']):
	if (clus_iden[2]<=z_cut):
		if clus_iden[3] >=20 and clus_iden[3]<lmda_cut:
			kappa_1_1.append(redmapper_data['kappa_qe_arr'][j])
		else:
			kappa_1_2.append(redmapper_data['kappa_qe_arr'][j])
	if clus_iden[2]>z_cut:
		if clus_iden[3] >=20 and clus_iden[3]<lmda_cut:
			kappa_2_1.append(redmapper_data['kappa_qe_arr'][j])
		else:
			kappa_2_2.append(redmapper_data['kappa_qe_arr'][j])

kappa_1_1 -= mean_field
kappa_1_2 -= mean_field
kappa_2_1 -= mean_field
kappa_2_2 -= mean_field
print len(kappa_1_1)
print len(kappa_1_2)
print len(kappa_2_1)
print len(kappa_2_2)
kappa_1_1 = np.asarray(kappa_1_1);kappa_1_2 = np.asarray(kappa_1_2); kappa_2_1 = np.asarray(kappa_2_1); kappa_2_2 = np.asarray(kappa_2_2)
subplot(221);imshow(np.mean(kappa_1_1,axis =0)[75:125,75:125]);colorbar()
subplot(222);imshow(np.mean(kappa_1_2,axis = 0)[75:125,75:125]);colorbar()
subplot(223);imshow(np.mean(kappa_2_1,axis =0)[75:125,75:125]);colorbar()
subplot(224);imshow(np.mean(kappa_2_2,axis =0)[75:125,75:125]);colorbar()
show()


