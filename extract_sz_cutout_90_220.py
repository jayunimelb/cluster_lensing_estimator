"""

To extract sz data cutouts

"""
import numpy as np, pickle,gzip
from scipy.io import readsav
import sys, modules
import sky_local
sims = modules.scl_cmb.simulations()
sz_dirctory = '/data15/tcrawfor/run1_clusters_fullsurvey'
sz_fields = readsav('%s/sptsz_field_info.sav'%(sz_dirctory))['field_info']
op_fname = 'data_sp/sptsz/data/spt/sz_cutouts.pkl.gz'
op_fnam_90 = 'data_sp/sptsz/data/spt/sz_cutouts_90GHz.pkl.gz'
op_fnam_220 = 'data_sp/sptsz/data/spt/sz_cutouts_220GHz.pkl.gz'
catalog_file = 'data_sp/sptsz/data/spt/catalog.pkl.gz'
catalog = pickle.load(gzip.open(catalog_file))
RA, DEC, XI, zvals,field_names = catalog['RA'], catalog['DEC'], catalog['XI'],catalog['REDSHIFT'],catalog['field_name']

field_exceptions = ['ra5h30dec-55','ra23h30dec-55']
cutouts = {}
cutouts_90,cutouts_220 ={},{}
boxsize = 100 # arcmins
sz_cutouts_90 = {}
sz_cutouts_220 = {}
desired_reso = 0.5 # arcmins
use_mask = 1

mask_file_name = 'SPTSZ_field_masks_psource.pkl.gz'
mask_file_name = 'data_sp/sptsz/data/spt/SPTSZ_field_masks_psource.pkl.gz'
mask_fld_dic = pickle.load(gzip.open(mask_file_name))['masks']
for i,(ra,dec) in enumerate(zip(RA,DEC)):
    sub_dir = field_names[i]
    data_file_name = '%s/%s/coadd_all_02Oct12.sav'%(sz_dirctory,sub_dir)
    if zvals[i] ==0 or None:
        continue
    mask = mask_fld_dic[field_names[i]]
    
    if field_names[i] in field_exceptions:
        sub_dir = '%s_2year'%(sub_dir)
        data_file_name = '%s/%s/coadd_all_02Oct.sav'%(sz_dirctory,sub_dir)

    
    try:
        map_data = readsav(data_file_name) 
    except:
        data_file_name = '%s/%s/coadd_all_02Oct.sav'%(sz_dirctory,sub_dir)
        map_data = readsav(data_file_name)

    map_150 = map_data['maps'][1]
    map_90 = map_data['maps'][0]
    map_220 = map_data['maps'][2]
    # multiply mask
    if use_mask:
        mask = mask.T
    map_150 = map_150*mask
    map_90, map_220  = map_90*mask, map_220*mask
    map_info= map_data['mapinfo']
    map_reso = map_info[0][4]
    ra0, dec0 = map_info[0][5],map_info[0][6]
    ra_dec_center = (ra0, dec0)
    ra_dec = (ra,dec)
#    from IPython import embed;embed()
    x_coord, y_coord = sky_local.ang2Pix(ra_dec, ra_dec_center, map_reso, map_150.shape, proj=0)
    
    cutout_size = boxsize/map_reso
    cutout = map_150[x_coord-(cutout_size/2.):x_coord+(cutout_size/2.),y_coord-(cutout_size/2.):y_coord+(cutout_size/2.)]
    cutout_90 = map_90[x_coord-(cutout_size/2.):x_coord+(cutout_size/2.),y_coord-(cutout_size/2.):y_coord+(cutout_size/2.)]
    cutout_220 =map_220[x_coord-(cutout_size/2.):x_coord+(cutout_size/2.),y_coord-(cutout_size/2.):y_coord+(cutout_size/2.)]
    
    cutout_ds = sims.downsample_map(cutout,int(desired_reso/map_reso))
    cutout_ds_90 = sims.downsample_map(cutout_90,int(desired_reso/map_reso))
    cutout_ds_220 = sims.downsample_map(cutout_220,int(desired_reso/map_reso))
    
    keyname = ra,dec,zvals[i],XI[i],1,field_names[i]
    cutouts[keyname] = np.asarray([cutout_ds])
    cutouts_90[keyname] = np.asarray([cutout_ds_90])
    cutouts_220[keyname]= np.asarray([cutout_ds_220])
sz_cutouts_90['cutouts'] = cutouts_90
sz_cutouts_220['cutouts'] = cutouts_220

sz_cutouts_90['resol_arcmins'] =sz_cutouts_220['resol_arcmins'] = desired_reso
sz_cutouts_90['map_source']= sz_cutouts_220['map_source'] = '/data15/tcrawfor/run1_clusters_fullsurvey/'
sz_cutouts_90['cluster_catalogue'] = sz_cutouts_220['cluster_catalogue'] = 'data_sp/sptsz/list_for_sanjay.txt'
sz_cutouts_90['cutout_size_arcmins'] =sz_cutouts_220['cutout_size_arcmins'] = boxsize
'''
if use_mask:
    op_fname = 'data_sp/sptsz/data/spt/ptsrc_masked_sz_cutouts.pkl.gz'
'''  
pickle.dump(sz_cutouts_90,gzip.open(op_fname_90,'w'))
pickle.dump(sz_cutouts_220,gzip.open(op_fname_90,'w'))
