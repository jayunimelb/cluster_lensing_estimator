def fn_rad_profile(xarr,yarr,zarr,ZARR,binsize=0.01,arcmins=True):

        if len(xarr)==0:
                # Calculate the indices from the image
                y,x = np.indices(ZARR.shape)
                center = np.array([(x.max()-x.min())/2.0, (y.max()-y.min())/2.0])
                r = np.hypot(x - center[0], y - center[1])
        else:
                #converting into polar coordinates
                r=np.asarray((xarr**2+yarr**2)**0.5)

        #r*=60. #degrees to arcmins

        #bin in radial direction
        maxbin=1. #degrees
	binsize = 0.5/60.
        if arcmins:
                binsize*=60. #arcmins
                maxbin*=60. #arcmins 

        binarr=np.arange(0,maxbin,binsize)

        radprf=np.zeros((0,3))
        b=0
        hit_count=[]
        for bin in binarr:
                ind=np.where((r>=bin) & (r<bin+binsize))[0]

                radprf.resize((b+1,3))
                radprf[b,0]=bin
                if len(ind)>0:
                        radprf[b,1]=np.average(zarr[ind])
                        radprf[b,2]=np.std(zarr[ind])

                hit_count.append(len(np.where(abs(zarr[ind])>0.)[0]))
                b+=1

        hit_count=np.asarray(hit_count)
        std_mean=np.sum(radprf[:,2]*hit_count)/np.sum(hit_count)

	errval = np.zeros(len(hit_count))
	errval[hit_count>0] =std_mean/(hit_count[hit_count>0])**0.5
        radprf[:,2]=errval

        return radprf


def fn_format_axis(ax,fx,fy,maxxloc=5):
        for label in ax.get_xticklabels():
                label.set_fontsize(fx)
        for label in ax.get_yticklabels():
                label.set_fontsize(fy)

        ax.xaxis.set_major_locator(MaxNLocator(maxxloc))

        return ax

import numpy as np, pickle, gzip, glob, sys, os
from pylab import *
import modules.radialProfile as radialProfile
from matplotlib.font_manager import fontManager, FontProperties


make_rad_plot = 1
isolated_clus = 0

########################################################
########################################################
boxsize = 100
dx = float(sys.argv[1])
#beamsize = float(sys.argv[2])
#### dx, dy = .25, .25 #arcmins
dy = dx

nx, ny = int(boxsize/dx), int(boxsize/dy)
boxsize = int(nx * dx)

minval, maxval = -boxsize/2/60.,  boxsize/2/60.
raarr = np.linspace(minval, maxval, nx)
decarr = np.linspace(minval, maxval, ny)
RA, DEC = np.meshgrid(raarr,decarr)

########################################################

'''
if len(sys.argv) == 3:
	#ipfolder = 'detections'
	#ipfiles = sorted(glob.glob('%s/*/*.pkl.gz' %(ipfolder)))
	ipfolder = sys.argv[2]
	ipfiles = sorted(glob.glob('%s/*.pkl.gz' %(ipfolder)))
elif len(sys.argv) == 4:
	ipfolder = sys.argv[3] 
	ipfiles = sorted(glob.glob('%s/*%s*.pkl.gz*' %(ipfolder,beamsize)))
'''
ipfolder = 'detections/boxsize_%s_%samresol' %(boxsize,dx)
ipfiles = sorted(glob.glob('%s/*.pkl.gz*' %(ipfolder)))

if not make_rad_plot:
	if len(ipfiles)>1:
		totrows, totcols = 1, 1
	else:
		totrows, totcols = 1, 1

	subplots_adjust(hspace=0.3)

else:
	ax = subplot(111)#, yscale='log')#, xscale = 'log')


reqd_box, dx  = 14., 0.5 #arcmins
minval,maxval = -reqd_box/2., reqd_box/2. 
reqd_box_pix = reqd_box/dx

kappa_true_FILE = '%s/kappa_m_2e+14_z_0.6.pkl.gz' %(ipfolder)
kappa_true = pickle.load(gzip.open(kappa_true_FILE,'r'))

x1,x2 = int(len(kappa_true)/2-reqd_box_pix/2), int(len(kappa_true)/2+reqd_box_pix/2)
y1,y2 = int(len(kappa_true)/2-reqd_box_pix/2), int(len(kappa_true)/2+reqd_box_pix/2)

kappa_true = kappa_true[y1:y2,x1:x2]
#imshow(kappa_true);colorbar();show();quit()
#kappa_true /= max(kappa_true.ravel())


###kappa_true =np.exp(-(((X)/5.)**2+((Y)/5.)**2)/2.); imshow(kappa_true);show();quit()

azradprof_true, azradprof_std = radialProfile.azimuthalAverage(kappa_true)

'''
x = y = np.arange(minval,maxval,dx)
X, Y = np.meshgrid(x,y)
xarr=X.ravel()
yarr=Y.ravel()
zarr=kappa_true.ravel()

azradprof = fn_rad_profile(xarr,yarr,zarr,kappa_true,binsize=0.01,arcmins=True)
#azradprof/=max(azradprof)
'''

theta_arr = np.arange(len(azradprof_true)) * dx
lab = 'true kappa: Mass = 2e14; z = 0.6'
if make_rad_plot:
	#ax = subplot(111, xscale = 'log')
	plot(theta_arr,azradprof_true, 'k', lw = 2., label = lab , ls = 'solid')
	#plot(azradprof[:,0],azradprof[:,1], 'k', lw = 3., label = lab , ls = 'solid')
	#show();quit()
sbpl = 1

#colorarr = ['k','r','orange','lime','brown','m','gold']
colorarrdic = {1.2:'r',2.5:'g',3.5:'b'}
fill_colorarrdic = {1.2:'Lightcoral',2.5:'lightgreen',3.5:'lightskyblue'}
alphadic = {1.2:.8,2.5:1.,3.5:1.}
#colorarrdic = {3.5:'b'}
for beam in sorted(colorarrdic.keys(), reverse=True):
	azradprof_for_beam = []
	for fcnt,ff in enumerate(ipfiles):

		if ff.find('KAPPA')>-1:
			continue

		if isolated_clus:
			if ff.find('non')>-1:
				continue
		else:
			if ff.find('non')==-1:
				continue


		if ff.find('stacked_%s_' %(beam))==-1:
			continue

		#print ff, fcnt
		dic = pickle.load(gzip.open(ff,'r'))
		tit = ff[ff.find('/')+1:ff.find('.pkl')]
		lab = ff[ff.rfind('/')+1:ff.find('.pkl')]
		tit = tit.replace('/stacked','\n/stacked')
	
		kappa = dic.values()[0] #using just the first key/values of the dic for now
		#kappa/=max(kappa.ravel())

		kappa = kappa[y1:y2,x1:x2]
		RA = RA[y1:y2,x1:x2]
		DEC = DEC[y1:y2,x1:x2]

		if not make_rad_plot:
			ax = subplot(totrows,totcols,sbpl)
			imshow(kappa,extent=[minval,maxval,minval,maxval])#;colorbar()
			grid(True,ls='solid')
			title(tit,fontsize=6)
			#show();quit()
			ax = fn_format_axis(ax,6,6)
			sbpl+=1
			continue

		'''
		colorval = colorarrdic[beam]
		zarr=kappa.ravel()
		azradprof = fn_rad_profile(xarr,yarr,zarr,kappa,binsize=0.01,arcmins=True)
		#plot(azradprof[:,0],azradprof[:,1], 'k', lw = 3., label = lab , ls = 'solid')
		errorbar(azradprof[:,0], azradprof[:,1], yerr = [azradprof[:,2], azradprof[:,2]], marker = 'o', label = lab, color = colorval, capsize=0., mec = 'None' , ls = 'None')
		break
		'''

		azradprof, azradprof_std = radialProfile.azimuthalAverage(kappa)
		azradprof_for_beam.append(azradprof)

		'''
		#azradprof/=max(azradprof)
		#colorval = colorarr[fcnt]
		startpos = ff.find('stacked_')+8
		beamval = float(ff[startpos:startpos+3])
		colorval = colorarrdic[beamval]
		plot(theta_arr, azradprof, marker = 'o', label = lab, color = colorval, mec = 'None')# , ls = 'None')
		'''

		'''
		xarr, yarr = RA.ravel(), DEC.ravel()
		zarr = kappa.ravel()
		ZARR = kappa
		#azradprof = fn_rad_profile(xarr,yarr,zarr,ZARR,binsize=.1,arcmins=true)
		#plot(azradprof[:,0], azradprof[:,1], marker = 'o')#, ls = 'None')
		'''

		#break

	colorval = colorarrdic[beam]
	fill_colorval = fill_colorarrdic[beam]
	alphaval = alphadic[beam]
	
	azradprof_for_beam = np.asarray(azradprof_for_beam)
	#plot(theta_arr, np.mean(azradprof_for_beam,axis = 0), marker = 'o', label = lab, color = colorval, mec = 'None')# , ls = 'None')
	yerr = np.std(azradprof_for_beam,axis = 0)

	plot(theta_arr, np.mean(azradprof_for_beam,axis = 0), label = lab, color = colorval, mec = 'None' , ls = 'dashed', lw = 3.)
        fill_between(theta_arr,azradprof_true - yerr, azradprof_true + yerr,color=fill_colorval, alpha = alphaval)#,edgecolor = colorval)#, alpha = alphaval)#, label = lab)

	'''
	if beam == 3.5:
	        fill_between(theta_arr,azradprof_true - yerr, azradprof_true + yerr,color='None',hatch = '///',edgecolor = colorval)#, label = lab)
	elif beam == 2.5:
	        fill_between(theta_arr,azradprof_true - yerr, azradprof_true + yerr,color='None',alpha=0.5,hatch = '+',edgecolor = colorval)#, label = lab)
	elif beam == 1.2:
	        fill_between(theta_arr,azradprof_true - yerr, azradprof_true + yerr,color='None',edgecolor = colorval)#, label = lab)
	'''

	'''
	errorbar(theta_arr[1:], np.mean(azradprof_for_beam,axis = 0)[1:], yerr = [yerr[1:], yerr[1:]], marker = 'o', label = lab, color = colorval, mec = 'None' , ls = 'None')
        fill_between(theta_arr,np.mean(azradprof_for_beam,axis = 0) - yerr, np.mean(azradprof_for_beam,axis = 0) + yerr,color=colorval,alpha=0.5)#,edgecolor = colorval, label = lab)
	'''

	'''
	errorbar(theta_arr[1:], azradprof_true[1:], yerr = [yerr[1:], yerr[1:]], marker = 'o', label = lab, color = colorval, mec = 'None' , ls = 'None')
        fill_between(theta_arr,azradprof_true - yerr, azradprof_true + yerr,color=colorval,alpha=0.5)#,edgecolor = colorval, label = lab)
	'''
        
#show();quit()

if make_rad_plot: 
	legend(loc = 'Best', fancybox = 1,prop=FontProperties(size=12))
	xlim(0.,6.5)

	xlabel('Radial distance [arcmins] -->', fontsize = 14)
	ylabel('kappa -->', fontsize = 14)
	
show();quit()

