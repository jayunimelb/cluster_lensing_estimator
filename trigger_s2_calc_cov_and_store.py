import numpy as np, os, sys, time, pickle, gzip

fname = sys.argv[1]
noofprocs = int(sys.argv[2])
#wait_time = float(sys.argv[3]) * 60.

noofsims = 500
#minM,maxM,delM = 1.,6.,0.25
#minc,maxc,delc =  2.6,3.6,0.2

minM,maxM,delM = 0.,7.,1.
minc,maxc,delc =  2.6,3.4,0.2

M_arr = np.arange(minM,maxM+delM,delM) * 1e14
c_arr = np.arange(minc,maxc+delc,delc)
c_arr = [3.0]

totalclus = int(fname.split('/')[-1].split('_')[4])
covdic_folder = '/'.join(fname.split('/')[0:-1])

#totalclus, noofsims = 10, 10
cmd_arr = []
for MM in M_arr:
	for cc in c_arr:

		covdic_name = '%s/COVS_%sclusters_%ssims_%smass_%sconc' %(covdic_folder, totalclus, noofsims, MM, cc)
		if os.path.exists(covdic_name):
			print '\n\n %s already exists \n\n' %(covdic_name)
			continue

		cmd = 'python s2_calc_cov_and_store.py %s %s %s %s &' %(fname,MM,cc,noofsims)
		cmd_arr.append(cmd)

		print len(cmd_arr)

		if len(cmd_arr)>= noofprocs or ( MM == M_arr[-1] and cc == c_arr[-1]):
			cmds = ' '.join(cmd_arr)
			print cmds#;quit()
			os.system(cmds)
			cmd_arr = []

			while not os.path.exists(covdic_name):
				pass
			'''
			start_time = time.time()
			while time.time() - start_time < wait_time:
				pass
			'''

