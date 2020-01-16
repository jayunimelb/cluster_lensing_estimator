import numpy as np, os, sys, glob

args = sys.argv[1:]
dx, no_lensing, prefix, clustype, totalclus, noofsims, paramfile, add_tSZ, max_jobs = args

no_lensing = int(no_lensing)
max_jobs = int(max_jobs)
noofsims = int(noofsims)
add_tSZ = int(add_tSZ)

#clustype = 'spt'#'whu'
spartan = 0
#noofsims = 25#500
#totalclus = 250#1000
#dx = 0.5
#add_tSZ = 0
use_data = 0
#paramfile = 'params_sptpol_tszfree.txt' #'params_for_S4.txt'
#prefix = 'tszfree'#'1am_beam'

#noise_arr = [1.0, 2.0, 3.0, 5.0, 7.0, 10.0, 15.0, 25.0, 30.0, 40.0, 50.0]
#noise_arr = [23.0]#, 5.0, 7.0, 15.0, 25.0, 40.0, 50.0]
#noise_arr = [55.0]
#noise_arr = [10.0]#, 'sptpol']
#noise_arr = [0.1, 0.25, 0.5, 0.7, 1., 3., 5., 7., 10.]
#noise_arr = [0.1,0.5,1.0]
#noise_arr = [50.0, 10.0]
noise_arr = [12.0]#7.0,16.0]
#max_jobs = 10
fno = 1
for expnoiselevel in noise_arr:

	expnoiselevel = str(expnoiselevel)

	if add_tSZ:
		tsz_str = 'tSZ'
	else:
		tsz_str = 'no_tSZ'

	if expnoiselevel.lower() == 'sptpol':
		noise_folder = 'sptpol'
	else:
		noise_folder = 'white_%s' %(expnoiselevel)


	if no_lensing:
		prefix = '%s_no_lensing' %(prefix)

	opfolder = 'data/sims/%s/%s/%s_%s_%s' %(clustype, prefix, dx, tsz_str, noise_folder)

	cmd_arr =[]
	searchstr = '%s/stacked*_%s_clusters*' %(opfolder, totalclus)
	inifiles = len(glob.glob(searchstr))

	totexec = 0
	for simcnt in range(noofsims):
		cmbrandomseedval = np.random.randint(1e6)
		cmd = 'python s1_sims_QE_tester.py %s %s %s %s %s %s %s %s %s %s &' %(totalclus, dx, cmbrandomseedval, expnoiselevel, add_tSZ, clustype, paramfile, prefix, no_lensing, use_data)
		cmd_arr.append(cmd)
		if len(cmd_arr)>=max_jobs or simcnt == noofsims-1:
			if not spartan:
				cmds = ' '.join(cmd_arr)
				totexec += len(cmd_arr)
				print totexec, len(cmd_arr), cmd_arr
				#quit()
				os.system(cmds)
			else:
				template = open('template_s1.sh','r')
				op_file = 's1_%s.sh' %(fno)
				opfile = open(op_file,'w')
				for lines in template:
				        opfile.writelines('%s\n' %lines.strip())

				for opline in cmd_arr:
				        opfile.writelines('%s\n' %opline.strip('&'))

				opfile.close()
				template.close()
				cmd = 'cd ..; sbatch batch_jobs/%s' %(op_file)
				os.system(cmd)
				fno +=1

			cmd_arr = []
			if not spartan:
				totfiles = len(glob.glob(searchstr)) - inifiles
				#print totfiles, totexec, searchstr
				while totfiles<totexec:
					totfiles = len(glob.glob(searchstr)) - inifiles
					pass
			#print totfiles, totexec, searchstr

		if len(glob.glob(searchstr)) >= noofsims:
			break
		
	
	
