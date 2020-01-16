import numpy as np
noise = np.arange(1)
mass = np.arange(1,7,0.1)
for nn in noise:
	outF = open('sys_checks.sh','w')
	outF.write('#!/bin/bash \n')
	outF.write('#SBATCH --nodes=1\n')
	outF.write('#SBATCH --ntasks=1\n')
	outF.write('#SBATCH --cpus-per-task=1\n')
	outF.write('#SBATCH --mem 40000\n')
	outF.write('#SBATCH --time=30:00:00\n')
	outF.write('module load Anaconda2/4.2.0\n\n')
	for mm in mass:
		line_to_be_written = 'python s1_2_sims_QE_tester_crossmaps.py -cmbrandomseedval %s -clustype des -paramfile plnk_sz_grad_maps.txt -prefix 2019038 -totalclus 1000 -use_data\
 0 -add_tSZ 0 -use_lgmca 0 -matched_filter 0 -cluster_mass 2.8 -noiselevels 3'%(np.random.randint(1e8))
 		outF.write(line_to_be_written)
 		outF.write('\n')

outF.close()