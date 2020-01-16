import numpy as np
mass = np.arange(1,5,0.2)
noise = np.arange(3,3.1,0.2)

for mm in mass:#nn in noise:
	outF = open('spartan_jobs/dummy_%s.sh'%(mm),'w')
	outF.write('#!/bin/bash \n')
	outF.write('#SBATCH --nodes=1\n')
	outF.write('#SBATCH --ntasks=1\n')
	outF.write('#SBATCH --cpus-per-task=1\n')
	outF.write('#SBATCH --mem 40000\n')
	outF.write('#SBATCH --time=60:00:00\n')
	outF.write('module load Anaconda2/4.2.0\n\n')
	for nn in noise:#mm in mass:
		line_to_be_written = 'python s1_2_sims_QE_tester_crossmaps.py -cmbrandomseedval 17819867 -clustype des -paramfile plnk_sz_grad_maps.txt -prefix test_mass_dependence -totalclus 1000 -use_data\
 0 -add_tSZ 1 -use_lgmca 0 -matched_filter 0 -no_lensing 1 -cluster_mass %s -noiselevels %s'%(mm,nn)
 		outF.write(line_to_be_written)
 		outF.write('\n')

outF.close()