import os
from ase.io import read
import numpy as np
from shutil import copyfile

phase2_results_path = '/projects/p30148/vasp_jobs/MOFs/oxidized_oms/results/'
final_folder_name = 'final_spe'
phase2_cleaned_results_path = '/projects/p30148/vasp_jobs/MOFs/oxidized_oms/results_cleaned/'

if not os.path.isdir(phase2_cleaned_results_path):
	os.makedirs(phase2_cleaned_results_path)
for folder in os.listdir(phase2_results_path):
	if not os.path.isdir(phase2_results_path+folder):
		continue
	mof_path = phase2_results_path+folder+'/'+final_folder_name+'/'
	if not os.path.isdir(mof_path):
		continue
	E = np.inf
	for subfolder in os.listdir(mof_path):
		mof_subpath = mof_path+subfolder+'/'
		mof_name = folder+'_'+subfolder
		mof = read(mof_subpath+'OUTCAR')
		E_temp = mof.get_potential_energy()
		if E_temp < E:
			E = E_temp
			best_mof = mof_subpath
			best_mof_name = mof_name
	copyfile(best_mof+'CONTCAR',phase2_cleaned_results_path+'POSCAR_'+best_mof_name)
