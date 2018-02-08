import os
import numpy as np
from shutil import copyfile
from settings import mofpath, basepath
from calc_swaps import update_calc
from error_handler import get_error_msgs, update_calc_after_errors, continue_mof
from magmom_handler import continue_magmoms, get_mag_indices
from compute_environ import choose_vasp_version
from ase.optimize import BFGSLineSearch
from writers import pprint
from metal_types import sblock_metals, poor_metals
from ase.io import read

def mof_run(mof,calc,cif_file,gpt_version,nprocs,calc_swaps):
#Get the optimized structure of the MOF
	copyfile(mofpath+cif_file,basepath+'working/'+cif_file)
	success = False
	calc = update_calc(calc,calc_swaps)
	mof.set_calculator(calc)
	try:
		mof.get_potential_energy()
		success = True
	except:
		old_error_len = 0
		refcode = cif_file.split('.cif')[0]
		if os.path.isfile('WAVECAR'):
			os.remove('WAVECAR')
		while True:
			old_swaps = calc_swaps
			errormsg = get_error_msgs('OUTCAR',refcode)
			calc, calc_swaps = update_calc_after_errors(calc,calc_swaps,errormsg)
			error_len = len(errormsg)
			if error_len == old_error_len:
				break
			print('Attempting to solve VASP issue(s): ',errormsg)
			mof = continue_mof()
			if 'vtst' in calc_swaps and 'vtst' not in old_swaps:
				choose_vasp_version(gpt_version,nprocs,calc_swaps,'vasp')
				static_swaps = ['nsw=0','ibrion=-1','iopt=0','isif=0','algo=All']
				calc = update_calc(calc,calc_swaps+static_swaps)
				choose_vasp_version(gpt_version,nprocs,calc_swaps+static_swaps,'vasp')
				mof.set_calculator(calc)
				mof.get_potential_energy()
				mof = read('OUTCAR')
				mof, abs_magmoms = continue_magmoms(mof,'INCAR')
			if mof != None:
				choose_vasp_version(gpt_version,nprocs,calc_swaps)
				calc = update_calc(calc,calc_swaps)
				mof.set_calculator(calc)
				try:
					mof.get_potential_energy()
				except:
					pass
			old_error_len = error_len
	if success == False:
		mof = None
	return mof, calc_swaps

def mof_bfgs_run(mof,calc,cif_file,calc_swaps,steps,fmax):
#Optimize with BFGSLineSearch
	copyfile(mofpath+cif_file,basepath+'working/'+cif_file)
	success = False
	calc = update_calc(calc,calc_swaps)
	mof.set_calculator(calc)
	try:
		dyn = BFGSLineSearch(mof,trajectory='opt.traj')
		dyn.run(fmax=fmax,steps=steps)
		success = True
	except:
		pprint('Original run failed. Trying to auto-solve issue.')
		old_error_len = 0
		refcode = cif_file.split('.cif')[0]
		if os.path.isfile('WAVECAR'):
			os.remove('WAVECAR')
		while True:
			errormsg = get_error_msgs('OUTCAR',refcode)
			calc, calc_swaps = update_calc_after_errors(calc,calc_swaps,errormsg)
			error_len = len(errormsg)
			if error_len == old_error_len:
				break
			print('VASP failed with error(s). Trying to auto-solve issue.',errormsg)
			mof = continue_mof()
			try:				
				mof.set_calculator(calc)
				dyn = BFGSLineSearch(mof,trajectory='opt.traj')
				dyn.run(fmax=fmax,steps=steps)
				success = True
			except:
				pass
			old_error_len = error_len
	if success == False:
		mof = None
	return mof, dyn, calc_swaps

def prep_next_run(acc_level,run_i,refcode,spin_level):
#Update counter and decide if next job should be skipped
	skip_spin2 = False
	success_path = basepath+'results/'+refcode+'/'+acc_level+'/'+spin_level
	incarpath = success_path+'/INCAR'
	outcarpath = success_path+'/OUTCAR'
	errorpath = basepath+'errors/'+refcode+'/'+acc_level+'/'+spin_level
	if os.path.exists(errorpath) == True:
		mof = None
	else:
		mof = read(outcarpath)
		if acc_level != 'scf_test':
			mof, abs_magmoms = continue_magmoms(mof,incarpath)
			mag_indices = get_mag_indices(mof)
			mag_nums = mof[mag_indices].get_atomic_numbers()
			if np.sum(abs_magmoms < 0.1) == len(abs_magmoms) or all(num in sblock_metals+poor_metals for num in mag_nums) == True:
				skip_spin2 = True
	run_i += 1
	return mof, run_i, skip_spin2