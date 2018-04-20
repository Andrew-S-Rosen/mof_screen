import os
import numpy as np
from shutil import copyfile
from ase.io import read
from ase.optimize import BFGSLineSearch
from pymofscreen.calc_swaps import update_calc
from pymofscreen.error_handler import get_niter, get_error_msgs, update_calc_after_errors, continue_mof
from pymofscreen.magmom_handler import continue_magmoms, get_mag_indices
from pymofscreen.compute_environ import choose_vasp_version
from pymofscreen.metal_types import spblock_metals, poor_metals

def mof_run(screener,mof,calc,gpt_version):
#Get the optimized structure of the MOF
	cif_file = screener.cif_file
	stdout_file = screener.stdout_file
	nprocs = screener.nprocs
	calc_swaps = screener.calc_swaps
	mofpath = screener.mofpath
	basepath = screener.basepath
	success = False
	copyfile(os.path.join(mofpath,cif_file),os.path.join(basepath,'working',cif_file))
	calc, calc_swaps = update_calc(calc,calc_swaps)
	mof.set_calculator(calc)
	try:
		mof.get_potential_energy()
		niter = get_niter('OUTCAR')
		if niter < mof.calc.int_params['nsw'] and mof.calc.converged != True:
			raise SystemError('VASP stopped but did not crash and burn')
		success = True
	except:
		if not os.path.isfile('STOPCAR'):
			old_error_len = 0
			refcode = cif_file.split('.cif')[0]
			restart_files = ['WAVECAR','CHGCAR']
			for file in restart_files:
				if os.path.isfile(file):
					os.remove(file)
			while True:
				errormsg = get_error_msgs('OUTCAR',refcode,stdout_file)
				calc, calc_swaps = update_calc_after_errors(calc,calc_swaps,errormsg)
				error_len = len(errormsg)
				if error_len == old_error_len:
					break
				mof = continue_mof()
				choose_vasp_version(gpt_version,nprocs,calc_swaps)
				mof.set_calculator(calc)
				try:
					mof.get_potential_energy()
					niter = get_niter('OUTCAR')
					if niter < mof.calc.int_params['nsw'] and mof.calc.converged != True:
						raise SystemError('VASP stopped but did not crash and burn')
					success = True
				except:
					pass
				old_error_len = error_len
	if success == False:
		mof = None

	return mof, calc_swaps

def mof_bfgs_run(screener,mof,calc,steps,fmax):
#Optimize with BFGSLineSearch
	cif_file = screener.cif_file
	calc_swaps = screener.calc_swaps
	stdout_file = screener.stdout_file
	mofpath = screener.mofpath
	basepath = screener.basepath
	copyfile(os.path.join(mofpath,cif_file),os.path.join(basepath,'working',cif_file))
	success = False
	calc, calc_swaps = update_calc(calc,calc_swaps)
	mof.set_calculator(calc)
	dyn = BFGSLineSearch(mof,trajectory='opt.traj')
	try:
		dyn.run(fmax=fmax,steps=steps)
		success = True
	except:
		if not os.path.isfile('STOPCAR'):
			old_error_len = 0
			refcode = cif_file.split('.cif')[0]
			restart_files = ['WAVECAR','CHGCAR']
			for file in restart_files:
				if os.path.isfile(file):
					os.remove(file)
			while True:
				errormsg = get_error_msgs('OUTCAR',refcode,stdout_file)
				calc, calc_swaps = update_calc_after_errors(calc,calc_swaps,errormsg)
				error_len = len(errormsg)
				if error_len == old_error_len:
					break
				mof = continue_mof()
				mof.set_calculator(calc)
				dyn = BFGSLineSearch(mof,trajectory='opt.traj')
				try:				
					dyn.run(fmax=fmax,steps=steps)
					success = True
				except:
					pass
				old_error_len = error_len
	if success == False:
		mof = None

	return mof, dyn, calc_swaps

def prep_next_run(screener):
#Update counter and decide if next job should be skipped
	acc_levels = screener.acc_levels
	acc_level = acc_levels[screener.run_i]
	refcode = screener.refcode
	spin_level = screener.spin_level
	basepath = screener.basepath

	skip_spin2 = False
	success_path = os.path.join(basepath,'results',refcode,acc_level,spin_level)
	incarpath = os.path.join(success_path,'INCAR')
	outcarpath = os.path.join(success_path,'OUTCAR')
	errorpath = os.path.join(basepath,'errors',refcode,acc_level,spin_level)
	if os.path.exists(errorpath) == True:
		mof = None
	else:
		mof = read(outcarpath)
		if acc_level != 'scf_test':
			mof, abs_magmoms = continue_magmoms(mof,incarpath)
			mag_indices = get_mag_indices(mof)
			mag_nums = mof[mag_indices].get_atomic_numbers()
			if np.sum(abs_magmoms < 0.1) == len(abs_magmoms) or all(num in spblock_metals+poor_metals for num in mag_nums) == True:
				skip_spin2 = True
	screener.run_i += 1
	print(screener.run_i)

	return mof, skip_spin2