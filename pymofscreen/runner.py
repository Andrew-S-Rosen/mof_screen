import os
import numpy as np
from shutil import copyfile
from ase.io import read
from ase.optimize import BFGSLineSearch
from pymofscreen.calc_swaps import update_calc, check_nprocs
from pymofscreen.error_handler import get_niter, get_error_msgs, update_calc_after_errors, continue_mof
from pymofscreen.magmom_handler import continue_magmoms, get_mag_indices
from pymofscreen.compute_environ import choose_vasp_version
from pymofscreen.metal_types import spblock_metals, poor_metals
from pymofscreen.writers import pprint
from pymofscreen.janitor import clean_files

def mof_run(workflow,mof,calc,kpts):
	"""
	Run an atoms.get_potential_energy() calculation
	Args:
		workflow (class): pymofscreen.screen_phases.worfklow class
		mof (ASE Atoms object): ASE Atoms object for MOF
		calc (dict): ASE Vasp calculator
		kpts (list of ints): k-point grid
	Returns:
		mof (ASE Atoms object): updated ASE Atoms object
		calc_swaps (list of strings): calc swaps
	"""

	nprocs = workflow.nprocs
	ppn = workflow.ppn
	spin_level = workflow.spin_level
	acc_level = workflow.acc_levels[workflow.run_i]
	calc_swaps = workflow.calc_swaps
	cif_file = workflow.cif_file
	stdout_file = workflow.stdout_file
	calc_swaps = workflow.calc_swaps
	mofpath = workflow.mofpath
	basepath = workflow.basepath
	gamma = workflow.kpts_dict['gamma']

	if sum(kpts) == 3:
		gpt_version = True
	else:
		gpt_version = False
	nprocs = check_nprocs(len(mof),nprocs,ppn)
	choose_vasp_version(gpt_version,nprocs)
	calc.input_params['kpts'] = kpts
	calc.input_params['gamma'] = gamma

	copyfile(os.path.join(mofpath,cif_file),os.path.join(basepath,'working',cif_file))
	calc, calc_swaps = update_calc(calc,calc_swaps)
	mof.set_calculator(calc)

	pprint('Running '+spin_level+', '+acc_level)
	success = False

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
			clean_files(restart_files)

			while True:

				errormsg = get_error_msgs('OUTCAR',refcode,stdout_file)
				calc, calc_swaps = update_calc_after_errors(calc,calc_swaps,errormsg)
				error_len = len(errormsg)
				if error_len == old_error_len:
					break

				mof = continue_mof()
				choose_vasp_version(gpt_version,nprocs)
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

def mof_bfgs_run(workflow,mof,calc,kpts,steps=100,fmax=0.05):
	"""
	Run ASE BFGSLineSearch calculation
	Args:
		workflow (class): pymofscreen.screen_phases.worfklow class
		mof (ASE Atoms object): ASE Atoms object for MOF
		calc (dict): ASE Vasp calculator
		kpts (list of ints): k-point grid
		steps (int): maximum number of steps
		fmax (int): force tolerance
	Returns:
		mof (ASE Atoms object): updated ASE Atoms object
		dyn (class): ASE dynamics class
		calc_swaps (list of strings): calc swaps
	"""

	spin_level = workflow.spin_level
	acc_level = workflow.acc_level
	nprocs = workflow.nprocs
	ppn = workflow.ppn
	cif_file = workflow.cif_file
	stdout_file = workflow.stdout_file
	calc_swaps = workflow.calc_swaps
	mofpath = workflow.mofpath
	basepath = workflow.basepath
	gamma = workflow.gamma

	if sum(kpts) == 3:
		gpt_version = True
	else:
		gpt_version = False

	nprocs = check_nprocs(len(mof),nprocs,ppn)
	choose_vasp_version(gpt_version,nprocs)
	calc.input_params['kpts'] = kpts
	calc.input_params['gamma'] = gamma
	
	copyfile(os.path.join(mofpath,cif_file),os.path.join(basepath,'working',cif_file))
	calc, calc_swaps = update_calc(calc,calc_swaps)
	mof.set_calculator(calc)
	dyn = BFGSLineSearch(mof,trajectory='opt.traj')

	pprint('Running '+spin_level+', '+acc_level)
	success = False

	try:
		dyn.run(fmax=fmax,steps=steps)
		success = True
	except:

		if not os.path.isfile('STOPCAR'):

			old_error_len = 0
			refcode = cif_file.split('.cif')[0]
			restart_files = ['WAVECAR','CHGCAR']
			clean_files(restart_files)

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

def prep_next_run(workflow):
	"""
	Prepare for the next run
	Args:
		workflow (class): pymofscreen.screen_phases.worfklow class
	Returns:
		mof (ASE Atoms object): updated ASE Atoms object
		skip_spin2 (bool): True if next spin state should be skipped
	"""

	acc_levels = workflow.acc_levels
	acc_level = acc_levels[workflow.run_i]
	refcode = workflow.refcode
	spin_level = workflow.spin_level
	basepath = workflow.basepath

	success_path = os.path.join(basepath,'results',refcode,acc_level,spin_level)
	incarpath = os.path.join(success_path,'INCAR')
	outcarpath = os.path.join(success_path,'OUTCAR')
	errorpath = os.path.join(basepath,'errors',refcode,acc_level,spin_level)
	skip_spin2 = False

	if os.path.exists(errorpath):
		mof = None
	else:
		mof = read(outcarpath)
		if acc_level != 'scf_test':
			mof, abs_magmoms = continue_magmoms(mof,incarpath)
			mag_indices = get_mag_indices(mof)
			mag_nums = mof[mag_indices].get_atomic_numbers()
			if np.sum(abs_magmoms < 0.1) == len(abs_magmoms) or all(num in spblock_metals+poor_metals for num in mag_nums) == True:
				skip_spin2 = True

	workflow.run_i += 1

	return mof, skip_spin2