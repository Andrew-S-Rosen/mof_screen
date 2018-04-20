import os
import numpy as np
from ase.io import read
from pymofscreen.compute_environ import get_nprocs, choose_vasp_version
from pymofscreen.writers import pprint, write_success, write_errors
from pymofscreen.kpts_handler import get_gpt_version
from pymofscreen.janitor import manage_restart_files
from pymofscreen.runner import mof_run, prep_next_run, mof_bfgs_run
from pymofscreen.cif_handler import cif_to_mof
from pymofscreen.magmom_handler import set_initial_magmoms, continue_magmoms, get_mag_indices
from pymofscreen.calculators import calcs_ads, calcs_vol
from pymofscreen.error_handler import get_warning_msgs

class workflows():

	def __init__(self,screener,defaults,cif_file,spin_level,prior_spin=None,vasp_files=['INCAR','POSCAR','KPOINTS','POTCAR','OUTCAR',
		'CONTCAR','CHGCAR','WAVECAR']):
		self.cif_file = cif_file
		self.acc_levels = screener.acc_levels
		self.spin_level = spin_level
		self.vasp_files = vasp_files
		self.calc_swaps = []
		self.run_i = 0
		self.refcode = cif_file.split('.cif')[0]
		self.defaults = defaults
		self.stdout_file = screener.stdout_file
		self.mofpath = screener.mofpath
		self.submit_script = screener.submit_script
		self.basepath = screener.basepath
		self.niggli = screener.niggli
		results_partial_paths = []
		error_outcar_partial_paths = []
		error_outcar_paths = []
		outcar_paths = []
		for acc_level in self.acc_levels:
			results_partial_paths.append(os.path.join(self.basepath,'results',self.refcode,acc_level))
			error_outcar_partial_paths.append(os.path.join(self.basepath,'errors',self.refcode,acc_level))
		for results_partial_path in results_partial_paths:
			outcar_paths.append(os.path.join(results_partial_path,spin_level,'OUTCAR'))
		for error_outcar_partial_path in error_outcar_partial_paths:
			error_outcar_paths.append(os.path.join(error_outcar_partial_path,spin_level,'OUTCAR'))
		self.outcar_paths = outcar_paths
		self.error_outcar_paths = error_outcar_paths
		self.nprocs, self.ppn = get_nprocs(self.submit_script)
		self.clean_files(vasp_files)
		if prior_spin is None:
			self.spin1_final_mof_path = None
		else:
			self.spin1_final_mof_path = os.path.join(results_partial_paths[-1],prior_spin,'OUTCAR')
		pprint('***STARTING '+self.refcode+': '+spin_level+'***')

	def clean_files(self,remove_files):
	#Clean files
		for file in remove_files:
			if os.path.isfile(file) == True:
				os.remove(file)
				
	def scf_test(self):
		refcode = self.refcode
		vasp_files = self.vasp_files
		acc_levels = self.acc_levels
		outcar_paths = self.outcar_paths
		error_outcar_paths = self.error_outcar_paths
		spin_level = self.spin_level
		cif_file = self.cif_file
		nprocs = self.nprocs
		ppn = self.ppn
		mofpath = self.mofpath
		spin1_final_mof_path = self.spin1_final_mof_path
		kpts_lo = self.defaults['kpts_lo']
		acc_level = acc_levels[self.run_i]
		niggli = self.niggli
		if os.path.isfile(outcar_paths[self.run_i]) != True and os.path.isfile(error_outcar_paths[self.run_i]) != True:
			if spin1_final_mof_path is None:
				mof = cif_to_mof(mofpath,cif_file,niggli)
			else:
				mof = read(spin1_final_mof_path)
			mof = set_initial_magmoms(mof,spin_level)
			gpt_version, nprocs = get_gpt_version(kpts_lo,len(mof),nprocs,ppn)
			choose_vasp_version(gpt_version,nprocs,self.calc_swaps)
			pprint('Running '+spin_level+', '+acc_level)
			mof, self.calc_swaps = mof_run(self,mof,calcs_ads(self.run_i),gpt_version)
			if mof != None:
				write_success(self,acc_level,vasp_files)
			else:
				pprint('^ VASP crashed')
				write_errors(refcode,spin_level,acc_level,vasp_files,cif_file)
		elif os.path.isfile(outcar_paths[self.run_i]) == True:
			pprint('COMPLETED: '+spin_level+', '+acc_level)
		mof, skip_spin2 = prep_next_run(self)
		if mof == None:
			pprint('Skipping rest because of errors')
			return None, None
		warnings = get_warning_msgs(outcar_paths[self.run_i-1])
		self.calc_swaps.extend(warnings)
		return mof

	def isif2_lowacc(self):
		vasp_files = self.vasp_files
		acc_levels = self.acc_levels
		outcar_paths = self.outcar_paths
		error_outcar_paths = self.error_outcar_paths
		spin_level = self.spin_level
		cif_file = self.cif_file
		nprocs = self.nprocs
		ppn = self.ppn
		mofpath = self.mofpath
		spin1_final_mof_path = self.spin1_final_mof_path
		kpts_lo = self.defaults['kpts_lo']
		acc_level = acc_levels[self.run_i]
		niggli = self.niggli
		if os.path.isfile(outcar_paths[self.run_i-1]) == True and os.path.isfile(outcar_paths[self.run_i]) != True and os.path.isfile(error_outcar_paths[self.run_i]) != True:
			if spin1_final_mof_path is None:
				mof = cif_to_mof(mofpath,cif_file,niggli)
			else:
				mof = read(spin1_final_mof_path)
			manage_restart_files(outcar_paths[self.run_i-1].split('OUTCAR')[0])
			mof = set_initial_magmoms(mof,spin_level)
			gpt_version, nprocs = get_gpt_version(kpts_lo,len(mof),nprocs,ppn)
			choose_vasp_version(gpt_version,nprocs,self.calc_swaps)
			pprint('Running '+spin_level+', '+acc_level)
			steps = 100
			fmax = 5.0
			mof, dyn, self.calc_swaps = mof_bfgs_run(self,mof,calcs_ads(self.run_i),steps,fmax)
			if mof != None and dyn and mof.calc.scf_converged == True:
				loop_i = 0
				converged = False
				self.clean_files(['opt.traj'])
				while mof != None and loop_i < 5 and converged == False and mof.calc.scf_converged == True:
					mof = read('OUTCAR')
					mof, abs_magmoms = continue_magmoms(mof,'INCAR')
					choose_vasp_version(gpt_version,nprocs,self.calc_swaps)
					mof, self.calc_swaps = mof_run(self,mof,calcs_ads(1.5),gpt_version)
					if mof == None:
						break
					converged = mof.calc.converged
					loop_i += 1
			if mof != None and mof.calc.scf_converged == True and mof.calc.converged == True:
				write_success(self,acc_level,vasp_files)
			else:
				write_errors(self,acc_level,vasp_files)
				if mof == None:
					pprint('^ VASP crashed')
				elif mof.calc.scf_converged == False:
					pprint('^ SCF did not converge')
				elif mof.calc.converged == False:
					pprint('^ Convergence not reached')
		elif os.path.isfile(outcar_paths[self.run_i]) == True:
			pprint('COMPLETED: '+spin_level+', '+acc_level)
		mof, skip_spin2 = prep_next_run(self)
		if mof == None:
			pprint('Skipping rest because of errors')
			return None
		if spin_level == 'spin2':
			mag_indices = get_mag_indices(mof)
			old_mof = read(spin1_final_mof_path)
			if np.sum(np.abs(mof.get_initial_magnetic_moments()[mag_indices] - old_mof.get_magnetic_moments()[mag_indices]) >= 0.05) == 0:
				pprint('Skipping rest because SPIN2 converged to SPIN1')
				return None
		return mof

	def isif2_medacc(self,mof):
		vasp_files = self.vasp_files
		acc_levels = self.acc_levels
		outcar_paths = self.outcar_paths
		error_outcar_paths = self.error_outcar_paths
		spin_level = self.spin_level
		nprocs = self.nprocs
		ppn = self.ppn
		kpts_lo = self.defaults['kpts_lo']
		kpts_hi = self.defaults['kpts_hi']
		acc_level = acc_levels[self.run_i]
		if os.path.isfile(outcar_paths[self.run_i-1]) == True and os.path.isfile(outcar_paths[self.run_i]) != True and os.path.isfile(error_outcar_paths[self.run_i]) != True:
			gpt_version, nprocs = get_gpt_version(kpts_hi,len(mof),nprocs,ppn)
			pprint('Running '+spin_level+', '+acc_level)
			if sum(kpts_lo) == 3 and sum(kpts_hi) > 3:
				self.clean_files(['CHGCAR','WAVECAR'])
			else:
				manage_restart_files(outcar_paths[self.run_i-1].split('OUTCAR')[0])
			choose_vasp_version(gpt_version,nprocs,self.calc_swaps)
			mof,self.calc_swaps = mof_run(self,mof,calcs_ads(self.run_i),gpt_version)
			if mof != None and mof.calc.scf_converged == True and mof.calc.converged == True:
				write_success(self,acc_level,vasp_files)
			else:
				write_errors(self,acc_level,vasp_files)
				if mof == None:
					pprint('^ VASP crashed')
				elif mof.calc.scf_converged == False:
					pprint('^ SCF did not converge')
				elif mof.calc.converged == False:
					pprint('^ Convergence not reached')
		elif os.path.isfile(outcar_paths[self.run_i]) == True:
			pprint('COMPLETED: '+spin_level+', '+acc_level)
		mof, skip_spin2 = prep_next_run(self)
		if mof == None:
			pprint('Skipping rest because of errors')
			return None
		return mof

	def isif2_final(self,mof):
		vasp_files = self.vasp_files
		acc_levels = self.acc_levels
		outcar_paths = self.outcar_paths
		error_outcar_paths = self.error_outcar_paths
		spin_level = self.spin_level
		nprocs = self.nprocs
		ppn = self.ppn
		kpts_hi = self.defaults['kpts_hi']
		acc_level = acc_levels[self.run_i]
		if os.path.isfile(outcar_paths[self.run_i-1]) == True and os.path.isfile(outcar_paths[self.run_i]) != True and os.path.isfile(error_outcar_paths[self.run_i]) != True:
			gpt_version, nprocs = get_gpt_version(kpts_hi,len(mof),nprocs,ppn)
			choose_vasp_version(gpt_version,nprocs,self.calc_swaps)
			manage_restart_files(outcar_paths[self.run_i-1].split('OUTCAR')[0])
			pprint('Running '+spin_level+', '+acc_level)
			mof,self.calc_swaps = mof_run(self,mof,calcs_ads(self.run_i),gpt_version)
			if mof != None and mof.calc.scf_converged == True and mof.calc.converged == True:
				if 'large_supercell' in self.calc_swaps:
					pprint('Running '+spin_level+', '+acc_level+' (LREAL=False)')
					self.calc_swaps.remove('large_supercell')
					mof = read('OUTCAR')
					mof, abs_magmoms = continue_magmoms(mof,'INCAR')
					mof, self.calc_swaps = mof_run(self,mof,calcs_ads(self.run_i),gpt_version)
					if mof != None and mof.calc.scf_converged == True and mof.calc.converged == True:
						write_success(self,acc_level,vasp_files)
					else:
						write_errors(self,acc_level,vasp_files)
						if mof == None:
							pprint('^ VASP crashed')
						elif mof.calc.scf_converged == False:
							pprint('^ SCF did not converge')
						elif mof.calc.converged == False:
							pprint('^ Convergence not reached')
				else:
					write_success(self,acc_level,vasp_files)
			else:
				write_errors(self,acc_level,vasp_files)
				if mof == None:
					pprint('^ VASP crashed')
				elif mof.calc.scf_converged == False:
					pprint('^ SCF did not converge')
				elif mof.calc.converged == False:
					pprint('^ Convergence not reached')
		elif os.path.isfile(outcar_paths[self.run_i]) == True:
			pprint('COMPLETED: '+spin_level+', '+acc_level)
		mof, skip_spin2 = prep_next_run(self)
		if mof == None:
			pprint('Skipping rest because of errors')
			return None
		return mof

	def final_spe(self,mof):
		vasp_files = self.vasp_files
		acc_levels = self.acc_levels
		outcar_paths = self.outcar_paths
		error_outcar_paths = self.error_outcar_paths
		spin_level = self.spin_level
		nprocs = self.nprocs
		ppn = self.ppn
		kpts_hi = self.defaults['kpts_hi']
		acc_level = acc_levels[self.run_i]
		calc_swaps = []
		if os.path.isfile(outcar_paths[self.run_i-1]) == True and os.path.isfile(outcar_paths[self.run_i]) != True and os.path.isfile(error_outcar_paths[self.run_i]) != True:
			gpt_version, nprocs = get_gpt_version(kpts_hi,len(mof),nprocs,ppn)
			choose_vasp_version(gpt_version,nprocs,calc_swaps)
			manage_restart_files(outcar_paths[self.run_i-1].split('OUTCAR')[0])
			pprint('Running '+spin_level+', '+acc_level)
			mof,calc_swaps = mof_run(self,mof,calcs_ads(self.run_i),gpt_version)
			if mof != None and mof.calc.scf_converged == True:
				write_success(self,acc_level,vasp_files)
			else:
				write_errors(self,acc_level,vasp_files)
				if mof == None:
					pprint('^ VASP crashed')
				elif mof.calc.scf_converged == False:
					pprint('^ SCF did not converge')
		elif os.path.isfile(outcar_paths[self.run_i]) == True:
			pprint('COMPLETED: '+spin_level+', '+acc_level)
		mof, skip_spin2 = prep_next_run(self)
		if mof == None:
			pprint('Skipping rest because of errors')
			return None
		return mof, skip_spin2