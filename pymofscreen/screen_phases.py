import os
import numpy as np
from ase.io import read
from pymofscreen.compute_environ import get_nprocs
from pymofscreen.writers import pprint, write_success, write_errors
from pymofscreen.janitor import manage_restart_files, clean_files
from pymofscreen.runner import mof_run, prep_next_run, prep_new_run, mof_bfgs_run
from pymofscreen.cif_handler import cif_to_mof
from pymofscreen.magmom_handler import set_initial_magmoms, continue_magmoms
from pymofscreen.error_handler import get_warning_msgs
from pymofscreen.default_calculators import defaults

class workflows():
	"""
	This class constructs a workflow for a given calculation stage
	"""

	def __init__(self,screener,cif_file,kpts_dict,spin_level,prior_spin=None,
		vasp_files=None):
		"""
		Initialize variables that should be used on all MOFs in a database
		Args:
			screener (class): pymofscreen.screen.screener class
			cif_file (string): name of CIF file
			kpts_dict (dict): dictionary containing kpoint and gamma information
			spin_level (string): name of spin level
			prior_spin (string): name of previous spin level (if applicable)
			vasp_files (list of strings): VASP files to save
		"""
		self.cif_file = cif_file
		self.kpts_dict = kpts_dict
		self.acc_levels = screener.acc_levels
		self.spin_level = spin_level
		if vasp_files is None:
			self.vasp_files = ['INCAR','POSCAR','KPOINTS','POTCAR','OUTCAR',
		'CONTCAR','CHGCAR','WAVECAR']
		self.calc_swaps = []
		self.run_i = 0
		self.refcode = cif_file.split('.cif')[0]
		self.stdout_file = screener.stdout_file
		self.mofpath = screener.mofpath
		self.submit_script = screener.submit_script
		self.basepath = screener.basepath
		self.niggli = screener.niggli
		self.calcs = screener.calcs

		self.nprocs, self.ppn = get_nprocs(self.submit_script)
		if defaults.get('ncore') is None and defaults.get('npar') is None:
			defaults['ncore'] = int(self.ppn/2.0)

		clean_files(self.vasp_files)

		results_partial_paths = []
		error_outcar_partial_paths = []
		error_outcar_paths = []
		outcar_paths = []
		for acc_level in self.acc_levels:
			results_partial_paths.append(os.path.join(self.basepath,'results',
				self.refcode,acc_level))
			error_outcar_partial_paths.append(os.path.join(self.basepath,
				'errors',self.refcode,acc_level))
		for results_partial_path in results_partial_paths:
			outcar_paths.append(os.path.join(results_partial_path,spin_level,
				'OUTCAR'))
		for error_outcar_partial_path in error_outcar_partial_paths:
			error_outcar_paths.append(os.path.join(error_outcar_partial_path,
				spin_level,'OUTCAR'))
		self.outcar_paths = outcar_paths
		self.error_outcar_paths = error_outcar_paths
		if prior_spin is None:
			self.spin1_final_mof_path = None
		else:
			self.spin1_final_mof_path = os.path.join(results_partial_paths[-1],
				prior_spin,'OUTCAR')
		pprint('***STARTING '+self.refcode+': '+spin_level+'***')
				
	def scf_test(self):
		"""
		Run SCF test job to check for errors
		Returns:
			scf_pass (bool): True if passed SCF test
		"""
		outcar_paths = self.outcar_paths
		error_outcar_paths = self.error_outcar_paths
		spin_level = self.spin_level
		cif_file = self.cif_file
		mofpath = self.mofpath
		spin1_final_mof_path = self.spin1_final_mof_path
		kpts_lo = self.kpts_dict['kpts_lo']
		acc_level = self.acc_levels[self.run_i]
		niggli = self.niggli
		calcs = self.calcs
		
		if not os.path.isfile(outcar_paths[self.run_i]) and not os.path.isfile(error_outcar_paths[self.run_i]):
			if spin1_final_mof_path is None:
				mof = cif_to_mof(os.path.join(mofpath,cif_file),niggli)
			else:
				mof = read(spin1_final_mof_path)
			mof = set_initial_magmoms(mof,spin_level)
			pprint('Running '+spin_level+', '+acc_level)
			mof, self.calc_swaps = mof_run(self,mof,calcs('scf_test'),kpts_lo)
			if mof != None:
				write_success(self)
			else:
				pprint('^ VASP crashed')
				write_errors(self,mof)
		elif os.path.isfile(outcar_paths[self.run_i]):
			pprint('COMPLETED: '+spin_level+', '+acc_level)
		mof = prep_next_run(self)
		if mof == None:
			pprint('Skipping rest because of errors')
			return False
		warnings = get_warning_msgs(outcar_paths[self.run_i-1])
		self.calc_swaps.extend(warnings)

		return True

	def isif2_lowacc(self):
		"""
		Run low accuracy ISIF2
		Returns:
			mof (ASE Atoms object): updated ASE Atoms object
		"""
		acc_levels = self.acc_levels
		outcar_paths = self.outcar_paths
		error_outcar_paths = self.error_outcar_paths
		spin_level = self.spin_level
		cif_file = self.cif_file
		mofpath = self.mofpath
		spin1_final_mof_path = self.spin1_final_mof_path
		kpts_lo = self.kpts_dict['kpts_lo']
		acc_level = acc_levels[self.run_i]
		niggli = self.niggli
		calcs = self.calcs

		if os.path.isfile(outcar_paths[self.run_i-1]) and not os.path.isfile(outcar_paths[self.run_i]) and not os.path.isfile(error_outcar_paths[self.run_i]):
			if spin1_final_mof_path is None:
				mof = cif_to_mof(os.path.join(mofpath,cif_file),niggli)
			else:
				mof = read(spin1_final_mof_path)
			manage_restart_files(outcar_paths[self.run_i-1].split('OUTCAR')[0])
			mof = set_initial_magmoms(mof,spin_level)
			fmax = 5.0
			pprint('Running '+spin_level+', '+acc_level)
			mof, dyn, self.calc_swaps = mof_bfgs_run(self,mof,calcs('ase_bfgs'),
				kpts_lo,fmax=fmax)
			if mof != None and dyn and mof.calc.scf_converged == True:
				loop_i = 0
				converged = False
				clean_files(['opt.traj'])
				while mof != None and loop_i < 4 and converged == False and mof.calc.scf_converged == True:
					if loop_i == 2 and 'fire' not in self.calc_swaps and 'zbrent' not in self.calc_swaps:
						self.calc_swaps.append('fire')
					mof = read('OUTCAR')
					mof = continue_magmoms(mof,'INCAR')
					mof, self.calc_swaps = mof_run(self,mof,
						calcs('isif2_lowacc'),kpts_lo)
					if mof == None:
						break
					converged = mof.calc.converged
					loop_i += 1
			if 'fire' in self.calc_swaps:
				self.calc_swaps.remove('fire')
			if mof != None and mof.calc.scf_converged == True and mof.calc.converged == True:
				write_success(self)
			else:
				write_errors(self,mof)
		elif os.path.isfile(outcar_paths[self.run_i]) == True:
			pprint('COMPLETED: '+spin_level+', '+acc_level)
		mof = prep_next_run(self)
		if mof == None:
			pprint('Skipping rest because of errors')
			return None

		return mof

	def isif2_medacc(self):
		"""
		Run medium accuracy ISIF2
		Returns:
			mof (ASE Atoms object): updated ASE Atoms object
		"""
		acc_levels = self.acc_levels
		outcar_paths = self.outcar_paths
		error_outcar_paths = self.error_outcar_paths
		spin_level = self.spin_level
		kpts_lo = self.kpts_dict['kpts_lo']
		kpts_hi = self.kpts_dict['kpts_hi']
		acc_level = acc_levels[self.run_i]
		calcs = self.calcs

		if os.path.isfile(outcar_paths[self.run_i-1]) and not os.path.isfile(outcar_paths[self.run_i]) and not os.path.isfile(error_outcar_paths[self.run_i]):
			mof = prep_new_run(self)
			if sum(kpts_lo) == 3 and sum(kpts_hi) > 3:
				clean_files(['CHGCAR','WAVECAR'])
			else:
				manage_restart_files(outcar_paths[self.run_i-1].split('OUTCAR')[0])
			pprint('Running '+spin_level+', '+acc_level)
			mof = prep_new_run()
			mof,self.calc_swaps = mof_run(self,mof,calcs('isif2_medacc'),kpts_hi)
			if mof != None and mof.calc.scf_converged == True and mof.calc.converged == True:
				write_success(self)
			else:
				write_errors(self,mof)
		elif os.path.isfile(outcar_paths[self.run_i]) == True:
			pprint('COMPLETED: '+spin_level+', '+acc_level)
		mof = prep_next_run(self)
		if mof == None:
			pprint('Skipping rest because of errors')
			return None

		return mof

	def isif2_highacc(self):
		"""
		Run high accuracy ISIF2
		Returns:
			mof (ASE Atoms object): updated ASE Atoms object
		"""
		acc_levels = self.acc_levels
		outcar_paths = self.outcar_paths
		error_outcar_paths = self.error_outcar_paths
		spin_level = self.spin_level
		kpts_hi = self.kpts_dict['kpts_hi']
		acc_level = acc_levels[self.run_i]
		calcs = self.calcs

		if os.path.isfile(outcar_paths[self.run_i-1]) and not os.path.isfile(outcar_paths[self.run_i]) and not os.path.isfile(error_outcar_paths[self.run_i]):
			mof = prep_new_run(self)
			manage_restart_files(outcar_paths[self.run_i-1].split('OUTCAR')[0])
			pprint('Running '+spin_level+', '+acc_level)
			mof,self.calc_swaps = mof_run(self,mof,calcs('isif2_highacc'),kpts_hi)
			if mof != None and mof.calc.scf_converged == True and mof.calc.converged == True:
				if 'large_supercell' in self.calc_swaps:
					self.calc_swaps.remove('large_supercell')
					mof = read('OUTCAR')
					mof = continue_magmoms(mof,'INCAR')
					mof, self.calc_swaps = mof_run(self,mof,calcs('isif2_highacc'),kpts_hi)
					if mof != None and mof.calc.scf_converged == True and mof.calc.converged == True:
						write_success(self)
					else:
						write_errors(self,mof)
				else:
					write_success(self)
			else:
				write_errors(self,mof)
		elif os.path.isfile(outcar_paths[self.run_i]) == True:
			pprint('COMPLETED: '+spin_level+', '+acc_level)
		mof = prep_next_run(self)
		if mof == None:
			pprint('Skipping rest because of errors')
			return None
			
		return mof

	def isif3_lowacc(self):
		"""
		Run low accuracy ISIF3
		Returns:
			mof (ASE Atoms object): updated ASE Atoms object
		"""
		acc_levels = self.acc_levels
		outcar_paths = self.outcar_paths
		error_outcar_paths = self.error_outcar_paths
		spin_level = self.spin_level
		kpts_lo = self.kpts_dict['kpts_lo']
		acc_level = acc_levels[self.run_i]
		calcs = self.calcs

		if os.path.isfile(outcar_paths[self.run_i-1]) and not os.path.isfile(outcar_paths[self.run_i]) and not os.path.isfile(error_outcar_paths[self.run_i]):
			mof = prep_new_run(self)
			converged = False
			loop_i = 0
			n_runs = 15
			manage_restart_files(outcar_paths[self.run_i-1].split('OUTCAR')[0])
			while converged == False and loop_i < n_runs:
				pprint('Running '+spin_level+', '+acc_level+': iteration '+str(loop_i)+'/'+str(n_runs-1))
				if loop_i == 10 and 'fire' not in self.calc_swaps and 'zbrent' not in self.calc_swaps:
					self.calc_swaps.append('fire')
				mof,self.calc_swaps = mof_run(self,mof,calcs('isif3_lowacc'),kpts_lo)
				if mof == None:
					break
				converged = mof.calc.converged
				mof = read('OUTCAR')
				mof = continue_magmoms(mof,'INCAR')
				loop_i += 1
			if 'fire' in self.calc_swaps:
				self.calc_swaps.remove('fire')
			if mof != None and converged == True:
				write_success(self)
			else:
				write_errors(self,mof)
		elif os.path.isfile(outcar_paths[self.run_i]) == True:
			pprint('COMPLETED: '+spin_level+', '+acc_level)
		mof = prep_next_run(self)
		if mof == None:
			pprint('Skipping rest because of errors')
			return None

		return mof

	def isif3_highacc(self):
		"""
		Run high accuracy ISIF3
		Returns:
			mof (ASE Atoms object): updated ASE Atoms object
		"""
		acc_levels = self.acc_levels
		outcar_paths = self.outcar_paths
		error_outcar_paths = self.error_outcar_paths
		spin_level = self.spin_level
		kpts_lo = self.kpts_dict['kpts_lo']
		kpts_hi = self.kpts_dict['kpts_hi']
		acc_level = acc_levels[self.run_i]
		calcs = self.calcs

		if os.path.isfile(outcar_paths[self.run_i-1]) and not os.path.isfile(outcar_paths[self.run_i]) and not os.path.isfile(error_outcar_paths[self.run_i]):
			mof = prep_new_run(self)
			converged = False
			loop_i = 0
			n_runs = 15
			V_diff = np.inf
			V_cut = 0.01
			V0 = mof.get_volume()
			if sum(kpts_lo) == 3 and sum(kpts_hi) > 3:
				clean_files(['CHGCAR','WAVECAR'])
			else:
				manage_restart_files(outcar_paths[self.run_i-1].split('OUTCAR')[0])
			while (converged == False or V_diff > V_cut) and loop_i < n_runs:
				pprint('Running '+spin_level+', '+acc_level+': iteration '+str(loop_i)+'/'+str(n_runs-1))
				if loop_i == 10 and 'fire' not in self.calc_swaps and 'zbrent' not in self.calc_swaps:
					self.calc_swaps.append('fire')
				mof,self.calc_swaps = mof_run(self,mof,calcs('isif3_highacc'),
					kpts_hi)
				if mof == None:
					break
				if loop_i > 0:
					converged = mof.calc.converged
				mof = read('OUTCAR')
				V = mof.get_volume()
				mof = continue_magmoms(mof,'INCAR')
				if loop_i > 0:
					V_diff = np.abs((V-V0))/V0
				V0 = V
				loop_i += 1
			if mof != None and converged == True and V_diff <= V_cut and 'large_supercell' in self.calc_swaps:
				pprint('Running '+spin_level+', '+acc_level+' (LREAL=False)')
				self.calc_swaps.append('nsw=100')
				self.calc_swaps.remove('large_supercell')
				mof,self.calc_swaps = mof_run(self,mof,calcs('isif3_highacc'),
					kpts_hi)
				self.calc_swaps.remove('nsw=100')
				if mof != None and mof.calc.converged == True:
					write_success(self)
				else:
					write_errors(self,mof)
			else:
				write_errors(self,mof)
				if mof is not None and V_diff > V_cut:
					pprint('^ Change in V of '+str(V_diff)+' percent')
			if 'fire' in self.calc_swaps:
				self.calc_swaps.remove('fire')
		elif os.path.isfile(outcar_paths[self.run_i]) == True:
			pprint('COMPLETED: '+spin_level+', '+acc_level)
		mof = prep_next_run(self)
		if mof == None:
			pprint('Skipping rest because of errors')
			return None

		return mof

	def final_spe(self):
		"""
		Run final single point
		Returns:
			mof (ASE Atoms object): updated ASE Atoms object
		"""
		acc_levels = self.acc_levels
		outcar_paths = self.outcar_paths
		error_outcar_paths = self.error_outcar_paths
		spin_level = self.spin_level
		kpts_hi = self.kpts_dict['kpts_hi']
		acc_level = acc_levels[self.run_i]
		calcs = self.calcs
		self.calc_swaps = []
		if os.path.isfile(outcar_paths[self.run_i-1]) and not os.path.isfile(outcar_paths[self.run_i]) and not os.path.isfile(error_outcar_paths[self.run_i]):
			mof = prep_new_run(self)
			manage_restart_files(outcar_paths[self.run_i-1].split('OUTCAR')[0])
			pprint('Running '+spin_level+', '+acc_level)
			mof,self.calc_swaps = mof_run(self,mof,calcs('final_spe'),kpts_hi)
			if mof != None and mof.calc.scf_converged == True:
				write_success(self)
			else:
				write_errors(self,mof)
		elif os.path.isfile(outcar_paths[self.run_i]) == True:
			pprint('COMPLETED: '+spin_level+', '+acc_level)
		mof = prep_next_run(self)
		if mof == None:
			pprint('Skipping rest because of errors')
			return None

		return mof