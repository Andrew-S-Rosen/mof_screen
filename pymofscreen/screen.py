import os
import numpy as np
import sys
from copy import deepcopy
from shutil import copyfile
from pymofscreen.writers import pprint
from pymofscreen.kpts_handler import get_kpts
from pymofscreen.screen_phases import workflows
from pymofscreen.janitor import prep_paths
from pymofscreen.default_calculators import calcs
from pymofscreen.magmom_handler import check_if_new_spin, check_if_skip_low_spin

class screener():
	"""
	This class constructs a high-throughput screening workflow
	"""

	def __init__(self,basepath,mofpath=None,kpts_path='Auto',kppas=None,
		submit_script=None,stdout_file=None):
		"""
		Initialize variables that should be used on all MOFs in a databasehttp://sites.northwestern.edu/
		Args:
			basepath (string): path to the base directory for the DFT screening
			mofpath (string): path to the directory containing the CIF files
			kpts_path (string): can be either 'Auto' for an automatic generation
			of the kpoints based on KPPAs or a string representing the path to a
			text file with all the kpoint information
			kppas (list of ints): KPPAs to use if kpts_path == 'Auto'
			submit_script (string): path to job submission script
			stdout_file (string): path to the stdout file
		"""

		#Setup default parameters
		self.mofpath = mofpath
		self.basepath = basepath
		pwd = os.getcwd()
		if submit_script is None:
			submit_script = os.path.join(pwd,'sub_screen.job')
		self.submit_script = submit_script
		if stdout_file is None:
			stdout_file = os.path.join(pwd,sys.argv[0].split('.py')[0]+'.out')
		self.stdout_file = stdout_file
		self.kpts_path = kpts_path
		if kppas is None:
			self.kppas = [100,1000]
		prep_paths(basepath)

	def run_screen(self,cif_file,mode,spin_levels=None,acc_levels=None,niggli=True,calcs=calcs):
		"""
		Run high-throughput ionic or volume relaxations
		Args:
			cif_file (string): name of CIF file
			mode (string): 'ionic' or 'volume'
			spin_levels (list of strings): spin states to consider
			acc_levels (list of strings): accuracy levels to consider
			niggli (bool): True/False if Niggli-reduction should be done
			calcs (function): function to call respective calculator
		Returns:
			best_mof (ASE Atoms objects): ASE Atoms object for optimized
			MOF given by cif_file (lowest energy spin state)
		"""

		#Setup default parameters
		basepath = self.basepath
		self.calcs = calcs
		self.niggli = niggli
		if mode == 'ionic':
			if acc_levels is None:
				acc_levels = ['scf_test','isif2_lowacc','isif2_medacc',
				'isif2_highacc','final_spe']
		elif mode == 'ionic_legacy':
			if acc_levels is None:
				acc_levels = ['scf_test','isif2_lowacc','isif2_medacc',
				'final','final_spe']
		elif mode == 'volume':
			if acc_levels is None:
				acc_levels = ['scf_test','isif2_lowacc','isif3_lowacc',
				'isif3_highacc','isif2_highacc','final_spe']
		elif mode == 'volume_legacy':
			if acc_levels is None:
				acc_levels = ['scf_test','isif2','isif3_lowacc',
				'isif3_highacc','final','final_spe']
		else:
			raise ValueError('Unsupported DFT screening mode')
		if 'scf_test' not in acc_levels:
			acc_levels = ['scf_test']+acc_levels
		self.acc_levels = acc_levels
		if spin_levels is None:
			spin_levels = ['spin1','spin2']
		self.spin_levels = spin_levels

		#Make sure MOF isn't running on other process
		refcode = cif_file.split('.cif')[0]
		working_cif_path = os.path.join(basepath,'working',refcode)

		if os.path.isfile(working_cif_path):
			pprint('SKIPPED: Running on another process')
			return None
		open(working_cif_path,'w').close()

		#Get the kpoints
		kpts_lo, gamma = get_kpts(self,cif_file,'low')
		kpts_hi, gamma = get_kpts(self,cif_file,'high')
		kpts_dict = {}
		kpts_dict['kpts_lo'] = kpts_lo
		kpts_dict['kpts_hi'] = kpts_hi
		kpts_dict['gamma'] = gamma

		#Initialize variables
		E = np.inf
		mof = None

		#for each spin level, optimize the structure
		for i, spin_level in enumerate(spin_levels):

			pprint('***STARTING '+refcode+': '+spin_level+'***')

			#Check if spin state should be skipped
			if spin_level != spin_levels[0]:
				prior_spin = spin_levels[i-1]
			else:
				prior_spin = None
			if i > 0:
				skip_low_spin = check_if_skip_low_spin(self,mof,refcode,prior_spin)
				if (prior_spin == 'spin1' or prior_spin == 'high_spin') and skip_low_spin:
					pprint('Skipping '+spin_level+' run')
					continue
			same_spin = False

			#Set up workflow object
			wf = workflows(self,cif_file,kpts_dict,spin_level,prior_spin)

			#for each accuracy level, optimize structure
			for acc_level in acc_levels:

				if acc_level == 'scf_test':
					scf_pass = wf.scf_test()
					if not scf_pass:
						os.remove(working_cif_path)
						return None

				elif acc_level == 'isif2_lowacc' or (acc_level == 'isif2' and mode == 'volume_legacy'):
					mof = wf.isif2_lowacc()
					if mof is None:
						os.remove(working_cif_path)
						return None

					if i > 0:
						is_new_spin = check_if_new_spin(self,mof,refcode,acc_level,spin_level)
						if not is_new_spin:
							same_spin = True
							break

				elif acc_level == 'isif2_medacc':
					mof = wf.isif2_medacc()
					if mof is None:
						os.remove(working_cif_path)
						return None

				elif acc_level == 'isif2_highacc' or (acc_level == 'final' and 'legacy' in mode):
					mof = wf.isif2_highacc()
					if mof is None:
						os.remove(working_cif_path)
						return None
				
				elif acc_level == 'isif3_lowacc':
					mof = wf.isif3_lowacc()
					if mof is None:
						os.remove(working_cif_path)
						return None

				elif acc_level == 'isif3_highacc':
					mof = wf.isif3_highacc()
					if mof is None:
						os.remove(working_cif_path)
						return None

				elif acc_level == 'final_spe':
					mof = wf.final_spe()
					if mof is None:
						os.remove(working_cif_path)
						return None

				else:
					raise ValueError('Unsupported accuracy level')
		
			#***********SAVE and CONTINUE***********
			if same_spin:
				continue
			E_temp = mof.get_potential_energy()
			if E_temp < E:
				best_mof = deepcopy(mof)

		os.remove(working_cif_path)
		return best_mof

	def run_ts_screen(self,name,initial_atoms,final_atoms,n_images=4,cif_file=None,spin_levels=None,acc_levels=None,calcs=calcs):
		"""
		Run high-throughput TS calculation
		Args:
			name (string): name of CIF file
			initial_atoms (ASE Atoms object): initial structure
			final_atoms (ASE Atoms object): final structure
			n_images (int): number of NEB images
			cif_file (string): name of CIF file to generate kpoints if
			set to 'Auto'
			spin_levels (list of strings): spin states to consider
			acc_levels (list of strings): accuracy levels to consider
			calcs (function): function to call respective calculator
		Returns:
			best_mof (ASE Atoms objects): ASE Atoms object for optimized
			MOF given by cif_file (lowest energy spin state)
		"""

		#Setup default parameters
		basepath = self.basepath
		self.niggli = False
		self.spin_levels = spin_levels
		self.calcs = calcs
		if spin_levels is None:
			spin_levels = ['spin1','spin2']
		self.spin_levels = spin_levels
		if acc_levels is None:
			acc_levels = ['scf_test','cineb_lowacc','dimer_lowacc','dimer_medacc',
			'dimer_highacc','final_spe']
		if 'cineb_lowacc' not in acc_levels:
			acc_levels = ['cineb_lowacc']+acc_levels
		if 'scf_test' not in acc_levels:
			acc_levels = ['scf_test']+acc_levels
		self.acc_levels = acc_levels
		kpts_path = self.kpts_path
		if kpts_path == 'Auto':
			if cif_file is None:
				raise ValueError('Specify a CIF file if not using automatic KPPA')
			elif cif_file.split('.cif')[0] == name:
				raise ValueError('Input name and name of CIF file must not be identical')

		#Ensure initial/final state have the same composition
		if initial_atoms.get_chemical_formula() != final_atoms.get_chemical_formula():
			pprint('SKIPPED: Atoms not identical between initial and final state')
			return None

		#Make sure MOF isn't running on other process
		working_cif_path = os.path.join(basepath,'working',name)
		if os.path.isfile(working_cif_path):
			pprint('SKIPPED: Running on another process')
			return None
		open(working_cif_path,'w').close()

		#Get the kpoints
		kpts_lo, gamma = get_kpts(self,name,'low')
		kpts_hi, gamma = get_kpts(self,name,'high')
		kpts_dict = {}
		kpts_dict['kpts_lo'] = kpts_lo
		kpts_dict['kpts_hi'] = kpts_hi
		kpts_dict['gamma'] = gamma

		#Initialize variables
		E = np.inf
		mof = None

		#for each spin level, optimize the structure
		for i, spin_level in enumerate(spin_levels):

			pprint('***STARTING '+name+': '+spin_level+'***')

			#Check if spin state should be skipped
			if spin_level != spin_levels[0]:
				prior_spin = spin_levels[i-1]
			else:
				prior_spin = None
			if i > 0:
				skip_low_spin = check_if_skip_low_spin(self,mof,name,prior_spin)
				if (prior_spin == 'spin1' or prior_spin == 'high_spin') and skip_low_spin:
					pprint('Skipping '+spin_level+' run')
					continue
			same_spin = False

			#Set up workflow object
			wf = workflows(self,name,kpts_dict,spin_level,prior_spin)

			#for each accuracy level, optimize structure
			for acc_level in acc_levels:

				if acc_level == 'scf_test':
					scf_pass = wf.scf_test(atoms_overwrite=initial_atoms,quick_test=True)
					if not scf_pass:
						os.remove(working_cif_path)
						return None					

				elif acc_level == 'cineb_lowacc' and i == 0:
					neb_conv = wf.cineb_lowacc(initial_atoms,final_atoms,n_images)
					if not neb_conv:
						os.remove(working_cif_path)
						return None

				elif acc_level == 'cineb_lowacc' and i > 0:
					wf.run_i += 1
					continue

				elif 'dimer' in acc_level:
					mof = wf.dimer()
					if mof is None:
						os.remove(working_cif_path)
						return None

				elif acc_level == 'final_spe':
					mof = wf.final_spe()
					if mof is None:
						os.remove(working_cif_path)
						return None
					result_path = os.path.join(basepath,'results',name)
					newmodecar = os.path.join(result_path,acc_levels[i-2],spin_level,'NEWMODECAR')
					newmodecar_spe = os.path.join(result_path,acc_level,spin_level,'NEWMODECAR')
					copyfile(newmodecar,newmodecar_spe)

				else:
					raise ValueError('Unsupported accuracy level')

				if acc_level == 'dimer_lowacc' and i > 0:
					is_new_spin = check_if_new_spin(self,mof,name,acc_level,spin_level)
					if not is_new_spin:
						same_spin = True
						break

			#***********SAVE and CONTINUE***********
			if same_spin:
				continue
			E_temp = mof.get_potential_energy()
			if E_temp < E:
				best_mof = deepcopy(mof)

		os.remove(working_cif_path)
		return best_mof