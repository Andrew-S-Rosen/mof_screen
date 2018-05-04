import os
import numpy as np
import sys
from copy import deepcopy
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

	def __init__(self,mofpath,basepath,kpts_path='Auto',kppas=None,niggli=True,
		submit_script='sub_screen.job',stdout_file=None):
		"""
		Initialize variables that should be used on all MOFs in a database
		Args:
			mofpath (string): path to the directory containing the CIF files
			basepath (string): path to the base directory for the DFT screening
			kpts_path (string): can be either 'Auto' for an automatic generation
			of the kpoints based on KPPAs or a string representing the path to a
			text file with all the kpoint information
			kppas (list of ints): KPPAs to use if kpts_path == 'Auto'
			niggli (bool): True if Niggli-reduction should be performed
			submit_script (string): name of job submission script
			stdout_file (string): name of the stdout file
		"""
		self.mofpath = mofpath
		self.basepath = basepath
		self.submit_script = submit_script
		self.stdout_file = stdout_file
		self.kpts_path = kpts_path
		self.niggli = niggli
		if kppas is None:
			self.kppas = [100,1000]
		if stdout_file is None:
			self.stdout_file = sys.argv[0].split('.py')[0]+'.out'
		prep_paths(basepath)

	def run_screen(self,cif_file,mode,spin_levels=None,acc_levels=None,calcs=calcs):
		"""
		Run high-throughput ionic or volume relaxations
		Args:
			cif_file (string): name of CIF file
			spin_levels (list of strings): spin states to consider
			acc_levels (list of strings): accuracy levels to consider
			calcs (function): function to call respective calculator
		Returns:
			best_mof (ASE Atoms objects): ASE Atoms object for optimized
			MOF given by cif_file (lowest energy spin state)
		"""

		basepath = self.basepath
		self.calcs = calcs
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
		working_cif_path = os.path.join(basepath,'working',cif_file)
		refcode = cif_file.split('.cif')[0]
		if os.path.isfile(working_cif_path) == True:
			pprint('SKIPPED: Running on another process')
			return None

		#Get the kpoints
		kpts_lo, gamma = get_kpts(self,cif_file,'low')
		kpts_hi, gamma = get_kpts(self,cif_file,'high')
		kpts_dict = {}
		kpts_dict['kpts_lo'] = kpts_lo
		kpts_dict['kpts_hi'] = kpts_hi
		kpts_dict['gamma'] = gamma
		E = np.inf

		#for each spin level, optimize the structure
		for i, spin_level in enumerate(spin_levels):

			#***********PREP FOR RUN***********
			if spin_level != spin_levels[0]:
				prior_spin = spin_levels[i-1]
			else:
				prior_spin = None
			same_spin = False

			wf = workflows(self,cif_file,kpts_dict,spin_level,prior_spin)
			for acc_level in acc_levels:

				if acc_level == 'scf_test':
					scf_pass = wf.scf_test()
					if not scf_pass:
						return None

				elif acc_level == 'isif2_lowacc' or (acc_level == 'isif2' and mode == 'volume_legacy'):
					mof = wf.isif2_lowacc()
					if mof is None:
						return None

					if i > 0:
						is_new_spin = check_if_new_spin(self,mof,refcode,acc_level,spin_level)
						if not is_new_spin:
							same_spin = True
							break

				elif acc_level == 'isif2_medacc':
					mof = wf.isif2_medacc()
					if mof is None:
						return None

				elif acc_level == 'isif2_highacc' or (acc_level == 'final' and 'legacy' in mode):
					mof = wf.isif2_highacc()
					if mof is None:
						return None
				
				elif acc_level == 'isif3_lowacc':
					mof = wf.isif3_lowacc()
					if mof is None:
						return None

				elif acc_level == 'isif3_highacc':
					mof = wf.isif3_highacc()
					if mof is None:
						return None

				elif acc_level == 'final_spe':
					mof = wf.final_spe(newmodecar=True)
					if mof is None:
						return None

				else:
					raise ValueError('Unsupported accuracy level')
		
			#***********SAVE and CONTINUE***********
			if same_spin == True:
				continue
			E_temp = mof.get_potential_energy()
			if E_temp < E:
				best_mof = deepcopy(mof)
			if len(spin_levels) > 1:
				skip_low_spin = check_if_skip_low_spin(self,mof,refcode,spin_levels[i])
				if (spin_level == 'spin1' or spin_level == 'high_spin') and skip_low_spin == True:
					pprint('Skipping '+spin_levels[i+1]+' run')
					continue

		return best_mof

	def run_ts_screen(self,name,initial_atoms,final_atoms,n_images=6,cif_file=None,spin_levels=None,acc_levels=None,calcs=calcs):
		"""
		Run high-throughput TS calculation
		Args:
			name (string): name of CIF file
			spin_levels (list of strings): spin states to consider
			acc_levels (list of strings): accuracy levels to consider
			calcs (function): function to call respective calculator
		Returns:
			best_mof (ASE Atoms objects): ASE Atoms object for optimized
			MOF given by cif_file (lowest energy spin state)
		"""

		basepath = self.basepath
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
			# elif cif_file.split('.cif')[0] == name:
			# 	raise ValueError('Input name and name of CIF file must not be identical')

		#Make sure MOF isn't running on other process
		working_cif_path = os.path.join(basepath,'working',name)
		if os.path.isfile(working_cif_path) == True:
			pprint('SKIPPED: Running on another process')
			return None

		#Get the kpoints
		kpts_lo, gamma = get_kpts(self,cif_file,'low')
		kpts_hi, gamma = get_kpts(self,cif_file,'high')
		kpts_dict = {}
		kpts_dict['kpts_lo'] = kpts_lo
		kpts_dict['kpts_hi'] = kpts_hi
		kpts_dict['gamma'] = gamma
		E = np.inf

		#for each spin level, optimize the structure
		for i, spin_level in enumerate(spin_levels):

			#***********PREP FOR RUN***********
			if spin_level != spin_levels[0]:
				prior_spin = spin_levels[i-1]
			else:
				prior_spin = None
			same_spin = False

			wf = workflows(self,name,kpts_dict,spin_level,prior_spin)
			for acc_level in acc_levels:

				if acc_level == 'scf_test':
					scf_pass = wf.scf_test(quick_test=True)
					if not scf_pass:
						return None					

				elif acc_level == 'cineb_lowacc' and i == 0:
					neb_conv = wf.cineb_lowacc(initial_atoms,final_atoms,n_images)
					if not neb_conv:
						return None

				elif acc_level == 'cineb_lowacc' and i > 0:
					wf.run_i += 1
					continue

				elif 'dimer' in acc_level:
					mof = wf.dimer()
					if mof is None:
						return None

				elif acc_level == 'final_spe':
					mof = wf.final_spe()
					if mof is None:
						return None
				else:
					raise ValueError('Unsupported accuracy level')

				if acc_level == 'dimer_lowacc' and i > 0:
					is_new_spin = check_if_new_spin(self,mof,name,acc_level,spin_level)
					if not is_new_spin:
						same_spin = True
						break

			#***********SAVE and CONTINUE***********
			if same_spin == True:
				continue
			E_temp = mof.get_potential_energy()
			if E_temp < E:
				best_mof = deepcopy(mof)
			if len(spin_levels) > 1:
				skip_low_spin = check_if_skip_low_spin(self,mof,name,spin_levels[i])
				if (spin_level == 'spin1' or spin_level == 'high_spin') and skip_low_spin == True:
					pprint('Skipping '+spin_levels[i+1]+' run')
					continue

		return best_mof