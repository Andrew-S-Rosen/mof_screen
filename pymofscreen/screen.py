import os
from pymofscreen.writers import pprint
from pymofscreen.kpts_handler import get_kpts
from pymofscreen.screen_phases import workflows
from pymofscreen.janitor import prep_paths
from pymofscreen.default_calculators import calcs
import sys
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
			basepath (string): path to the base
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
			stdout_file = sys.argv[0].split('.py')[0]+'.out'
		prep_paths(basepath)

	def run_screen(self,cif_file,mode,spin_levels=None,acc_levels=None,calcs=calcs):
		"""
		Run high-throughput ionic relaxations
		Args:
			cif_file (string): name of CIF file
			spin_levels (list of strings): spin states to consider
			acc_levels (list of strings): accuracy levels to consider
			calcs (function): function to call respective calculator
		Returns:
			mofs (list of ASE Atoms objects): ASE Atoms objects for optimized
			MOF given by cif_file for each spin_level
		"""

		basepath = self.basepath
		self.spin_levels = spin_levels
		self.calcs = calcs
		if mode == 'ionic':
			if acc_levels is None:
				acc_levels = ['scf_test','isif2_lowacc','isif2_medacc',
				'isif2_highacc','final_spe']
		elif mode == 'volume':
			if acc_levels is None:
				acc_levels = ['scf_test','isif2_lowacc','isif3_lowacc',
				'isif3_highacc','isif2_highacc','final_spe']
		else:
			raise ValueError('Unsupported DFT screening mode')
		self.acc_levels = acc_levels
		if spin_levels is None:
			spin_levels = ['spin1','spin2']

		#Make sure MOF isn't running on other process
		working_cif_path = os.path.join(basepath,'working',cif_file)
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
		mofs = []

		#for each spin level, optimize the structure
		for i, spin_level in enumerate(spin_levels):

			#***********PREP FOR RUN***********
			if spin_level != spin_levels[0]:
				prior_spin = spin_levels[i-1]
			else:
				prior_spin = None
			wf = workflows(self,cif_file,kpts_dict,spin_level,prior_spin)

			for acc_level in acc_levels:

				if acc_level == 'scf_test':
					scf_pass = wf.scf_test()
					if not scf_pass:
						return None

				elif acc_level == 'isif2_lowacc':
					if mode == 'volume':
						acc_level = 'isif2' #for legacy reasons
					mof = wf.isif2_lowacc()
					if mof is None:
						return None

				elif acc_level == 'isif2_medacc':
					mof = wf.isif2_medacc(mof)
					if mof is None:
						return None

				elif acc_level == 'isif2_highacc':
					acc_level = 'final' #for legacy reasons
					mof = wf.isif2_highacc(mof)
					if mof is None:
						return None
				
				elif acc_level == 'isif3_lowacc':
					mof = wf.isif3_lowacc(mof)
					if mof is None:
						return None

				elif acc_level == 'isif3_highacc':
					mof = wf.isif3_highacc(mof)
					if mof is None:
						return None

				elif acc_level == 'final_spe':
					mof, skip_spin2 = wf.final_spe(mof)
					if mof is None:
						return None

				else:
					raise ValueError('Unsupported accuracy level')

			#***********SAVE and CONTINUE***********
			mofs.append(mof)
			if skip_spin2 == True:
				pprint('Skipping '+spin_levels[i+1]+' run')
				return mofs

		return mofs