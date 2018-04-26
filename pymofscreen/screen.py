import os
from pymofscreen.writers import pprint
from pymofscreen.kpts_handler import get_kpts
from pymofscreen.screen_phases import workflows
from pymofscreen.janitor import prep_paths
from pymofscreen.default_calculators import calcs_ads

class screener():
	"""
	This class constructs a high-throughput screening workflow
	"""

	def __init__(self,mofpath,basepath,ads_species=None,kpts_path='Auto',kppas=[100,1000],niggli=True,
		submit_script='sub_screen.job',stdout_file='driver.out'):
		"""
		Initialize variables that should be used on all MOFs in a database
		Args:
			mofpath (string): path to the directory containing the CIF files
			basepath (string): path to the base
			ads_species (string): (optional) string representing adsorbate species (e.g. 'CH4')
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
		self.ads_species = ads_species
		self.submit_script = submit_script
		self.stdout_file = stdout_file
		self.kpts_path = kpts_path
		self.niggli = niggli
		self.kppas = kppas
		prep_paths(basepath)

	def run_ads_screen(self,cif_file,spin_levels=['spin1','spin2'],acc_levels=['scf_test','isif2_lowacc',
		'isif2_medacc','final','final_spe'],calcs_ads=calcs_ads):
		"""
		Run high-throughput ionic relaxations
		Args:
			cif_file (string): name of CIF file
			spin_levels (list of strings): spin states to consider
			acc_levels (list of strings): accuracy levels to consider
			calcs_ads (function): function to call respective calculator
		Returns:
			mofs (list of ASE Atoms objects): ASE Atoms objects for optimized MOF given by
			cif_file for each spin_level
		"""

		basepath = self.basepath
		self.spin_levels = spin_levels
		self.acc_levels = acc_levels
		self.calcs_ads = calcs_ads

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

			#***********SCF TEST************
			scf_pass = wf.scf_test()
			if not scf_pass:
				return None

			#***********ISIF 2 (lowacc)************
			mof = wf.isif2_lowacc()
			if mof is None:
				return None

			#***********ISIF 2 (medacc)************
			mof = wf.isif2_medacc(mof)
			if mof is None:
				return None

			#***********ISIF 2 (final)************
			mof = wf.isif2_final(mof)
			if mof is None:
				return None
				
			#***********Final SPE************
			mof, skip_spin2 = wf.final_spe(mof)
			if mof is None:
				return None

			#***********SAVE and CONTINUE***********
			mofs.append(mof)
			if skip_spin2 == True:
				pprint('Skipping '+spin_levels[i+1]+' run')
				return mofs

		return mofs